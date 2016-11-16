#!/bin/bash

##
## NAME: Post_Medusa_Processing_V1.0.sh - Performs various post processing steps on Medusa results. 
##
## Created by Austin Christofferson on 06/03/2016.
## Copyright 2016 Translational Genomics Research Institute. All rights reserved.
##
## Warning     - This script must be run from the results directory on a box that 
##               has qsub (merckx).
##
## Description - Depending on options used this script will loop through 
##               new "*ps*" time stamped results directories renaming them, copying 
##               files over to pnap for reduing snpEFF with updated databases and 
##               tool versions as well as starting the advanced cna scripts.
##               After all pbs jobs on pnap are finished the new data will
##               automatically copy back to isilon.
##
## Usage = Post_Medusa_Processing_V1.0.sh [-help] [-rp] [-b y] [-l file.txt]
##
## Options:
##      -r      Rename "*_ps*" timestamped directories
##      -R      IN-Development: Needs directories to first go through -p
##                              Restart the post processing directories that have a failed or in 
##				progress stamp
##      -p      IN-Development: Folders must have the "*_ps*" timestamp from Medusa with one 
##				underscore seperating the study from the Patient ID (MMRF_1234_ps201605130955) 
##				or the -l option must be used with only one underscore in the name.
##				same as -r and -b [y/yes] and performes re-annotate with SnpEFF, Star, Kallisto and cna_manual
##      -b      Flag [y/yes|n/no] to keep or remove fastqc and seurat/*REVERSE 
##      -l      IN-Development: -l <file_of_directory_names.txt> Uses a list of directories to operate on.
##				Depending on the other options you use [-p|-r|-R] this could have unintended
##                              side effects.
##	-c	Used with the -R and -l option to only restart Copy Number analysis. 
##				
##
####################################################################
####################################################################


if [ "$*" == "--help" ] || [ "$*" == "-h" ] || [ "$*" == "" ]
then
        grep "^##" $0
	echo
        exit 0
fi

# Set default settings

#{{{

RENAME=0
RESTART=0
RMBLOAT=no
POSTPROCESSING=0
CNAONLY=0

while getopts ":rRb:l:pc" opt
do
        case $opt in
                r)
			RENAME=1
                        FIND() {
                                ls -d *_ps*
                                } ;;
                R) 
			RESTART=1
                        FIND() {
                                find . -maxdepth 2 -name "*Fail" -o -name "*In_Progress" | awk -F'/' '{ print $2 }' | sort | uniq
                                } ;;
                p)
			RESTART=0
			RMBLOAT=yes
			RENAME=1
			POSTPROCESSING=1
                        FIND() {
                                ls -d *_ps*
                                } ;;
                b)
                        RMBLOAT=$OPTARG
                        if [ $RMBLOAT == "y" ] || [ $RMBLOAT == "yes" ]
                        then
                                RMBLOAT=yes
                        elif [ $RMBLOAT == "n" ] || [ $RMBLOAT == "no" ]
                        then
                                RMBLOAT=no
                        else
                                echo "$RMBLOAT is not a valid argument for -b "
                                echo "Please use one of the following [ y/yes/n/no ] "
                                exit 1
                        fi ;;
		c)
			CNAONLY=1 ;;
                l)
                        DIRLIST=$OPTARG
                        if [ ! -f $DIRLIST ] ; then
                                echo
                                echo $DIRLIST is not a valid file path.
                                echo Please check the path to your input list.
                                exit 1
                        else
                                FIND() {
                                        cat $DIRLIST
                                        }
                        fi ;;
                :)
                        echo "Option -$OPTARG requires an argument." >&2
                        exit 1 ;;
                \?)
                echo "Invalid option: -$OPTARG" >&2 ;;
        esac
done

#}}}

# Post Processing Variables

#{{{

EXOMEPBS=/home/tgenref/pecan/post_central_pipe_processing/post_medusa_V1.0/jobScripts/rsync_Exome_Files.pbs
GENOMEPBS=/home/tgenref/pecan/post_central_pipe_processing/post_medusa_V1.0/jobScripts/rsync_Genome_Files.pbs
RNAPBS=/home/tgenref/pecan/post_central_pipe_processing/post_medusa_V1.0/jobScripts/rsync_RNA_Files.pbs
CONSTANTS=/home/tgenref/pecan/post_central_pipe_processing/post_medusa_V1.0/constants.txt
STARTDIR=`pwd`

#}}}

# Make pbs out directory if does not exist

#{{{

if [ ! -d "${HOME}/STD_OUT" ]
then
        mkdir ${HOME}/STD_OUT
fi

#}}}

# Remove tracking list if it exists

#{{{

if [ -f ~/Post_processing_projects_to_be_restarted.txt ]
then
	rm ~/Post_processing_projects_to_be_restarted.txt
fi

#}}}

# Start loop on patient find

for line in `FIND`
do
	PATIENT_NAME=`echo ${line} | awk -F'_' '{ OFS = "_" ; print $1,$2 }'`

	cd ${STARTDIR}	
	cd ${line}

	# Check if project has completed copying back to isilon
	
	#{{{

	if [ ! -f "copyDone.txt" ]
	then
		cd ..
                #Message
                echo
                echo "---------------------------------"
                echo "WARNING: Transfere is not complete"
                echo "${line} WILL NOT BE RENAMED!!"
		continue
	fi

	#}}}	

	# Remove Post processing directory on /scratch if performing -R or -p

	#{{{

	if [ ${RESTART} == 1 ] || [ ${POSTPROCESSING} == 1 ]
	then
		ssh ${USER}@pnap-data1.tgen.org "if [ -d /scratch/${USER}/post_medusa_processing/${PATIENT_NAME} ] ; then rm -rf /scratch/${USER}/post_medusa_processing/${PATIENT_NAME} ; fi"
        	if [ $? -ne 0 ]
        	then
        		echo
        	        echo Warning!!! ssh failed to remove the pre-existing Post Processing directory for on scratch.
        	        echo /scratch/${USER}/post_medusa_processing/${PATIENT_NAME} 
        	        echo
        	        echo Would you like to try to auto remove it again or remove it manualy [ auto/manual ]
        	        read RESPONSE
	
        	        if [ ${RESPONSE} == "auto" ]
        	        then
        	        	ssh ${USER}@pnap-data1.tgen.org "if [ -d /scratch/${USER}/post_medusa_processing/${PATIENT_NAME} ] ; then rm -rf /scratch/${USER}/post_medusa_processing/${PATIENT_NAME} ; fi"
	
        	                if [ $? -ne 0 ]
        	                then
        	                	echo
        	                	echo Auto rm post processing directory failed again. Moving on. Please check ${PATIENT_NAME}!! 
        	                fi
        	        elif [ ${RESPONSE} == "manual" ]
        	        then
				echo
				echo Press enter after you have deleted the folder on scratch.
				read PRESS
        	        else
				echo Response was not auto or manual. exiting script.
				exit 1
        	        fi
        	fi
	fi

	#}}}

	# Remove all progress tags and analysis diretories 	

	#{{{

	if [ ${RESTART} == 1 ]
	then
		# Remove cna_manual directory
		if [ -d cna_manual_2016 ] 
		then
			rm -rf cna_manual_2016
		fi

		# Remove CNA_Manual progress tag
		if [ -f CNA_Manual_In_Progress ]
		then
			rm CNA_Manual_In_Progress
		fi
		
		if [ -f CNA_Manual_Fail ]
		then
			rm CNA_Manual_Fail 
		fi
		
		if [ -f CNA_Manual_Complete ]
		then
			rm CNA_Manual_Complete
		fi	

		# Remove SnpEFF progress tags
		
		if [ ${CNAONLY} = 0 ]
		then
			find . -name SnpEFF_ANN_Complete -exec rm {} \;
			find . -name SnpEFF_ANN_In_Progress -exec rm {} \;
			find . -name SnpEFF_ANN_Fail -exec rm {} \;
		fi		

		# Remove salmon and Kallisto tags
		if [ ${CNAONLY} = 0 ]
                then
			find . -name Salmon_Complete -exec rm {} \;
			find . -name Kallisto_Complete -exec rm {} \;
			find . -name Salmon_In_Progress -exec rm {} \;
			find . -name Kallisto_In_Progress -exec rm {} \;
			find . -name Salmon_Fail -exec rm {} \;
			find . -name Kallisto_Fail -exec rm {} \;
		fi
	fi

	#}}}

	# Remove bloat folders that keats lab does not use or is now redundant

	#{{{	

	if [ ${RMBLOAT} == "yes" ]
	then
		# Remove all fastqc directories and files
                find . -name *_fastqc -type d -exec rm -rf {} +

                # Remove all seurat REVERSE dir
                if [ "` ls seurat/ |grep REVERSE |wc -l`" -gt 0 ]
                then
                        rm -rf seurat/*REVERSE
                fi
	fi

	#}}}	

	# Rename Directories to make compatible with KBase import

	#{{{

	if [ ${RENAME} = 1 ]
	then
		cd ..		
		if [ -d ${PATIENT_NAME} ]
                then
                        echo
                        echo "---------------------------------"
                        echo "Warning: ${PATIENT_NAME} already exists"
                        echo
                        DATE="`date +%Y%m%d-%H:%M`"
                        echo "Adding time stamp to older folder"
                        echo "Old folder name=${PATIENT_NAME}"
                        echo "New folder name=${PATIENT_NAME}_TS${DATE}"
                        mv ${PATIENT_NAME} ${PATIENT_NAME}_ts${DATE}

                        #Rename files
                        mv ${line} ${PATIENT_NAME}

                        #Message
                        echo
                        echo "---------------------------------"
                        echo "RENAMING: ${line}"
                        echo Old Name=${line}
                        echo New Name=${PATIENT_NAME}
                else
                        #Rename files
                        mv ${line} ${PATIENT_NAME}

                        #Message
                        echo
                        echo "---------------------------------"
                        echo "RENAMING: ${line}"
                        echo Old Name=${line}
                        echo New Name=${PATIENT_NAME}
                fi
		cd ${PATIENT_NAME}
	fi

	#}}}

	# Start PostProcessing or Restart if flag was set

	if [ ${RESTART} = 1 ] || [ ${POSTPROCESSING} = 1 ]
	then
		# Check for any DNAPAIR= lines and start appropriate snpEff and CNA jobs
		
		DNAPAIRCOUNT="`grep "DNAPAIR=" ${PATIENT_NAME}.config | wc -l`"

		if [ ${DNAPAIRCOUNT} = 0 ]
	        then
	        	echo ${PATIENT_NAME} has no DNA pair lines. Skipping snpEff and CNA launcher.
	                touch CNA_Manual_Complete
			
			if [ ${CNAONLY} = 0 ]
			then
	                	touch SnpEFF_ANN_Complete
			fi
	        else
	                # Find out if there are any exome DNA pair lines. If none then skip snpEff and CNA launcher

			#{{{

	                EXOMECOUNT=0
	                EXOMEPAIRS=()
	                EXPECTEDPAIRS=()

	                for PAIRTEST in `grep "DNAPAIR=" ${PATIENT_NAME}.config`
	                do
	                	NORM=`echo ${PAIRTEST} | awk -F"[=,]" '{print $2}'`
	                        TUMR=`echo ${PAIRTEST} | awk -F"[=,]" '{print $3}'`
	
	                        NORMTYPE=`grep "SAMPLE=.*${NORM}" ${PATIENT_NAME}.config | awk -F',' '{ print $3 }'`
	                        TUMRTYPE=`grep "SAMPLE=.*${TUMR}" ${PATIENT_NAME}.config | awk -F',' '{ print $3 }'`

	                        if [ ${NORMTYPE} == "Exome" ] && [ ${TUMRTYPE} == "Exome" ]
	                        then
	                        	((EXOMECOUNT+=1))
	                                EXOMEPAIRS+=("${NORM}-${TUMR}")
	                                EXPECTEDPAIRS+=("${NORM}-${TUMR}_exo")
	                        fi
	                done

			#}}}

	                # Make array of exome pairs and their corresponding genome

			#{{{

	                for COMBINED in `echo ${EXOMEPAIRS[@]} | sed 's/\s/\n/g'`
	                do
        	        	TUMORSAMPLE="`echo ${COMBINED} | cut -d- -f2 `"
        	                TUMORSHORT="`echo ${TUMORSAMPLE} | cut -d_ -f1,2,3,4,5`"

        	                for PAIRTEST in `grep "DNAPAIR=" ${PATIENT_NAME}.config | grep ${TUMORSHORT}`
        	                do
                	        	NORM=`echo ${PAIRTEST} | awk -F"[=,]" '{print $2}'`
        	                        TUMR=`echo ${PAIRTEST} | awk -F"[=,]" '{print $3}'`

                	                NORMTYPE=`grep "SAMPLE=.*${NORM}" ${PATIENT_NAME}.config | awk -F',' '{ print $3 }'`
                	                TUMRTYPE=`grep "SAMPLE=.*${TUMR}" ${PATIENT_NAME}.config | awk -F',' '{ print $3 }'`

                        	        if [ ${NORMTYPE} = "Genome" ] && [ ${TUMRTYPE} = "Genome" ]
                	                then
                        	        	EXPECTEDPAIRS+=("${NORM}-${TUMR}_filt2012")
						EXPECTEDPAIRS+=("${NORM}-${TUMR}_filt2016")
                        	                EXPECTEDPAIRS+=("${NORM}-${TUMR}_unfi")
                        	        fi
                        	done
                	done		
			
			#}}}

			if [ ${EXOMECOUNT} = 0 ]
                	then
                		echo ${PATIENT_NAME} has no Exome DNA pair lines. Skipping snpEff and CNA launcher.
                	        touch CNA_Manual_Complete
			
				if [ ${CNAONLY} = 0 ]
				then
                	        	touch SnpEFF_ANN_Complete
				fi
                	else
                	        mkdir cna_manual_2016
                	        touch CNA_Manual_In_Progress
                	        
				if [ ${CNAONLY} = 0 ]
				then
					touch SnpEFF_ANN_In_Progress
				fi
		
				# Tracking for merged vcfs that have already gone through re-annotation or corrupt file.
                      
				#{{{
 
			 	MERGEDVCFCOMPLETED=0
                	        for PAIR in `echo ${EXOMEPAIRS[@]} | sed 's/\s/\n/g'`
                        	do
                	                NORMALSAMPLE="`echo ${PAIR} | cut -d- -f1`"
                        	        TUMORSAMPLE="`echo ${PAIR} | cut -d- -f2`"
                        	        EXOMEPAIR="${NORMALSAMPLE}-${TUMORSAMPLE}"

                        	        # Test if their are Original merged vcfs and if they have already already gone through the SnpEFF ANN anotation
                        	        
					if [ ${CNAONLY} = 0 ]
					then
						if [ -f vcfMerger/${EXOMEPAIR}/${EXOMEPAIR}.merged.all*Original ]
						then
							ORIGTEST=`ls vcfMerger/${EXOMEPAIR}/${EXOMEPAIR}*final.vcf.Original | wc -l `
                        	  		else
							ORIGTEST=0      
						fi
			
						if [ ${ORIGTEST} = 2 ]
                        	        	then
                                		        for MERGEDVCFORIGINAL in `ls vcfMerger/${EXOMEPAIR}/${EXOMEPAIR}*final.vcf.Original `
                        	        	        do
                                		                BASENAME=`basename ${MERGEDVCFORIGINAL} .Original`
                                		                rsync ${MERGEDVCFORIGINAL} vcfMerger/${EXOMEPAIR}/${BASENAME}
                                		        done
                                		else
                                		        ANNTEST=`grep "##INFO=<ID=ANN" vcfMerger/${EXOMEPAIR}/${EXOMEPAIR}*final.vcf | wc -l`
                                		        EFFTEST=`grep "##INFO=<ID=EFF" vcfMerger/${EXOMEPAIR}/${EXOMEPAIR}*final.vcf | wc -l`

                                		        if [ ${ANNTEST} = 2 ]
                                		        then
                                	        	        echo
                                	        	        echo Warning!!! The vcf for ${PATIENT_NAME} has already been gone through re-annotation and does not have an Original.
                                	        	        echo Skipping all steps for this patient.
                                	        	        echo
                                	        	        echo You should restart ${PATIENT_NAME} from KBase!!
                                	        	        echo
                                	        	        MERGEDVCFCOMPLETED=1
                                	        	        continue

                                	        	elif [ ${EFFTEST} = 2 ]
							then
								echo Original ${EXOMEPAIR} VCF passed.	
							else
								echo
								echo Warning!!! The vcf for ${PATIENT_NAME} does not have an Original or a re-annotated version.
								echo Skipping all steps for this patient.
                                	        	        echo
                                	        	        echo You should restart ${PATIENT_NAME} from KBase!!
                                	        	        echo
                                	        	        MERGEDVCFCOMPLETED=1
                                	        	        continue
                                	        	fi
                                		fi
					fi
                        	done

				#}}}

				# If one of the exome pairs has already gone through re-annotation then exit the the loop for this patient.		

				#{{{
	
				if [ ${MERGEDVCFCOMPLETED} = 1 ]
				then
					echo ${PATIENT_NAME} >> ~/Post_processing_projects_to_be_restarted.txt
					continue
				fi

				#}}}
	
				# loop through each DNA pair line Capturing needed variables and Start snpEff and CNA launcher for exomes       

                        	EXOMEPAIRSLIST=`echo ${EXOMEPAIRS[@]} | sed 's/\s/@/g'`
                        	EXPECTEDPAIRSLIST=`echo ${EXPECTEDPAIRS[@]} | sed 's/\s/@/g'`

                        	for PAIR in `echo ${EXOMEPAIRS[@]} | sed 's/\s/\n/g'`
                        	do
                        		echo $PAIR
                        	        LIBTYPE="Exome"
                        	        NORMALSAMPLE="`echo ${PAIR} | cut -d- -f1`"
                        	        TUMORSAMPLE="`echo ${PAIR} | cut -d- -f2`"
                        	        EXOMEPAIR="${NORMALSAMPLE}-${TUMORSAMPLE}"
                        	        VCF=${EXOMEPAIR}.HC_All.snpEff.vcf
                        	        NORMALDAT=${NORMALSAMPLE}.proj.md.jr.bam.clc.cln.dat
                        	        TUMORDAT=${TUMORSAMPLE}.proj.md.jr.bam.clc.cln.dat
                        	        TUMORSHORT="`echo ${TUMORSAMPLE} | cut -d_ -f1,2,3,4,5`"
                        	        NASSAY="`echo ${NORMALSAMPLE} | cut -d_ -f7`"
                        	        TASSAY="`echo ${TUMORSAMPLE} | cut -d_ -f7`"
                        	        CNAEXOMETARGET=`grep "${NASSAY}_CNABEDP=" ${CONSTANTS} | cut -d= -f2 | tr -d '\n'`

					if [[ "$TASSAY" == "TSE61" ]] ; then
                				bedFile="/home/tgenref/pipeline_v0.3/annotations/exome_capture/illumina_truseq/TruSeq_exome_targeted_regions_b37_padded.bed"
        				elif [[ "$TASSAY" == *S5U ]] || [[ "$kitName" == *S5X ]] ; then
                				bedFile="/home/tgenref/pipeline_v0.3/ensembl70/Ensembl_v70_hs37d5_exonic_coordinates_touched_v5UTR_padded25.bed"
        				elif [[ "$TASSAY" == *STX ]] ; then
                				bedFile="/home/tgenref/pipeline_v0.4/annotations/exome_capture/strexome/Strexome_targets_intersect_sorted_padded100.bed"
        				elif [[ "$TASSAY" == *SCR ]] ; then
                				bedFile="/home/tgenref/pecan/annotations/exome_capture/agilent_clinical_research_exome/Agilent_Clinical_Research_Exome_hs37d5_PaddedTargets_intersect_sorted_padded100.bed"
        				elif [[ "$TASSAY" == *S2X ]] ; then
                				bedFile="/home/tgenref/pipeline_v0.4/annotations/exome_capture/Agilent_V2_hs37d5/Agilent_V2_hs37d5_Targets_intersect_sorted_padded100.bed"
        				elif [[ "$TASSAY" == *STL ]] ; then
                				bedFile="/home/tgenref/pecan/annotations/exome_capture/strexome_lite/Strexome_Lite_Targets_intersect_sorted_padded100.bed"
        				elif [[ "$TASSAY" == *S1X ]] ; then
                				bedFile="/home/tgenref/pecan/annotations/exome_capture/agilent_SureSelectV1/SureSelectV1_hs37d5_PaddedTargets_intersect_sorted_padded100.bed"
        				elif [[ "$TASSAY" == *S6X ]] ; then
                				bedFile="/home/tgenref/pecan/annotations/exome_capture/agilent_v6_noUTR/Agilent_V6_noUTR_hs37d5_Targets_intersect_sorted_padded100.bed"
        				elif [[ "$TASSAY" == *SXP ]] ; then
                				bedFile="/home/tgenref/pecan/annotations/exome_capture/Agilent_SureSelect_V6R2_plusUTR/Agilent_SureSelect_V6R2_plusUTR_hs37d5_GRCh37.74_PaddedTargets_intersect_sorted_padded100.bed"
        				elif [[ "$TASSAY" == *S4X ]] ; then
                				bedFile="/home/tgenref/pecan/annotations/exome_capture/agilent_v4_noUTR/Agilent_V4_noUTR_hs37d5_Targets_intersect_sorted_padded100.bed"
        				elif [[ "$TASSAY" == *E62 ]] ; then
                				bedFile="/home/tgenref/pecan/annotations/exome_capture/illumina_nextera_expanded/NexteraExpandedExome_hs37d5_Targets_PicardPadded100.bed"
        				fi
	
					mkdir cna_manual_2016/${EXOMEPAIR}_exo
        	                        touch cna_manual_2016/${EXOMEPAIR}_exo/CNA_Manual_In_Progress
        	                        
					if [ ${CNAONLY} = 0 ]
					then
						touch vcfMerger/${EXOMEPAIR}/SnpEFF_ANN_In_Progress
					fi

        	                        qsub -v CNAONLY="${CNAONLY}",BEDFILE=${bedFile},PBSNAME="${EXOMEPAIR}",STARTDIR="${STARTDIR}",CNAEXOMETARGET="${CNAEXOMETARGET}",NASSAY="${NASSAY}",TASSAY="${TASSAY}",LIBTYPE="${LIBTYPE}",USER="${USER}",PATIENT_NAME="${PATIENT_NAME}",NORMALSAMPLE="${NORMALSAMPLE}",TUMORSAMPLE="${TUMORSAMPLE}",VCF="${VCF}",NORMALDAT="${NORMALDAT}",TUMORDAT="${TUMORDAT}",EXOMEPAIR="${EXOMEPAIR}",EXOMEPAIRS=${EXOMEPAIRSLIST},EXPECTEDPAIRS=${EXPECTEDPAIRSLIST} ${EXOMEPBS}

        	                        if [ $? -ne 0 ]
        	                        then
        	                        	rm CNA_Manual_In_Progress
						rm cna_manual_2016/${EXOMEPAIR}_exo/CNA_Manual_In_Progress
						echo Failed to start rsync qsub job for ${EXOMEPAIR} >> CNA_Manual_Fail
						echo Failed to start rsync qsub job for ${EXOMEPAIR} >> cna_manual_2016/${EXOMEPAIR}_exo/CNA_Manual_Fail

						if [ ${CNAONLY} = 0 ]
						then
        	                                	rm SnpEFF_ANN_In_Progress
        	                                	rm vcfMerger/${EXOMEPAIR}/SnpEFF_ANN_In_Progress
        	                                	echo Failed to start rsync qsub job for ${EXOMEPAIR} >> SnpEFF_ANN_Fail
        	                                	echo Failed to start rsync qsub job for ${EXOMEPAIR} >> vcfMerger/${EXOMEPAIR}/SnpEFF_ANN_Fail
        	                        	fi
					fi

					# Start LI CNA if available

        	                        LICOUNT=0
        	                        LIPAIRS=()
	
	                                for PAIRTEST in `grep "DNAPAIR=" ${PATIENT_NAME}.config | grep ${TUMORSHORT}`
	                                do
	                                	NORM=`echo ${PAIRTEST} | awk -F"[=,]" '{print $2}'`
	                                        TUMR=`echo ${PAIRTEST} | awk -F"[=,]" '{print $3}'`

	                                        NORMTYPE=`grep "SAMPLE=.*${NORM}" ${PATIENT_NAME}.config | awk -F',' '{ print $3 }'`
	                                        TUMRTYPE=`grep "SAMPLE=.*${TUMR}" ${PATIENT_NAME}.config | awk -F',' '{ print $3 }'`

	                                        if [ ${NORMTYPE} = "Genome" ] && [ ${TUMRTYPE} = "Genome" ]
	                                        then
	                                        	((LICOUNT+=1))
	                                                LIPAIRS+=("${NORM}-${TUMR}")
	                                        fi
	                                done

	                                if [ ${LICOUNT} -eq 1 ]
	                                then
	                                	for GPAIR in `echo ${LIPAIRS[@]} | sed 's/\s/\n/g'`
	                                        do
	                                        	LIBTYPE="Genome"
	                                                NORMALSAMPLEG="`echo ${GPAIR} | cut -d- -f1`"
	                                                TUMORSAMPLEG="`echo ${GPAIR} | cut -d- -f2`"
	                                                NORMALDAT=${NORMALSAMPLEG}.proj.md.jr.bam.clc.cln.dat
	                                                TUMORDAT=${TUMORSAMPLEG}.proj.md.jr.bam.clc.cln.dat
	                                                GENOMEPAIR="${NORMALSAMPLEG}-${TUMORSAMPLEG}"
	                                                NASSAY="`echo ${NORMALSAMPLEG} | cut -d_ -f7`"
	                                                TASSAY="`echo ${TUMORSAMPLEG} | cut -d_ -f7`"
	
        	                                        mkdir cna_manual_2016/${GENOMEPAIR}_filt2012
        	                                        touch cna_manual_2016/${GENOMEPAIR}_filt2012/CNA_Manual_In_Progress
							mkdir cna_manual_2016/${GENOMEPAIR}_filt2016
                                                        touch cna_manual_2016/${GENOMEPAIR}_filt2016/CNA_Manual_In_Progress
        	                                        mkdir cna_manual_2016/${GENOMEPAIR}_unfi
        	                                        touch cna_manual_2016/${GENOMEPAIR}_unfi/CNA_Manual_In_Progress
	
        	                                        qsub -v PBSNAME=${GENOMEPAIR},BEDFILE=${bedFile},STARTDIR="${STARTDIR}",GENOMEPAIR=${GENOMEPAIR},EXOMEPAIR=${EXOMEPAIR},NORMALSAMPLEG="${NORMALSAMPLEG}",TUMORSAMPLEG="${TUMORSAMPLEG}",NASSAY="${NASSAY}",TASSAY="${TASSAY}",LIBTYPE="${LIBTYPE}",PATIENT_NAME="${PATIENT_NAME}",NORMALSAMPLE="${NORMALSAMPLE}",TUMORSAMPLE="${TUMORSAMPLE}",VCF="${VCF}",NORMALDAT="${NORMALDAT}",TUMORDAT="${TUMORDAT}",EXPECTEDPAIRS=${EXPECTEDPAIRSLIST} ${GENOMEPBS}

        	                                        if [ $? -ne 0 ]
        	                                        then
        	                                        	rm CNA_Manual_In_Progress
        	                                                rm cna_manual_2016/${GENOMEPAIR}_filt2012/CNA_Manual_In_Progress
								rm cna_manual_2016/${GENOMEPAIR}_filt2016/CNA_Manual_In_Progress
        	                                                rm cna_manual_2016/${GENOMEPAIR}_unfi/CNA_Manual_In_Progress
        	                                                echo Failed to start rsync qsub job for ${GENOMEPAIR} >> CNA_Manual_Fail
        	                                                echo Failed to start rsync qsub job for ${GENOMEPAIR} >> cna_manual_2016/${GENOMEPAIR}_filt2012/CNA_Manual_Fail
								echo Failed to start rsync qsub job for ${GENOMEPAIR} >> cna_manual_2016/${GENOMEPAIR}_filt2016/CNA_Manual_Fail
        	                                                echo Failed to start rsync qsub job for ${GENOMEPAIR} >> cna_manual_2016/${GENOMEPAIR}_unfi/CNA_Manual_Fail
        	                                        fi
        	                                done
        	                        elif [ ${LICOUNT} -ne 0 ]
        	                        then
        	                        	echo There is more than one LI for the Tumor Isolation.
        	                                echo "There is more than one LI for the Tumor Isolation." >> CNA_Manual_Fail
        	                                echo "There is more than one LI for the Tumor Isolation. ${TUMORSHORT}" | mailx -s "LI_cna_manual_failed" ${USER}@tgen.org 
        	                        fi
        	                done
        	        fi
        	fi
	
	
		# Start Salmon and Kallisto if RNA is available

		if [ ${CNAONLY} = 0 ]
		then

        	RNACOUNT="`grep "SAMPLE=.*,RNA," ${PATIENT_NAME}.config | wc -l`"
	
	        if [ ${RNACOUNT} = 0 ]
	        then
	        	echo ${PATIENT_NAME} has no RNA.
	                touch Salmon_Complete
	                touch Kallisto_Complete
		else
	        	touch Salmon_In_Progress
	                touch Kallisto_In_Progress
	
	                # Make list of all RNA samples for final copy back check

	                RNACHECK=()

	                for RNALIST in `grep "SAMPLE=.*,RNA," ${PATIENT_NAME}.config`
	                do
	        	        ASSAY="`echo ${RNALIST} | awk -F"[=,]" '{print $2}'`"
	                        RNASAMPLEID="`echo ${RNALIST} | awk -F"[=,]" '{print $3}'`"
	                        RNACHECK+=("${ASSAY}/${RNASAMPLEID}/${RNASAMPLEID}")
	                done

	                RNACHECKLIST=`echo ${RNACHECK[@]} | sed 's/\s/@/g'`

	                echo $RNACHECKLIST

	                for RNASAMPLE in `grep "SAMPLE=.*,RNA," ${PATIENT_NAME}.config`
	                do
	         		LIBTYPE="RNA"
	                        ASSAY="`echo ${RNASAMPLE} | awk -F"[=,]" '{print $2}'`"
	                       	STUDY="`echo ${PATIENT_NAME} | awk -F'_' '{print $1}'`"
	                       	RNASAMPLEID="`echo ${RNASAMPLE} | awk -F"[=,]" '{print $3}'`"
				
				if [ ${RESTART} == 1 ]
				then
					if [ -d ${ASSAY}/${RNASAMPLEID}/${RNASAMPLEID}.salmonDir ]
					then
				       		rm -rf ${ASSAY}/${RNASAMPLEID}/${RNASAMPLEID}.salmonDir	                       
					fi

					if [ -d ${ASSAY}/${RNASAMPLEID}/${RNASAMPLEID}.kallistoDir ]
					then
						rm -rf ${ASSAY}/${RNASAMPLEID}/${RNASAMPLEID}.kallistoDir
					fi
 				fi

			        mkdir -p ${ASSAY}/${RNASAMPLEID}/${RNASAMPLEID}.salmonDir/{ensembl74_cDNA,ensembl74_GTF}
	                        mkdir -p ${ASSAY}/${RNASAMPLEID}/${RNASAMPLEID}.kallistoDir/{ensembl74_cDNA,ensembl74_GTF}
	                        touch ${ASSAY}/${RNASAMPLEID}/${RNASAMPLEID}.salmonDir/{ensembl74_cDNA,ensembl74_GTF}/Salmon_In_Progress
	                        touch ${ASSAY}/${RNASAMPLEID}/${RNASAMPLEID}.kallistoDir/{ensembl74_cDNA,ensembl74_GTF}/Kallisto_In_Progress

	                        # Will need code to get fastq directory and library type of RNA library for salmon      

	                        FASTQDIR="`grep "${STUDY}_STUDY=" ${CONSTANTS} | cut -d= -f2 | tr -d '\n'`"
	                        RNATYPE="`grep "${ASSAY}_SALMONlibType=" ${CONSTANTS} | cut -d= -f2 | tr -d '\n'`"

	                        qsub -v PBSNAME=${RNASAMPLEID},STARTDIR="${STARTDIR}",RNATYPE="${RNATYPE}",USER="${USER}",ASSAY="${ASSAY}",FASTQDIR="${FASTQDIR}",RNASAMPLEID="${RNASAMPLEID}",PATIENT_NAME="${PATIENT_NAME}",LIBTYPE="${LIBTYPE}",RNACHECK=${RNACHECKLIST} ${RNAPBS}

	                        if [ $? -ne 0 ]
	                        then
	                		rm Salmon_In_Progress
	                                rm Kallisto_In_Progress
	                                rm ${ASSAY}/${RNASAMPLEID}/${RNASAMPLEID}.salmonDir/{ensembl74_cDNA,ensembl74_GTF}/Salmon_In_Progress
	                                rm ${ASSAY}/${RNASAMPLEID}/${RNASAMPLEID}.kallistoDir/{ensembl74_cDNA,ensembl74_GTF}/Kallisto_In_Progress
	                                echo Failed to start rsync qsub job for ${RNASAMPLEID} >> Salmon_Fail
	                                echo Failed to start rsync qsub job for ${RNASAMPLEID} >> Kallisto_Fail
	                                echo Failed to start rsync qsub job for ${RNASAMPLEID} | tee -a ${ASSAY}/${RNASAMPLEID}/${RNASAMPLEID}.salmonDir/{ensembl74_cDNA,ensembl74_GTF}/Salmon_Fail
	                                echo Failed to start rsync qsub job for ${RNASAMPLEID} | tee -a ${ASSAY}/${RNASAMPLEID}/${RNASAMPLEID}.kallistoDir/{ensembl74_cDNA,ensembl74_GTF}/Kallisto_Fail
	                        fi
	                done
		fi
		fi
	fi
done





