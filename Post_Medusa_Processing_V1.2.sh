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
## Description - Depending on options used, this script will loop through 
##               new "*ps*" time stamped results directories renaming them, copying 
##               files over to dback for reduing snpEFF with updated databases and 
##               tool versions as well as starting the advanced cna scripts.
##               After all pbs jobs on dback are finished the new data will
##               automatically copy back to isilon.
##
## Usage = Post_Medusa_Processing_V1.0.sh [--help] [-rp] [-b y] [-l file.txt]
##
## Options:
##	-r	Rename "*_ps*" timestamped directories
##	-R	IN-Development: Needs directories to first go through -p
##					 Restart the post processing directories that have a failed or in 
##					 progress stamp
##	-p	IN-Development: Folders must have the "*_ps*" timestamp from Medusa with one 
##					 underscore seperating the study from the Patient ID (MMRF_1234_ps201605130955) 
##					 or the -l option must be used with only one underscore in the name.
##					 same as -r and -b [y/yes] and performes re-annotate with SnpEFF, Star, Kallisto and cna_manual
##	-b	Flag [y/yes|n/no] to keep or remove fastqc and seurat/*REVERSE 
##	-l	IN-Development: -l <file_of_directory_names.txt> Uses a list of directories to operate on.
##					 Depending on the other options you use [-p|-r|-R] this could have unintended
##					 side effects.
##	-d	Number [1|2|3] Select which data mover to use for rsyncing files to dback. Default is 3.
##	-D	Number [123|12|23|13] Cycle the data mover used after each patient.
##	-c	Used with the -R and -l option to only restart Copy Number analysis. 
##	-m	Used with the -R and -l option to only restart vcfMerger using the pegasus pipe protocol.
##	-t	Used with the -R and -l option to only restart TopHat Fusion.
##	-f	Used with the -R and -l option to only restart digar/fusion validator
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
MERGERONLY=0
TOPHATONLY=0
FUSIONONLY=0
DATAMOVER=dback-data3.tgen.org
DATAMOVERIP=10.48.73.18
CYCLEDM=No

while getopts ":rRb:l:pcmtd:D:f" opt
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
		m)
			MERGERONLY=1 ;;
		t)
			TOPHATONLY=1 ;;
		d)
			DATAMOVER=$OPTARG
			if [ ${DATAMOVER} == "1" ] || [ ${DATAMOVER} == "2" ] || [ ${DATAMOVER} == "3" ] 
			then
				if [ ${DATAMOVER} == "1" ]
				then
					DATAMOVERIP=10.48.73.16
				elif [ ${DATAMOVER} == "2" ]
				then
					DATAMOVERIP=10.48.73.17
				elif [ ${DATAMOVER} == "3" ]
				then
					DATAMOVERIP=10.48.73.18
				fi
				DATAMOVER=dback-data${DATAMOVER}.tgen.org
			else
				echo "${DATAMOVER} is not a valid data mover."
				echo "Please use one of the following [ 1/2/3 ]"
				exit 1
			fi ;;
		D)
			CYCLEDATAMOVER=$OPTARG
			CYCLEDM="Yes"
			if [ ${CYCLEDATAMOVER} == "123" ] || [ ${CYCLEDATAMOVER} == "12" ] || [ ${CYCLEDATAMOVER} == "23" ] || [ ${CYCLEDATAMOVER} == "13" ]
			then
				DMCOUNT=`echo $CYCLEDATAMOVER | wc -c | awk '{ print $1 - 1 }'`
			else
				echo "${CYCLEDATAMOVER} is not a valid cycle strategy for the data movers."
				echo "Please use one of the following [ 123/12/23/13 ]"
				exit 1
			fi ;;
		f)
			FUSIONONLY=1 ;;
		l)
			DIRLIST=$OPTARG
			if [ ! -f $DIRLIST ]
			then
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

# Test if more than one ONLY option is used

ONLYSUM=`echo $(($MERGERONLY + $CNAONLY + $TOPHATONLY + $FUSIONONLY))`

if [ ${ONLYSUM} -gt 1 ] 
then
	echo
	echo You can only use one ONLY option at a time.
	echo Please choose which one you want to use and try again.
	exit 1
fi

#}}}   

# Post Processing Variables

#{{{

BASEDIR=$(dirname $(readlink -f "$0"))

EXOMEPBS=${BASEDIR}/jobScripts/rsync_Exome_Files.pbs
GENOMEPBS=${BASEDIR}/jobScripts/rsync_Genome_Files.pbs
RNAPBS=${BASEDIR}/jobScripts/rsync_RNA_Files.pbs
CONSTANTS=${BASEDIR}/constants.txt
STARTDIR=`pwd`

#set Data mover cycles if used

if [ ${CYCLEDM} == "Yes" ]
then
	if [ ${DMCOUNT} == "3" ]
	then
		DATAMOVERNUMBER=3
	else
		DATAMOVERNUMBER=2
	fi
fi


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

	# Set datamover to use if cycling

	#{{{
	if [ ${CYCLEDM} == "Yes" ]
	then
		if [ ${DMCOUNT} == "3" ]
		then
			if [ ${DATAMOVERNUMBER} == "3" ]
			then
				DATAMOVERNUMBER=1
				DATAMOVER=dback-data${DATAMOVERNUMBER}.tgen.org
				DATAMOVERIP=10.48.73.16
			else
				DATAMOVERNUMBER=`echo $DATAMOVERNUMBER | awk '{ print $1 + 1 }'`
				DATAMOVERNUMBER2=`echo $CYCLEDATAMOVER | cut -c$DATAMOVERNUMBER`
				if [ $DATAMOVERNUMBER2 == "2" ]
				then
					DATAMOVER=dback-data${DATAMOVERNUMBER2}.tgen.org
					DATAMOVERIP=10.48.73.17
				else
					DATAMOVER=dback-data${DATAMOVERNUMBER2}.tgen.org
					DATAMOVERIP=10.48.73.18
				fi
				
			fi
		else
			if [ ${DATAMOVERNUMBER} == "2" ]
			then
				DATAMOVERNUMBER=1
				DATAMOVERNUMBER2=`echo $CYCLEDATAMOVER | cut -c$DATAMOVERNUMBER`
				if [ $DATAMOVERNUMBER2 == "1" ]
				then
					DATAMOVER=dback-data${DATAMOVERNUMBER2}.tgen.org
					DATAMOVERIP=10.48.73.16
				elif [ $DATAMOVERNUMBER2 == "2" ]
				then
					DATAMOVER=dback-data${DATAMOVERNUMBER2}.tgen.org
					DATAMOVERIP=10.48.73.17
				fi
		
			else
				DATAMOVERNUMBER=2
				DATAMOVERNUMBER2=`echo $CYCLEDATAMOVER | cut -c$DATAMOVERNUMBER`
				if [ $DATAMOVERNUMBER2 == "2" ]
				then
					DATAMOVER=dback-data${DATAMOVERNUMBER2}.tgen.org
					DATAMOVERIP=10.48.73.17
				elif [ $DATAMOVERNUMBER2 == "3" ]
				then
					DATAMOVER=dback-data${DATAMOVERNUMBER2}.tgen.org
					DATAMOVERIP=10.48.73.18
				fi
			fi
		fi
	fi
	#}}}


	# Remove Post processing directory on /scratch if performing -R or -p

	#{{{

	if [ ${RESTART} == 1 ] || [ ${POSTPROCESSING} == 1 ]
	then
		ssh ${USER}@${DATAMOVER} "if [ -d /scratch/${USER}/post_medusa_processing/${PATIENT_NAME} ] ; then rm -rf /scratch/${USER}/post_medusa_processing/${PATIENT_NAME} ; fi"
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
				ssh ${USER}@${DATAMOVER} "if [ -d /scratch/${USER}/post_medusa_processing/${PATIENT_NAME} ] ; then rm -rf /scratch/${USER}/post_medusa_processing/${PATIENT_NAME} ; fi"
	
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

	if [ ${RESTART} == 1 ] && [ ${MERGERONLY} == 0 ] && [ ${TOPHATONLY} == 0 ] && [ ${FUSIONONLY} == 0 ]
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
			if [ -d vcfMerger_pegasus ]
			then
				rm -rf vcfMerger_pegasus
			fi
			
			if [ -f SnpEFF_ANN_Complete ]
			then
				rm SnpEFF_ANN_Complete
			fi 

			if [ -f SnpEFF_ANN_In_Progress ]
			then
				rm SnpEFF_ANN_In_Progress
			fi

			if [ -f SnpEFF_ANN_Fail ]
			then
				rm SnpEFF_ANN_Fail
			fi
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
		
		# Remove delly directory and progress tags

		if [ ${CNAONLY} = 0 ]
		then
			if [ -d delly ]
			then
				rm -rf delly
			fi
		
			if [ -f Delly_In_Progress ]
			then
				rm -rf Delly_In_Progress
			fi
	
			if [ -f Delly_Fail ]
			then
				rm -rf Delly_Fail
			fi
			
			if [ -f Delly_Complete ]
			then
				rm -rf Delly_Complete
			fi
		fi

		# Remove Digar/validation directories
		if [ ${CNAONLY} = 0 ]
		then
			if [ -d fusionValidator ]
			then
				rm -rf fusionValidator
			fi

			if [ -f FusionValidator_In_Progress ]
			then
				rm FusionValidator_In_Progress
			fi

			if [ -f FusionValidator_Fail ]
			then
				rm FusionValidator_Fail
			fi

			if [ -f FusionValidator_Complete ]
			then
				rm FusionValidator_Complete
			fi
		fi
	fi

	# Remove vcf merger results if only redueing merger
	if [ ${MERGERONLY} == 1 ]
	then
		if [ -d vcfMerger_pegasus ]
		then
			rm -rf vcfMerger_pegasus
		fi

		if [ -f SnpEFF_ANN_Complete ]
		then
			rm SnpEFF_ANN_Complete
		fi

		if [ -f SnpEFF_ANN_In_Progress ]
		then
			rm SnpEFF_ANN_In_Progress
		fi

		if [ -f SnpEFF_ANN_Fail ]
		then
          	rm SnpEFF_ANN_Fail
		fi
	fi

	# Restart digar/FusionValidator if only
	if [ ${FUSIONONLY} == 1 ]
	then
		if [ -d fusionValidator ]
		then
			rm -rf fusionValidator
		fi

		if [ -f FusionValidator_Complete ]
		then
			rm FusionValidator_Complete
		fi

		if [ -f FusionValidator_In_Progress ]
		then
			rm FusionValidator_In_Progress
		fi

		if [ -f FusionValidator_Fail ]
		then
          	rm FusionValidator_Fail
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
	
		# Remove alCount if it exists
		if [ -d alCount  ]
		then
			rm -rf alCount
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
		if [ ${TOPHATONLY} = 0 ] && [ ${FUSIONONLY} = 0 ]
		then
		# Check for any DNAPAIR= lines and start appropriate snpEff and CNA jobs
		
		DNAPAIRCOUNT="`grep "DNAPAIR=" ${PATIENT_NAME}.config | wc -l`"

		if [ ${DNAPAIRCOUNT} = 0 ] && [ ${MERGERONLY} == 0 ]
		then
			echo ${PATIENT_NAME} has no DNA pair lines. Skipping snpEff and CNA launcher.
			touch CNA_Manual_Complete
			
			if [ ${CNAONLY} = 0 ]
			then
				touch SnpEFF_ANN_Complete
				touch Delly_Complete
			fi
		elif [ ${DNAPAIRCOUNT} = 0 ] && [ ${MERGERONLY} == 1 ]
		then
			echo ${PATIENT_NAME} has no DNA pair lines. Skipping snpEff launcher.
			touch SnpEFF_ANN_Complete
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

			if [ ${EXOMECOUNT} = 0 ] && [ ${MERGERONLY} = 0 ]
			then
				echo ${PATIENT_NAME} has no Exome DNA pair lines. Skipping snpEff and CNA launcher.
				touch CNA_Manual_Complete
			
				if [ ${CNAONLY} = 0 ]
				then
					touch SnpEFF_ANN_Complete
				fi
			elif [ ${EXOMECOUNT} = 0 ] && [ ${MERGERONLY} == 1 ]
			then
				echo ${PATIENT_NAME} has no Exome DNA pair lines. Skipping snpEff launcher.
				touch SnpEFF_ANN_Complete
			else
				if [ ${MERGERONLY} = 0 ]
				then
					mkdir cna_manual_2016
					touch CNA_Manual_In_Progress
				fi
                	        
				if [ ${CNAONLY} = 0 ] 
				then
					mkdir vcfMerger_pegasus
					touch SnpEFF_ANN_In_Progress
				fi
		
				# loop through each DNA pair line Capturing needed variables and Start snpEff and CNA launcher for exomes       

				EXOMEPAIRSLIST=`echo ${EXOMEPAIRS[@]} | sed 's/\s/@/g'`
				EXPECTEDPAIRSLIST=`echo ${EXPECTEDPAIRS[@]} | sed 's/\s/@/g'`

				for PAIR in `echo ${EXOMEPAIRS[@]} | sed 's/\s/\n/g'`
				do
					# {{{

					echo $PAIR
					LIBTYPE="Exome"
					NORMALSAMPLE="`echo ${PAIR} | cut -d- -f1`"
					TUMORSAMPLE="`echo ${PAIR} | cut -d- -f2`"
					EXOMEPAIR="${NORMALSAMPLE}-${TUMORSAMPLE}"
					VCF=${EXOMEPAIR}.HC_All.snpEff.vcf
					NORMALDAT=${NORMALSAMPLE}.proj.md.jr.bam.clc.cln.dat
					TUMORDAT=${TUMORSAMPLE}.proj.md.jr.bam.clc.cln.dat
					TUMORSHORT="`echo ${TUMORSAMPLE} | cut -d_ -f1,2,3,4,5`"
					TUMORSPECIMEN="`echo ${TUMORSAMPLE} | cut -d_ -f1,2,3,4`"

					#############!!!!!!!!!!!!!!!!!!!!!!!!
					# MMRF specific
					RNASAMPLE="`grep TRIPLET4ALLELECOUNT ${PATIENT_NAME}.config | awk -F',' -v TUMORSAMPLE=${TUMORSAMPLE} '$2 == TUMORSAMPLE { print $3 }' | grep TSMRU`"
					############!!!!!!!!!!!!!!!!!!!!!!!!                       	        

					NASSAY="`grep "SAMPLE=" ${PATIENT_NAME}.config | grep ${NORMALSAMPLE} | cut -d, -f1 | cut -d= -f2`"
					TASSAY="`grep "SAMPLE=" ${PATIENT_NAME}.config | grep ${TUMORSAMPLE} | cut -d, -f1 | cut -d= -f2`"
					CNAEXOMETARGET=`grep "${NASSAY}_CNABEDP=" ${CONSTANTS} | cut -d= -f2 | tr -d '\n'`
					
					if [ -z ${RNASAMPLE} ]
					then
						RNAFLAG=NO
						RNAASSAY=NO
					else
						RNAFLAG=YES
						RNAASSAY="`grep "SAMPLE=" ${PATIENT_NAME}.config | grep $RNASAMPLE | awk -F '[=,]' '{ print $2 }' `"
					fi					

					if [[ "$TASSAY" == "TSE61" ]]
					then
						bedFile="/home/tgenref/homo_sapiens/grch37_hg19/capture_targets/illumina_truseq/TruSeq_exome_targeted_regions_b37_padded.bed"
					elif [[ "$TASSAY" == *S5U ]] || [[ "$kitName" == *S5X ]]
					then
						#bedFile="/home/tgenref/pipeline_v0.3/ensembl70/Ensembl_v70_hs37d5_exonic_coordinates_touched_v5UTR_padded25.bed"
						bedFile="/home/tgenref/homo_sapiens/grch37_hg19/hs37d5_tgen/gene_model/ensembl_v74/tool_specific_resources/vcfmerger/agilent_sureselect_v5_plusUTR/Agilent_SureSelect_V5_plusUTR_hs37d5_GRCh37.74_PaddedTargets_intersect_sorted_padded100.bed"
					elif [[ "$TASSAY" == *STX ]]
					then
						bedFile="/home/tgenref/homo_sapiens/grch37_hg19/capture_targets/agilent_sureselect_v5_strexome/Strexome_targets_intersect_sorted_padded100.bed"
					elif [[ "$TASSAY" == *SCR ]]
					then
						bedFile="/home/tgenref/homo_sapiens/grch37_hg19/capture_targets/agilent_sureselect_cre_v1/Agilent_Clinical_Research_Exome_hs37d5_PaddedTargets_intersect_sorted_padded100.bed"
						# missing after the move to new tgenref. also not in kbase
					#elif [[ "$TASSAY" == *S2X ]] ; then
						#bedFile="/home/tgenref/homo_sapiens/grch37_hg19//Agilent_V2_hs37d5_Targets_intersect_sorted_padded100.bed"
					elif [[ "$TASSAY" == *STL ]]
					then
						bedFile="/home/tgenref/homo_sapiens/grch37_hg19/capture_targets/agilent_custom_strexome_lite/Strexome_Lite_Targets_intersect_sorted_padded100.bed"
					elif [[ "$TASSAY" == *S1X ]]
					then
						bedFile="/home/tgenref/homo_sapiens/grch37_hg19/capture_targets/agilent_sureselect_v1_NA/SureSelectV1_hs37d5_PaddedTargets_intersect_sorted_padded100.bed"
					elif [[ "$TASSAY" == *S6X ]]
					then
						bedFile="/home/tgenref/homo_sapiens/grch37_hg19/capture_targets/agilent_sureselect_v6r2_noUTR/Agilent_SureSelect_V6R2_noUTR_hs37d5_GRCh37.74_PaddedTargets_intersect_sorted_padded100.bed"
					elif [[ "$TASSAY" == *SXP ]]
					then
						bedFile="/home/tgenref/homo_sapiens/grch37_hg19/capture_targets/agilent_sureselect_v6r2_plusUTR/Agilent_SureSelect_V6R2_plusUTR_hs37d5_GRCh37.74_PaddedTargets_intersect_sorted_padded100.bed"
					elif [[ "$TASSAY" == *S4X ]]
					then
						bedFile="/home/tgenref/homo_sapiens/grch37_hg19/capture_targets/agilent_sureselect_v4_noUTR/Agilent_V4_noUTR_hs37d5_Targets_intersect_sorted_padded100.bed"
					elif [[ "$TASSAY" == *E62 ]]
					then
						bedFile="/home/tgenref/homo_sapiens/grch37_hg19/capture_targets/illumina_nextera_expanded/NexteraExpandedExome_hs37d5_Targets_PicardPadded100.bed"
					fi
					
					if [ ${MERGERONLY} = 0 ]	
					then
						mkdir cna_manual_2016/${EXOMEPAIR}_exo
						touch cna_manual_2016/${EXOMEPAIR}_exo/CNA_Manual_In_Progress
					fi
					
					if [ ${CNAONLY} = 0 ]
					then
						mkdir vcfMerger_pegasus/${EXOMEPAIR}
						touch vcfMerger_pegasus/${EXOMEPAIR}/SnpEFF_ANN_In_Progress
					fi

					qsub -v DATAMOVERIP=${DATAMOVERIP},DATAMOVER=${DATAMOVER},BASEDIR="${BASEDIR}",RNAFLAG="${RNAFLAG}",RNASAMPLE="${RNASAMPLE}",RNAASSAY="${RNAASSAY}",CNAONLY="${CNAONLY}",MERGERONLY="${MERGERONLY}",BEDFILE=${bedFile},PBSNAME="${EXOMEPAIR}",STARTDIR="${STARTDIR}",CNAEXOMETARGET="${CNAEXOMETARGET}",NASSAY="${NASSAY}",TASSAY="${TASSAY}",LIBTYPE="${LIBTYPE}",USER="${USER}",PATIENT_NAME="${PATIENT_NAME}",NORMALSAMPLE="${NORMALSAMPLE}",TUMORSAMPLE="${TUMORSAMPLE}",VCF="${VCF}",NORMALDAT="${NORMALDAT}",TUMORDAT="${TUMORDAT}",EXOMEPAIR="${EXOMEPAIR}",EXOMEPAIRS=${EXOMEPAIRSLIST},EXPECTEDPAIRS=${EXPECTEDPAIRSLIST} ${EXOMEPBS}

					if [ $? -ne 0 ]
					then
						if [ ${MERGERONLY} == 0 ]
						then
							rm CNA_Manual_In_Progress
							rm cna_manual_2016/${EXOMEPAIR}_exo/CNA_Manual_In_Progress
							echo Failed to start rsync qsub job for ${EXOMEPAIR} >> CNA_Manual_Fail
							echo Failed to start rsync qsub job for ${EXOMEPAIR} >> cna_manual_2016/${EXOMEPAIR}_exo/CNA_Manual_Fail
						fi
						if [ ${CNAONLY} = 0 ]
						then
							rm SnpEFF_ANN_In_Progress
							rm vcfMerger_pegasus/${EXOMEPAIR}/SnpEFF_ANN_In_Progress
							echo Failed to start rsync qsub job for ${EXOMEPAIR} >> SnpEFF_ANN_Fail
							echo Failed to start rsync qsub job for ${EXOMEPAIR} >> vcfMerger_pegasus/${EXOMEPAIR}/SnpEFF_ANN_Fail
						fi
					fi

					#}}}

					# Start LI CNA if available

					 #{{{

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

					if [ ${LICOUNT} -eq 1 ] && [ ${MERGERONLY} == 0 ]
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
							DELLY="No"
	
							mkdir cna_manual_2016/${GENOMEPAIR}_filt2012
							touch cna_manual_2016/${GENOMEPAIR}_filt2012/CNA_Manual_In_Progress
							mkdir cna_manual_2016/${GENOMEPAIR}_filt2016
							touch cna_manual_2016/${GENOMEPAIR}_filt2016/CNA_Manual_In_Progress
							mkdir cna_manual_2016/${GENOMEPAIR}_unfi
							touch cna_manual_2016/${GENOMEPAIR}_unfi/CNA_Manual_In_Progress
	
							qsub -v DATAMOVERIP=${DATAMOVERIP},DATAMOVER=${DATAMOVER},BASEDIR="${BASEDIR}",PBSNAME=${GENOMEPAIR},DELLY=${DELLY},BEDFILE=${bedFile},STARTDIR="${STARTDIR}",GENOMEPAIR=${GENOMEPAIR},EXOMEPAIR=${EXOMEPAIR},NORMALSAMPLEG="${NORMALSAMPLEG}",TUMORSAMPLEG="${TUMORSAMPLEG}",NASSAY="${NASSAY}",TASSAY="${TASSAY}",LIBTYPE="${LIBTYPE}",PATIENT_NAME="${PATIENT_NAME}",NORMALSAMPLE="${NORMALSAMPLE}",TUMORSAMPLE="${TUMORSAMPLE}",VCF="${VCF}",NORMALDAT="${NORMALDAT}",TUMORDAT="${TUMORDAT}",EXPECTEDPAIRS=${EXPECTEDPAIRSLIST} ${GENOMEPBS}

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
					elif [ ${LICOUNT} -ne 0 ] && [ ${MERGERONLY} == 0 ]
					then
						echo There is more than one LI for the Tumor Isolation.
						echo "There is more than one LI for the Tumor Isolation." >> CNA_Manual_Fail
						echo "There is more than one LI for the Tumor Isolation. ${TUMORSHORT}" | mailx -s "LI_cna_manual_failed" ${USER}@tgen.org 
					fi

					#}}} 

				done
			fi

			# Start Delly

			#{{{

			if [ ${CNAONLY} = 0 ] && [ ${MERGERONLY} == 0 ]
			then
				LICOUNT=0
				LIPAIRS=()
				EXPECTEDPAIRS=()

				for PAIRTEST in `grep "DNAPAIR=" ${PATIENT_NAME}.config`
				do
					NORM=`echo ${PAIRTEST} | awk -F"[=,]" '{print $2}'`
					TUMR=`echo ${PAIRTEST} | awk -F"[=,]" '{print $3}'`

					NORMTYPE=`grep "SAMPLE=.*${NORM}" ${PATIENT_NAME}.config | awk -F',' '{ print $3 }'`
					TUMRTYPE=`grep "SAMPLE=.*${TUMR}" ${PATIENT_NAME}.config | awk -F',' '{ print $3 }'`

					if [ ${NORMTYPE} = "Genome" ] && [ ${TUMRTYPE} = "Genome" ]
					then
						((LICOUNT+=1))
						LIPAIRS+=("${NORM}-${TUMR}")
						EXPECTEDPAIRS+=("${NORM}-${TUMR}")
					fi
				done

				if [ ${LICOUNT} = 0 ]
				then
					echo ${PATIENT_NAME} has no LI DNA pair lines. Skipping delly.
					touch Delly_Complete
				else
					mkdir delly
					touch Delly_In_Progress	
	
					EXPECTEDPAIRSLIST=`echo ${EXPECTEDPAIRS[@]} | sed 's/\s/@/g'`

					for GPAIR in `echo ${LIPAIRS[@]} | sed 's/\s/\n/g'`
					do
						NORMALSAMPLEG="`echo ${GPAIR} | cut -d- -f1`"
						TUMORSAMPLEG="`echo ${GPAIR} | cut -d- -f2`"
						GENOMEPAIR="${NORMALSAMPLEG}-${TUMORSAMPLEG}"
						NASSAY="`echo ${NORMALSAMPLEG} | cut -d_ -f7`"
						TASSAY="`echo ${TUMORSAMPLEG} | cut -d_ -f7`"
						DELLY="Yes"

						mkdir delly/${GENOMEPAIR}
						touch delly/${GENOMEPAIR}/Delly_In_Progress
						
						qsub -v DATAMOVERIP=${DATAMOVERIP},DATAMOVER=${DATAMOVER},BASEDIR="${BASEDIR}",PBSNAME=${GENOMEPAIR},DELLY=${DELLY},STARTDIR="${STARTDIR}",NASSAY="${NASSAY}",TASSAY="${TASSAY}",GENOMEPAIR=${GENOMEPAIR},NORMALSAMPLEG="${NORMALSAMPLEG}",TUMORSAMPLEG="${TUMORSAMPLEG}",PATIENT_NAME="${PATIENT_NAME}",EXPECTEDPAIRS=${EXPECTEDPAIRSLIST} ${GENOMEPBS}

						if [ $? -ne 0 ]
						then
							rm Delly_In_Progress
							rm delly/${GENOMEPAIR}/Delly_In_Progress
							echo Failed to start rsync qsub job for ${GENOMEPAIR} >> Delly_Fail
							echo Failed to start rsync qsub job for ${GENOMEPAIR} >> delly/${GENOMEPAIR}/Delly_Fail
						fi
					done
				fi
			fi
			#}}}

		fi
		fi
	
		# Start Salmon and Kallisto, Fusion Validator TopHat Fusion if RNA is available
		
		#{{{
		
		if [ ${CNAONLY} = 0 ] && [ ${MERGERONLY} = 0 ]
		then
			RNACOUNT="`grep "SAMPLE=.*,RNA," ${PATIENT_NAME}.config | wc -l`"
	
			if [ ${RNACOUNT} = 0 ]
			then
				echo ${PATIENT_NAME} has no RNA.
				touch Salmon_Complete
				touch Kallisto_Complete
				touch FusionValidator_Complete
			else
				if [ ${TOPHATONLY} = 0 ] && [ ${FUSIONONLY} = 0 ]
				then
					touch Salmon_In_Progress
					touch Kallisto_In_Progress
					touch FusionValidator_In_Progress
				fi

				if [ ${FUSIONONLY} = 1 ]
				then
					touch FusionValidator_In_Progress
				fi

				# Make list of all RNA samples for final copy back check

				RNACHECK=()
				RNACHECK2=()

				for RNALIST in `grep "SAMPLE=.*,RNA," ${PATIENT_NAME}.config`
				do
					ASSAY="`echo ${RNALIST} | awk -F"[=,]" '{print $2}'`"
					RNASAMPLEID="`echo ${RNALIST} | awk -F"[=,]" '{print $3}'`"
					RNACHECK+=("${ASSAY}/${RNASAMPLEID}/${RNASAMPLEID}")
					RNACHECK2+=("fusionValidator/${RNASAMPLEID}")
				done

				RNACHECKLIST=`echo ${RNACHECK[@]} | sed 's/\s/@/g'`
				RNACHECKLIST2=`echo ${RNACHECK2[@]} | sed 's/\s/@/g'`

				for RNASAMPLE in `grep "SAMPLE=.*,RNA," ${PATIENT_NAME}.config`
				do
					LIBTYPE="RNA"
					ASSAY="`echo ${RNASAMPLE} | awk -F"[=,]" '{print $2}'`"
					STUDY="`echo ${PATIENT_NAME} | awk -F'_' '{print $1}'`"
					RNASAMPLEID="`echo ${RNASAMPLE} | awk -F"[=,]" '{print $3}'`"
					RNATUMORSHORT="`echo ${RNASAMPLEID} | cut -d_ -f1,2,3,4,5`"
					
					if [ ${RESTART} == 1 ] && [ ${TOPHATONLY} = 0 ] && [ ${FUSIONONLY} = 0 ]
					then
						if [ -d ${ASSAY}/${RNASAMPLEID}/${RNASAMPLEID}.salmonDir ]
						then
							rm -rf ${ASSAY}/${RNASAMPLEID}/${RNASAMPLEID}.salmonDir	                       
						fi

						if [ -d ${ASSAY}/${RNASAMPLEID}/${RNASAMPLEID}.kallistoDir ]
						then
							rm -rf ${ASSAY}/${RNASAMPLEID}/${RNASAMPLEID}.kallistoDir
						fi
						
					elif [ ${RESTART} == 1 ] && [ ${TOPHATONLY} = 1 ]
					then
						if [ -d ${ASSAY}/${RNASAMPLEID}/${RNASAMPLEID}.topHatFusionDir ]
						then
							# change this once confirmed that results are good
							rm -rf ${ASSAY}/${RNASAMPLEID}/${RNASAMPLEID}.topHatFusionDir
							#mv ${ASSAY}/${RNASAMPLEID}/${RNASAMPLEID}.topHatFusionDir ${ASSAY}/${RNASAMPLEID}/${RNASAMPLEID}.topHatFusionDir_old
						fi
					fi

					if [ ${TOPHATONLY} = 0 ] && [ ${FUSIONONLY} = 0 ]
					then
						mkdir -p ${ASSAY}/${RNASAMPLEID}/${RNASAMPLEID}.salmonDir/{ensembl74_cDNA,ensembl74_GTF,ensembl74_GTF_V7.2}
						mkdir -p ${ASSAY}/${RNASAMPLEID}/${RNASAMPLEID}.kallistoDir/{ensembl74_cDNA,ensembl74_GTF}
						mkdir -p fusionValidator/${RNASAMPLEID}
						touch ${ASSAY}/${RNASAMPLEID}/${RNASAMPLEID}.salmonDir/{ensembl74_cDNA,ensembl74_GTF,ensembl74_GTF_V7.2}/Salmon_In_Progress
						touch ${ASSAY}/${RNASAMPLEID}/${RNASAMPLEID}.kallistoDir/{ensembl74_cDNA,ensembl74_GTF}/Kallisto_In_Progress
						touch fusionValidator/${RNASAMPLEID}/FusionValidator_In_Progress
					elif [ ${FUSIONONLY} = 1 ]
					then
						mkdir -p fusionValidator/${RNASAMPLEID}
						touch fusionValidator/${RNASAMPLEID}/FusionValidator_In_Progress
					else
						mkdir -p ${ASSAY}/${RNASAMPLEID}/${RNASAMPLEID}.topHatFusionDir/tophatfusion_out
					fi

					# Get matching genomes for fusion validation

					MATCHEDWGSCOUNT=0
					unset NORM
					unset TUMR
					unset NORMTYPE
					unset TUMRTYPE
					unset NORMALSAMPLEG
					unset TUMORSAMPLEG
					unset NASSAY
					unset TASSAY


					for MATCHEDWGS in `grep "DNAPAIR=" ${PATIENT_NAME}.config | grep ${RNATUMORSHORT}`
					do
						NORM=`echo ${MATCHEDWGS} | awk -F"[=,]" '{print $2}'`
						TUMR=`echo ${MATCHEDWGS} | awk -F"[=,]" '{print $3}'`

						NORMTYPE=`grep "SAMPLE=.*${NORM}" ${PATIENT_NAME}.config | awk -F',' '{ print $3 }'`
						TUMRTYPE=`grep "SAMPLE=.*${TUMR}" ${PATIENT_NAME}.config | awk -F',' '{ print $3 }'`

						if [ ${NORMTYPE} = "Genome" ] && [ ${TUMRTYPE} = "Genome" ]
						then
							((MATCHEDWGSCOUNT+=1))
							
							if [ ${MATCHEDWGSCOUNT} -lt 2 ]	
							then
								NORMALSAMPLEG=`echo ${MATCHEDWGS} | awk -F"[=,]" '{print $2}'`
								NASSAY="`echo ${NORMALSAMPLEG} | cut -d_ -f7`"
								TUMORSAMPLEG=`echo ${MATCHEDWGS} | awk -F"[=,]" '{print $3}'`
								TASSAY="`echo ${TUMORSAMPLEG} | cut -d_ -f7`"
							else
								echo Warning!!!! ${RNASAMPLEID} pairs to more than one WGS pair.
								echo 	This is not expected and will cause many problems downstream.
								echo 	Exiting script...
								exit 1
							fi
						fi

					done

					if [ -z $TUMORSAMPLEG ]
					then
						NORMALSAMPLEG=NotAvailable
						TUMORSAMPLEG=NotAvailable
					fi

					if [ -z $NORMALSAMPLEG ]
					then
						NORMALSAMPLEG=NotAvailable
						TUMORSAMPLEG=NotAvailable
					fi

					# Will need code to get fastq directory and library type of RNA library for salmon      

					FASTQDIR="`grep "${STUDY}_STUDY=" ${CONSTANTS} | cut -d= -f2 | tr -d '\n'`"
					RNATYPE="`grep "${ASSAY}_SALMONlibType=" ${CONSTANTS} | cut -d= -f2 | tr -d '\n'`"
	
					qsub -v DATAMOVERIP=${DATAMOVERIP},DATAMOVER=${DATAMOVER},NASSAY=${NASSAY},TASSAY=${TASSAY},TUMORSAMPLEG=${TUMORSAMPLEG},NORMALSAMPLEG=${NORMALSAMPLEG},FUSIONONLY="${FUSIONONLY}",TOPHATONLY="${TOPHATONLY}",BASEDIR="${BASEDIR}",PBSNAME=${RNASAMPLEID},STARTDIR="${STARTDIR}",RNATYPE="${RNATYPE}",USER="${USER}",ASSAY="${ASSAY}",FASTQDIR="${FASTQDIR}",RNASAMPLEID="${RNASAMPLEID}",PATIENT_NAME="${PATIENT_NAME}",LIBTYPE="${LIBTYPE}",RNACHECK=${RNACHECKLIST},RNACHECK2=${RNACHECKLIST2} ${RNAPBS}

					STATUS=$?

					if [ $STATUS -ne 0 ] && [ ${TOPHATONLY} = 0 ]
					then
						rm Salmon_In_Progress
						rm Kallisto_In_Progress
						rm ${ASSAY}/${RNASAMPLEID}/${RNASAMPLEID}.salmonDir/{ensembl74_cDNA,ensembl74_GTF,ensembl74_GTF_V7.2}/Salmon_In_Progress
						rm ${ASSAY}/${RNASAMPLEID}/${RNASAMPLEID}.kallistoDir/{ensembl74_cDNA,ensembl74_GTF}/Kallisto_In_Progress
						echo Failed to start rsync qsub job for ${RNASAMPLEID} >> Salmon_Fail
						echo Failed to start rsync qsub job for ${RNASAMPLEID} >> Kallisto_Fail
						echo Failed to start rsync qsub job for ${RNASAMPLEID} | tee -a ${ASSAY}/${RNASAMPLEID}/${RNASAMPLEID}.salmonDir/{ensembl74_cDNA,ensembl74_GTF,ensembl74_GTF_V7.2}/Salmon_Fail
						echo Failed to start rsync qsub job for ${RNASAMPLEID} | tee -a ${ASSAY}/${RNASAMPLEID}/${RNASAMPLEID}.kallistoDir/{ensembl74_cDNA,ensembl74_GTF}/Kallisto_Fail
					elif [ $STATUS -ne 0 ] && [ ${TOPHATONLY} = 1 ]
					then
						echo Failed to start TopHat Fusion rsync qsub job for ${RNASAMPLEID} >> TopHatFusion_Fail
					fi
				done
			fi
		fi
		#}}}
	
	fi
done






