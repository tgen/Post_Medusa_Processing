#!/bin/bash

##
## Post_Medusa_Processing_V1.0.sh
##
## Created by Austin Christofferson on 04/06/2016.
## Copyright 2016 Translational Genomics Research Institute. All rights reserved.
##
## Warning     - This script must be run from the results directory on a box that 
##		 has qsub (merckx) and removes the following folders.
##			- fastqc
##			- seurat/*REVERSE
##
## Description - This script will loop through new "ps*" time stamped results 
##		 directories renaming them, copying files over to pnap
##               for reduing snpEFF with updated databases and tool versions
##		 as well as starting the advanced cna scripts.
##               After all pbs jobs on pnap are finished the new data will
##               automatically copy back to isilon.
##
## Usage = Post_Medusa_Processing_V1.0.sh 
##
####################################################################
####################################################################

if [ ! -z "$1" ]
then
	if [ $1 == "--help" ] || [ $1 == "-h" ]
	then
	        grep "^##" $0
	        exit 0
	fi
fi

EXOMEPBS=/home/tgenref/pecan/post_central_pipe_processing/post_medusa_V1.0/jobScripts/rsync_Exome_Files.pbs
GENOMEPBS=/home/tgenref/pecan/post_central_pipe_processing/post_medusa_V1.0/jobScripts/rsync_Genome_Files.pbs
RNAPBS=/home/tgenref/pecan/post_central_pipe_processing/post_medusa_V1.0/jobScripts/rsync_RNA_Files.pbs
CONSTANTS=/home/tgenref/pecan/post_central_pipe_processing/post_medusa_V1.0/constants.txt
STARTDIR=`pwd`

# Make pbs out directory if does not exist
if [ ! -d "${HOME}/STD_OUT" ]
then
	mkdir ${HOME}/STD_OUT
fi


# Loop through all files in current directory ending with "*_ps*"
for line in `ls -d *_ps*`
do
	cd ${STARTDIR}

	# Check if copy is complete DO NOT RENAME IF NOT COMPLETE, THAT WOULD BREAK THE PIPELINE
	cd ${line}
	
	if [ -f "copyDone.txt" ]
	then
		# Remove all fastqc directories and files
		find . -name *_fastqc -type d -exec rm -rf {} + 
		
		# Remove all seurat REVERSE dir
		if [ "` ls seurat/ |grep REVERSE |wc -l`" -gt 0 ] 
		then
			rm -rf seurat/*REVERSE
		fi
			
		cd ..
		#Capture the patient name from each file as a variable
		PATIENT_NAME=`echo ${line} | awk -F'_' '{ OFS = "_" ; print $1,$2 }'`
	
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
		# Remove scratch Post processing directory if it exists
		ssh ${USER}@pnap-data3.tgen.org "if [ -d /scratch/${USER}/post_medusa_processing/${PATIENT_NAME} ] ; then rm -rf /scratch/${USER}/post_medusa_processing/${PATIENT_NAME} ; fi"
		if [ $? -ne 0 ]
		then
			echo
			echo Warning!!! ssh failed to remove the pre-existing Post Processing directory for on scratch.
			echo /scratch/${USER}/post_medusa_processing/${PATIENT_NAME} 
			echo
			echo Would you like to try to auto remove it again or remove it manualy [ auto/manual ]
			read RESPONSE

			if [ ${RESPONSE} = "auto" ]
			then
				ssh ${USER}@pnap-data3.tgen.org "if [ -d /scratch/${USER}/post_medusa_processing/${PATIENT_NAME} ] ; then rm -rf /scratch/${USER}/post_medusa_processing/${PATIENT_NAME} ; fi"

				if [ $? -ne 0 ]
                		then
					echo
					echo Auto rm post processing directory failed again. Moving on. Please check ${PATIENT_NAME}!! 
				fi
			elif [ ${RESPONSE} = "manual" ]
			then
				
			else

			fi 
		fi
	
		# Check for any DNAPAIR= lines and start appropriate snpEff and CNA jobs
		
		cd ${PATIENT_NAME}

		DNAPAIRCOUNT="`grep "DNAPAIR=" ${PATIENT_NAME}.config | wc -l`" 	
		
		if [ ${DNAPAIRCOUNT} = 0 ]
		then
			echo ${PATIENT_NAME} has no DNA pair lines. Skipping snpEff and CNA launcher.
			touch CNA_Manual_Complete
			touch SnpEFF_ANN_Complete
		else
			# Find out if there are any exome DNA pair lines. If none then skip snpEff and CNA launcher

			EXOMECOUNT=0
			EXOMEPAIRS=()
			EXPECTEDPAIRS=()			

			for PAIRTEST in `grep "DNAPAIR=" ${PATIENT_NAME}.config`
			do
				NORM=`echo ${PAIRTEST} | awk -F"[=,]" '{print $2}'`
				TUMR=`echo ${PAIRTEST} | awk -F"[=,]" '{print $3}'`
				
				NORMTYPE=`grep "SAMPLE=.*${NORM}" ${PATIENT_NAME}.config | awk -F',' '{ print $3 }'` 
				TUMRTYPE=`grep "SAMPLE=.*${TUMR}" ${PATIENT_NAME}.config | awk -F',' '{ print $3 }'` 
				
				if [ ${NORMTYPE} = "Exome" ] && [ ${TUMRTYPE} = "Exome" ]
				then
					((EXOMECOUNT+=1))
					EXOMEPAIRS+=("${NORM}-${TUMR}")
					EXPECTEDPAIRS+=("${NORM}-${TUMR}_exo")	
				fi
			done
			
			# Make array of exome pairs and their corresponding genome
			
			for COMBINED in `echo ${EXOMEPAIRS[@]} | sed 's/\s/\n/g'`
			do
				TUMORSAMPLE="`echo ${COMBINED} | cut -d- -f2 `"
				TUMORSHORT="`echo ${TUMORSAMPLE} | cut -d_ -f1,2,3,4,5,6`"
				
				for PAIRTEST in `grep "DNAPAIR=" ${PATIENT_NAME}.config | grep ${TUMORSHORT}`
                                do
                                	NORM=`echo ${PAIRTEST} | awk -F"[=,]" '{print $2}'`
                                        TUMR=`echo ${PAIRTEST} | awk -F"[=,]" '{print $3}'`

                                        NORMTYPE=`grep "SAMPLE=.*${NORM}" ${PATIENT_NAME}.config | awk -F',' '{ print $3 }'`
                                        TUMRTYPE=`grep "SAMPLE=.*${TUMR}" ${PATIENT_NAME}.config | awk -F',' '{ print $3 }'`

                                        if [ ${NORMTYPE} = "Genome" ] && [ ${TUMRTYPE} = "Genome" ]
                                        then
                                        	EXPECTEDPAIRS+=("${NORM}-${TUMR}_filt")
						EXPECTEDPAIRS+=("${NORM}-${TUMR}_unfi")
                                	fi
                        	done
			done
					
	
			if [ ${EXOMECOUNT} = 0 ]
			then
				echo ${PATIENT_NAME} has no Exome DNA pair lines. Skipping snpEff and CNA launcher.
				touch CNA_Manual_Complete
				touch SnpEFF_ANN_Complete
			else
				mkdir cna_manual
				touch CNA_Manual_In_Progress
				touch SnpEFF_ANN_In_Progress
				
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
					TUMORSHORT="`echo ${TUMORSAMPLE} | cut -d_ -f1,2,3,4,5,6`"
					NASSAY="`echo ${NORMALSAMPLE} | cut -d_ -f7`"
					TASSAY="`echo ${TUMORSAMPLE} | cut -d_ -f7`"
					CNAEXOMETARGET=`grep "${NASSAY}_CNABEDP=" ${CONSTANTS} | cut -d= -f2 | tr -d '\n'`
					
					mkdir cna_manual/${EXOMEPAIR}_exo
					touch cna_manual/${EXOMEPAIR}_exo/CNA_Manual_In_Progress
					touch vcfMerger/${EXOMEPAIR}/SnpEFF_ANN_In_Progress
					
					qsub -v PBSNAME="${EXOMEPAIR}",STARTDIR="${STARTDIR}",CNAEXOMETARGET="${CNAEXOMETARGET}",NASSAY="${NASSAY}",TASSAY="${TASSAY}",LIBTYPE="${LIBTYPE}",USER="${USER}",PATIENT_NAME="${PATIENT_NAME}",NORMALSAMPLE="${NORMALSAMPLE}",TUMORSAMPLE="${TUMORSAMPLE}",VCF="${VCF}",NORMALDAT="${NORMALDAT}",TUMORDAT="${TUMORDAT}",EXOMEPAIR="${EXOMEPAIR}",EXOMEPAIRS=${EXOMEPAIRSLIST},EXPECTEDPAIRS=${EXPECTEDPAIRSLIST} ${EXOMEPBS}
					
					if [ $? -ne 0 ]
                                        then
                                        	rm CNA_Manual_In_Progress
						rm SnpEFF_ANN_In_Progress
						rm cna_manual/${EXOMEPAIR}_exo/CNA_Manual_In_Progress
						rm vcfMerger/${EXOMEPAIR}/SnpEFF_ANN_In_Progress
						echo Failed to start rsync qsub job for ${EXOMEPAIR} >> CNA_Manual_Fail
						echo Failed to start rsync qsub job for ${EXOMEPAIR} >> SnpEFF_ANN_Fail
                                                echo Failed to start rsync qsub job for ${EXOMEPAIR} >> cna_manual/${EXOMEPAIR}_exo/CNA_Manual_Fail
                                                echo Failed to start rsync qsub job for ${EXOMEPAIR} >> vcfMerger/${EXOMEPAIR}/SnpEFF_ANN_Fail
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
							
							mkdir cna_manual/${GENOMEPAIR}_filt
							touch cna_manual/${GENOMEPAIR}_filt/CNA_Manual_In_Progress
							mkdir cna_manual/${GENOMEPAIR}_unfi
							touch cna_manual/${GENOMEPAIR}_unfi/CNA_Manual_In_Progress
					
							qsub -v PBSNAME=${GENOMEPAIR},STARTDIR="${STARTDIR}",GENOMEPAIR=${GENOMEPAIR},EXOMEPAIR=${EXOMEPAIR},NORMALSAMPLEG="${NORMALSAMPLEG}",TUMORSAMPLEG="${TUMORSAMPLEG}",NASSAY="${NASSAY}",TASSAY="${TASSAY}",LIBTYPE="${LIBTYPE}",PATIENT_NAME="${PATIENT_NAME}",NORMALSAMPLE="${NORMALSAMPLE}",TUMORSAMPLE="${TUMORSAMPLE}",VCF="${VCF}",NORMALDAT="${NORMALDAT}",TUMORDAT="${TUMORDAT}",EXPECTEDPAIRS=${EXPECTEDPAIRSLIST} ${GENOMEPBS}
							
							if [ $? -ne 0 ]
                                			then
                                        			rm CNA_Manual_In_Progress
                                        			rm cna_manual/${GENOMEPAIR}_filt/CNA_Manual_In_Progress
								rm cna_manual/${GENOMEPAIR}_unfi/CNA_Manual_In_Progress
                                        			echo Failed to start rsync qsub job for ${GENOMEPAIR} >> CNA_Manual_Fail
                                        			echo Failed to start rsync qsub job for ${GENOMEPAIR} >> cna_manual/${GENOMEPAIR}_filt/CNA_Manual_Fail
								echo Failed to start rsync qsub job for ${GENOMEPAIR} >> cna_manual/${GENOMEPAIR}_unfi/CNA_Manual_Fail
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

	else
		cd ..
		#Message
		echo
		echo "---------------------------------"
		echo "WARNING: Transfere is not complete"
		echo "${line} WILL NOT BE RENAMED!!"
	fi
done
