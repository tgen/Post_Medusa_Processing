#!/bin/sh


##
# Help Menu
##

if [ "$1" = "--help" ]
then
	echo
	echo
	echo "This Script will start the analysis of a folder of exome BAM files by haplotypecaller"
	echo "Requires a configured parameters file (*.ini)"
	echo
	echo "Usage: Start_HaplotypeCaller_gVCF_Pipeline.sh <Parameters.ini>"
	echo
	echo
	exit 1
elif [ -z $1 ]
then
	echo
	echo
	echo "You did not provide an PARAMETERS.ini file.... exiting now"
	echo
	echo
	exit 1
elif [ -n $1 ]
then
	echo
	echo
	echo Starting Process
	echo
	echo
fi

##
## This script will take a folder of exome BAM and BAI files as input:
##          - It will create a subfolder for each BAM and BAI pair
##	    - It will submit a qsub job for each BAM to generate the gVCF
##      NOTE: You will need to copy the gVCF to the correct location


# Read in the ini file
. $1

#Copy ini file to processing folder
cp $1 INPUT.ini

# Output setup
echo
echo
echo Enviroment Setup
echo Study: ${STUDY}
echo GATK version: ${GATKPATH}
echo
echo

# Create the subfolders and launch gVCF creation process
for line in `ls *.bam`
do

	#INPUT_BAM=
	SAMPLE=`basename ${line} ".bwa.final.bam"`
	echo "-------------------------------------------"
	echo Processing ${SAMPLE}

	#Make sample directory
	mkdir ${SAMPLE}
	#Move BAM and BAI into the folder
	mv ${SAMPLE}.bwa.final.bai ${SAMPLE}/
	mv ${SAMPLE}.bwa.final.bam ${SAMPLE}/
	#Move to folder and start job
	cd ${SAMPLE}
	#Copy ini filde to processing folder
	cp $1 INPUT.ini
	qsub -v INPUT_BAM=$line,SAMPLE=$SAMPLE /home/jkeats/toolkit_jkeats/haplotypeCaller_gVCF_Pipeline/start_at_gVCF/gVCF_Make_HaplotypeCaller_gVCF_Process.pbs
	cd ..

done

echo
echo
echo "Finished Launching all Jobs"
echo
echo
