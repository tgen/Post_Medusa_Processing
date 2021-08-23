#!/bin/sh


##
# Help Menu
##

if [ "$1" = "--help" ]
then
	echo
	echo
	echo "This Script will start the joint calling of a folder of gVCF files by haplotypecaller"
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
## This script will take a folder of gVCF VCF and IDX files as input:
##          - It will separate sets of X defined in PARAMETERS.ini into sub-batch folders
##	    - It will start a combineGVCF process for all sub-batch gVCFs
##	    - After the final sub-batch gVCF is created it will start a GenotypeGVCF process
##	    - After the GenotypeGVCF process finishes it will start a VQSR process
##


# Read in the ini file
. $1

#Copy ini file to processing folder
cp $1 INPUT.ini

# Output setup
echo
echo
echo Enviroment Setup
echo Study: ${STUDY}
echo gVCF Folder Location: ${GVCF_PATH}
echo GATK version: ${GATKPATH}
echo Batch Size: ${BATCH_SIZE}
echo
echo

# Create the subfolders and move sets of defined batch size into each folder

SUBBATCH=1
COUNT=0

for line in `ls ${GVCF_PATH}/*.vcf`
do

	## Make batch directory
	# Test if batch directory exits
	if [ -e subBatch_${SUBBATCH} ]
	then
		# Do Nothing
		echo Folder Already Exists
	else
		# Make the needed directory
		echo Making Subfolder
		mkdir subBatch_${SUBBATCH}
	fi

	## Move files into directory
	# Create variable so you can easily move VCF and IDX
    SAMPLE=`basename ${line} ".g.vcf"`
    # Move files
    mv ${GVCF_PATH}/${SAMPLE}.g.vcf subBatch_${SUBBATCH}/
	mv ${GVCF_PATH}/${SAMPLE}.g.vcf.idx subBatch_${SUBBATCH}/

	## Increment batch trackers
	COUNT=`expr $COUNT + 1`
	# Test count to determine if you need to increase the batch size
	if [ $COUNT -lt $BATCH_SIZE ]
	then
		# Do Nothing
		echo Still in current batch
	else
		# Need to increment the batch and reset count
		COUNT=0
		SUBBATCH=`expr $SUBBATCH + 1`
	fi

done

## Now walk through each sub-batch directory and start the combine gVCF step

for folder in `ls -d subBatch_*`
do

cd ${folder}
#Copy ini filde to processing folder
cp $1 INPUT.ini

# Start combine gVCF job
qsub /home/jkeats/toolkit_jkeats/haplotypeCaller_gVCF_Pipeline/start_at_gVCF/gVCF_Step1_CombineGVCFs.pbs

cd ..

done
