#!/bin/bash
#SBATCH --job-name=Combine_gVCFs
#SBATCH -N 1
#SBATCH -t 0-48:00
#SBATCH -n 1
#SBATCH --cpus-per-task 16
#SBATCH --mem-per-cpu=3000


######################################
# Needed variables
cd $OUT

. INPUT.ini

FOLDER=${PWD##*/}
OUTPUT_NAME=${STUDY}_${FOLDER}

#####################################################################################################


#Step1 - GENERATE THE LIST OF INPUT VCF FILES
ls *.vcf > temp_samples.txt

#Step2 - Pull off first record
head -n1 temp_samples.txt > temp_first.txt
FIRST_SAMPLE=`cat temp_first.txt`
INPUT_VARIANT_LIST="--variant ${FIRST_SAMPLE} "

#Step 3 - Pull off second to last records
awk 'NR>1' temp_samples.txt > temp_others.txt

#Step4 - Generate the input VCF list
for line in `cat temp_others.txt`
do

echo Adding sample ${line}
INPUT_VARIANT_LIST="${INPUT_VARIANT_LIST}--variant ${line} "

done

echo
echo "-------------------------------------------"
echo "COMPLETE VARIANT LIST:"
echo ${INPUT_VARIANT_LIST}

#Progress Marker
touch ../InProgress_${FOLDER}_CombineGVCFs
	
java -jar -Xmx32g ${GATKPATH}/GenomeAnalysisTK.jar \
	-T CombineGVCFs \
	-R ${REF} \
	${INPUT_VARIANT_LIST}\
	-o ${OUTPUT_NAME}.g.vcf

#Progress Marker
if [ "$?" = "0" ]
then
	touch ../Completed_${FOLDER}_CombineGVCFs
	rm ../InProgress_${FOLDER}_CombineGVCFs
else
	touch ../Fail_${FOLDER}_CombineGVCFs
	rm ../InProgress_${FOLDER}_CombineGVCFs
	exit 1
fi
	
rm temp_samples.txt
rm temp_first.txt
rm temp_others.txt

########
##  Test if next step can be started
#######

# Backup one directory
cd ..
## Test if all process are complete, if so  start CombineGVCFs
# Get a count of completed haplotypeCaller runs
FINISHED=`ls Completed_*_CombineGVCFs | wc -l`
# Get a count of the Study folders in current directory
FOLDERS=`ls -d subBatch_* | wc -l`

echo Found ${FOLDERS} folders and ${FINISHED} completed haplotypecaller runs

if [ $FINISHED -eq $FOLDERS ]
then
    # Start CombineGVCFs
	out=`pwd`
    sbatch --export ALL,OUT=${out} ${SCRIPTS_DIR}/gVCF_Step2_GenotypeSubBatchGVCFs.sh
else
    echo Other jobs still running, exiting now
fi
