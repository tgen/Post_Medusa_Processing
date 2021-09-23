#!/bin/bash
#SBATCH --job-name=GenotypeSubBatch_gVCFs
#SBATCH -N 1
#SBATCH -t 0-48:00
#SBATCH -n 1
#SBATCH --cpus-per-task 16
#SBATCH --mem-per-cpu=3000


######################################
# Needed variables
cd $OUT

. INPUT.ini

OUTPUT_CALLS=${STUDY}_Exome_Raw.vcf

#####################################################################################################

# This script is called from the original starting folder that had all the gVCFs originally

#Step1 - GENERATE THE LIST OF INPUT SUB-BATCH COMBINED gVCF FILES
ls -d subBatch_* > temp_samples.txt

#Step2 - Pull off first record
head -n1 temp_samples.txt > temp_first.txt
FIRST_SAMPLE=`cat temp_first.txt`
INPUT_VARIANT_LIST="--variant ${FIRST_SAMPLE}/${STUDY}_${FIRST_SAMPLE}.g.vcf "

#Step 3 - Pull off second to last records
awk 'NR>1' temp_samples.txt > temp_others.txt

#Step4 - Generate the input bam list
for line in `cat temp_others.txt`
do

echo Adding sample ${line}
INPUT_VARIANT_LIST="${INPUT_VARIANT_LIST}--variant ${line}/${STUDY}_${line}.g.vcf "

done

echo
echo "-------------------------------------------"
echo "COMPLETE VARIANT LIST:"
echo ${INPUT_VARIANT_LIST}

#Progress Marker
touch InProgress_Final_GenotypeSubBatchGVCFs

java -jar -Xmx4g ${GATKPATH}/GenomeAnalysisTK.jar \
	-T GenotypeGVCFs \
	-nt 8 \
	-R ${REF} \
	${INPUT_VARIANT_LIST}\
	--dbsnp ${DBSNP} \
	--out ${OUTPUT_CALLS}

#Progress Marker
if [ "$?" = "0" ]
then
	touch Completed_Final_GenotypeSubBatchGVCFs
	rm InProgress_Final_GenotypeSubBatchGVCFs
else
	touch Fail_Final_GenotypeSubBatchGVCFs
	rm InProgress_Final_GenotypeSubBatchGVCFs
	exit 1
fi
	
rm temp_samples.txt
rm temp_first.txt
rm temp_others.txt

########
## If we got to this point we can start the VQSR step
#######

# Start VQSR
out=`pwd`
sbatch --export ALL,OUT=${out} ${SCRIPTS_DIR}/gVCF_Step3_VQSR_snpANDindel.sh



