#!/bin/bash
#SBATCH --job-name=Haplotyper_Test
#SBATCH -N 1
#SBATCH -t 0-24:00
#SBATCH -n 1
#SBATCH --cpus-per-task 2
#SBATCH --mem-per-cpu=4000

cd $OUT
# Read in variables
. INPUT.ini

########
##   Operational Code
########


#Progress Marker
touch ../InProgress_${SAMPLE}_HaplotyeCaller

java -jar -Xmx4g ${GATKPATH}/GenomeAnalysisTK.jar \
	-T HaplotypeCaller \
	-R ${REF} \
	-I ${INPUT_BAM} \
	--dbsnp ${DBSNP} \
	--emitRefConfidence GVCF \
	--min_base_quality_score 10 \
	--annotation HaplotypeScore \
	--out ${SAMPLE}.g.vcf \
	-L ${TARGET_INTERVALS}

# -nct X , was removed after testing showed it did not decrease processing time

#Progress Marker
if [ "$?" = "0" ]
then
	touch ../Completed_${SAMPLE}_HaplotyeCaller
	rm ../InProgress_${SAMPLE}_HaplotyeCaller
else
	touch ../FAIL_${SAMPLE}_HaplotyeCaller
	rm ../InProgress_${SAMPLE}_HaplotyeCaller
	exit 1
fi

echo "Finished Generating per sample gVCF"
