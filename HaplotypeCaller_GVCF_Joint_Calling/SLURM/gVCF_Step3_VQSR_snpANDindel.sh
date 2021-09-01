#!/bin/bash
#SBATCH --job-name=VQSR_Test
#SBATCH -N 1
#SBATCH -t 0-48:00
#SBATCH -n 1
#SBATCH --cpus-per-task 16
#SBATCH --mem-per-cpu=3000

#Load R so images are created
module load R/3.2.1

cd $OUT
# Load .ini file with variables
. INPUT.ini

#****************************************************************************************************************************************
###Resources:   
#https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantrecalibration_ApplyRecalibration.php
#https://www.broadinstitute.org/gatk/guide/article?id=1259
#****************************************************************************************************************************************


INPUT_VCF=${STUDY}_Exome_Raw.vcf
INPUT_VCF_BASE=${STUDY}_Exome_Raw

#### Build the SNP recalibration model

#Progress Marker
touch InProgress_SNP_CalculateRecal

java -jar -Xmx4g ${GATKPATH}/GenomeAnalysisTK.jar \
	-T VariantRecalibrator \
	-nt 8 \
	-R ${REF} \
	-input ${INPUT_VCF} \
	-resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${HAPMAP} \
	-resource:omni,known=false,training=true,truth=true,prior=12.0 ${OMNI} \
	-resource:1000G,known=false,training=true,truth=false,prior=10.0 ${KG} \
	-resource:dbsnp138,known=false,training=false,truth=false,prior=7.0 ${DBSNP138} \
	-resource:dbsnp129,known=true,training=false,truth=false,prior=3.0 ${DBSNP129} \
	--maxGaussians 6 \
	-mode SNP \
	-an QD \
	-an FS \
	-an SOR \
	-an MQ \
	-an MQRankSum \
	-an ReadPosRankSum \
	-an InbreedingCoeff \
	-allPoly -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.7 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.2 -tranche 99.1 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 92.5 -tranche 90.0 \
	-recalFile ${INPUT_VCF_BASE}.recalibrate_SNP.recal \
	-tranchesFile ${INPUT_VCF_BASE}.recalibrate_SNP.tranches \
	-rscriptFile ${INPUT_VCF_BASE}_SNP.plots.R

## REMOVED FOR PRODUCTION, 
## Because (https://www.broadinstitute.org/gatk/guide/article?id=1259)
## For exomes, we do not recommend using DP for variant recalibration (see below for details of why).
# -an DP \

## REMOVED FOR TESTING, NOT CALCULATED DUE TO TEST COHORT SIZE, REQUIRES AT LEAST 10 UNRELATED INDIVIDUALS
# -an InbreedingCoeff \

#Error Capture
if [ "$?" = "0" ]
then
        touch Completed_SNP_CalculateRecal
        rm InProgress_SNP_CalculateRecal
else
        touch FAIL_SNP_CalculateRecal
        rm InProgress_SNP_CalculateRecal
        exit 1
fi

#### Apply the desired level of recalibration to the SNPs in the call set

#Progress Marker
touch InProgress_SNP_ApplyRecal

java -jar -Xmx4g ${GATKPATH}/GenomeAnalysisTK.jar \
	-T ApplyRecalibration \
	-nt 8 \
	-R ${REF} \
	-input ${INPUT_VCF} \
	-mode SNP \
	--ts_filter_level 99.0 \
	-recalFile ${INPUT_VCF_BASE}.recalibrate_SNP.recal \
	-tranchesFile ${INPUT_VCF_BASE}.recalibrate_SNP.tranches \
	-o ${INPUT_VCF_BASE}_rcSNP.vcf 

## Removed, so all variants survive and can be post filtered
#--excludeFiltered \

#Error Capture
if [ "$?" = "0" ]
then
        touch Completed_SNP_ApplyRecal
        rm InProgress_SNP_ApplyRecal
else
        touch FAIL_SNP_ApplyRecal
        rm InProgress_SNP_ApplyRecal
        exit 1
fi
 
#### Build the Indel recalibration model

#Progress Marker
touch InProgress_INDEL_CalculateRecal

java -jar -Xmx16g ${GATKPATH}/GenomeAnalysisTK.jar \
	-T VariantRecalibrator \
	-nt 8 \
	-R ${REF} \
	-resource:mills,known=false,training=true,truth=true,prior=12.0 ${Mills_1000g} \
	-resource:dbsnp138,known=true,training=false,truth=false,prior=2.0 ${DBSNP138} \
	-resource:axiomPoly,known=false,training=true,truth=false,prior=10.0 ${AXIOM} \
	-input ${INPUT_VCF_BASE}_rcSNP.vcf \
	-mode INDEL \
	--maxGaussians 6 \
	-an QD \
	-an FS \
	-an SOR \
	-an MQRankSum \
	-an ReadPosRankSum \
	-an InbreedingCoeff \
	-allPoly -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
	-recalFile ${INPUT_VCF_BASE}.recalibrate_INDEL.recal \
	-tranchesFile ${INPUT_VCF_BASE}.recalibrate_INDEL.tranches \
	-rscriptFile ${INPUT_VCF_BASE}_INDEL.plots.R  


## REMOVED FOR PRODUCTION, 
## Because (https://www.broadinstitute.org/gatk/guide/article?id=1259)
## For exomes, we do not recommend using DP for variant recalibration (see below for details of why).
# -an DP \

## REMOVED FOR TESTING, NOT CALCULATED DUE TO TEST COHORT SIZE, REQUIRES AT LEAST 10 UNRELATED INDIVIDUALS
# -an InbreedingCoeff \



#Error Capture
if [ "$?" = "0" ]
then
        touch Completed_INDEL_CalculateRecal
        rm InProgress_INDEL_CalculateRecal
else
        touch FAIL_INDEL_CalculateRecal
        rm InProgress_INDEL_CalculateRecal
        exit 1
fi

#### Final Recalibrated vcf

#Progress Marker
touch InProgress_INDEL_ApplyRecal

java -jar -Xmx4g ${GATKPATH}/GenomeAnalysisTK.jar \
	-T ApplyRecalibration \
	-nt 8 \
	-R ${REF} \
	-input ${INPUT_VCF_BASE}_rcSNP.vcf \
	-mode INDEL \
	--ts_filter_level 95.0 \
	-recalFile ${INPUT_VCF_BASE}.recalibrate_INDEL.recal \
	-tranchesFile ${INPUT_VCF_BASE}.recalibrate_INDEL.tranches \
	-o ${INPUT_VCF_BASE}_rcSNP_rcINDEL.vcf 

## Removed, so all variants survive and can be post filtered
#--excludeFiltered \

#Error Capture
if [ "$?" = "0" ]
then
        touch Completed_INDEL_ApplyRecal
        rm InProgress_INDEL_ApplyRecal
else
        touch FAIL_INDEL_ApplyRecal
        rm InProgress_INDEL_ApplyRecal
        exit 1
fi

## SHOULD ADD ANNOTATION STEPS
