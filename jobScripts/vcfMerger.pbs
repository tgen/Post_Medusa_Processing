#!/bin/bash
#SBATCH --job-name="vcfMerger"
#SBATCH --time=24:00:00 
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --cpus-per-task 8
#SBATCH --mail-user=${USER}@tgen.org

cd ${DIR}

set -o pipefail

# Does not exist. need to use newer version
#module load BEDTools/2.14.0
module load BEDTools/2.26.0

# Does not exist. need to use newer version
#module load R/3.1.1
module load R/3.2.2

module load java/jdk-1.8.0_91

PATH=$PATH:/home/tgenref/binaries/vt/
export PATH

export PERL5LIB=/home/achristofferson/local/vcftools_v0.1.14/src/perl/

# Tool and Data Base Paths

#{{{ 
PICARDPATH=/packages/picard-tools/1.128
EXAC=/home/tgenref/homo_sapiens/grch37_hg19/public_databases/exac/r0.3.0/ExAC.r0.3.sites.vep_fixed_norm.vcf
COSMICC=/home/tgenref/homo_sapiens/grch37_hg19/public_databases/cosmic/v74/CosmicCodingMuts_with_contigs.norm.vcf
COSMICNC=/home/tgenref/homo_sapiens/grch37_hg19/public_databases/cosmic/v74/CosmicNonCodingVariants_with_contigs_norm.vcf
DBNSFP=/home/tgenref/homo_sapiens/grch37_hg19/public_databases/dbnsfp/v2.9/dbNSFP2.9.txt.gz
SNPEFFPATH=/home/achristofferson/local/snpEff_4.2_2015-12-05
SNPSIFT=/home/achristofferson/local/snpEff_4.2_2015-12-05/SnpSift.jar
SNPEFF=/home/achristofferson/local/snpEff_4.2_2015-12-05/snpEff.jar
SAMTOOLS=/packages/samtools/1.4.1/bin
VARSCAN=/home/tgenref/binaries/varscan/varscan-2.3.7
REF=/home/tgenref/homo_sapiens/grch37_hg19/hs37d5_tgen/genome_reference/hs37d5_plusRibo_plusOncoViruses_plusERCC.fa
DICT=/home/tgenref/homo_sapiens/grch37_hg19/hs37d5_tgen/genome_reference/hs37d5_plusRibo_plusOncoViruses_plusERCC.dict
KG=/home/tgenref/homo_sapiens/grch37_hg19/public_databases/1000G/ALL.ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_no_samples_OLD_VARIANT_norm.vcf
NHLBI=/home/tgenref/homo_sapiens/grch37_hg19/public_databases/nhlbi/esp6500si_v2_ssa137/ESP6500SI-V2-SSA137.GRCh38-liftover.all.snps_indels_with_contigs_norm.vcf
INDELS=/home/tgenref/homo_sapiens/grch37_hg19/public_databases/1000G/Mills_and_1000G_gold_standard.indels.b37.vcf
GATK=/home/achristofferson/local/GenomeAnalysisTK-3.8-0-ge9d806836
VCFMERGER=/home/tgenref/binaries/vcfMerger/mergeVcf/normalization/vcfMerger/pecan.merge.3vcfs.main.sh
VCFMERGER_DIR=/home/tgenref/binaries/vcfMerger/mergeVcf/normalization/vcfMerger
VCFSORTER=/home/tgenref/binaries/vcfMerger/vcfMerger/vcfsorter.pl
RNA_VCF_HEADER=${BASEDIR}/RNA_VCF_HEADER.vcf
POST_MERGE_VENN=/home/tgenref/binaries/vcfMerger/mergeVcf/normalization/vcfMerger/pecan.Venn_postMMRF_specific_filtering.sh
DBSNP=/home/tgenref/homo_sapiens/grch37_hg19/public_databases/dbsnp/b147/dbsnp_147.b37_norm.vcf
DBVERSION=GRCh37.74

VCFVALIDATOR=/home/achristofferson/local/vcftools_v0.1.14/bin/vcf-validator
VCFVALIDATORPATH=/home/achristofferson/local/vcftools_v0.1.14/bin
PATH=$VCFVALIDATORPATH:$PATH
BCFTOOLS=/home/achristofferson/local/bcftools-1.6/bcftools
BCFTOOLSPATH=/home/achristofferson/local/bcftools-1.6
PATH=$BCFTOOLSPATH:$PATH
export PATH

HEADER=${BASEDIR}/vcf_RNA_INFO_header.hdr

#}}}

FAILED=0

# Filter the seurat vcf

cat ${EXOMEPAIR}.seurat.vcf | \
	java -jar -Xmx20g ${SNPSIFT} filter "( TYPE='somatic_SNV' )" > ${EXOMEPAIR}_seurat_snv.vcf

if [ $? -ne 0 ]
then
        FAILED=1
	echo
	echo !!!!!!!!!!!!  Failed Filter the seurat vcf !!!!!!!!!!
	echo
fi

# Filter the seurat INDELS, use bed if no matched normal

cat ${EXOMEPAIR}.seurat.vcf | \
	java -jar -Xmx20g ${SNPSIFT} filter "(( TYPE='somatic_deletion' ) | ( TYPE='somatic_insertion' ))" > ${EXOMEPAIR}_seurat_indel.vcf

if [ $? -ne 0 ]
then
        FAILED=1
	echo
        echo !!!!!!!!!!!!  Failed Filter the seurat INDELS !!!!!!!!!!
        echo
fi

# Filter the INDELs from SEURAT

cat ${EXOMEPAIR}_seurat_snv.vcf | \
	java -jar -Xmx20g ${SNPSIFT} filter "(( QUAL >= 15 ) & ( DP1 >= 10 ) & ( DP2 >= 10 ) & ( AR1 <= 0.02 ) & ( AR2 >= 0.05 ))" > ${EXOMEPAIR}.seurat_snv_filt.vcf

if [ $? -ne 0 ]
then
        FAILED=1
	echo
        echo !!!!!!!!!!!!  Failed Filter the seurat snv by QUAL !!!!!!!!!!
        echo
fi

cat ${EXOMEPAIR}_seurat_indel.vcf | \
	java -jar -Xmx20g ${SNPSIFT} filter "( QUAL >= 25)" > ${EXOMEPAIR}.seurat_indel_filt.vcf

if [ $? -ne 0 ]
then
        FAILED=1
	echo
        echo !!!!!!!!!!!!  Failed Filter the seurat INDELS by QUAL !!!!!!!!!!
        echo
fi

# Filter the mutect vcf
# Check the Mutect sample order, leverage our standardized naming to find if its a C or T (STUDY_PATIENT_VISIT_SOURCE_FRACTIONincrement_ASSAYCODE_LIBRARY)
# Need to test the first charcter of the "FRACTIONincrement" is a C or T

FIRST_GENOTYPE_COLUMN=`grep -v "##" ${EXOMEPAIR}_MuTect_All.vcf | grep "#" | cut -f10`

if [ $? -ne 0 ]
then
        FAILED=1
	echo
        echo !!!!!!!!!!!!  Failed grep header !!!!!!!!!!
        echo
fi

if [ "${FIRST_GENOTYPE_COLUMN}" == "${TUMOR}" ]
then
	# Filter the MUTECT calls
        cat ${EXOMEPAIR}_MuTect_All.vcf | java -jar -Xmx20g ${SNPSIFT} filter "(( FILTER = 'PASS') & ( GEN[1].FA <= 0.02 ) & ( GEN[0].FA >= 0.05 ))" > ${EXOMEPAIR}.mutect_snv_filt.vcf

	if [ $? -ne 0 ]
	then
        	FAILED=1
		echo
        	echo !!!!!!!!!!!!  Failed Filter the mutect tumor section !!!!!!!!!!
        	echo
	fi
elif [ "${FIRST_GENOTYPE_COLUMN}" == "${CONTROL}" ]
then
        # Reorder the genotype columns and then filter
        awk '{ FS = "\t" ; OFS = "\t" ; print $1, $2, $3, $4, $5, $6, $7, $8, $9, $11, $10}' ${EXOMEPAIR}_MuTect_All.vcf | \
		java -jar ${SNPSIFT} filter "(( FILTER = 'PASS') & ( GEN[1].FA <= 0.02 ) & ( GEN[0].FA >= 0.05 ))" > ${EXOMEPAIR}.mutect_snv_filt.vcf

	if [ $? -ne 0 ]
	then
	        FAILED=1
		echo
        	echo !!!!!!!!!!!!  Failed Filter the mutect control section !!!!!!!!!!
        	echo
	fi
else
        #This should not happen
        echo "ERROR - The mutect vcf did not contain the tumor or the control listed first.  Something is wrong here."
	FAILED=1
	echo
        echo !!!!!!!!!!!!  Failed Filter the mutect bad !!!!!!!!!!
        echo
fi

sed -i 's/\t*$//g' ${EXOMEPAIR}.mutect_snv_filt.vcf

if [ $? -ne 0 ]
then
        FAILED=1
	echo
        echo !!!!!!!!!!!!  Failed sed part !!!!!!!!!!
        echo
fi

# Merge vcf from the 3 callers

${VCFMERGER} --dirscript ${VCFMERGER_DIR} \
	--seusnv ${EXOMEPAIR}.seurat_snv_filt.vcf \
	--seuindel ${EXOMEPAIR}.seurat_indel_filt.vcf \
	--slksnv ${EXOMEPAIR}.strelka.passed.somatic.snvs.vcf \
	--slkindel ${EXOMEPAIR}.strelka.passed.somatic.indels.vcf \
	--mtcsnv ${EXOMEPAIR}.mutect_snv_filt.vcf \
	--refgenfa ${REF} \
	--force \
	--outprefix ${EXOMEPAIR}

if [ $? -ne 0 ]
then
        FAILED=1
	echo
        echo !!!!!!!!!!!!  Failed vcfMerger !!!!!!!!!!
        echo
fi
	
# Sort the Merged VCF for GATK compatibility

${VCFSORTER} ${DICT} ${EXOMEPAIR}.merge.norm.vcf > ${EXOMEPAIR}.merge.sort.vcf

if [ $? -ne 0 ]
then
        FAILED=1
	echo
        echo !!!!!!!!!!!!  Failed sort merged vcf !!!!!!!!!!
        echo
fi

# Remove unwanted INFO keys

java -jar ${SNPSIFT} rmInfo ${EXOMEPAIR}.merge.sort.vcf \
	SEURAT_DNA_ALT_ALLELE_FORWARD_FRACTION \
	SEURAT_DNA_ALT_ALLELE_FORWARD \
	SEURAT_DNA_ALT_ALLELE_REVERSE_FRACTION \
	SEURAT_DNA_ALT_ALLELE_REVERSE \
        SEURAT_DNA_ALT_ALLELE_TOTAL_FRACTION \
        SEURAT_DNA_ALT_ALLELE_TOTAL \
        SEURAT_DNA_REF_ALLELE_FORWARD \
        SEURAT_DNA_REF_ALLELE_REVERSE \
        SEURAT_DNA_REF_ALLELE_TOTAL > ${EXOMEPAIR}.merge.sort.clean.vcf

if [ $? -ne 0 ]
then
        FAILED=1
	echo
        echo !!!!!!!!!!!!  Failed remove info keys  !!!!!!!!!!
        echo
fi

# Filter to target regions

cat ${EXOMEPAIR}.merge.sort.clean.vcf | java -jar ${SNPSIFT} intervals ${BEDFILE} > ${EXOMEPAIR}.merge.sort.clean.f2t.vcf

if [ $? -ne 0 ]
then
        FAILED=1
	echo
        echo !!!!!!!!!!!!  Failed filter to targets  !!!!!!!!!!
        echo
fi

# Annotate the merged vcf using GATK

java -Xmx24g -jar ${GATK}/GenomeAnalysisTK.jar -R ${REF} \
	-T VariantAnnotator \
	-nt 4 \
	-o ${EXOMEPAIR}.merge.sort.clean.f2t.ann.vcf \
	-U ALLOW_SEQ_DICT_INCOMPATIBILITY \
	--variant ${EXOMEPAIR}.merge.sort.clean.f2t.vcf \
	--dbsnp ${DBSNP} \
	--comp:EXAC ${EXAC} \
	--comp:NHLBI ${NHLBI} \
	--comp:1000G ${KG} \
	--comp:COSMIC_NC ${COSMICNC} \
	--comp:COSMIC_C ${COSMICC} \
	--disable_auto_index_creation_and_locking_when_reading_rods

if [ $? -ne 0 ]
then
	FAILED=1
	echo
	echo !!!!!!!!!!!!  Failed GATK annotation  !!!!!!!!!!
	echo
fi

# Determine if you need to run RNA allele counts or not 

if [[ "${RNAFLAG}" == YES ]]
then
	# Extract the variant positions to a list to provide to samtools for mpileup
	grep -v "#" ${EXOMEPAIR}.merge.sort.clean.f2t.ann.vcf | awk '{print $1"-"$2"-"$3"-"$4"-"$5}' > ${EXOMEPAIR}_positions.txt

	# This creates a file with each variant postion "chr-position", 
	# to test indel calling we will provide the 25bp upstream and downstream for pileup generation

	for line in `cat ${EXOMEPAIR}_positions.txt`
	do
		CHR=`echo ${line} | cut -d- -f1`
		POSITION=`echo ${line} | cut -d- -f2`

		awk -v var1="$CHR" -v var2="$POSITION" 'BEGIN{ OFS = "\t" ; for (i = 25; i >= 1; i-- ) { print var1"-"var2-i } ; { print var1"-"var2 } ; for (i = 1; i <= 25; i++ )  { print var1"-"var2+i }}' >> ${EXOMEPAIR}_positions_expanded.txt

	done

	# Since variants might be within 25 bp of one another ensure the positions file for mpileup only contains unique positions
	sort -n ${EXOMEPAIR}_positions_expanded.txt | uniq | awk '{gsub("-","\t",$0); print;}' > ${EXOMEPAIR}_positions_expanded_unique.txt

	# Generate samtools pileup
	# -B to turn of BAQ, no idea the effect on STAR BAM
	# -d to prevent under counting in high coverage regions that might exist in RNAseq
	# -Q to prevent exclusion of reads based on base quality - set to 5, no good reason other than 0 seemed to low and default of 13 to high
	${SAMTOOLS}/samtools mpileup -B -d10000000 -Q5 -f ${REF} -l ${EXOMEPAIR}_positions_expanded_unique.txt ${RNASAMPLE}.starAligned.final.bam > ${EXOMEPAIR}.RNA.pileup
	
	if [ $? -ne 0 ]
	then
        	FAILED=1
		echo
        	echo !!!!!!!!!!!!  Failed mpileup RNAflag yes  !!!!!!!!!!
        	echo
	fi

	# Create vcf file header
	cat ${RNA_VCF_HEADER} > ${EXOMEPAIR}.RNA.calls.vcf

	for line in `cat ${EXOMEPAIR}_positions.txt`
	do
		CHR=`echo ${line} | cut -d- -f1`
		POSITION=`echo ${line} | cut -d- -f2`
		REF=`echo ${line} | cut -d- -f4`
		ALT=`echo ${line} | cut -d- -f5`
		# Determine if the position is in the pileup
		LINE_COUNT=`awk -v var1="$CHR" -v var2="$POSITION" '{if($1 == var1 && $2 == var2) print $0}' ${EXOMEPAIR}.RNA.pileup | wc -l`	
		if [ ${LINE_COUNT} -eq 0 ]
		then
			# This means the position is not in the pileup at all so the REF and ALT counts are then 0
			echo -e ${CHR}"\t"${POSITION}"\t"".""\t"${REF}"\t"${ALT}"\t"".""\t"".""\t""RNA_REF_COUNT=0;RNA_ALT_COUNT=0;RNA_ALT_FREQ=0.00" >> ${EXOMEPAIR}.RNA.calls.vcf
		else
			# This means the line is in the pileup, it does not tell if its just because of spanning reads due to splicing
			echo -e ${CHR}"-"${POSITION}"-.-"${REF}"-"${ALT} >> ${EXOMEPAIR}.RNA.calls.ToQueryPileup
		fi
	done

	# Tabulate pileup with varscan - using VERY minimal criteria
	java -jar ${VARSCAN}/VarScan.v2.3.7.jar mpileup2cns ${EXOMEPAIR}.RNA.pileup \
		--min-coverage 1 \
		--min-reads2 1 \
        	--min-var-freq 0.02 \
        	--min-freq-for-hom 0.9 \
        	--min-avg-qual 5 \
        	--strand-filter 0 \
        	--output-vcf > ${EXOMEPAIR}.RNA.calls.varscan.vcf

	if [ $? -ne 0 ]
        then
                FAILED=1
		echo
        	echo !!!!!!!!!!!!  Failed tabulate pileup  !!!!!!!!!!
        	echo
        fi

	# Convert the VCF to at table
	java -jar ${SNPSIFT} extractFields ${EXOMEPAIR}.RNA.calls.varscan.vcf CHROM POS REF ALT GEN[0].SDP GEN[0].RD GEN[0].AD > ${EXOMEPAIR}.RNA.calls.varscan.pre.table
	
	if [ $? -ne 0 ]
        then
                FAILED=1
		echo
        	echo !!!!!!!!!!!!  Failed rna vcf to table !!!!!!!!!!
        	echo
        fi

	cut -f7 ${EXOMEPAIR}.RNA.calls.varscan.pre.table | sed -e 's/^$/./g' > ${EXOMEPAIR}.RNA.calls.varscan.pre2.table 
	cut -f6 ${EXOMEPAIR}.RNA.calls.varscan.pre.table | sed -e 's/^$/./g' > ${EXOMEPAIR}.RNA.calls.varscan.pre3.table
	cut -f5 ${EXOMEPAIR}.RNA.calls.varscan.pre.table | sed -e 's/^$/./g' > ${EXOMEPAIR}.RNA.calls.varscan.pre4.table
	cut -f4 ${EXOMEPAIR}.RNA.calls.varscan.pre.table | sed -e 's/^$/./g' > ${EXOMEPAIR}.RNA.calls.varscan.pre5.table
	cut -f3 ${EXOMEPAIR}.RNA.calls.varscan.pre.table | sed -e 's/^$/./g' > ${EXOMEPAIR}.RNA.calls.varscan.pre6.table
	cut -f2 ${EXOMEPAIR}.RNA.calls.varscan.pre.table | sed -e 's/^$/./g' > ${EXOMEPAIR}.RNA.calls.varscan.pre7.table
	cut -f1 ${EXOMEPAIR}.RNA.calls.varscan.pre.table | sed -e 's/^$/./g' > ${EXOMEPAIR}.RNA.calls.varscan.pre8.table

	paste ${EXOMEPAIR}.RNA.calls.varscan.pre8.table ${EXOMEPAIR}.RNA.calls.varscan.pre7.table ${EXOMEPAIR}.RNA.calls.varscan.pre6.table ${EXOMEPAIR}.RNA.calls.varscan.pre5.table ${EXOMEPAIR}.RNA.calls.varscan.pre4.table ${EXOMEPAIR}.RNA.calls.varscan.pre3.table ${EXOMEPAIR}.RNA.calls.varscan.pre2.table > ${EXOMEPAIR}.RNA.calls.varscan.table

	# Query each variant that was in the pileup for its varscan tabulation result
	for line in `cat ${EXOMEPAIR}.RNA.calls.ToQueryPileup`
	do
		CHR=`echo ${line} | cut -d- -f1`
		POSITION=`echo ${line} | cut -d- -f2`
		REF=`echo ${line} | cut -d- -f4`
		ALT=`echo ${line} | cut -d- -f5`
		
		# Determine if the position is in the tabulated varscan results
		LINE_COUNT=`awk -v var1="$CHR" -v var2="$POSITION" '{if($1 == var1 && $2 == var2) print $0}' ${EXOMEPAIR}.RNA.calls.varscan.table | wc -l`
		if [ ${LINE_COUNT} -eq 0 ]
		then
			# This means the position was in the pileup but did not make it to the varscan tabulated results #### NEED A REASON FOR THIS!!! ####
                	# These are added to the RNA calls as 0:0:0
			echo -e ${CHR}"\t"${POSITION}"\t"".""\t"${REF}"\t"${ALT}"\t"".""\t"".""\t""RNA_REF_COUNT=0;RNA_ALT_COUNT=0;RNA_ALT_FREQ=0.00" >> ${EXOMEPAIR}.RNA.calls.vcf
						
		elif [ ${LINE_COUNT} -eq 1 ]
		then
			# This means the line is in the VARSCAN TABULATION RESULTS
			# Now we need to confirm the alternate is actually present as expected and get the values to add to the vcf
                	# This gets the call line to a single
			awk -v var1="$CHR" -v var2="$POSITION" '{if($1 == var1 && $2 == var2) print $0}' ${EXOMEPAIR}.RNA.calls.varscan.table > temp.RNA.calls
			ALT_READS=`awk -v var1="$CHR" -v var2="$POSITION" -v var3="$REF" -v var4="$ALT" '{if($1 == var1 && $2 == var2 && $7 >= 1) print $0}' temp.RNA.calls | wc -l`
		
			if [ ${ALT_READS} -eq 1 ]
                	then		
				# Found expected variant allele
                        	REF_COUNT=`awk '{print $6}' temp.RNA.calls`
                        	ALT_COUNT=`awk '{print $7}' temp.RNA.calls`
                        	# Determine if the REF_COUNT variable is a "."
                        	if [ "${REF_COUNT}" == "." ]
                        	then
                        	        # Reset variable to 0
                        	        REF_COUNT=0
                        	        echo "Found Reference allele count was empty, forcing to 0"
                        	else
                        	        # Varialbe is good
                        	        echo "Found ${REF_COUNT} bases"
                        	        # Calculate the ratio
                        	fi

				ALT_FREQ=`echo "${REF_COUNT} ${ALT_COUNT}" | awk '{print $2/($1+$2)}'`
				# Print new lines to output
				echo -e ${CHR}"\t"${POSITION}"\t"".""\t"${REF}"\t"${ALT}"\t"".""\t"".""\t""RNA_REF_COUNT="${REF_COUNT}";RNA_ALT_COUNT="${ALT_COUNT}";RNA_ALT_FREQ="${ALT_FREQ} >> ${EXOMEPAIR}.RNA.calls.vcf

			elif [ ${ALT_READS} -eq 0 ]
			then
				# No alt reads found at all, just print ref counts and set alt to 0 ### NOT SO EASY ACTUALLY
				REF_COUNT=`awk '{print $6}' temp.RNA.calls`
				ALT_COUNT=0
				ALT_FREQ=0.00
			
				# Determine if the REF_COUNT variable is a "."
				if [ "${REF_COUNT}" == "." ]
                        	then
                                	# Reset variable to 0
                                	REF_COUNT=0
                                	echo "Found Reference allele count was empty, forcing to 0"
                                else
                                	# Varialbe is good
                                	echo "Found ${REF_COUNT} bases"
                        	fi	

				# Print new lines to output
				echo -e ${CHR}"\t"${POSITION}"\t"".""\t"${REF}"\t"${ALT}"\t"".""\t"".""\t""RNA_REF_COUNT="${REF_COUNT}";RNA_ALT_COUNT="${ALT_COUNT}";RNA_ALT_FREQ="${ALT_FREQ} >> ${EXOMEPAIR}.RNA.calls.vcf

			else
				# Unexpected
				echo "BAD THINGS ARE OCCURING.... chr${CHR}:${POSITION}"
			fi
		else
			# THIS IS UNEXPECTED
			echo "ERROR - The same variant is in the tabulated results more than once"
		fi
	done

	# Sort RNA calls VCF
	
	${VCFSORTER} ${DICT} ${EXOMEPAIR}.RNA.calls.vcf > ${EXOMEPAIR}.RNA.calls.sorted.vcf

	if [ $? -ne 0 ]
        then
                FAILED=1
		echo
        	echo !!!!!!!!!!!!  Failed sort rna vcf calls  !!!!!!!!!!
        	echo
        fi

	# Add RNA values to merged calls

	java -jar ${SNPSIFT} annotate ${EXOMEPAIR}.RNA.calls.sorted.vcf ${EXOMEPAIR}.merge.sort.clean.f2t.ann.vcf > ${EXOMEPAIR}.merge.sort.clean.f2t.ann.rna.vcf

	if [ $? -ne 0 ]
        then
                FAILED=1
		echo
        	echo !!!!!!!!!!!!  Failed  add rna to merged vcf !!!!!!!!!!
        	echo
        fi

	
else
	# CHange the name of the vcf to match the RNA one
	cp ${EXOMEPAIR}.merge.sort.clean.f2t.ann.vcf ${EXOMEPAIR}.merge.sort.clean.f2t.ann.rna.vcf

fi


###########

#{{{

# vcf file names

if [ "`grep -f ${HEADER} ${EXOMEPAIR}.merge.sort.clean.f2t.ann.rna.vcf | wc -l`" = 3  ]
then
	ARNA="Yes"
else
	ARNA="No"
fi


#############################
#
# Update All Transcripts VCF
#
#############################

# Remove info and obsolete header lines  

java -jar ${SNPSIFT} rmInfo ${EXOMEPAIR}.merge.sort.clean.f2t.ann.rna.vcf INDEL | \
	grep -v "##INFO=<ID=INDEL" | \
	grep -v "##FILTER=<ID=REJECT" > ${EXOMEPAIR}.merge.sort.clean.f2t.ann.rna.clean.vcf

if [ $? -ne 0 ]
then
	FAILED=1
	echo
        echo !!!!!!!!!!!!  Failed snpsift rminfo !!!!!!!!!!
        echo
fi

# Add dbNSFP annotations

java -jar ${SNPSIFT} dbnsfp \
	-v \
	-a \
	-db ${DBNSFP} \
	-f CADD_raw,CADD_raw_rankscore,CADD_phred,Interpro_domain,Polyphen2_HVAR_pred,GERP++_NR,GERP++_RS,LRT_score,MutationTaster_score,MutationAssessor_score,FATHMM_score,Polyphen2_HVAR_score,SIFT_score,Polyphen2_HDIV_score,MetaSVM_score,MetaSVM_rankscore,MetaSVM_pred,MetaLR_score,MetaLR_rankscore,MetaLR_pred,Reliability_index \
	${EXOMEPAIR}.merge.sort.clean.f2t.ann.rna.clean.vcf > ${EXOMEPAIR}.merge.sort.clean.f2t.ann.rna.clean.dbnsfp.vcf

if [ $? -ne 0 ]
then
	FAILED=1
	echo
        echo !!!!!!!!!!!!  Failed dbNSFP annotations  !!!!!!!!!!
        echo
fi

# Add SnpEff ANN annotations (default with version 4.2)

java -Xmx14g -jar ${SNPEFF} -c ${SNPEFFPATH}/snpEff.config -lof GRCh37.74 ${EXOMEPAIR}.merge.sort.clean.f2t.ann.rna.clean.dbnsfp.vcf | \
        java -jar ${SNPSIFT} varType - | \
        grep -v "##INFO=<ID=HOM,Number=A,Type=Flag,Description=\"Variant is homozygous\">" | \
        grep -v "##INFO=<ID=HET,Number=A,Type=Flag,Description=\"Variant is heterozygous\">" > ${EXOMEPAIR}.merge.sort.clean.f2t.ann.rna.clean.dbnsfp.se74lof.vcf

if [ $? -ne 0 ]
then
	FAILED=1
	echo
        echo !!!!!!!!!!!!  Failed snpeff annotations  !!!!!!!!!!
        echo
fi

java -Xmx14g -jar ${SNPEFF} -canon -c ${SNPEFFPATH}/snpEff.config -lof GRCh37.74 ${EXOMEPAIR}.merge.sort.clean.f2t.ann.rna.clean.dbnsfp.vcf | \
        java -jar ${SNPSIFT} varType - | \
        grep -v "##INFO=<ID=HOM,Number=A,Type=Flag,Description=\"Variant is homozygous\">" | \
        grep -v "##INFO=<ID=HET,Number=A,Type=Flag,Description=\"Variant is heterozygous\">" > ${EXOMEPAIR}.merge.sort.clean.f2t.ann.rna.clean.dbnsfp.se74lofcan.vcf

if [ $? -ne 0 ]
then
        FAILED=1
	echo
        echo !!!!!!!!!!!!  Failed canon anno  !!!!!!!!!!
        echo
fi

# Make final call list venn

${POST_MERGE_VENN} --vcf ${EXOMEPAIR}.merge.sort.clean.f2t.ann.rna.clean.dbnsfp.se74lofcan.vcf --dirscript ${VCFMERGER_DIR}/ --dir-snpsift ${SNPEFFPATH} --outprefix ${EXOMEPAIR}_finalVenn --maintitle ${EXOMEPAIR} --

if [ $? -ne 0 ]
then
        FAILED=1
	echo
        echo !!!!!!!!!!!!  Failed final call list venn  !!!!!!!!!!
        echo
fi

# Fix header lines that are missing, have spaces, dbNSFP errors, RNA_ALT_FREQ float problem.

if [ ${ARNA} = "No"  ]
then
	${BCFTOOLS} annotate -h ${HEADER} ${EXOMEPAIR}.merge.sort.clean.f2t.ann.rna.clean.dbnsfp.se74lof.vcf | \
	sed 's/'\'' '\"'/'\'\"'/g' | \
	sed 's/##FORMAT=<ID=AD,Number=R/##FORMAT=<ID=AD,Number=./g' | \
	sed 's|##FILTER=<ID=PASS>|##FILTER=<ID=PASS,Description="Accept as a confident somatic mutation">|g' | \
	sed 's|##FILTER=<ID=PASS,Description="All filters passed">|##FILTER=<ID=PASS,Description="Accept as a confident somatic mutation">|g' | \
	sed 's/dbNSFP_CADD_raw,Number=A/dbNSFP_CADD_raw,Number=./g' | \
	sed 's/dbNSFP_CADD_raw_rankscore,Number=A/dbNSFP_CADD_raw_rankscore,Number=./g' | \
	sed 's/dbNSFP_CADD_phred,Number=A/dbNSFP_CADD_phred,Number=./g' | \
	sed 's/dbNSFP_Polyphen2_HVAR_pred,Number=A/dbNSFP_Polyphen2_HVAR_pred,Number=./g' | \
	sed 's/dbNSFP_GERP___NR,Number=A/dbNSFP_GERP___NR,Number=./g' | \
	sed 's/dbNSFP_GERP___RS,Number=A/dbNSFP_GERP___RS,Number=./g' | \
	sed 's/dbNSFP_LRT_score,Number=A/dbNSFP_LRT_score,Number=./g' | \
	sed 's/dbNSFP_MutationTaster_score,Number=A/dbNSFP_MutationTaster_score,Number=./g' | \
	sed 's/dbNSFP_MutationAssessor_score,Number=A/dbNSFP_MutationAssessor_score,Number=./g' | \
	sed 's/dbNSFP_MetaLR_rankscore,Number=A/dbNSFP_MetaLR_rankscore,Number=./g' | \
	sed 's/dbNSFP_MetaSVM_rankscore,Number=A/dbNSFP_MetaSVM_rankscore,Number=./g' | \
	sed 's/dbNSFP_MetaSVM_pred,Number=A/dbNSFP_MetaSVM_pred,Number=./g' | \
	sed 's/dbNSFP_MetaSVM_score,Number=A/dbNSFP_MetaSVM_score,Number=./g' | \
	sed 's/dbNSFP_MetaLR_pred,Number=A/dbNSFP_MetaLR_pred,Number=./g' | \
	sed 's/dbNSFP_MetaLR_score,Number=A/dbNSFP_MetaLR_score,Number=./g' | \
	sed 's/dbNSFP_Reliability_index,Number=A/dbNSFP_Reliability_index,Number=./g' | \
	sed 's/dbNSFP_Polyphen2_HVAR_score,Number=A/dbNSFP_Polyphen2_HVAR_score,Number=./g' | \
        sed 's/dbNSFP_SIFT_score,Number=A/dbNSFP_SIFT_score,Number=./g' | \
        sed 's/dbNSFP_Polyphen2_HVAR_pred,Number=A/dbNSFP_Polyphen2_HVAR_pred,Number=./g' | \
        sed 's/dbNSFP_LRT_score,Number=A/dbNSFP_LRT_score,Number=./g' | \
        sed 's/dbNSFP_FATHMM_score,Number=A/dbNSFP_FATHMM_score,Number=./g' | \
        sed 's/dbNSFP_Polyphen2_HDIV_score,Number=A/dbNSFP_Polyphen2_HDIV_score,Number=./g' | \
        sed 's/dbNSFP_Interpro_domain,Number=A/dbNSFP_Interpro_domain,Number=./g' | \
	sed 's/##INFO=<ID=SNP,Number=A,Type=Flag,Description="Variant is a SNP">/##INFO=<ID=SNP,Number=0,Type=Flag,Description="Variant is a SNP">/g' | \
	sed 's/##INFO=<ID=MNP,Number=A,Type=Flag,Description="Variant is a MNP">/##INFO=<ID=MNP,Number=0,Type=Flag,Description="Variant is a MNP">/g' | \
	sed 's/##INFO=<ID=INS,Number=A,Type=Flag,Description="Variant is a INS">/##INFO=<ID=INS,Number=0,Type=Flag,Description="Variant is a INS">/g' | \
	sed 's/##INFO=<ID=DEL,Number=A,Type=Flag,Description="Variant is a DEL">/##INFO=<ID=DEL,Number=0,Type=Flag,Description="Variant is a DEL">/g' | \
	sed 's/##INFO=<ID=MIXED,Number=A,Type=Flag,Description="Variant is a MIXED">/##INFO=<ID=MIXED,Number=0,Type=Flag,Description="Variant is a MIXED">/g' > ${EXOMEPAIR}.merge.sort.clean.f2t.ann.rna.clean.dbnsfp.se74lof.fix.vcf 

	if [ $? -ne 0 ]
	then
		FAILED=1
		echo
        	echo !!!!!!!!!!!!  Failed ARNA no vcf fix  !!!!!!!!!!
        	echo
	else
		VALIDATOR="`${VCFVALIDATOR} ${EXOMEPAIR}.merge.sort.clean.f2t.ann.rna.clean.dbnsfp.se74lof.fix.vcf 2>&1 | wc -l`"
                if [ ${VALIDATOR} = "0" ]
                then
                        mv ${EXOMEPAIR}.merge.sort.clean.f2t.ann.rna.clean.dbnsfp.se74lof.fix.vcf ${EXOMEPAIR}.merged.allTranscripts.rna.final.vcf
                else
			FAILED=1
			echo
        		echo !!!!!!!!!!!!  Failed validator ARNA no on all vcf  !!!!!!!!!!
        		echo
                fi	
	fi

	${BCFTOOLS} annotate -h ${HEADER} ${EXOMEPAIR}.merge.sort.clean.f2t.ann.rna.clean.dbnsfp.se74lofcan.vcf | \
        sed 's/'\'' '\"'/'\'\"'/g' | \
		sed 's/##FORMAT=<ID=AD,Number=R/##FORMAT=<ID=AD,Number=./g' | \
	sed 's|##FILTER=<ID=PASS>|##FILTER=<ID=PASS,Description="Accept as a confident somatic mutation">|g' | \
        sed 's|##FILTER=<ID=PASS,Description="All filters passed">|##FILTER=<ID=PASS,Description="Accept as a confident somatic mutation">|g' | \
        sed 's/dbNSFP_CADD_raw,Number=A/dbNSFP_CADD_raw,Number=./g' | \
        sed 's/dbNSFP_CADD_raw_rankscore,Number=A/dbNSFP_CADD_raw_rankscore,Number=./g' | \
        sed 's/dbNSFP_CADD_phred,Number=A/dbNSFP_CADD_phred,Number=./g' | \
        sed 's/dbNSFP_Polyphen2_HVAR_pred,Number=A/dbNSFP_Polyphen2_HVAR_pred,Number=./g' | \
        sed 's/dbNSFP_GERP___NR,Number=A/dbNSFP_GERP___NR,Number=./g' | \
        sed 's/dbNSFP_GERP___RS,Number=A/dbNSFP_GERP___RS,Number=./g' | \
        sed 's/dbNSFP_LRT_score,Number=A/dbNSFP_LRT_score,Number=./g' | \
        sed 's/dbNSFP_MutationTaster_score,Number=A/dbNSFP_MutationTaster_score,Number=./g' | \
        sed 's/dbNSFP_MutationAssessor_score,Number=A/dbNSFP_MutationAssessor_score,Number=./g' | \
        sed 's/dbNSFP_MetaLR_rankscore,Number=A/dbNSFP_MetaLR_rankscore,Number=./g' | \
        sed 's/dbNSFP_MetaSVM_rankscore,Number=A/dbNSFP_MetaSVM_rankscore,Number=./g' | \
        sed 's/dbNSFP_MetaSVM_pred,Number=A/dbNSFP_MetaSVM_pred,Number=./g' | \
        sed 's/dbNSFP_MetaSVM_score,Number=A/dbNSFP_MetaSVM_score,Number=./g' | \
        sed 's/dbNSFP_MetaLR_pred,Number=A/dbNSFP_MetaLR_pred,Number=./g' | \
        sed 's/dbNSFP_MetaLR_score,Number=A/dbNSFP_MetaLR_score,Number=./g' | \
        sed 's/dbNSFP_Reliability_index,Number=A/dbNSFP_Reliability_index,Number=./g' | \
        sed 's/dbNSFP_Polyphen2_HVAR_score,Number=A/dbNSFP_Polyphen2_HVAR_score,Number=./g' | \
        sed 's/dbNSFP_SIFT_score,Number=A/dbNSFP_SIFT_score,Number=./g' | \
        sed 's/dbNSFP_Polyphen2_HVAR_pred,Number=A/dbNSFP_Polyphen2_HVAR_pred,Number=./g' | \
        sed 's/dbNSFP_LRT_score,Number=A/dbNSFP_LRT_score,Number=./g' | \
        sed 's/dbNSFP_FATHMM_score,Number=A/dbNSFP_FATHMM_score,Number=./g' | \
        sed 's/dbNSFP_Polyphen2_HDIV_score,Number=A/dbNSFP_Polyphen2_HDIV_score,Number=./g' | \
        sed 's/dbNSFP_Interpro_domain,Number=A/dbNSFP_Interpro_domain,Number=./g' | \
	sed 's/##INFO=<ID=SNP,Number=A,Type=Flag,Description="Variant is a SNP">/##INFO=<ID=SNP,Number=0,Type=Flag,Description="Variant is a SNP">/g' | \
        sed 's/##INFO=<ID=MNP,Number=A,Type=Flag,Description="Variant is a MNP">/##INFO=<ID=MNP,Number=0,Type=Flag,Description="Variant is a MNP">/g' | \
        sed 's/##INFO=<ID=INS,Number=A,Type=Flag,Description="Variant is a INS">/##INFO=<ID=INS,Number=0,Type=Flag,Description="Variant is a INS">/g' | \
        sed 's/##INFO=<ID=DEL,Number=A,Type=Flag,Description="Variant is a DEL">/##INFO=<ID=DEL,Number=0,Type=Flag,Description="Variant is a DEL">/g' | \
        sed 's/##INFO=<ID=MIXED,Number=A,Type=Flag,Description="Variant is a MIXED">/##INFO=<ID=MIXED,Number=0,Type=Flag,Description="Variant is a MIXED">/g' > ${EXOMEPAIR}.merge.sort.clean.f2t.ann.rna.clean.dbnsfp.se74lofcan.fix.vcf

        if [ $? -ne 0 ]
        then
                FAILED=1
		echo
                echo !!!!!!!!!!!!  Failed ARNA no vcf fix on canon  !!!!!!!!!!
                echo
        else
                VALIDATOR="`${VCFVALIDATOR} ${EXOMEPAIR}.merge.sort.clean.f2t.ann.rna.clean.dbnsfp.se74lofcan.fix.vcf 2>&1 | wc -l`"
                if [ ${VALIDATOR} = "0" ]
                then
                        mv ${EXOMEPAIR}.merge.sort.clean.f2t.ann.rna.clean.dbnsfp.se74lofcan.fix.vcf ${EXOMEPAIR}.merged.canonicalOnly.rna.final.vcf
                else
                        FAILED=1
			echo
                        echo !!!!!!!!!!!!  Failed validator ARNA no on canon vcf  !!!!!!!!!!
                        echo
                fi
        fi
else
	sed 's/'\'' '\"'/'\'\"'/g' ${EXOMEPAIR}.merge.sort.clean.f2t.ann.rna.clean.dbnsfp.se74lof.vcf | \
	sed 's|##FILTER=<ID=PASS>|##FILTER=<ID=PASS,Description="Accept as a confident somatic mutation">|g' | \
	sed 's/dbNSFP_CADD_raw,Number=A/dbNSFP_CADD_raw,Number=./g' | \
		sed 's/##FORMAT=<ID=AD,Number=R/##FORMAT=<ID=AD,Number=./g' | \
        sed 's/dbNSFP_CADD_raw_rankscore,Number=A/dbNSFP_CADD_raw_rankscore,Number=./g' | \
        sed 's/dbNSFP_CADD_phred,Number=A/dbNSFP_CADD_phred,Number=./g' | \
        sed 's/dbNSFP_Polyphen2_HVAR_pred,Number=A/dbNSFP_Polyphen2_HVAR_pred,Number=./g' | \
        sed 's/dbNSFP_GERP___NR,Number=A/dbNSFP_GERP___NR,Number=./g' | \
        sed 's/dbNSFP_GERP___RS,Number=A/dbNSFP_GERP___RS,Number=./g' | \
        sed 's/dbNSFP_LRT_score,Number=A/dbNSFP_LRT_score,Number=./g' | \
        sed 's/dbNSFP_MutationTaster_score,Number=A/dbNSFP_MutationTaster_score,Number=./g' | \
        sed 's/dbNSFP_MutationAssessor_score,Number=A/dbNSFP_MutationAssessor_score,Number=./g' | \
        sed 's/dbNSFP_MetaLR_rankscore,Number=A/dbNSFP_MetaLR_rankscore,Number=./g' | \
        sed 's/dbNSFP_MetaSVM_rankscore,Number=A/dbNSFP_MetaSVM_rankscore,Number=./g' | \
        sed 's/dbNSFP_MetaSVM_pred,Number=A/dbNSFP_MetaSVM_pred,Number=./g' | \
        sed 's/dbNSFP_MetaSVM_score,Number=A/dbNSFP_MetaSVM_score,Number=./g' | \
        sed 's/dbNSFP_MetaLR_pred,Number=A/dbNSFP_MetaLR_pred,Number=./g' | \
        sed 's/dbNSFP_MetaLR_score,Number=A/dbNSFP_MetaLR_score,Number=./g' | \
        sed 's/dbNSFP_Reliability_index,Number=A/dbNSFP_Reliability_index,Number=./g' | \
        sed 's/dbNSFP_Polyphen2_HVAR_score,Number=A/dbNSFP_Polyphen2_HVAR_score,Number=./g' | \
        sed 's/dbNSFP_SIFT_score,Number=A/dbNSFP_SIFT_score,Number=./g' | \
        sed 's/dbNSFP_Polyphen2_HVAR_pred,Number=A/dbNSFP_Polyphen2_HVAR_pred,Number=./g' | \
        sed 's/dbNSFP_LRT_score,Number=A/dbNSFP_LRT_score,Number=./g' | \
        sed 's/dbNSFP_FATHMM_score,Number=A/dbNSFP_FATHMM_score,Number=./g' | \
        sed 's/dbNSFP_Polyphen2_HDIV_score,Number=A/dbNSFP_Polyphen2_HDIV_score,Number=./g' | \
        sed 's/dbNSFP_Interpro_domain,Number=A/dbNSFP_Interpro_domain,Number=./g' | \
	sed 's/##INFO=<ID=SNP,Number=A,Type=Flag,Description="Variant is a SNP">/##INFO=<ID=SNP,Number=0,Type=Flag,Description="Variant is a SNP">/g' | \
        sed 's/##INFO=<ID=MNP,Number=A,Type=Flag,Description="Variant is a MNP">/##INFO=<ID=MNP,Number=0,Type=Flag,Description="Variant is a MNP">/g' | \
        sed 's/##INFO=<ID=INS,Number=A,Type=Flag,Description="Variant is a INS">/##INFO=<ID=INS,Number=0,Type=Flag,Description="Variant is a INS">/g' | \
        sed 's/##INFO=<ID=DEL,Number=A,Type=Flag,Description="Variant is a DEL">/##INFO=<ID=DEL,Number=0,Type=Flag,Description="Variant is a DEL">/g' | \
        sed 's/##INFO=<ID=MIXED,Number=A,Type=Flag,Description="Variant is a MIXED">/##INFO=<ID=MIXED,Number=0,Type=Flag,Description="Variant is a MIXED">/g' > ${EXOMEPAIR}.merge.sort.clean.f2t.ann.rna.clean.dbnsfp.se74lof.fix.vcf

	if [ $? -ne 0 ]
	then
		FAILED=1
		echo
		echo !!!!!!!!!!!!  Failed ARNA yes vcf fix  !!!!!!!!!!
		echo
	else
		VALIDATOR="`${VCFVALIDATOR} ${EXOMEPAIR}.merge.sort.clean.f2t.ann.rna.clean.dbnsfp.se74lof.fix.vcf 2>&1 | wc -l`"
		if [ ${VALIDATOR} = "0" ]
		then
			mv ${EXOMEPAIR}.merge.sort.clean.f2t.ann.rna.clean.dbnsfp.se74lof.fix.vcf ${EXOMEPAIR}.merged.allTranscripts.rna.final.vcf
		else
			FAILED=1
			echo
			echo !!!!!!!!!!!!  Failed validator ARNA yes  !!!!!!!!!!
			echo
		fi        
	fi

	sed 's/'\'' '\"'/'\'\"'/g' ${EXOMEPAIR}.merge.sort.clean.f2t.ann.rna.clean.dbnsfp.se74lofcan.vcf | \
        sed 's|##FILTER=<ID=PASS>|##FILTER=<ID=PASS,Description="Accept as a confident somatic mutation">|g' | \
        sed 's/dbNSFP_CADD_raw,Number=A/dbNSFP_CADD_raw,Number=./g' | \
		sed 's/##FORMAT=<ID=AD,Number=R/##FORMAT=<ID=AD,Number=./g' | \
        sed 's/dbNSFP_CADD_raw_rankscore,Number=A/dbNSFP_CADD_raw_rankscore,Number=./g' | \
        sed 's/dbNSFP_CADD_phred,Number=A/dbNSFP_CADD_phred,Number=./g' | \
        sed 's/dbNSFP_Polyphen2_HVAR_pred,Number=A/dbNSFP_Polyphen2_HVAR_pred,Number=./g' | \
        sed 's/dbNSFP_GERP___NR,Number=A/dbNSFP_GERP___NR,Number=./g' | \
        sed 's/dbNSFP_GERP___RS,Number=A/dbNSFP_GERP___RS,Number=./g' | \
        sed 's/dbNSFP_LRT_score,Number=A/dbNSFP_LRT_score,Number=./g' | \
        sed 's/dbNSFP_MutationTaster_score,Number=A/dbNSFP_MutationTaster_score,Number=./g' | \
        sed 's/dbNSFP_MutationAssessor_score,Number=A/dbNSFP_MutationAssessor_score,Number=./g' | \
        sed 's/dbNSFP_MetaLR_rankscore,Number=A/dbNSFP_MetaLR_rankscore,Number=./g' | \
        sed 's/dbNSFP_MetaSVM_rankscore,Number=A/dbNSFP_MetaSVM_rankscore,Number=./g' | \
        sed 's/dbNSFP_MetaSVM_pred,Number=A/dbNSFP_MetaSVM_pred,Number=./g' | \
        sed 's/dbNSFP_MetaSVM_score,Number=A/dbNSFP_MetaSVM_score,Number=./g' | \
        sed 's/dbNSFP_MetaLR_pred,Number=A/dbNSFP_MetaLR_pred,Number=./g' | \
        sed 's/dbNSFP_MetaLR_score,Number=A/dbNSFP_MetaLR_score,Number=./g' | \
        sed 's/dbNSFP_Reliability_index,Number=A/dbNSFP_Reliability_index,Number=./g' | \
        sed 's/dbNSFP_Polyphen2_HVAR_score,Number=A/dbNSFP_Polyphen2_HVAR_score,Number=./g' | \
        sed 's/dbNSFP_SIFT_score,Number=A/dbNSFP_SIFT_score,Number=./g' | \
        sed 's/dbNSFP_Polyphen2_HVAR_pred,Number=A/dbNSFP_Polyphen2_HVAR_pred,Number=./g' | \
        sed 's/dbNSFP_LRT_score,Number=A/dbNSFP_LRT_score,Number=./g' | \
        sed 's/dbNSFP_FATHMM_score,Number=A/dbNSFP_FATHMM_score,Number=./g' | \
        sed 's/dbNSFP_Polyphen2_HDIV_score,Number=A/dbNSFP_Polyphen2_HDIV_score,Number=./g' | \
        sed 's/dbNSFP_Interpro_domain,Number=A/dbNSFP_Interpro_domain,Number=./g' | \
	sed 's/##INFO=<ID=SNP,Number=A,Type=Flag,Description="Variant is a SNP">/##INFO=<ID=SNP,Number=0,Type=Flag,Description="Variant is a SNP">/g' | \
        sed 's/##INFO=<ID=MNP,Number=A,Type=Flag,Description="Variant is a MNP">/##INFO=<ID=MNP,Number=0,Type=Flag,Description="Variant is a MNP">/g' | \
        sed 's/##INFO=<ID=INS,Number=A,Type=Flag,Description="Variant is a INS">/##INFO=<ID=INS,Number=0,Type=Flag,Description="Variant is a INS">/g' | \
        sed 's/##INFO=<ID=DEL,Number=A,Type=Flag,Description="Variant is a DEL">/##INFO=<ID=DEL,Number=0,Type=Flag,Description="Variant is a DEL">/g' | \
        sed 's/##INFO=<ID=MIXED,Number=A,Type=Flag,Description="Variant is a MIXED">/##INFO=<ID=MIXED,Number=0,Type=Flag,Description="Variant is a MIXED">/g' > ${EXOMEPAIR}.merge.sort.clean.f2t.ann.rna.clean.dbnsfp.se74lofcan.fix.vcf

        if [ $? -ne 0 ]
        then
                FAILED=1
		echo
                echo !!!!!!!!!!!!  Failed ARNA yes vcf fix on canon  !!!!!!!!!!
                echo
        else
                VALIDATOR="`${VCFVALIDATOR} ${EXOMEPAIR}.merge.sort.clean.f2t.ann.rna.clean.dbnsfp.se74lofcan.fix.vcf 2>&1 | wc -l`"
                if [ ${VALIDATOR} = "0" ]
                then
                        mv ${EXOMEPAIR}.merge.sort.clean.f2t.ann.rna.clean.dbnsfp.se74lofcan.fix.vcf ${EXOMEPAIR}.merged.canonicalOnly.rna.final.vcf
                else
                        FAILED=1
			echo
                        echo !!!!!!!!!!!!  Failed validator ARNA yes on canon  !!!!!!!!!!
                        echo
                fi
        fi
fi

fails=0

if [ ${FAILED} = 0 ]
then
	ssh ${USER}@${DATAMOVERIP} "rsync -r ${DIR}/Venns ${STARTDIR}/${PATIENT_NAME}/vcfMerger_pegasus/${EXOMEPAIR}/"

	if [ $? -ne 0 ]
        then
                fails=1
        fi
	
	ssh ${USER}@${DATAMOVERIP} "rsync ${DIR}/${EXOMEPAIR}.merge.vcf ${STARTDIR}/${PATIENT_NAME}/vcfMerger_pegasus/${EXOMEPAIR}/"

	if [ $? -ne 0 ]
        then
                fails=1
        fi

	for extToCopy in tsv png
	do
		for fileToSave in `find ${DIR} -maxdepth 1 -name "*$extToCopy"`
		do
			ssh ${USER}@${DATAMOVERIP} "rsync $fileToSave ${STARTDIR}/${PATIENT_NAME}/vcfMerger_pegasus/${EXOMEPAIR}/"
			
			if [ $? -ne 0 ]
                        then
                                fails=1
                        fi
		done
	done

	ssh ${USER}@${DATAMOVERIP} "rsync ${DIR}/${EXOMEPAIR}.merged.canonicalOnly.rna.final.vcf ${STARTDIR}/${PATIENT_NAME}/vcfMerger_pegasus/${EXOMEPAIR}/"
	
	if [ $? = 0 ] && [ $fails -eq 0 ] 
	then
		ssh ${USER}@${DATAMOVERIP} "rsync ${DIR}/${EXOMEPAIR}.merged.allTranscripts.rna.final.vcf ${STARTDIR}/${PATIENT_NAME}/vcfMerger_pegasus/${EXOMEPAIR}/"

		if [ $? = 0 ] && [ $fails -eq 0 ]
		then
			ssh ${USER}@${DATAMOVERIP} "rm ${STARTDIR}/${PATIENT_NAME}/vcfMerger_pegasus/${EXOMEPAIR}/SnpEFF_ANN_In_Progress ; \
				touch ${STARTDIR}/${PATIENT_NAME}/vcfMerger_pegasus/${EXOMEPAIR}/SnpEFF_ANN_Complete ; \
				COMPLETE=0 ; \
				for XPAIR in \`echo ${EXOMEPAIRS} | sed 's/@/\n/g'\`
				do 
					if [ ! -f ${STARTDIR}/${PATIENT_NAME}/vcfMerger_pegasus/\${XPAIR}/SnpEFF_ANN_Complete ] 
					then 
						COMPLETE=1  
					fi  
				done ; \
				if [ \${COMPLETE} = 0 ]
				then
					rm ${STARTDIR}/${PATIENT_NAME}/SnpEFF_ANN_In_Progress
					touch ${STARTDIR}/${PATIENT_NAME}/SnpEFF_ANN_Complete
				fi"
		else
			ssh ${USER}@${DATAMOVERIP} "rm ${STARTDIR}/${PATIENT_NAME}/vcfMerger_pegasus/${EXOMEPAIR}/SnpEFF_ANN_In_Progress ; \
				echo Failed to rsync the allTranscripts merged vcf back to isilon for ${EXOMEPAIR} >> ${STARTDIR}/${PATIENT_NAME}/vcfMerger_pegasus/${EXOMEPAIR}/SnpEFF_ANN_Fail ; \
				echo Failed to rsync the allTranscripts merged vcf back to isilon for ${EXOMEPAIR} | mailx -s "Post Medusa Processing Failed" ${USER}@tgen.org "	
		fi
	else
		ssh ${USER}@${DATAMOVERIP} "rm ${STARTDIR}/${PATIENT_NAME}/vcfMerger_pegasus/${EXOMEPAIR}/SnpEFF_ANN_In_Progress ; \
                	echo Failed to rsync the canonical merged vcf back to isilon for ${EXOMEPAIR} >> ${STARTDIR}/${PATIENT_NAME}/vcfMerger_pegasus/${EXOMEPAIR}/SnpEFF_ANN_Fail ; \
                        echo Failed to rsync the canonical merged vcf back to isilon for ${EXOMEPAIR} | mailx -s "Post Medusa Processing Failed" ${USER}@tgen.org "
	
	fi
else
	ssh ${USER}@${DATAMOVERIP} "rm ${STARTDIR}/${PATIENT_NAME}/vcfMerger_pegasus/${EXOMEPAIR}/SnpEFF_ANN_In_Progress ; \
        	echo SnpEFF failed somewere in the process for ${EXOMEPAIR} >> ${STARTDIR}/${PATIENT_NAME}/vcfMerger_pegasus/${EXOMEPAIR}/SnpEFF_ANN_Fail ; \
                echo SnpEFF failed somewere in the process for ${EXOMEPAIR} | mailx -s "Post Medusa Processing Failed" ${USER}@tgen.org "

fi

 # }}} 


