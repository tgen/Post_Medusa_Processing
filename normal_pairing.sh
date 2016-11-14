#!/bin/bash

CONSTANTS=/home/tgenref/pecan/post_central_pipe_processing/post_medusa_V1.0/constants.txt
STARTDIR=`pwd`
USER=`whoami`
EXOMEPBS=/home/tgenref/pecan/post_central_pipe_processing/post_medusa_V1.0/jobScripts/ngs_cna2015_Exome_normal_pairing_V1.pbs

for PAIR in `cat $1`
do
	cd ${STARTDIR}
	LIBTYPE="Exome"
        NORMALSAMPLE="`echo ${PAIR} | cut -d- -f1`"
        TUMORSAMPLE="`echo ${PAIR} | cut -d- -f2`"
	PATIENT_NAME="`echo ${TUMORSAMPLE} | cut -d_ -f1-2`"
        EXOMEPAIR="${NORMALSAMPLE}-${TUMORSAMPLE}"
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

	mkdir ${PATIENT_NAME}
	cd ${PATIENT_NAME}

	qsub -q overflow -v OFILE=${EXOMEPAIR},BEDFILE=${bedFile},CNAEXOMETARGET="${CNAEXOMETARGET}",NASSAY="${NASSAY}",TASSAY="${TASSAY}",LIBTYPE="${LIBTYPE}",USER="${USER}",PATIENT_NAME="${PATIENT_NAME}",NORMALSAMPLE="${NORMALSAMPLE}",TUMORSAMPLE="${TUMORSAMPLE}",NORMALDAT="${NORMALDAT}",TUMORDAT="${TUMORDAT}",EXOMEPAIR="${EXOMEPAIR}" ${EXOMEPBS}

done
