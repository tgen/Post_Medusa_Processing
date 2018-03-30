#!/bin/bash

for line in `ls -d */ |cut -d/ -f1`
do 
	cd $line
	
	if [ ! -f rename_stats_complete.txt ]
	then
	
	for row in `grep "FQ=" ${line}.config | rev | cut -d/ -f1 | rev | awk -F'_' '{ OFS = "_" ; print $1,$2,$3,$4,$5,$6,$7,$8 }' | sort | uniq`
	do 
		SAMPLE_NAME=`echo $row | awk -F'_' '{ OFS = "_" ; print $1,$2,$3,$4,$5,$6,$7 }'` 
		SAMPLE_ID=`echo $row | awk -F'_' '{ OFS = "_" ; print $8 }'`
		ASSAY=`echo $row | awk -F'_' '{ OFS = "_" ; print $7 }'`

		sed "s/${SAMPLE_NAME}	${ASSAY}/${row}	${SAMPLE_ID}/g" stats/Summary_allStats.txt > stats/Summary_allStats1.txt
		
		mv stats/Summary_allStats1.txt stats/Summary_allStats.txt		
		for file in `ls stats/${SAMPLE_NAME}*`
		do 
			ENDING=`echo $file | cut -d. -f2-`
			mv ${file} stats/${row}.${ENDING}
		done
	done	
	
	touch rename_stats_complete.txt
	
	fi
	
	cd ..
done
