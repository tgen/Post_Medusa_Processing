#!/bin/sh

#####################################
##
##   Define required variables
##
#####################################

# Folder to copy constitutional exome files
DESTINATION_DIR=/scratch/MMRF/haplotypeCaller_Results/gVCFs

#####################################
##
##   Parameterized Script
##
#####################################

echo
echo
echo Starting Process
echo
echo

# Loop through each folder in current directory looking for constitutional exome BAM files
for line in `ls -d */`
do
    echo
    echo "################################################"
    echo "Processing folder ${line}"

    # Move into folder
    cd ${line}

    # Check if gVCF creation is complete
    if [ -f Exome_gVCF_Created ]
    then
        echo "gVCF Already Created..."
        cd ..
    else
        # Determine if a constitutional exome exists
        echo "Checking folder for constitutional Exome files:"

	    # Ensure the Exome BAM list file is not present
	    if [ -e "EXOME_BAM_LIST" ]
	    then
		    echo ERROR - Found Exome_BAM_LIST file
		    echo THIS IS NOT EXPECTED, EXITING
		    exit 1
	    fi

	    # Test which folders are present if found execute find and add BAM files to Exome list
        if [ -e "TSE61" ]
        then
                echo Found TSE61 Folder
                find TSE61 -name "*.bam" > Exome_BAM_LIST
        else
                echo No TSE61 Folder
        fi

        if [ -e "KAS5U" ]
        then
                echo Found KAS5U Folder
                find KAS5U -name "*.bam" >> Exome_BAM_LIST
        else
                echo No KAS5U Folder
        fi

        if [ -e "KBS5U" ]
        then
                echo Found KBS5U Folder
                find KBS5U -name "*.bam" >> Exome_BAM_LIST
        else
                echo No KBS5U Folder
        fi

        if [ -e "KHS5U" ]
        then
                echo Found KHS5U Folder
                find KHS5U -name "*.bam" >> Exome_BAM_LIST
        else
                echo No KHS5U Folder
        fi

        # Test if any DNA BAM Files were found if not don't do all the subsequent steps
	    if [ -e "Exome_BAM_LIST" ]
	    then
            echo
            echo Found the following Exome BAM files:
            cat Exome_BAM_LIST
            echo
            echo Determining which file is the single expected constitutional exome

            # Find the single expected constitutional BAM file
            CONSTITUTIONAL_COUNT=0
            for row in `cat Exome_BAM_LIST`
            do
                echo "--------------------------------------------"
                # Capture folder names and BAM name
                ASSAY_FOLDER=`echo ${row} | cut -d/ -f1`
                EXOME_FOLDER=`echo ${row} | cut -d/ -f2`
                BAM_NAME=`echo ${row} | cut -d/ -f3`
                # Capture SubGroup Flag C=Constitutional/Normal/Germline and T=Tumor
                SUBGROUP=`echo ${row} | cut -d_ -f6 | cut -c1`
                if [ $SUBGROUP = C ]
                then
                    echo Found Constitutional Exome ${BAM_NAME}
                    if [ $CONSTITUTIONAL_COUNT -ne 0 ]
                    then
                        echo ERROR - Found more than one Constitutional Exome
                        echo exiting...
                        exit 1
                    else
                        CONSTITUTIONAL_COUNT=1
                        echo ${row} > CON_EXOME
                    fi
                elif [ $SUBGROUP = T ]
                then
                    echo Found Tumor Exome ${BAM_NAME}
                else
                    echo Unexpected Event, exiting....
                    echo Unrecognized SUBGROUP Variable - ${SUBGROUP}
                    exit 1
                fi
            done

            # Determine if a Constitutional Exome was Found
            if [ -e "CON_EXOME" ]
            then
                # Extract needed information from path for rsync processes
                echo
                echo STARTING RSYNC PROCESSES
                echo
                ASSAY_FOLDER=`cat CON_EXOME | cut -d/ -f1`
                EXOME_FOLDER=`cat CON_EXOME | cut -d/ -f2`
                BAM_NAME=`cat CON_EXOME | cut -d/ -f3`
                BAM_BASENAME=`basename ${BAM_NAME} ".bam"`

                # Copy the BAI file
                echo "--------------------------------------------"
                echo
                # Create Progress Marker
	            touch InProgress_rsync_ConExomeBAI

	            rsync -hu --progress ${ASSAY_FOLDER}/${EXOME_FOLDER}/${BAM_BASENAME}.bai ${DESTINATION_DIR}

	            # Error Capture
	            if [ "$?" = "0" ]
	            then
	                echo "BAI rsync complete"
	                echo
		            rm InProgress_rsync_ConExomeBAI
	            else
	                echo "FAIL --- BAI rsync --- FAIL"
		            touch FAIL_rsync_ConExomeBAI
		            rm InProgress_rsync_ConExomeBAI
		            exit 1
	            fi

                # Copy the BAM file
                echo "--------------------------------------------"
                echo
                # Create Progress Marker
	            touch InProgress_rsync_ConExomeBAM

	            rsync -hu --progress ${ASSAY_FOLDER}/${EXOME_FOLDER}/${BAM_BASENAME}.bam ${DESTINATION_DIR}

	            # Error Capture
	            if [ "$?" = "0" ]
	            then
	                echo "BAM rsync complete"
	                echo
		            rm InProgress_rsync_ConExomeBAM
	            else
	                echo "FAIL --- BAM rsync --- FAIL"
		            touch FAIL_rsync_ConExomeBAM
		            rm InProgress_rsync_ConExomeBAM
		            exit 1
	            fi

	            # Clean-up the temp constitutional exome file list
	            rm CON_EXOME

            else
                echo
                echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
                echo
                echo WARNING - Found Exomes but no Constitutional Exome
                echo
                echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
                echo

            fi
            # Clean-up temp exome BAM List file
            rm Exome_BAM_LIST
            # Back-up to the starting directory
            cd ..
        else
            echo No Exome files found
            # Clean-up temp exome BAM List file
            rm Exome_BAM_LIST
            # Back-up to the starting directory
            cd ..
        fi
    fi
done

echo
echo
echo All Processes Complete
echo
echo
