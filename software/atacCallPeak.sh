#!/usr/bin/env bash


# =============================================================
# set some global variables, such as genome directories and Bowtie2 indexes
# this variables need to be modified by the users
BOWTIE2_INDEXES="/data2/Epi/likai/genome/human/UCSC/hg19/bowtie2_index/hg19"
GENOME_SIZE="/data2/Epi/likai/genome/human/UCSC/hg19/hg19.chrom.sizes"
THREAD_NUM=5
GENOME_NAME="hg19"
ATAC_SCRIPTS_PATH="/data2/Epi/likai/ATAC_annotation_Project/ATAC_raw_tech/scripts"
# =============================================================

PEAK_DIR=$(echo $1 | sed -r 's/\/+/\//g' | sed -r 's/\/$//g')
if [[ "$#" != 1 ]]
then
    echo
    echo -e "\033[0;31mUsage: bash atacCallPeak.sh [PEAK_DIR]\033[0m"
    echo
    exit 1
fi


function callPeak {
    BED_FILE=$1
    OUT_PREFIX=$(echo ${BED_FILE} | sed 's/\.bed//g')
    PEAK_LOG="${OUT_PREFIX}.peak.log"
    echo -e "  \033[0;32mStart calling peak: $1\033[0m"
    macs2 callpeak -t ${BED_FILE} -g hs --nomodel --extsize 75 -n ${OUT_PREFIX} -f BED --keep-dup all --call-summits > ${PEAK_LOG} 2>&1
}


function multiCallPeak {
    mkfifo tmp
    exec 5<>tmp
    rm -rf tmp

    for((i=1;i<=$THREAD_NUM;i++))
    do
        echo -ne "\n" 1>&5;
    done

    for BED_FILE in $(find . -name '*.bed')
    do
	    read -u5
	    {
	        callPeak ${BED_FILE}
	        echo -ne "\n" 1>&5
	    }&
    done
    wait
    exec 5>&-
    echo -e "\033[0;32mFinished calling peaks\033[0m"
}

# call peaks with macs2
source activate /data/Epi/software/miniconda3/envs/atac
cd ${PEAK_DIR}
multiCallPeak

LOG_FILES=$(find . -name "*.log")
if [[ ${#LOG_FILES} -gt 0 ]]
then
    for file in ${LOG_FILES}
    do
        mv ${file} "../log"
    done
fi
