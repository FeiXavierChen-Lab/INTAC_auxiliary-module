#!/bin/bash



## align
echo -e "\n
***************************
Align begins at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)
***************************"
mkdir -p ${align_exp_dir}
mkdir -p ${alignexp_log_dir}
mkdir -p ${align_myc_dir}
mkdir -p ${alignmyc_log_dir}
mkdir -p ${exp_bam_rmdup}
mkdir -p ${rmdup_exp_log}
mkdir -p ${spike_bam_rmdup}
mkdir -p ${rmdup_spike_log}

for fq1 in `ls ${trimmedFastq_dir}/*R1.fq.gz`
do
    fq2=${fq1/R1.fq.gz/R2.fq.gz}
    sample="$(basename ${fq1/_trimmed_R1.fq.gz/})"
    echo "Generating file: ${align_exp_dir}/${sample}_${exp_info}.sort.rmdup.bam;"

    # -t:time -q:fastq(default) -L:seedlength no-mixed&no-discordant:for paired reads
    if [ ! -s ${align_exp_dir}/${sample}_${exp_info}.sort.bam.stat ]
    then
    bowtie2 -p 60 -t -q -N 1 -L 25 -X 700 --no-mixed --no-discordant -x ${GENOME_EXP} -1 $fq1 -2 $fq2 \
        2>${alignexp_log_dir}/${sample}_${exp_info}.log | samtools sort -@ 60 -O BAM -o ${align_exp_dir}/${sample}_${exp_info}.sort.bam;
    samtools flagstat ${align_exp_dir}/${sample}_${exp_info}.sort.bam  > ${align_exp_dir}/${sample}_${exp_info}.sort.bam.stat
    fi
    # align to myc
    if [ ! -s ${align_myc_dir}/${sample}_myc.sort.bam ]
    then
    bowtie2 -p 40 -t -q -N 1 -L 25 -X 700 --no-mixed --no-discordant -x ${GENOME_MYC} -1 $fq1 -2 $fq2 \
        2>${alignmyc_log_dir}/${sample}_myc.log | samtools sort -@ 40 -O BAM -o ${align_myc_dir}/${sample}_myc.sort.bam;
    fi
done


# Align to spikein  genome
for fq1 in `ls ${trimmedFastq_dir}/*R1.fq.gz`
do
    fq2=${fq1/R1.fq.gz/R2.fq.gz}
    sample="$(basename ${fq1/_trimmed_R1.fq.gz/})"
    if [ ! -s ${align_exp_dir}/${sample}_THP1.sort.filtered.bam ]
    then
        echo "Generating file: ${align_exp_dir}/${sample}_THP1.sort.bam;"
        # no-overlap: not concordant when mates overlap at all
        # --no-dovetail
        # to avoid possible cross-mapping of the experimental genome to that of the carry-over E. coli DNA that is used for calibration.
        bowtie2 -p 60 -t -q -N 1 -L 25 -X 700 --no-mixed --no-discordant -x ${GENOME_human} -1 $fq1 -2 $fq2 \
            2>${alignexp_log_dir}/${sample}_THP1.log |samtools view -bS -f 2 -q 10 | samtools sort -@ 35 -O BAM -o ${align_exp_dir}/${sample}_THP1.sort.filtered.bam;
        echo "logs in ${alignexp_log_dir}/${sample}_THP1.log"
    fi
done

###   rmdup
for fq1 in `ls ${trimmedFastq_dir}/*R1.fq.gz`
do
    sample="$(basename ${fq1/_trimmed_R1.fq.gz/})"
    if [ ! -s ${exp_bam_rmdup}/${sample}_THP1.sort.filtered.rmdup.bam ]
    then
        samtools flagstat ${align_exp_dir}/${sample}_THP1.sort.filtered.bam  > ${align_exp_dir}/${sample}_THP1.sort.bam.stat

        picard MarkDuplicates REMOVE_DUPLICATES=True INPUT=${align_exp_dir}/${sample}_THP1.sort.filtered.bam \
                    output=${exp_bam_rmdup}/${sample}_THP1.sort.filtered.rmdup.bam METRICS_FILE=${exp_bam_rmdup}/${sample}_THP1.metrics \
            2>${rmdup_exp_log}/${sample}_THP1.dup.log;
        samtools flagstat ${exp_bam_rmdup}/${sample}_THP1.sort.filtered.rmdup.bam >${rmdup_exp_log}/${sample}_THP1.sort.rmdup.bam.stat
    fi
done


nohup_number=0
for fq1 in `ls ${trimmedFastq_dir}/*R1.fq.gz`
do
    fq2=${fq1/R1.fq.gz/R2.fq.gz}
    sample="$(basename ${fq1/_trimmed_R1.fq.gz/})"
    if [ ! -s ${exp_bam_rmdup}/${sample}_${exp_info}.sort.rmdup.bam ]
    then
        picard MarkDuplicates REMOVE_DUPLICATES=True INPUT=${align_exp_dir}/${sample}_${exp_info}.sort.bam \
                output=${exp_bam_rmdup}/${sample}_${exp_info}.sort.rmdup.bam METRICS_FILE=${exp_bam_rmdup}/${sample}_${exp_info}.metrics \
            2>${rmdup_exp_log}/${sample}_${exp_info}.dup.log
    fi
done
wait

