#!/usr/bin/bash

for fq1 in `ls ${trimmedFastq_dir}/*R1.fq.gz`
do
    fq2=${fq1/R1.fq.gz/R2.fq.gz}
    sample="$(basename ${fq1/_trimmed_R1.fq.gz/})"
    if [ ! -s ${exp_bam_rmdup}/${sample}_${exp_info}.Shift.sort.bam ]
    then
        # -f2:PROPER_PAIR
        samtools flagstat ${exp_bam_rmdup}/${sample}_${exp_info}.sort.rmdup.bam >${rmdup_exp_log}/${sample}_${exp_info}.sort.rmdup.bam.stat
        samtools view -bS -f 2 -q 10 ${exp_bam_rmdup}/${sample}_${exp_info}.sort.rmdup.bam \
        | samtools sort -@ 60 -o ${exp_bam_rmdup}/${sample}_${exp_info}.sort.rmdup.q10.bam
        samtools flagstat ${exp_bam_rmdup}/${sample}_${exp_info}.sort.rmdup.q10.bam > ${rmdup_exp_log}/${sample}_${exp_info}.sort.rmdup.q10.stat
        samtools index -@ 60 ${exp_bam_rmdup}/${sample}_${exp_info}.sort.rmdup.q10.bam
        alignmentSieve --numberOfProcessors 60 --ATACshift \
            -b ${exp_bam_rmdup}/${sample}_${exp_info}.sort.rmdup.q10.bam \
            -o ${exp_bam_rmdup}/${sample}_${exp_info}.Shift.bam
        samtools sort -@ 60 ${exp_bam_rmdup}/${sample}_${exp_info}.Shift.bam -o ${exp_bam_rmdup}/${sample}_${exp_info}.Shift.sort.bam
        samtools index -@ 60 ${exp_bam_rmdup}/${sample}_${exp_info}.Shift.sort.bam
    fi

    ## size distribution
    if [ ! -s ${rmdup_exp_log}/${sample}_fragmentLen.txt ]
    then
        gatk=/share/home/Coconut/Software/gatk-4.2.0.0/gatk
        $gatk CollectInsertSizeMetrics -H ${rmdup_exp_log}/${sample}_InsertSize.pdf \
            -I ${exp_bam_rmdup}/${sample}_${exp_info}.sort.rmdup.q10.bam -O ${rmdup_exp_log}/${sample}_InsertSize.txt

        # Extract the 9th column from the alignment sam file which is the fragment length
        samtools view -F 0x04 ${exp_bam_rmdup}/${sample}_${exp_info}.Shift.bam | \
        awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | \
        awk -v OFS="\t" '{print $2, $1/2}' >${rmdup_exp_log}/${sample}_fragmentLen.txt
    fi
done
