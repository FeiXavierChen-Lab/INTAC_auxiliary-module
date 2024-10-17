#!/bin/bash



###call peak
echo -e "\n
***************************
Peak calling begins at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)
***************************"
mkdir -p ${peak_calling_dir}

nohup_number=0
for sample in `cat ${scalefactor_file} |cut -f 1`
do
    bam_file=${exp_bam_rmdup}/${sample}*.Shift.sort.bam
    if [ ! -s ${peak_calling_dir}/${sample}_peaks.narrowPeak ]
    then
    macs2 callpeak -f BAMPE --keep-dup all --nolambda --bdg -g hs -n ${sample} -p 1e-5 \
    -t $bam_file --outdir ${peak_calling_dir} 2>${peak_calling_dir}/${sample}_macs2.log &
    nohup_number=`echo $nohup_number+1 | bc`

    if [[ $nohup_number -eq 60 ]]
    then
        wait
        echo "waiting..."
        nohup_number=0
    fi
fi
done
wait


for sample in `cat ${scalefactor_file} |cut -f 1`
do
    if [ ! -s ${peak_calling_dir}/${sample}_peaks.final.narrowPeak ]
    then
    bedtools intersect -a ${peak_calling_dir}/${sample}_peaks.narrowPeak -b $blacklist -f 0.25 -v >${peak_calling_dir}/${sample}_peaks.final.narrowPeak &
    fi
done
