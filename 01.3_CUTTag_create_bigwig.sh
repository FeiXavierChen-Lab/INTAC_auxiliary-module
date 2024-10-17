#!/usr/bin/bash

### track
### Making RPKM-normalized bigWig files with full-length reads without spike-in ###
echo -e "\n
***************************
Track begins at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)
***************************"
mkdir -p ${bw_fulllength_dir}
mkdir -p ${bw_fulllength_dir_EcoliSpike}

ls ${exp_bam_rmdup}/*.Shift.sort.bam |while read id;
do
    sample=$(basename ${id%%.Shift.sort.bam})
    echo "Bigwig for ${sample}"
    if [ ! -s "${bw_fulllength_dir}/${sample}_rpkm.bw" ]
    then
        bamCoverage -p 60 \
        --skipNonCoveredRegions \
        --normalizeUsing RPKM \
        --binSize 1 \
        --blackListFileName $blacklist \
        -b $id \
        -o ${bw_fulllength_dir}/${sample}_rpkm.bw
    fi
done

ls ${exp_bam_rmdup}/*.Shift.sort.bam |while read id;
do
    sample=$(basename ${id%%_${exp_info}.Shift.sort.bam})
    SCALEFACTOR=`grep $sample ${scalefactor_file} | awk '{print $NF}'`
    echo "$sample $SCALEFACTOR"
    if [ ! -s "${bw_fulllength_dir_EcoliSpike}/${sample}_THP1Spike.bw" ]
    then
        bamCoverage -p 60 \
        --skipNonCoveredRegions \
        --binSize 1 \
        --scaleFactor $SCALEFACTOR \
        --blackListFileName $blacklist \
	--normalizeUsing None \
        -b $id \
        -o ${bw_fulllength_dir_EcoliSpike}/${sample}_THP1Spike.bw
    fi
done

