#!/bin/bash

path_dir=""
p_logs="$path_dir/logs"
p_rawdata="$path_dir/00_rawdata"
p_rawqc="$path_dir/00-1_rawdataqc"
p_cleanLog="${path_dir}/logs/01-0_cleandata"
input_path=""
sampleinfo="${path_dir}/sampleinfo.txt"
p_clean="$path_dir/01_0_trimmed_cleandata"
sampleinfo1="${path_dir}/sampleinfo.txt"
trimmed_clean_log=${path_dir}/logs/01_0_trimmed_cleandata
spkFactorDir=${path_dir}/03_spkFactor
geneCountsDir=${path_dir}/05_geneCounts
geneCountsLogDir=${path_dir}/logs/featureCounts

#
##step 0.1
#####change file name#####
echo -e "\n***************************\nRenaming files at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"

mkdir -p ${p_rawdata}

cat $sampleinfo| while read id;
do
arr=($id)
sample1=${arr[0]}
sample2=${arr[1]}
if [ ! -s ${p_rawdata}/${sample2}_R1.fastq.gz ]
then
    
    fq1=$(ls ${input_path}/*R1.fastq.gz|grep "$sample1")
    fq2=$(ls ${input_path}/*R2.fastq.gz|grep "$sample1")
    ln -s $fq1 ${p_rawdata}/${sample2}_R1.fastq.gz
    ln -s $fq2 ${p_rawdata}/${sample2}_R2.fastq.gz
fi

done

#step 0.2 rawdata fastqc
mkdir -p ${p_rawqc}

fastqc ${p_rawdata}/*.gz -t 10 -o ${p_rawqc}

wait
multiqc ${p_rawqc}/ -n trimmed_multiqc -o ${p_rawqc}




#step 1.1 trim_galore

echo -e "\n
***************************
Trim_galore begins at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)
***************************"
mkdir -p ${p_logs}
mkdir -p ${p_clean}
mkdir -p ${trimmed_clean_log}


nohup_number=0
for fq1 in `ls ${p_rawdata}/*R1.fastq.gz`
do
fq2=${fq1/R1.fastq.gz/R2.fastq.gz}
    if [ ! -s ${p_clean}/"$(basename ${fq1/.fastq.gz/_val_1.fq.gz})" ]
    then
        echo "Generating file: ${p_clean}/"$(basename ${fq1/.fastq.gz/_val_1.fq.gz})";"
     trim_galore -q 25 --phred33 --length 36 -e 0.1 --stringency 4 --paired -o ${p_clean} $fq1 $fq2 \
    > ${trimmed_clean_log}/"$(basename ${fq1/_R1.fq.gz/_trimmed.log})" 2>&1 &
    
    nohup_number=`echo $nohup_number+1 | bc`
    fi
    if [[ $nohup_number -eq 28 ]]
    then
        echo "waiting..."
        wait
        nohup_number=0
    fi
done

wait
##
#step1.2 trimmeddata fastqc
p_trimmedqc=${path_dir}/01_1_trimmeddataqc
mkdir -p ${p_trimmedqc}
fastqc ${p_clean}/*.gz -t 8 -o ${p_trimmedqc}

multiqc ${p_trimmedqc}/ -n trimmed_multiqc -o ${p_trimmedqc}


##step2.1
#####align###
echo -e "\n***************************\nalign begins at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"
##
p_align_exp=$path_dir/02_1_aligneddata_exp
p_align_spk=$path_dir/02_2_aligneddata_spk
mkdir ${p_align_exp}
mkdir ${p_align_spk}

cd $path_dir/02_1_aligneddata
#
index=/share/home/Blueberry/reference/index/star/hg19_star_index
mm10index=/share/home/Blueberry/reference/index/star/mm10_star_index
#

cat $sampleinfo | while read id;do
    arr=($id)
    i=${arr[1]}

    if [ ! -s ${p_align_exp}/${i}_hg19.rmdup.bam ]
    then
   cleanfq1=`ls $path_dir/01_0_trimmed_cleandata/${i}*R1*gz`;
   cleanfq2=`ls $path_dir/01_0_trimmed_cleandata/${i}*R2*gz`;
    mkdir -p ${path_dir}/logs/02-1_align/${i}
    STAR --runThreadN 18 --genomeDir $index \
        --readFilesIn $cleanfq1 $cleanfq2 \
        --readFilesCommand  zcat \
        --outSAMtype BAM SortedByCoordinate --twopassMode Basic --outFilterMismatchNmax 2 --outSJfilterReads Unique \
        --quantMode GeneCounts --outFileNamePrefix ${p_align_exp}/${i}_hg19.
    STAR --runThreadN 32 --genomeDir $mm10index \
        --readFilesIn $cleanfq1 $cleanfq2 \
        --readFilesCommand  zcat \
        --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 36 --twopassMode Basic --outFilterMismatchNmax 2 --outSJfilterReads Unique \
        --quantMode GeneCounts --outFileNamePrefix ${p_align_spk}/${i}_mm10.
    

####step2.2
####dedup##
p_filter_exp=$path_dir/02_3_filterdata_exp
p_filter_spk=$path_dir/02_4_filterdata_spk
mkdir ${p_filter_exp}
mkdir ${p_filter_spk}
mkdir $path_dir/logs/exp_filter
mkdir $path_dir/logs/spki_filter

echo "$cleanfq1"
   input=`ls ${p_align_exp}/${i}_hg19*bam`
   picard MarkDuplicates REMOVE_DUPLICATES=True INPUT=$input output=${p_filter_exp}/${i}_hg19.rmdup.bam METRICS_FILE=${p_filter_exp}/${i}_hg19.metrics 2>$path_dir/logs/exp_filter/${i}_hg19.dup.log;
   input=`ls ${p_align_spk}/${i}_mm10*bam`
   picard MarkDuplicates REMOVE_DUPLICATES=True INPUT=$input output=${p_filter_spk}/${i}_mm10.rmdup.bam METRICS_FILE=${p_filter_spk}/${i}_mm10.metrics 2>$path_dir/logs/spki_filter/${i}_mm10.dup.log;
   samtools index -@ 25 ${p_filter_exp}/${i}_hg19.rmdup.bam
   samtools index -@ 25 ${p_filter_spk}/${i}_mm10.rmdup.bam
   samtools flagstat ${p_filter_exp}/${i}_hg19.rmdup.bam > $path_dir/logs/exp_filter/${i}_hg19.rmdup.stat
   samtools flagstat ${p_filter_spk}/${i}_mm10.rmdup.bam > $path_dir/logs/spki_filter/${i}_mm10.rmdup.stat
   fi
done
#####step3:alculate spike-in normalization factor
p_filter_exp=$path_dir/02_3_filterdata_exp
p_filter_spk=$path_dir/02_4_filterdata_spk
p_align_exp=$path_dir/02_1_aligneddata_exp
p_align_spk=$path_dir/02_2_aligneddata_spk
p_trimmedqc=${path_dir}/01_1_trimmeddataqc
spkFactorDir=${path_diir}/03_spkFactor
if [[ ! -d $spkFactorDir ]];then
               mkdir -p $spkFactorDir
               echo -e "sample_named id;\tall_reads\thg19_reads\thg19_mapping_ratio\thg19_qc_reads\thg19_qc_ratio\tmm10_qc_reads\tmm10_qc_ratio\tmm10_qc_ratio_all_qc\tscale_factor" > ${spkFactorDir}/spkScaleFactor.txt
                basename -a $(ls ${path_dir}/02_1_aligneddata_exp/*_hg19.Log.final.out) | while read file;do
                sampleName=${file/_hg19.Log.final.out/}

                allReads=`grep "Number of input reads" ${path_dir}/02_1_aligneddata_exp/${file} | awk '{print $NF*2}'`
                hg19Reads=`grep -e "Number.* mapped" -e "mapped.*number" ${path_dir}/02_1_aligneddata_exp/${file} | cut -f 2 | awk '{sum+=$0};END{print sum*2}'`
                hg19MapRatio=`grep " mapped.*%" ${path_dir}/02_1_aligneddata_exp/${file} | cut -f 2 | awk '{sum+=$0};END{print sum}'`

        
                hg19QcReads=`cat ${path_dir}/logs/exp_filter/${sampleName}_hg19.rmdup.stat | grep "total (QC-passed reads" | cut -d " " -f 1`
                hg19QcRatio=`printf "%.2f\n" $(echo "100*${hg19QcReads}/${allReads}" | bc -l)`
                       
                mm10QcReads=`cat ${path_dir}/logs/spki_filter/${sampleName}_mm10.rmdup.stat | grep "total (QC-passed reads" | cut -d " " -f 1`
                mm10QcRatio=`printf "%.2f\n" $(echo "100*${mm10QcReads}/${allReads}" | bc -l)`

                scaleFactor=`echo "1000000/${hg19QcReads}" | bc -l`

                allQcReads=`echo "${hg19QcReads}+${mm10QcReads}" | bc -l`
                mm10QcRatio4AllQc=`printf "%.2f\n" $(echo "100*${mm10QcReads}/${allQcReads}" | bc -l)`

                echo -e ${sampleName}"\t"${allReads}"\t"${hg19Reads}"\t"${hg19MapRatio}"\t"${hg19QcReads}"\t"${hg19QcRatio}"\t"${mm10QcReads}"\t"${mm10QcRatio}"\t"${mm10QcRatio4AllQc}"\t"${scaleFactor} >> ${spkFactorDir}/spkScaleFactor.txt
               done
       fi


## track plus minus
mkdir -p $path_dir/06-1_trackStrand
mkdir -p $path_dir/logs/06-1_trackStrand
scalefactor_file=${path_dir}/03_spkFactor/spkScaleFactor.txt
cat ${scalefactor_file} | sed '1d' | while read id;
do
arr=($id)
sample=${arr[0]}
bam_file=${p_filter_spk}/${sample}_mm10.rmdup.bam
scalefactor=${arr[9]}
         bamCoverage \
        --bam ${bam_file} \
        --skipNonCoveredRegions \
        --blackListFileName /share/home/Blueberry/reference/annotation/encode/blacklist/mm10-blacklist.v2.bed \
        --outFileName $path_dir/06-1_trackStrand/${sample}_plus.bw \
        --binSize 1 \
	--scaleFactor ${scalefactor} \
        --numberOfProcessors 20 \
        --filterRNAstrand forward \
        --normalizeUsing None
         bamCoverage \
        --bam ${bam_file} \
        --skipNonCoveredRegions \
        --blackListFileName /share/home/Blueberry/reference/annotation/encode/blacklist/mm10-blacklist.v2.bed \
        --outFileName $path_dir/06-1_trackStrand/${sample}_minus.bw \
        --binSize 1 \
        --scaleFactor ${scalefactor} \
        --numberOfProcessors 20 \
        --filterRNAstrand reverse \
        --normalizeUsing None
           bamCoverage \
        --bam ${bam_file} \
        --skipNonCoveredRegions \
        --blackListFileName /share/home/Blueberry/reference/annotation/encode/blacklist/mm10-blacklist.v2.bed \
        --outFileName $path_dir/06-1_trackStrand/${sample}_rmdup-Q7_dr_rev_minus.bw \
        --binSize 1 \
        --scaleFactor -${scalefactor} \
        --numberOfProcessors 20 \
        --filterRNAstrand reverse \
        --normalizeUsing None


done



#step5.1:gene quantify by run-featurecounts.R
mkdir -p $path_dir/05_1_runfeaturecounts
for rmdup_bam in $path_dir/02_4_filterdata_spk/*rmdup.bam
do
    rmdup_bam_base=$(basename $rmdup_bam)
    Rscript /chenfeilab/Pomelo/scripts/RNAseq.run-featurecounts.R \
        -b ${rmdup_bam} \
        -g /share/home/Blueberry/reference/annotation/gencode/mm10.gencode_chr.gtf \
        -o $path_dir/05_1_runfeaturecounts/exon_counts/${rmdup_bam_base/_mm10.rmdup.bam/}
done





















































