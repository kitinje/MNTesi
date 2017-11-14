# Tools used
# 1. FastQC (pre-alignmed)
# 2. AfterQC (https://github.com/OpenGene/AfterQC) (pre-alignmed)
# 3. STAR (Alignment)
# 4. 
# 4. Bamtools stats (post-alignment)
# 5. RSeqQC (post-alignment)
# 6. Picard Tools

module load bowtie2
module load samtools
module load r
module load autoload python/2.7.9
module load bwa

fastq_DIR=/pico/scratch/userexternal/mfratell/JPT/TJ/fastq/fastq
AfterQC_DIR=/pico/scratch/userexternal/mfratell/SW/AfterQC
FastqScreen_DIR=/pico/scratch/userexternal/mfratell/SW/fastq_screen_v0.11.2
Bamtools_DIR=/pico/scratch/userexternal/mfratell/SW/bamtools/bin
RSeQC_DIR=/pico/home/userexternal/mfratell/.local/bin/
picard_DIR=/cineca/prod/applications/picard/2.3.0/binary/bin/
bwa_DIR=/cineca/prod/applications/bwa/0.7.15/intel--cs-xe-2015--binary/bin/
deeptools_DIR=/pico/home/userexternal/mfratell/.local/bin/
qualimap_DIR=/pico/scratch/userexternal/mfratell/SW/qualimap_v2.2.1/

refFile=/pico/scratch/userexternal/mfratell/GENOME/hg38_RefSeq.bed
srefFile=/pico/scratch/userexternal/mfratell/GENOME/hg38_RefSeq_top1000.bed
rRNAref=/pico/scratch/userexternal/mfratell/GENOME/hg38_rRNA.bed
refFlat=/pico/scratch/userexternal/mfratell/GENOME/hg38_refFlat.txt
GTF=/pico/scratch/userexternal/mfratell/GENOME/GENCODE/v27/gencode.v27.annotation.gtf

cd $fastq_DIR
python $AfterQC_DIR/after.py --qc_only --read1_flag fastq -1 

ls | grep ".fastq" > fastqlist.txt
for i in $(cat fastqlist.txt); do
            echo item: $i
            python $AfterQC_DIR/after.py --qc_only --read1_flag fastq -1 $i 
done


ls *.bam > bamlist.txt
$Bamtools_DIR/bamtools stats -list bamlist.txt

for i in $(cat bamlist.txt); do
            #echo item: $i
            newName="${i/_Aligned.sortedByCoord.out.bam/_resampled_01.bam}"
            #echo $newName
            samtools view -@ 20 -b -s 0.1 $i >$newName
done



# DO THIS FOR PAIRED-END
for i in $(cat bamlist.txt); do
            #echo item: $i
            pdf="${i/_Aligned.sortedByCoord.out.bam/_report.html}"
            temp="${i/_Aligned.sortedByCoord.out.bam/}"
            dirOUT=./qmapQC/$temp
            echo $dirOUT
            $qualimap_DIR/qualimap bamqc -bam $i --feature-file $GTF -nt 20 -outdir $dirOUT  -outfile $pdf -outformat HTML -gd HUMAN
done

for i in $(cat bamlist.txt); do
            #echo item: $i
            pdf="${i/_Aligned.sortedByCoord.out.bam/_report.html}"
            temp="${i/_Aligned.sortedByCoord.out.bam/}"
            dirOUT=./qmapQC2/$temp
            echo $dirOUT
            $qualimap_DIR/qualimap RNASeq -bam $i -gtf $GTF -outdir $dirOUT  -outfile $pdf -outformat HTML --java-mem-size=100G
done

#samtools view -@ 20 -b -s 0.1 11-VCaP-PAR_S21_Aligned.sortedByCoord.out.bam > 11-VCaP-PAR_S21_resampled_01.bam

mkdir resampled
mv *resampled* ./resampled

ls *resampled*.bam > rebamlist.txt

for i in $(cat rebamlist.txt); do
            samtools index $i
            echo $i
done

pip install RSeQC --user



$RSeQC_DIR/infer_experiment.py -r $refFile -i 10-VCaP-1916-RNT5_S13.fastq

for i in $(cat rebamlist.txt); do
            echo item: $i
            prefix=$i
            $RSeQC_DIR/geneBody_coverage.py -i $i -r $srefFile -o $prefix
done

for i in $(cat rebamlist.txt); do
            echo item: $i
            newName="${i/_resampled_01.bam/_read_distribution.txt}"
            #echo new: $newName
            $RSeQC_DIR/read_distribution.py -i $i -r $refFile > $newName
done

for i in $(cat rebamlist.txt); do
            echo item: $i
            newName="${i/_resampled_01.bam/_read_duplication.txt}"
            #echo new: $newName
            $RSeQC_DIR/read_duplication.py -i $i -o $newName
done

for i in $(cat rebamlist.txt); do
            echo item: $i
            ribName="${i/_resampled_01.bam/_rRNA.bam}"
            newName="${i/_resampled_01.bam/_ribosomalrRNA.txt}"
            $RSeQC_DIR/split_bam.py -i $i -r $rRNAref -o $ribName > $newName
done

# for i in $(cat rebamlist.txt); do
#             echo item: $i
#             newName="${i/_resampled_01.bam/_bamstat.txt}"
#             $RSeQC_DIR/bam_stat.py -i $i > $newName
# done

for i in $(cat rebamlist.txt); do
            echo item: $i
            newName="${i/_resampled_01.bam/_junctionsat.txt}"
            $RSeQC_DIR/junction_saturation.py -i $i -r $refFile -o $newName
done

# ONLY FOR PAIRED-END
for i in $(cat bamlist.txt); do
            echo item: $i
            rnaName="${i/_Aligned.sortedByCoord.out.bam/_RNA_Metrics.txt}"
            pdfName="${i/_Aligned.sortedByCoord.out.bam/_insertsize.pdf}"
            java -XX:ParallelGCThreads=20 -jar $picard_DIR/picard.jar CollectInsertSizeMetrics I=$i  O=$txtName H=$pdfName M=0.05
done
