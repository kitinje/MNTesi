#!/bin/bash
. $MODULESHOME/init/bash

# load CINECA pre-installed modules
module load profile/advanced
module load autoload star

# define VARIABLES
ncpus=20
trimDir=./trimmedFq
UnmappedDir=./unmapped
LogDir=./logs
bamDir=./BAM
bam2Dir=./BAM2
tabDir=./TAB
tab2Dir=./TAB2
wigDir=./wig
refgenomeFA=/pico/scratch/userexternal/mfratell/GENOME/GENCODE/v26/GRCh38.p10.genome.fa
GTF=/pico/scratch/userexternal/mfratell/GENOME/GENCODE/v26/gencode.v26.annotation.gtf
GENOME_DIR=/pico/scratch/userexternal/mfratell/GENOME/GENCODE/v26/star_GRCh38.p10; 

printf "\033c" # clean screen
FW_count=$(ls $trimDir/*.fastq | grep "_R1" | wc -l)
RV_count=$(ls $trimDir/*.fastq | grep "_R2" | wc -l)

if [ $FW_count != $RV_count ] 
      then
         printf "\n"
         echo "ERROR: Number of FW files does not match number of RV files"
         printf "\n"
         exit 1
fi

echo "FASTQ PAIRS IDENTIFIED:"
for i in `seq 1 $FW_count`;
   do
      cmd='FW_input=$(ls -lh $trimDir/*.fastq | grep "_R1" |  awk '"'{print \$9}'"' | awk '"'NR==$i'"')'
      eval $cmd
      cmd='RV_input=$(ls -lh $trimDir/*.fastq | grep "_R2" |  awk '"'{print \$9}'"' | awk '"'NR==$i'"')'
      eval $cmd
      printf "\n"
      printf "$FW_input -- $RV_input"
   done 
printf "\n"   
printf "\n" 
echo "IDENTIFIED $FW_count PAIRS OF FASTQ FILES, PLEASE CHECK IF PAIRS ARE CORRECTLY MATCHED!!"
printf "\n"
read -p "Press enter to confirm or CTRL+C to abort.."

for i in `seq 1 $FW_count`;
   do
      cmd='FW_input=$(ls -lh $trimDir/*.fastq | grep "_R1" | awk '"'{print \$9}'"' | awk '"'NR==$i'"')'
      eval $cmd
      cmd='RV_input=$(ls -lh $trimDir/*.fastq | grep "_R2" | awk '"'{print \$9}'"' | awk '"'NR==$i'"')'
      eval $cmd
      Prefix=$(echo $FW_input | sed 's/.\/trimmedFq\///g' | sed 's/_trimmed.fastq//g' | sed 's/_R1//g')
      Prefix=$(echo $Prefix"_")
      printf "\n"
      echo "Performing sequence aligment of PAIR $i of $FW_count:"
      STAR --runThreadN $ncpus --genomeDir $GENOME_DIR  --outFileNamePrefix $Prefix --readFilesIn $FW_input $RV_input --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic --outWigType bedGraph  --outWigStrand Stranded --outWigNorm RPM
   done 
 mkdir $bamDir
 mkdir $bam2Dir
 mkdir $tabDir
 mkdir $tab2Dir
 mkdir $wigDir
 mkdir $LogDir   
 mkdir $UnmappedDir
 mv *.out.bg $wigDir
 mv *Aligned.sortedByCoord.out.bam $bamDir
 mv *Aligned.toTranscriptome.out.bam $bam2Dir
 mv *Unmapped* $UnmappedDir
 mv *Log*.out $LogDir
 mv *ReadsPerGene.out.tab $tabDir
 mv *SJ.out.tab $tab2Dir
 rm -r $(ls | grep "STARgenome")
 rm -r $(ls | grep "STARpass1")
 
printf "\n"
echo "STAR alignment completed"
printf "\n"
# End of bash script
# Marco Bolis 2017
