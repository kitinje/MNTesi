#!/bin/bash
. $MODULESHOME/init/bash

# load CINECA pre-installed modules
module load profile/advanced
module load autoload cutadapt

# define VARIABLES
ncpus=20
fastqDir=./fastq
trimDir=./trimmedFq
ADAPTER_R1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
ADAPTER_R2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
MIN_LENGTH=18
ERROR_RATE=0.2
OVERLAP=5

#cutadapt -a ${ADAPTER_R1} -A ${ADAPTER_R2} --minimum-length $MIN_LENGTH --error-rate $ERROR_RATE --overlap=$OVERLAP -o $trimmed_R1_fq --paired-output $trimmed_R2_fq $R1_fq $R2_fq

mkdir $trimDir

printf "\033c" # clean screen
FW_count=$(ls $fastqDir/*.fastq | grep ".fastq" | grep "_R1" | wc -l)
RV_count=$(ls $fastqDir/*.fastq | grep ".fastq" | grep "_R2" | wc -l)

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
      cmd='FW_input=$(ls -lh $fastqDir/ |  grep ".fastq" | grep "_R1" | awk '"'{print \$9}'"' | awk '"'NR==$i'"')'
      eval $cmd
      cmd='RV_input=$(ls -lh $fastqDir/ |  grep ".fastq" | grep "_R2" | awk '"'{print \$9}'"' | awk '"'NR==$i'"')'
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
      cmd='FW_input=$(ls -lh $fastqDir/ |  grep ".fastq" | grep "_R1" | awk '"'{print \$9}'"' | awk '"'NR==$i'"')'
      eval $cmd
      cmd='RV_input=$(ls -lh $fastqDir/ |  grep ".fastq" | grep "_R2" | awk '"'{print \$9}'"' | awk '"'NR==$i'"')'
      eval $cmd
      #$(echo $FW_input | sed 's/.\/trimmedFq\///g' | sed 's/_trimmed_R1.fastq//g')
      FW_output=$(echo $FW_input | sed 's/.fastq//g' )
      FW_output=$(echo $FW_output"_trimmed.fastq")
      RV_output=$(echo $RV_input | sed 's/.fastq//g' )
      RV_output=$(echo $RV_output"_trimmed.fastq")
      printf "\n"
      echo "Performing adapter trimming in PAIR $i of $FW_count:"
      cmd="cutadapt -a ${ADAPTER_R1} -A ${ADAPTER_R2} --minimum-length $MIN_LENGTH --error-rate $ERROR_RATE --overlap=$OVERLAP -o $trimDir/$FW_output --paired-output $trimDir/$RV_output $fastqDir/$FW_input $fastqDir/$RV_input"
      eval $cmd
   done 
printf "\n"
echo "Adapter-trimming completed."
printf "\n"
# End of bash script
# Marco Bolis 2017
