#!/bin/bash
. $MODULESHOME/init/bash

# load CINECA pre-installed modules
module load profile/advanced
module load fastqc

# define VARIABLES
ncpus=20
fastqDir=./fastq
qcDir=./QC
mkdir $qcDir

printf "\033c" # clean screen
fastq_count=$(ls $fastqDir/*.fastq | wc -l)
echo "IDENTIFIED $fastq_count FASTQ FILES:"
ls -lh $fastqDir | grep "fastq" | awk '{print $9}'
printf "\n"
read -p "Press enter to confirm or CTRL+C to abort.."
printf "\n"
for i in `seq 1 $fastq_count`;
   do
      printf "\n"
      echo "Performing quality trimming of file $i of $fastq_count:"
      cmd="ls $fastqDir | grep \"fastq\" | awk '"'NR=='$i''"' "
      current=$(eval $cmd)
      fastqc -t $ncpus -o $qcDir $fastqDir/$current
   done 
printf "\n"
echo "Initial QC completed."
printf "\n"
# End of bash script
# Marco Bolis 2017
