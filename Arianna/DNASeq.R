# ti colleghi 

mount -t cifs //santabarbara/egarattini /mnt -o "username=egarattlab,password=nextseq500"
#carichi l'unità

ssh mfratell@login.pico.cineca.it
cd /pico/scratch/userexternal/mfratell/
mkdir CORNELIA
exit

rsync -avz -progress //santabarbara/egarattini mfratell@login.pico.cineca.it:/pico/scratch/userexternal/mfratell/CORNELIA/

oppure scp 

scp -r //santabarbara/egarattini mfratell@login.pico.cineca.it:/pico/scratch/userexternal/mfratell/CORNELIA/

bcl2fastq

/pico/scratch/userexternal/mfratell/NGS_RUNS/bcl2fastq/bin/./bcl2fastq --sample-sheet /pico/scratch/userexternal/mfratell/NGS_RUNS/CELEGANS_08MAR16/samplesheet_celegans.csv --runfolder-dir /pico/scratch/userexternal/mfratell/NGS_RUNS/160308_NB501314_0005_AHW5YNBGXX --no-bgzf-compression --no-lane-splitting


es 170919_NB501314_0055_AHNCVWBGX3

cd 170919_NB501314_0055_AHNCVWBGX3/Data/Intensities/Basecalls
ls -lh

gunzip *.gz

/pico/scratch/userexternal/mfratell/CORNELIA/170919_NB501314_0055_AHNCVWBGX3/Data/Intensities/BaseCalls

qsub -q parallel -I -l select=1:ncpus=1,walltime=02:00:00 -A ELIX2_prj9

sposta fastq
mv ./170919_NB501314_0055_AHNCVWBGX3/Data/Intensities/BaseCalls/*.fastq ./fastq/

# PER FARE QC
qsub -q parallel -I -l select=1:ncpus=20:mem=120G,walltime=24:00:00 -A ELIX2_prj9

cd $CINECA_SCRATCH/CORNELIA/analisi

# load CINECA pre-installed modules
module load profile/advanced
module load fastqc

# define VARIABLES
ncpus=20
fastqDir=./fastq
qcDir=./QC #ricordati di crearla
mkdir $qcDir

fastqc -t $ncpus -o $qcDir ./fastq/MADRE_S2_R1_001.fastq
fastqc -t $ncpus -o $qcDir ./fastq/MADRE_S2_R2_001.fastq

fastqc -t $ncpus -o $qcDir ./fastq/PADRE_S3_R1_001.fastq
fastqc -t $ncpus -o $qcDir ./fastq/PADRE_S3_R2_001.fastq

fastqc -t $ncpus -o $qcDir ./fastq/SOFIA_S1_R1_001.fastq
fastqc -t $ncpus -o $qcDir ./fastq/SOFIA_S1_R2_001.fastq

multiqc ./QC -o ./mQC

#allinea con bowtie2
module load bowtie2
Fastq_R1=/pico/scratch/userexternal/mfratell/CORNELIA/analisi/fastq/SOFIA_S1_R1_001.fastq
Fastq_R2=/pico/scratch/userexternal/mfratell/CORNELIA/analisi/fastq/SOFIA_S1_R2_001.fastq
BOWTIE2_INDEX=/pico/scratch/userexternal/mfratell/GENOME/BOWTIE2/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index
bowtie2 -x $BOWTIE2_INDEX -1 $Fastq_R1 -2 $Fastq_R2 --local -p 20 -S /pico/scratch/userexternal/mfratell/CORNELIA/analisi/sam/sofia.sam   

#allinea con bwa
#scarichi 
module load bwa
Fastq_R1=/pico/scratch/userexternal/mfratell/CORNELIA/analisi/fastq/MADRE_S2_R1_001.fastq
Fastq_R2=/pico/scratch/userexternal/mfratell/CORNELIA/analisi/fastq/MADRE_S2_R2_001.fastq
hg38_fasta=/pico/scratch/userexternal/mfratell/GENOME/iGENOMES/Homo_sapiens/hg38/Sequence/BWAIndex/genome.fa

bwa mem -M -R '@RG\tID:foo\tSM:bar' -t 20 $hg38_fasta $Fastq_R1 $Fastq_R2 > madre_bwa.sam

module load autoload samtools
samtools view -bS -@ 20 madre_bwa.sam | samtools sort -@ 20 -o madre_bwa_sorted.bam
samtools view -bS -@ 20 padre_bwa.sam | samtools sort -@ 20 -o padre_bwa_sorted.bam
samtools view -bS -@ 20 sofia_bwa.sam | samtools sort -@ 20 -o sofia_bwa_sorted.bam

samtools view -bS -@ 20 sofia.sam | samtools sort -@ 20 -o sofia_sorted.bam
samtools index sofia_sorted.bam

samtools view -@ 20 -b sofia_sorted.bam "chr5:36544862-37398309" > nipbl.bam

samtools sort -@ 20 nipbl.bam -o nipbl_sorted.bam
samtools index nipbl_sorted.bam

scp mfratell@login.pico.cineca.it:/pico/scratch/userexternal/mfratell/CORNELIA/analisi/sam/nipbl_sorted* ./

PER VISUALIZZARLI 
IGV

https://software.broadinstitute.org/software/igv/download




