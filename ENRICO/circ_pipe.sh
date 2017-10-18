#############################################################
##### TOPHAT & STAR ALIGNERS PIPELINE IN CIRCEXPLORER2 #####
#############################################################

# Reference Genome : GRCh38.p7
# Annotations : Gencode v25
# Aligner : Tophat-Fusion (paired-end) // Note: Tophat is the only pe output supported by CircExplorer
# CircRNA algorithm : circExplorer2

# Server
ssh -t mfratell@login.pico.cineca.it

# Job
qsub -q bigmem -I -l select=1:ncpus=20:mem=510G,walltime=24:00:00 -A ELIX2_prj9

# Cartella principale
cd /pico/scratch/userexternal/mfratell/ENRICO/circRNA

# PREFIX Ã¨ il nome della linea seguito da _ATRA o _CONTROL
PREFIX=CAL851_ATRA

# PREFIX List
CAL851_ATRA      MB157_ATRA        MDAMB361_ATRA
CAL851_CONTROL   MB157_CONTROL     MDAMB361_CONTROL
HCC1419_ATRA     MDAMB157_ATRA     MDAMB436_ATRA
HCC1419_CONTROL  MDAMB157_CONTROL  MDAMB436_CONTROL
HCC1599_ATRA     MDAMB231_ATRA     SKBR3_ATRA
HCC1599_CONTROL  MDAMB231_CONTROL  SKBR3_CONTROL

# Variabili
TRIM_DIR=/pico/scratch/userexternal/mfratell/ENRICO/circRNA/trimmedfq_new/merged_trimmedfq
BASE_DIR=/pico/scratch/userexternal/mfratell/ENRICO/circRNA
STAR_DIR=/pico/scratch/userexternal/mfratell/ENRICO/circRNA/star_new
TOPHAT_DIR=/pico/scratch/userexternal/mfratell/ENRICO/circRNA/tophat
CIRCBASE_DIR=/pico/scratch/userexternal/mfratell/ENRICO/circRNA/circ_out
GENOME_DIR=/pico/scratch/userexternal/mfratell/GENOME/GENCODE/v25/star_GRCh38.p5 
CIRC_PATH=/pico/home/userexternal/mfratell/.local/bin
LOG_DIR=/pico/scratch/userexternal/mfratell/circRNA/log
OUT_DIR=$TOPHAT_DIR/$PREFIX
STAR_PATH=/pico/scratch/userexternal/mfratell/SW/STAR-2.5.1b/bin/Linux_x86_64_static
ADAPTER_R1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
ADAPTER_R2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
refgenomeFA=/pico/scratch/userexternal/mfratell/GENOME/GENCODE/v25/GRCh38.p7.genome.fa
refGTF=/pico/scratch/userexternal/mfratell/GENOME/GENCODE/v25/gencode.v25.annotation.gtf
refFlat=/pico/scratch/userexternal/mfratell/GENOME/GENCODE/v25/gencode.v25.refFlat.txt
ncpus=20
GRCh38_bowtie1_index=/pico/scratch/userexternal/mfratell/GENOME/bowtie1_gencodeV25/GRCh38.p7

# 0. Modules
module load profile/advanced
module load autoload cutadapt
module load autoload samtools/0.1.19
module load fastqc
module load python/2.7.9
module load bowtie

# 1. Tophat-Fusion Pipe
trimmed_R1_fq=$TRIM_DIR/$PREFIX"_R1"  
trimmed_R2_fq=$TRIM_DIR/$PREFIX"_R2"
mkdir $TOPHAT_DIR/$PREFIX

/pico/scratch/userexternal/mfratell/ENRICO/SW/tophat_2.1.0/tophat2 -o $OUT_DIR -p $ncpus --fusion-search --keep-fasta-order --bowtie1 --no-coverage-search ${GRCh38_bowtie1_index} $trimmed_R1_fq,$trimmed_R2_fq

# 1.1 CIRCExplorer on Tophat
mkdir $CIRCBASE_DIR/$PREFIX
mkdir $CIRCBASE_DIR/$PREFIX/tophat

TopHatFus=$TOPHAT_DIR/$PREFIX/"accepted_hits.bam"
$CIRC_PATH/CIRCexplorer2 parse -o $CIRCBASE_DIR/$PREFIX/tophat -t TopHat-Fusion $TopHatFus 
$CIRC_PATH/CIRCexplorer2 annotate -r $refFlat -g $refgenomeFA $CIRCBASE_DIR/$PREFIX/tophat

# 2. STAR Pipe
trimmed_R1_fq=$TRIM_DIR/$PREFIX"_R1"
trimmed_R2_fq=$TRIM_DIR/$PREFIX"_R2"
mkdir $BASE_DIR/$PREFIX
star_prefix=$STAR_DIR/$PREFIX"_"
$STAR_PATH/STAR --runThreadN $ncpus --genomeDir $GENOME_DIR --sjdbGTFfile $refGTF --sjdbOverhang 100 --readFilesIn $trimmed_R1_fq $trimmed_R2_fq --chimSegmentMin 10 --outFileNamePrefix $star_prefix 

# 2.1 CIRCExplorer on STAR
mkdir $CIRCBASE_DIR
mkdir $CIRCBASE_DIR/$PREFIX
mkdir $CIRCBASE_DIR/$PREFIX/star
chim_jun=$STAR_DIR/$PREFIX"_Chimeric.out.junction"

$CIRC_PATH/CIRCexplorer2 parse -t STAR $chim_jun -o $CIRCBASE_DIR/$PREFIX/star > $LOG_DIR/$PREFIX"_"CIRCexplorer2_parse.log
$CIRC_PATH/CIRCexplorer2 annotate --low-confidence  -r $refFlat -g $refgenomeFA $CIRCBASE_DIR/$PREFIX/star > $LOG_DIR/$PREFIX"_"CIRCexplorer2_annotate.log

# END
