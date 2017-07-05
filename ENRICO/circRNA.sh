#############################################################
##### CIRCULAR-RNA ANNOTATION/CHARACTERIZATION PIPELINE #####
#############################################################

# Reference Genome : GRCh38.p7
# Annotations : Gencode v25
# Aligner : Tophat-Fusion (paired-end) // Note: Tophat is the only pe output supported by CircExplorer
# CircRNA algorithm : circExplorer2


PREFIX=HCC1599_ATRA_REP1
R1_fq=/pico/scratch/userexternal/mfratell/circRNA/fastq/HCC1599_ATRA_REP2_S5_R1_001.fastq
R2_fq=/pico/scratch/userexternal/mfratell/circRNA/fastq/HCC1599_ATRA_REP2_S5_R2_001.fastq

BASE_DIR=/pico/scratch/userexternal/mfratell/circRNA
TRIM_DIR=/pico/scratch/userexternal/mfratell/circRNA/trimmedfq
FQ_DIR=/pico/scratch/userexternal/mfratell/circRNA/fastq
QC_DIR=/pico/scratch/userexternal/mfratell/circRNA/fastQC
tQC_DIR=/pico/scratch/userexternal/mfratell/circRNA/tfastQC
ADAPTER_R1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
ADAPTER_R2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
refgenomeFA=/pico/scratch/userexternal/mfratell/GENOME/GENCODE/v25/GRCh38.p7.genome.fa
refGTF=/pico/scratch/userexternal/mfratell/GENOME/GENCODE/v25/gencode.v25.annotation.gtf
refFlat=/pico/scratch/userexternal/mfratell/GENOME/GENCODE/v25/gencode.v25.refFlat.txt 
STAR_PATH=/pico/scratch/userexternal/mfratell/SW/STAR-2.5.1b/bin/Linux_x86_64_static
GENOME_DIR=/pico/scratch/userexternal/mfratell/GENOME/GENCODE/v25/star_GRCh38.p5 
CIRC_PATH=/pico/home/userexternal/mfratell/.local/bin
LOG_DIR=/pico/scratch/userexternal/mfratell/circRNA/log
CIRCBASE_DIR=/pico/scratch/userexternal/mfratell/circRNA/circ_out
STAR_DIR=/pico/scratch/userexternal/mfratell/circRNA/star
TOPHAT_DIR=/pico/scratch/userexternal/mfratell/circRNA/tophat
GRCh38_bowtie1_index=/pico/scratch/userexternal/mfratell/GENOME/bowtie1_gencodeV25/GRCh38.p7
ncpus=20


# 0. Load Modules (CINECA)
module load profile/advanced
module load autoload cutadapt
module load autoload samtools/0.1.19
module load fastqc
module load python/2.7.9 # per circexplorer2



# 1. QUALITY CONTROL
mkdir $QC_DIR
fastqc -t 20 -o $QC_DIR $R1_fq
fastqc -t 20 -o $QC_DIR $R2_fq

# 2. ADAPTER TRIMMING (TRUSEQ TOTALRNA STRANDED)
mkdir $TRIM_DIR
MIN_LENGTH=18
ERROR_RATE=0.2
OVERLAP=5
trimmed_R1_fq=$TRIM_DIR/$PREFIX"_trimmed_R1.fastq"
trimmed_R2_fq=$TRIM_DIR/$PREFIX"_trimmed_R2.fastq"
cutadapt -a ${ADAPTER_R1} -A ${ADAPTER_R2} --minimum-length $MIN_LENGTH --error-rate $ERROR_RATE --overlap=$OVERLAP -o $trimmed_R1_fq --paired-output $trimmed_R2_fq $R1_fq $R2_fq

# 3. QUALITY CONTROL AFTER TRIMMING
mkdir $tQC_DIR
fastqc -t 20 -o $tQC_DIR $trimmed_R1_fq
fastqc -t 20 -o $tQC_DIR $trimmed_R2_fq

# 4A. STAR PIPE
trimmed_R1_fq=$TRIM_DIR/$PREFIX"_trimmed_R1.fastq"
trimmed_R2_fq=$TRIM_DIR/$PREFIX"_trimmed_R2.fastq"
mkdir $BASE_DIR/$PREFIX
mkdir $STAR_DIR
star_prefix=$STAR_DIR/$PREFIX"_"
$STAR_PATH/STAR --runThreadN $ncpus --genomeDir $GENOME_DIR --sjdbGTFfile $refGTF --sjdbOverhang 100 --readFilesIn $trimmed_R1_fq $trimmed_R2_fq --chimSegmentMin 10 --outFileNamePrefix $star_prefix 

mkdir $CIRCBASE_DIR
mkdir $CIRCBASE_DIR/$PREFIX
mkdir $CIRCBASE_DIR/$PREFIX/star
chim_jun=$STAR_DIR/$PREFIX"_Chimeric.out.junction"

$CIRC_PATH/CIRCexplorer2 parse -t STAR $chim_jun -o $CIRCBASE_DIR/$PREFIX/star > $LOG_DIR/$PREFIX"_"CIRCexplorer2_parse.log
$CIRC_PATH/CIRCexplorer2 annotate --low-confidence  -r $refFlat -g $refgenomeFA $CIRCBASE_DIR/$PREFIX/star > $LOG_DIR/$PREFIX"_"CIRCexplorer2_annotate.log


# 4B. TopHat-Fusion Pipe
mkdir $TOPHAT_DIR
mkdir $TOPHAT_DIR/$PREFIX
OUT_DIR=$TOPHAT_DIR/$PREFIX
tophat2 -o $OUT_DIR -p $ncpus --fusion-search --keep-fasta-order --bowtie1 --no-coverage-search ${GRCh38_bowtie1_index} $trimmed_R1_fq,$trimmed_R2_fq
mkdir $CIRCBASE_DIR
mkdir $CIRCBASE_DIR/$PREFIX
mkdir $CIRCBASE_DIR/$PREFIX/tophat
TopHatFus=$TOPHAT_DIR/$PREFIX"/accepted_hits.bam"
$CIRC_PATH/CIRCexplorer2 parse -o $CIRCBASE_DIR/$PREFIX/tophat -t TopHat-Fusion $TopHatFus 
$CIRC_PATH/CIRCexplorer2 annotate --low-confidence  -r $refFlat -g $refgenomeFA $CIRCBASE_DIR/$PREFIX/tophat


