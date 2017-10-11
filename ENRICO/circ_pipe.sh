#############################################################
##### TOPHAT & STAR ALIGNERS PIPELINE IN CIRCEXPLORER2 #####
#############################################################

# Reference Genome : GRCh38.p7
# Annotations : Gencode v25
# Aligner : Tophat-Fusion (paired-end) // Note: Tophat is the only pe output supported by CircExplorer
# CircRNA algorithm : circExplorer2

# Collegati al server
ssh -t mfratell@login.pico.cineca.it

# Fai partire tutta la pipeline da un job col comando
qsub -q rcm_visual -I -l select=1:ncpus=20:mem=120G,walltime=24:00:00 -A ELIX2_prj9
# se non va bigmem usa rcm_visual o parallel o cambia ELIX â€”> ELIX_prj7, se vuoi vedere questi dati digita â€œsaldo -bâ€

# Vai nella cartella principale. Per vedere come funziona un comando utilizzato nella piattaforma digita â€œhistory | grep â€œNOME_COMANDOâ€â€
cd /pico/scratch/userexternal/mfratell/ENRICO/circRNA

#Variabili mie. Per vedere il valore di una variabile digita â€œecho $NOME_VARIABILEâ€
mkdir /pico/scratch/userexternal/mfratell/ENRICO/circRNA/star_new
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

# PREFIX Ã¨ il nome della linea seguito da _ATRA o _CONTROL
PREFIX=CAL851_ATRA

# 0. Load Modules (CINECA)
module load profile/advanced
module load autoload cutadapt
module load autoload samtools/0.1.19
module load fastqc
module load python/2.7.9
module load bowtie

# 4C. Tophat-Fusion Pipe Corrected
# Just once
mkdir $TOPHAT_DIR
# - - - - - - - - - - - - - -
trimmed_R1_fq=$TRIM_DIR/$PREFIX"_R1"  
trimmed_R2_fq=$TRIM_DIR/$PREFIX"_R2"
mkdir $TOPHAT_DIR/$PREFIX

/pico/scratch/userexternal/mfratell/ENRICO/SW/tophat_2.1.0/tophat2 -o $OUT_DIR -p $ncpus --fusion-search --keep-fasta-order --bowtie1 --no-coverage-search ${GRCh38_bowtie1_index} $trimmed_R1_fq,$trimmed_R2_fq
mkdir $CIRCBASE_DIR
mkdir $CIRCBASE_DIR/$PREFIX
mkdir $CIRCBASE_DIR/$PREFIX/tophat
TopHatFus=$TOPHAT_DIR/$PREFIX/"accepted_hits.bam"
$CIRC_PATH/CIRCexplorer2 parse -o $CIRCBASE_DIR/$PREFIX/tophat -t TopHat-Fusion $TopHatFus 
$CIRC_PATH/CIRCexplorer2 annotate --low-confidence  -r $refFlat -g $refgenomeFA $CIRCBASE_DIR/$PREFIX/tophat


# 4D. STAR PIPE NEW
# Just once
mkdir $STAR_DIR
# - - - - - - - - - - - - - - - - - -
trimmed_R1_fq=$TRIM_DIR/$PREFIX"_R1"
trimmed_R2_fq=$TRIM_DIR/$PREFIX"_R2"
mkdir $BASE_DIR/$PREFIX
star_prefix=$STAR_DIR/$PREFIX"_"
$STAR_PATH/STAR --runThreadN $ncpus --genomeDir $GENOME_DIR --sjdbGTFfile $refGTF --sjdbOverhang 100 --readFilesIn $trimmed_R1_fq $trimmed_R2_fq --chimSegmentMin 10 --outFileNamePrefix $star_prefix 

mkdir $CIRCBASE_DIR
mkdir $CIRCBASE_DIR/$PREFIX
mkdir $CIRCBASE_DIR/$PREFIX/star
chim_jun=$STAR_DIR/$PREFIX"_Chimeric.out.junction"

$CIRC_PATH/CIRCexplorer2 parse -t STAR $chim_jun -o $CIRCBASE_DIR/$PREFIX/star > $LOG_DIR/$PREFIX"_"CIRCexplorer2_parse.log
$CIRC_PATH/CIRCexplorer2 annotate --low-confidence  -r $refFlat -g $refgenomeFA $CIRCBASE_DIR/$PREFIX/star > $LOG_DIR/$PREFIX"_"CIRCexplorer2_annotate.log
