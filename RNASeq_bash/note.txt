note:
vi qc.sh
premi A
ESC + duepunti + wq!
fai partire un job

cd /pico/scratch/userexternal/mfratell/MARTINA/RNASeq
cd $CINECA_SCRATCH/MARTINA/RNASeq

qsub -q parallel -I -l select=1:ncpus=20:mem=120G,walltime=24:00:00 -A ELIX_prj7
qsub -q rcm_visual -I -l select=1:ncpus=20:mem=120G -A ELIX_prj7
