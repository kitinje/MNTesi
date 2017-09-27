
qsub -q parallel -I -l select=1:ncpus=20:mem=120G,walltime=12:00:00 -A ELIX2_prj9
cd /pico/scratch/userexternal/mfratell/GIADA/ageing
module load autoload r
R

quit()
exit #uscire dal job


scp mfratell@login.pico.cineca.it:/pico/scratch/userexternal/mfratell/GIADA/ageing/WKSpace27.RData ./

scp /Users/Giada/Desktop/WKSpace27.RData mfratell@login.pico.cineca.it:/pico/scratch/userexternal/mfratell/GIADA/ageing/


