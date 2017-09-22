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




