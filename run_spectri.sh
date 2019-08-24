#!/bin/bash
Xvnc :2 -screen 2 1600x1200x16 & ## pour créer un serveur X virtuel sur la machine
export DISPLAY=localhost:2.2 ## pour rediriger la sortie vers ce serveur
nohup nice  idl -vm=rt_spectribase_manage.sav > out_rt.txt & ## pour lancer le job (nice permet de mettre une priorité inférieur au job, non nécessaire)
a=100
for (( i=0 ; 16 - $i ; i++ ))
	do 
	for (( j=0 ; 12-$j ; j++ ))
		do
		x=$((i*a))
		y=$((j*a))
		xdotool mousemove --screen 2 $x $y click  1 
	done
done
##rm -rf /tmp/.X2-lock 
