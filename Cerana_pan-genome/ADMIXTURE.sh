#!/usr/local/bin/bash
for k in {1..10};do
    #nohup /mnt/cheng/software/admixture_linux-1.3.0/admixture -j10 -C 0.01 --cv admixture.ped ${k} >admixture.log${k}.out &
    nohup /mnt/cheng/software/admixture_linux-1.3.0/admixture -j20 --cv QC.bed ${k} > log${k}.out &

done
