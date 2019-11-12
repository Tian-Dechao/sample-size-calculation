#!/bin/bash
source /hive/dechaot/results/sample_size/code/config
N=20
t=4
#N=$1
#t=$2
rep=1000
#o=$dirRes/methylome_"$N"_"$t".txt
#opdf=$dirRes/methylome_"$N"_"$t".pdf
o=$dirRes/methylome_"$N"_"$t"_lmm.txt
opdf=$dirRes/methylome_"$N"_"$t"_lmm.pdf
## part1 generate data
#rm -f $o
## parameters for parallel 
#Np=250
#for i in {50..500..25}; do
#    ((k=i%Np)); ((k++==0)) && wait
##    (j=$(bc<<<"scale=5; $i/1000");$rscript $dirCode/sample_power.R methylo $N $t $j $rep >> $o) &
##    (j=$(bc<<<"scale=5; $i/1000");$rscript $dirCode/sample_power.R methylo frideman $N $t $j $rep) 
#    (j=$(bc<<<"scale=5; $i/1000");$rscript $dirCode/sample_power.R methylo  lmm $N $t $j $rep >> $o) &
#done
#wait
$rscript $dirCode/sample_power.R methyloFig $o $opdf --xlim 0.33
scp $opdf imac:/Users/Dechao/Documents/sample-size/
exit
# part2 create figure 
$rscript $dirCode/sample_power.R methyloFig $o $opdf
#$rscript $dirCode/sample_power.R methyloFig $o $opdf --xlim 0.33
scp $opdf imac:/Users/Dechao/Downloads/
exit
exit
#$rscript $dirCode/sample_power.R methylo $N $t 0.1 $rep
#$rscript $dirCode/sample_power.R methylo $N $t 0.4 $rep

