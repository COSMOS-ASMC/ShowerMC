#!/bin/bash
ZA=(0 1 1 1 2 2 3);

for A in `seq  6`; do
   Z=${ZA[$A]}   
   sed  -e  "s/AAA/$A/" -e "s/ZZZ/$Z/"  dpmjet.inpSave > dpmjet.inp
   ./a.out > K-${A}.dat
done

for A in `seq  8 2 32`; do
   let Z=$A/2
   sed  -e  "s/AAA/$A/" -e "s/ZZZ/$Z/"  dpmjet.inpSave > dpmjet.inp
   ./a.out > K-${A}.dat
done

for A in `seq  36  4 76`; do
   let Z=$A/2; let Z-=2;
   sed  -e  "s/AAA/$A/" -e "s/ZZZ/$Z/"  dpmjet.inpSave > dpmjet.inp
   ./a.out > K-${A}.dat
done


for A in `seq  80 10 140`; do
   let Z=($A-10)/2;
   sed  -e  "s/AAA/$A/" -e "s/ZZZ/$Z/"  dpmjet.inpSave > dpmjet.inp
   ./a.out > K-${A}.dat
done


for A in `seq  150 10 210`; do
   let Z=$A/3+12;
   sed  -e  "s/AAA/$A/" -e "s/ZZZ/$Z/"  dpmjet.inpSave > dpmjet.inp
   ./a.out > K-${A}.dat
done

