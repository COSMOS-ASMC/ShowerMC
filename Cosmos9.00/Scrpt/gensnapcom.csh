#!/bin/csh 
#  generate snapshot command for a given number of frmaes
if( $#argv != 2 ) then
  echo  './gensnapcom.csh  fisrtframe# #offrames'
  exit
endif
awk 'END {for(i=0; i<nf; i++) printf( "(snapshot focus ts%5.5d.ppm ppmscreen)\n",n+i)}'  nf=$2 n=$1 /dev/null
