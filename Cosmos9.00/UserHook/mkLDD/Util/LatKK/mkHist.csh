#!/bin/csh
#  from output by getR...csh, this make histograms
#  of 2pir rho at R=50,30,20,10,7,5 m.u
#  input. file name.

set file=$1:r

awk '$1==50. {print log($2)*0.4343}' ${file}.data |histo -8 0.04 > ${file}at50.hist

awk '$1==30. {print log($2)*0.4343}' ${file}.data |histo -8 0.04 > ${file}at30.hist

awk '$1==20. {print log($2)*0.4343}' ${file}.data |histo -7 0.02 > ${file}at20.hist


awk '$1==10. {print log($2)*0.4343}' ${file}.data |histo -7 0.02 > ${file}at10.hist


awk '$1==7. {print log($2)*0.4343}' ${file}.data |histo -6 0.01 > ${file}at7.hist

awk '$1==5. {print log($2)*0.4343}' ${file}.data |histo -6 0.01 > ${file}at5.hist



