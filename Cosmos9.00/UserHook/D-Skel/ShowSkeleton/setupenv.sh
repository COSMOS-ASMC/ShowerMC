#!/bin/sh 
###  don't touch next  line.   test is needed for sge job.
test  $# -eq  0  &&    source ../confirm.inc

#------------ probably need not touch next line
TOPDIR=$COSMOSTOP/UserHook/DisParaForTA2
export TOPDIR
#   dir to store execution command (ssh)  or submit command (sge)
EXECDIR=$TOPDIR/ShowSkeleton/Exec
export EXECDIR
#     id  emmbeded in the command name as mentioned above. must not start with 
#    number for sge  jobs.
EXECID=p20E999-40
export EXECID
#   
#         at which depth do you take histogram; the index must not exceed
#         really existing  depth index; otherwise, some 0 data will appear
#         in the hist data (Not fatal though).
#        give a list depth index where you want to take histograms'
#        
#  2000  3000  4000  4250  4500  4750  5000  5250  5500  5750 
#  6000  6250  6500  6750  7000  7250  7500  7750  8000  8250
#  8500  8750  9000  9250  9500  9750 10000 10250 10500 10750
# 11000 11250 11500 11750 12000 14000

HISTDEP='11 13 15 17 19 21 23 25  27 31 35/'
export HISTDEP
#         at which depth do you output individual ptcl info.
#         give  that depth index 
#  
INDIVDEP='1 2 3 4  5 6  7 8  9 10  11 12  13 14 15 16  17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35/' 
export INDIVDEP
#   OUTPUT:   which histogram do you take.
#    first one is not for histogram, but for indiviual info.
#    (output x,y or only r) first one is not used in this prog.
#      1     2        3           4
#    recxy, tklat, tkelosslat, tkrespec,
#        5         6          7       8
#    tkrzspec, tkzfspec,  tkrfspec, tkefspec,
#          9         10         11         12      13      14 
#    tkertspec, tkrtspec,  tkrezspec, tkrzfspec, tkrefspec, tkarspec
#
OUTPUT='t f t t f f  f f f f f f f t/' 
export OUTPUT
#     where param001 etc exist
PARAMDIR=$TOPDIR/ShowSkeleton/ParamDir
export PARAMDIR
#     where to put data from each host
OUTDIR=$TOPDIR/ShowSkeleton/OutDir
#OUTDIR=/tmp/kasahara2
export OUTDIR
#     where to put error message
ERRDIR=$TOPDIR/ShowSkeleton/ErrDir
export ERRDIR

#  don't touch below.
#  if used from script, we skip the next lines.

if [ $# -eq  0 ] ; then
    confirm  $PARAMDIR
    confirm  $OUTDIR
    confirm  $ERRDIR
    confirm  $EXECDIR
fi
