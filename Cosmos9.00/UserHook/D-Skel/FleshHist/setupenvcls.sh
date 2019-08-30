#!/bin/bash 
###  don't touch next  line.   test is needed for sge job.
test  $# -eq  0  &&    source ../confirm.inc


source ../setupenvcls.sh
TOPDIR=$COSMOSTOP/UserHook/$ARENA
export TOPDIR
#     id  emmbeded in the command name as mentioned above. must not start with 
#    number for sge  jobs.
EXECID=qgsp17cos1.000
export EXECID

HISTDEP='29 /'
export HISTDEP
#         at which depth do you output individual ptcl info.
#         give  that depth index 
#  
INDIVDEP='13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39  /' 
export INDIVDEP
#   OUTPUT:   which histogram do you take.
#  
# 1  time histo in each web sector
# 2  lateral in each fai bin   
#
OUTPUT=' t t /'
export OUTPUT
# to see E<500keV(=standard KEminObs) contribution to 
#   dE/dx, make next yes else no . (if yes, KEminObs must< 500keV)
SeeLowdE=no
export SeeLowdE

#     where to put data from each host
#OUTDIR=$TOPDIR/Assemble/OutDir
OUTDIR=/tmp/$USER
export OUTDIR
# observe weighted particle (by thinning) as 1 individual
# particle with given weight or n partciles with weight
# 1.  (n is an integer proportinal to  weight. )
#------------ probably need not touch next line
KeepWeight=yes
# even if yes is given, ThinSampling is F==> set no
temp=`awk ' ($1=="ThinSampling" || $1=="THINSAMPLING")  && $3~"F" {print "no"}' $TOPDIR/FleshHist/Sparam` 
if [ x$temp =  "xno" ] ; then
 KeepWeight=no
fi
export KeepWeight

#   
#         at which depth do you take histogram;
#  histogram can be output in ascii or binary format
#  One format can be converted into another format.
#  However, if it is bianry, it cannot be read on
#  a different host than on the creator.
#  BINW=1 for ascii write =2 for binary write
#  For assembling, we must use binary ouput
BINW=2
export BINW
#  We record a maximum of ~7500 particles of each particle
#  type in a web sector (7500 is a default) by selecting
#  particles  falling on the sector randomly. 
#  However, it is difficult to estimate how many
#  particles fall on the web sector,  so is the probabilty 
#  of accepting them for each web sector.  Therefore,
#  we accept much larger number of particles at first and
#  if the resultant number exceed 7500, we drop some of them
#  randaomly.  The LIMIT next is the first target number
#  of particles to be accepted.
LIMIT="20000 20000 20000 20000"
export LIMIT
#   dir to store execution command (ssh)  or submit command (sge)
EXECDIR=$TOPDIR/FleshHist/Exec
export EXECDIR
#     where param001 etc exist
PARAMDIR=$TOPDIR/FleshHist/ParamDir
export PARAMDIR
#     where to put error message
ERRDIR=$TOPDIR/FleshHist/ErrDir
export ERRDIR

#  don't touch below.
#  if used from script, we skip the next lines.

if [ $# -eq  0 ] ; then
    confirm  $PARAMDIR
    confirm  $OUTDIR
    confirm  $ERRDIR
    confirm  $EXECDIR
 if [ -e Sparam ]; then
   primaryf=` awk '$1=="PRIMARYFILE" {i1=index($0, qm) ; i2=index(substr($0, i1+1), " "); print substr($0, i1+1,i2-1)}' qm="'" Sparam ` 
  echo $primaryf " seems to be the primary file used at SkelFlesh"
  echo "now it is copied to this directory"
  cp ${COSMOSTOP}/UserHook/SkelFlesh/$primaryf ./
 else
    echo "There is no Sparam file which is a standard copy of "
    echo " parameter file used in skeleflesh process"
    echo " **make a copy as Sparam and change Job='newskel' into"
    echo " Job='newflesh'"
    exit
 fi
fi
