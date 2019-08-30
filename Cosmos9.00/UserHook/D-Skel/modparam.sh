#!/bin/sh -f

#
#   input: 1)  param file given by  SkeletonFile in SkelFlesh
#              This should be fed to this script
#          2)  Smash/setupenv.sh
#          3)  ./setupenv.sh
#   output:  param files for smashed skeleton files
#
if [ $# != 1 ]; then
	echo "Usage: ../modparm.sh  paramfile(pecified by SkeletonFile in the skeleton making parameter)"
	echo "Typically ../modparm.sh Sparm"
	exit
fi
if [ ! -f $1 ]; then
    echo "The parameter file $1 not exist" 
    exit
fi

if [ ! -f ../Smash/setupenv.sh ];  then
    echo "../Smash/setupenv.sh not  exists"
    exit
fi
if [ ! -f ./setupenv.sh ];  then
    echo "./setupenv.sh  not exists"
    exit
fi

##  change stdin
if [ $MCPU = $NCPU ]; then
    hostx=$HOSTLIST
    cpux=$NCPU
else
    hostx=$THINHOSTLIST
    cpux=`expr $MCPU+$MARGIN`
fi

exec 3<&0  <$hostx
cpu=0
while read num host
do 
####  echo "num=" $num "host=" $host
  test $num = "#" && continue
   cpu=`expr $cpu + 1`
   numb=`echo $num |  awk '{printf("%5.5d"), $1}' `
###  echo "numb = $numb"
   userhookc="'$SMSKELDIR/$SKELNAME$numb'"
   paramfile="'$PARAMDIR/param$numb'"
   awk -f ../modify.awk  UserHookc=$userhookc paramfile=$paramfile $1 > $PARAMDIR/param$numb
done
if [ $cpu -ne $cpux ];  then
    echo "no. of cpu is inconsistnet; NCPU=$NCPU and"
    echo "  MCPU=$MCPU  in ../Smash/setupenv.sh"
    echo "but those in  $HOSTLIST is $cpu or"
    echo "those in  $THINHOSTLIST is $cpu"
    exit
fi
## restore stdin; close 3
exec 0<&3  3<&-
echo "parameter files have been created in $PARAMDIR"
























































