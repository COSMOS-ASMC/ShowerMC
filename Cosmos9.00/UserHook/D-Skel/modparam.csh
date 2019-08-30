#!/bin/csh -f
#
#   This is used by FleshBasic/execflesh command.
#
#   input: 1)  param file given by  SkeletonFile in SkelFlesh
#              This should be fed to this script
#          2)  Smash/setupenv
#          3)  ./setupenv
#   output:  param files for smashed skeleton files
#
if ( $#argv != 1 ) then
	echo "Usage: ../modparm  paramfile(pecified by SkeletonFile in the skeleton making parameter)
	echo "Typically ../modparm Sparm"
	exit
endif
if ( ! -f $1 ) then
    echo "The parameter file $1 not exist" 
    exit
endif

if ( ! -f ../Smash/setupenv )  then
    echo "../Smash/setupenv not  exists"
    exit
endif
if ( ! -f ./setupenv )  then
    echo "./setupenv  not exists"
    exit
endif


@  cpu =  0
foreach f( "`cat $HOSTLIST`" )
	@ cpu++;
	set sharp = ` echo $f | awk '{print $1}' `;
	if( $sharp == '#' ) continue ;
	set numb = ` echo $f |  awk '{printf("%5.5d"), $1}' `;
	set userhookc = "'$SKELDIR/$SKELNAME$numb'";
	set paramfile = "'$PARAMDIR/param$numb'";
	awk -f ../modify.awk  UserHookc=$userhookc paramfile=$paramfile $1 > $PARAMDIR/param$numb
end
if ( $cpu != $NCPU ) then
    echo "no. of cpu is inconsistnet; NCPU=$NCPU but count=$cpu"
endif

echo "parameter files have been created in $PARAMDIR"
























































