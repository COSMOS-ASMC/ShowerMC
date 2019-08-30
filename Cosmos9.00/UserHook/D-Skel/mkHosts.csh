#!/bin/csh 

set vcpu=`awk '$1~"NCPU=" {print substr($1,6);exit}' Smash/setupenvcls.sh`
#  echo "# of vertual cpus is $vcpu"


@ num=0
while ( $num < $vcpu )
    foreach f ("`cat allHosts`")
        set col1=`    echo $f | awk '{print $1}'`
	if( $col1 == "#" )  continue
	if( $col1 == "" ) continue
	set ncpu=`    echo $f | awk '{print $3}'`
	set host=`    echo $f | awk '{print $2}'`
	set cpupw=`   echo $f | awk '{print $4}'`


	@ c=0 
	while ( $c < $ncpu ) 
	    @ num++
	    if( $num > $vcpu ) exit
	    echo $num  $host $cpupw
	    @ c++
	end
    end
end
