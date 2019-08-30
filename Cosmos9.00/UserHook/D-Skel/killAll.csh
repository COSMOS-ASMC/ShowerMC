#!/bin/csh -f
#  kill all jobs of which job name coicides with the given name
#  host given in arg. is examined
if ( $#argv  != 2 ) then
    echo "Usage  ./killAll.csh jobname  Hostfile"
    exit
endif
echo $1 " is the job name (it may be a part of the job name)"
echo "Enter y, if it is ok"
set yesno=$<
if( $yesno  != "y") exit

sed s/jobname/$1/ killcom.csh > ~/akillcom.csh
chmod +x ~/akillcom.csh
foreach host(`awk '{print $2}'  $2`)
    echo $host " is being  examined"
    ssh $host ~/akillcom.csh
    sleep 1
end

