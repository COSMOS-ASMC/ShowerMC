#!/bin/sh
#   sge
#$ -S /bin/sh 
COSMOSTOP=~/Cosmos7.51
export COSMOSTOP
cd $COSMOSTOP/UserHook/mkLDD
source  /home/Users/sgeadmin/SGE/default/common/settings.sh 
source  /Loft1/Intel/ifc/bin/ifortvars.sh

###   $COSMOSTOP,  ~/Cosmos will not work in the next
#$  -e /home/Users/kasahara/CosmosData/NewLDD/qgsjet2/p/1.00E18/cos0.700/ErrDir/$JOB_ID.$TASK_ID.$HOSTNAME.err
#$  -o  /dev/null
 ./cosmosPCLinuxIFC < param
#  ./cosmosPCLinuxIFC10.1.018athlon  < param
# ./cosmosPCLinuxIFC10.1.022athlon < param
# ./cosmosPCLinuxIFC10.1.022Xeon < param

