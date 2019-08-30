#!/bin/sh
#   sge
#$ -S /bin/sh 
COSMOSTOP=~/Cosmos7.51
export COSMOSTOP
cd $COSMOSTOP/UserHook/mkLDD
source  /home/Users/sgeadmin/AMD/default/common/settings.sh 
source  /Loft1/Intel/ifc/bin/ifortvars.sh

#$  -e /home/Users/kasahara/CosmosData/NewLDD/qgsjet2/p/1.00E18/cos0.700/ErrDir/$JOB_ID.$TASK_ID.$HOSTNAME.err
#$  -o  /dev/null
 ./cosmosPCLinuxIFC < param
