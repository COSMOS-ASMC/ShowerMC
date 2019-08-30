#!/bin/bash                                                                       
#$  -S /bin/bash
COSMOSTOP=/nfs/RAID/home/kasahara/Develop/Cosmos7.635
export COSMOSTOP
#$ -o /nfs/RAID//home/kasahara/Cosmos/UserHook/modifyX/DataM/p17L-$JOB_ID-$HOSTNAME.dat
#$ -e /nfs/RAID/home/kasahara/Cosmos/UserHook/modifyX/ErrM/p17L-$JOB_ID-$HOSTNAME.err


cd /nfs/RAID/home/kasahara/Cosmos/UserHook/modifyX
# echo  start cosmos
#./cosmosPCLinuxIFC64 < paramSpecial
./cosmosPCLinuxIFC64 < param
# pwd
