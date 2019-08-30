#!/bin/sh
# icr 32
#source /tibet/work/kasahara/intel/ifc/ifortvars.sh ia32
#source /tibet/work/kasahara/intel/icc/ifortvars.sh ia32
# icr 64
source  /tibet/work/kasahara/intel/ifc/ifortvars.sh intel64
source  /tibet/work/kasahara/intel/icc/iccvars.sh intel64

#    tasim 32 bit
#  source  /Loft1/Intel/ifc/bin/ifortvars.sh
#  source  /Loft1/Intel/icc/bin/iccvars.sh
#   waseda 32 bit
#  source  /opt/intel/ifc/bin/ifortvars.sh
#  source  /opt/intel/icc/iccvars_ia32.sh
#    waseda intel 64 bit
#  source   /opt/intel/ifce/ifortvars_intel64.sh
#  source  /opt/intel/icce/iccvars_intel64.sh


export COSMOSTOP=~/Cosmos
export COSMOSINC=$COSMOSTOP/cosmos
export PATH=$COSMOSTOP/Scrpt:$PATH
# export GXX_ROOT=/usr/lib/gcc/i486-linux-gnu/4.3.3/
if [ ! "${HOST}" ];then
    export HOST=$HOSTNAME
fi

#   how to use. ********************************
#  0)  edit setup.sh (this file) so that compliler (ifort and icc) are
#       propery accessed.  Fix COSMOSTOP
#     If there is a change in the source file names and/or header files
#     reflect it in flist and/or hlist
#     
#  
#  1)  bash (if you shell is not bash)
#  2)  source setup.sh
#  3)  bash compile.sh flist hlist 0
#       (for 32 bit)  or
#      bash compile64.sh flist hlsit 0
#   Then you will get lib src include directories
#  4)  make executable by
#     ifort -O0 -Llib chook.o -o cosmosPCLinuxIFC -lcosmos 
#   or (next is for 64, though above is also ok for 64bit)
#     ifort -O0 -Llib chook.o -o cosmosPCLinuxIFC64 -lcosmos
#  next msg may come out depending on icc
#icc: command line remark #10010: option '-prof-genx' is deprecated
#  and will be removed in a future release. See '-help deprecated'


