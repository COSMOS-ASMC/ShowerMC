#!/bin/bash
#  remove those object files in libcosmos.a that were created in
#       each atmosphere  model; they should be listed in "atmoslibList"
#       in Cosmos/Atmosphere/
# echo arch is $ARCH

(cd $COSMOSTOP/Atmosphere/;
 source $COSMOSTOP/Scrpt/setarch.sh;
 for f in `cat atmoslibList`;do
#     echo arch is $ARCH
	      ##	  echo $f
     ar -d $COSMOSTOP/lib/$ARCH/libcosmos.a $f
 done
)
 
