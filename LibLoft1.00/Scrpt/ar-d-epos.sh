#!/bin/bash
#  remove those object files in libloft.a that were created in
#       each EPOS model; they should be listed in "eposlibList"
#       in Interface of each EPOS 
# echo arch is $ARCH

(cd $LIBLOFT/Had/Import/Hidden/EPOS/;
 for d in epos*; do
     if [ -d $d ]; then
	 # $d should be epos-lhc-3400 etc
	 (cd $d; (
	     cd Interface; 
	     for f in `cat eposlibList`;do
	      ##	  echo $f
		 ar -d $LIBLOFT/lib/$ARCH/libloft.a $f
	     done)
	 )
     fi
 done
)
