#!/bin/bash
#  remove those object files in libloft.a that were created in
#       each SIBYLL model; they should be listed in "sibyllLibList"
#       in "Interface" of each SIBYLL


(cd $LIBLOFT/Had/Import/Hidden/Sibyll/;
 for d in sibyll*; do
     if [ -d $d ]; then
	 # $d should be epos-lhc-3400 etc
	 (cd $d; (
	     cd Interface; 
	     for f in `cat sibyllLibList`;do
	      ##	  echo $f
		 ar -d $LIBLOFT/lib/$ARCH/libloft.a $f
	     done)
	 )
     fi
 done
)
