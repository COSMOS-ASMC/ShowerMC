#!/bin/bash
# echo arch is $ARCH
(cd $COSMOSTOP/Import/Hidden/Sibyll/sibyll2.1/Interface;
      for f in `cat sibyllLibList`;do
##	  echo $f
	  ar -d $COSMOSTOP/lib/$ARCH/libcosmos.a $f
      done
)
