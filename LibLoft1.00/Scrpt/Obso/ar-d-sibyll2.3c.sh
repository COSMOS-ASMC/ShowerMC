#!/bin/bash
# echo arch is  $ARCH
(cd $COSMOSTOP/Import/Hidden/Sibyll/sibyll2.3c/Interface;
      for f in `cat sibyllLibList`;do
##	  echo $f
	  ar -d $COSMOSTOP/lib/$ARCH/libcosmos.a $f
      done
)
