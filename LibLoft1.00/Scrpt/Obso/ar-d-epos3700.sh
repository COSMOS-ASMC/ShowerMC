#!/bin/bash
# echo arch is  $ARCH
(cd $LIBLOFT/Had/Import/Hidden/EPOS/epos-lhc-v3700/Interface;
      for f in `cat eposlibList`;do
##	  echo $f
	  ar -d $LIBLOFT/lib/$ARCH/libloft.a $f
      done
)
