#!/bin/bash
(cd $LIBLOFT/Had/Import/Hidden/EPOS/epos199/Interface;
      for f in `cat eposlibList`;do
	  ar -d $LIBLOFT/lib/$ARCH/libloft.a $f
      done
)
