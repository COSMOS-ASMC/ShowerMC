#!/bin/sh -f
export ARCH=`awk '$1=="ARCH" && $2=="=" {print $3}' $EPICSTOP/site.config`



