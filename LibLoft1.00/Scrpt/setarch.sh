#!/bin/bash
export ARCH=`awk '$1=="ARCH" && $2=="=" {print $3}' $LIBLOFT/site.config`

