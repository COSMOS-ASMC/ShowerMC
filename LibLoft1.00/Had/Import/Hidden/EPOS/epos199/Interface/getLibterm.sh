#!/bin/bash
#  to get list of files corresponding to complied  files
#  for EPOS.  (.o files)
#  This is used to change EPOS (3400<->3700)
awk '$0 ~ ".o"' ../Makefile | grep epos | grep LIB | awk -F "(" '{print $3}' | awk -F ")" '{print $1}' | grep ".o" | sort | uniq  > eposlibList

awk '$0 ~ ".o"'  Makefile | grep LIB | awk -F "(" '{print $3}' | awk -F ")" '{print $1}' | grep ".o" |sort |uniq >> eposlibList

