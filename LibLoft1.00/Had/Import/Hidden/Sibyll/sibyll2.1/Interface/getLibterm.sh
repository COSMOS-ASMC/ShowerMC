#!/bin/bash
#  to get list of files corresponding to complied  files
#  for sibyll.  (.o files)

awk '$0 ~ ".o"' ../Makefile | grep sibyll | grep LIB | awk -F "(" '{print $3}' | awk -F ")" '{print $1}' | grep ".o" | sort | uniq  > sibyllLibList

awk '$0 ~ ".o"'  Makefile | grep LIB | awk -F "(" '{print $3}' | awk -F ")" '{print $1}' | grep ".o" |sort |uniq >> sibyllLibList

