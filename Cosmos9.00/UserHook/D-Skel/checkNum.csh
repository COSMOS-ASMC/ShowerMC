#!/bin/csh -f
sort -n Hosts | awk '{print $1}' | uniq -d
