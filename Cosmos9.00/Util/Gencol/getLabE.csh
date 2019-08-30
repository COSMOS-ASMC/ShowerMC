#!/bin/csh
awk -f getLabE1.awk $1 > data$$
awk -f getLabE2.awk data$$
rm -f data$$

