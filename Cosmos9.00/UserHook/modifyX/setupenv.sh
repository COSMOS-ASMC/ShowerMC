#!/bin/sh
###  don't touch next  line.   
test  $# -eq  0  &&    source ./confirm.inc



if [ $# -eq  0 ] ; then
    confirm  OutDir
    confirm  ErrDir
    confirm  Seed
fi
