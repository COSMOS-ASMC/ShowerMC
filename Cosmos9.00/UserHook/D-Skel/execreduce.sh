#!/bin/sh -f
if [  $# != 1 ];  then
    echo 'Usage ../execreduce.sh ssh'
    echo 'This must be invoked from ReduceEachSize'
    exit
else
#    if [ $1 != "sge" ] && [ $1 != "ssh" ]; then
    if  [ $1 != "ssh" ]; then
	echo "$1 for the argument is invalid"
	echo "At  present ssh  job submission is possilbe"
	exit
    fi
fi

##   next argument is for avoiding QA in the script
source ./doreduce $0
#  complete name of $FLESHDIR is needed
REDUCEDIR=`pwd`
export REDUCEDIR
echo $REDUCEDIR

source  ../execreduce_all.sh  $HOSTLIST $1
