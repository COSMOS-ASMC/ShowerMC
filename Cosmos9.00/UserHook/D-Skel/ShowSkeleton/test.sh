#!/bin/sh -f
function canIrm {
    echo dolloar 1 is $1
    some="`ls $1/`"
    echo $some
    if [ -n "$some" ] ; then
	echo "some in "
    fi
}
function quit {
    exit
}
function e {
    echo $1
}
e Hello
e World
canIrm $PARAMDIR
echo $0

