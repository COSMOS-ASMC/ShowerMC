#!/bin/sh  -f
if [  $#  != 2  ]; then
	echo 'Usage: ./execreduce_all.sh HostList_path  ssh'
	exit
fi
how=$2

#       save current stdin and change it to $1
exec 3<&0  <$1

mkdir -p $EXECDIR

while read numb host cpupow
 do
    test $numb = "#" && continue
    test -z $numb  && continue
    numb=` echo $numb |  awk '{printf("%5.5d"), $1}' `
    exechost=$host
#    echo $numb  $exechost
#   actual exec command.
    if [ $how = "ssh" ]; then
	sed "/^[^#]/s/__/$numb/g"  execSSHtemplate.sh  > exec.sh
	chmod +x exec.sh
	echo "command used at host $exechost (with  $numb ) is"
        awk 'substr($1,1,1)!="#"' ./exec.sh
        mv ./exec.sh $EXECDIR/$EXECID"-"$numb".sh"
###              next background spec. is need otherwise stdin is finished
	ssh $exechost  $EXECDIR/$EXECID"-"$numb".sh" &
	sleep 1
    else
	echo "you have to use ssh"
	exit
    fi
 done
#   restore stdin; close 3
 exec 0<&3 3<&-
