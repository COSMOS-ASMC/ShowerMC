#!/bin/sh  
if [  $#  != 2  ]; then
	echo 'Usage: ./execflesh_all.sh HostList_path  [sge|ssh]'
	exit
fi
how=$2

#       save current stdin and change it to $1
exec 3<&0  <$1

mkdir -p $EXECDIR
count=0
while read numb host cpupow
 do
    test $numb = "#" && continue
    test -z $numb  && continue
    numb=` echo $numb |  awk '{printf("%5.5d"), $1}' `
    exechost=$host
    count=`expr $count + 1`

##    echo $numb  $exechost
#   actual exec command.
    if [ $how = "ssh" ]; then
	sed "/^[^#]/s/__/$numb/g"  execSSHtemplate.sh  > exec.sh
	chmod +x exec.sh
	echo "command used at host $exechost (with  $numb ) is"
        awk 'substr($1,1,1)!="#"' ./exec.sh
        mv ./exec.sh $EXECDIR/$EXECID"-"$numb".sh"
###              next background spec. is need otherwise stdin is finished
#	ssh $exechost  $FLESHDIR/exec.sh &
	ssh $exechost  $EXECDIR/$EXECID"-"$numb".sh" &
	if [ $[ $count - $count/10*10] -eq 0 ] ; then
	    sleep 2
	fi
    elif [ $how = "sge" ]; then
	sed -e "s|__|$numb|g"  -e "s|ERRDIR|$ERRDIR|g" -e "s|OUTDIR|$OUTDIR|g" \
	    -e "s|EXECID|$EXECID|g" execSGEtemplate.sh > exec.sh
       echo "command used for cpu $numb  is"
       awk 'substr($1,1,1)!="#"' ./exec.sh
       chmod +x exec.sh
       mv ./exec.sh $EXECDIR/$EXECID"-"$numb".sh"
       qsub   $EXECDIR/$EXECID"-"$numb".sh"
       if [ $[ $count - $count/10*10] -eq 0 ] ; then
	   sleep 2
       fi
    else
	echo "you have to use sge or ssh"
	exit
    fi
    if [ $count -eq $NOOFCPU ]; then
	break
    fi
 done
 if [ $count -lt $NOOFCPU ]; then
     echo "# of cpus in " $1 "=" $count "<  requested cpus="  $NOOFCPU
 fi

#   restore stdin; close 3
 exec 0<&3 3<&-
