#!/bin/sh 
if [ $# != 2 ]; then
	echo 'Usage: ./execflesh_one.sh HostList_path [sge|ssh]'
	exit
fi
how=$2
## echo "how=" $how
while [ 1 ]
do
 echo  "Give a list of numbers existing  in $1 for flesh."
 read  numbers
 echo "You gave "
 echo $numbers
 echo "OK?  Enter y if yes"
 read yesno

 test "x$yesno" != "xy" && continue
 nn=0
 for  numb  in $numbers 
 do
##   echo $numb
   if [ "x$how"  = "xssh" ]; then
#        find every host in the host field where the number
#        is given in the number field. numb will be eg.  003 004 
       host=` awk '$1==numb {print  $2}' numb=$numb  $1 `
       test -z "$host" && nn=`expr $nn + 1 `
   else
       exist=` awk '$1==numb {print  $1}' numb=$numb  $1 `
       test -z "$exist" && nn=`expr $nn + 1 `
   fi
 done
 if [ $nn -ne  0 ] ; then
     echo "$nn data is invalid (not found in $1)"
     test "x$how" = "xssh" && echo "or no corresponding host in $1"
     continue
 else
     break
 fi
done

echo "all numbers  have been verified ; they are in $1"

mkdir -p  $EXECDIR

for num  in  $numbers
do
  numb=` echo $num | awk '{printf("%5.5d"), $1}' `
#   numb is now eg.  014 003 etc
  if [ $how = "ssh" ] ; then
      exec 3<&0  <$1
      while read numinp exechost cpupow
      do
	test -z $numinp && continue
	test -z $exechost && continue
	if [ $numinp -eq $num ]; then
	    break
	fi
      done
#  restore stdin and close 3
      exec 0<&3 3<&-
  fi
  if [ $how = "sge" ] ; then
    sed -e "s|__|$numb|g"  -e "s|ERRDIR|$ERRDIR|g" -e "s|OUTDIR|$OUTDIR|g" \
	-e "s|EXECID|$EXECID|g" execSGEtemplate.sh > exec.sh
  else
     sed -e "s|__|$numb|g"   execSSHtemplate.sh > exec.sh 
  fi
  chmod +x exec.sh
  echo "command used for number $num  is :"
  awk 'substr($1,1,1) != "#"' ./exec.sh
  mv ./exec.sh $EXECDIR/$EXECID"-"$numb".sh"
  if [ $how = "sge" ]; then
     qsub  $EXECDIR/$EXECID"-"$numb".sh"
     sleep 1
  elif [ $how = "ssh" ]; then
      echo "exec host is $exechost"
#           next one must be b.g. (I don't know why?)
      ssh   $exechost  $EXECDIR/$EXECID"-"$numb".sh"  &
      sleep 1
  else
     echo "You have to give ssh or sge"
     exit 1
  fi
done



