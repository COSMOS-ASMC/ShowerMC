function confirm {
#  first we test $1 is /tmp... or not (because current host may not have /tmp/$USER then
#  if [ ! -d $1 ]; will succeed )

x=` echo  $1 | awk '$0 ~ "^/tmp"{print "tmp"}' `
if [ "x$x" != "xtmp" ] ; then
    if [ ! -d $1 ] ; then
#           create it
	if [ -f $1 ]; then
	   echo "$1 is file ; not directory "
        else
   	  mkdir -p $1	
        fi
    else
      some="`ls $1/`"
      if [ -n "$some" ] ; then
	echo "some files in $1"
        num=0
        while [ $num -le 0 ] || [ $num -gt 4 ] 
       	do
    	   echo "1--remove all in $1 (normal)"
	   echo "2--remove files with specified file extensions"
           echo "3--keep all"
           echo "4--keep all and quit" 
           echo "Select number"
	   read num
           test -z $num && num=0
       done	
       if [ $num -eq 1 ] ; then
#                      rm $1/*  dose no work here (why?)
	   (cd $1;rm -fr `ls`)
       elif [ $num -eq 2 ]; then
           echo "Enter file extensions like .hyb .dat"
	   read extens
           ../rmOnlySpecExt.csh $1  "$extens"
       elif [ $num -eq 3 ]; then
          :
       else
          exit
       fi
    fi
  fi		
else
#           /tmp/
    echo '================================'
    echo "$1 seems to be in each host. "
    echo "If some files remain there, you may have "
    echo "difficulty when gathering  files you will create "
    num=0
    while [ $num -le 0 ] || [ $num -gt 3 ]
    do
      echo "1--If any in $1, delete all. "
      echo "   Or If there is no $1, create it (normal)" 
      echo "2--Keep all in $1;  you should be sure it exists in all hosts."
      echo "3--Keep all in $1 and quit"
      echo "Select number"
      read num
      test -z $num && num=0
    done
    if [ $num -eq 1 ]; then
	echo "In the later fleshing step, will you use ssh for job submission ?"
	echo "Enter y, if so"
	read yesno
	if [ x$yesno = "xy" ]; then
#              inspect only specified hosts  
	  source ../Smash/setupenv.sh $0
	  if [ $MCPU = $NCPU ]; then
	      ../rmAllInTmp.csh ../Hosts $1
          else
    	      ../rmAllInTmp.csh ../ThinHosts $1
	  fi
	else		
   	  ../rmAllInTmp.csh ../allHosts $1
        fi
    elif [ $num -eq 2 ]; then
      :
    else		
      exit
    fi
fi
}
