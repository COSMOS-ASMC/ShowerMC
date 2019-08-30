#!/bin/bash
# to change /usr/local/bin/tcsh at the top of each script into 
# /bin/tcsh
for f in *; do
    if [ $f == "changeUserLocalBinTcsh.sh" ]; then
	echo $f is exception
    else
	x=`grep -c "/usr/local/bin/tcsh" $f`
	if [  $x -gt  0 ]; then
	    sed "s/\/usr\/local//" $f > temp
	    mv $f $f-save
	    mv temp $f
	    chmod +x $f
	fi
    fi
done
echo "if ok, remove *-save"
