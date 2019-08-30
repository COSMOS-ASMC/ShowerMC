#!/bin/bash
if [ $# != 1 ]; then
   cat <<EOF
     Usage: elsepa.sh  path-to-dir
     Here "path-to-dir" is the path to the directory
     where Elsepa data is contained.

     e.g   elsepa.sh  ~/CosmosLoft/Elsepa

     This is to make a link so that $COSMOSTOP/Data/Elsepa
     can be used as if "path-to-dir"
EOF
   exit   
fi
dirpath=$1

if [ ! -d $dirpath ]; then
    # $1 is invalid
    echo Your input:  $dirpath non existent or not directory
    exit
elif [ ! -f $dirpath/Z01/e-/dcs_2p000e02.dat ]; then
    # no data is contained
    echo There is no valid data in $dirpath
    exit
fi

if [ -L $COSMOSTOP/Data/Elsepa ]; then
    if [ -f $COSMOSTOP/Data/Elsepa/Z01/e-/dcs_2p000e02.dat ]; then
	echo A valid link already exists as
	ls -l $COSMOSTOP/Data/Elsepa 
	echo "Do you want to update this link ?"
	echo Enter yes, if so
	read ans
	if [ y$ans == "yyes" ]; then
	    rm -f $COSMOSTOP/Data/Elsepa
	    (cd  $COSMOSTOP/Data/; ln -s $dirpath Elsepa)
	    echo New link
	    ls -l $COSMOSTOP/Data/Elsepa
	    echo established.
	    exit
	else
	    echo So far nothing changed
	    exit
	fi
    else
	echo Elsepa link exists but no valid data inside
	echo "Do you want update the link ?"
	echo Enter yes if so.
	read ans
	if [ y$ans == "yyes" ]; then
	    rm -f $COSMOSTOP/Data/Elsepa
	    (cd  $COSMOSTOP/Data/; ln -s $dirpath Elsepa)
	    echo New link
	    ls -l $COSMOSTOP/Data/Elsepa
	    echo established.
	    exit
	else
	    echo So far nothing changed
	    exit
	fi
    fi
elif [ -d $COSMOSTOP/Data/Elsepa ]; then
    echo Non link directory  \"$COSMOSTOP/Data/Elsepa\"  exists
    echo Be sure if you want to remove this directory
    echo and make a new link.  If so manually remove the
    echo directory and again issue elsepa.sh
    echo So far nothing changed.
    exit
elif [ -f $COSMOSTOP/Data/Elsepa ]; then
    echo \"$COSMOSTOP/Data/Elsepa\" is a file.
    echo Be sure if you want to remove the file
    echo and make a new link.  If so manually remove the
    echo file and again issue elsepa.sh
    echo So far nothing changed.
    exit
else
    (cd $COSMOSTOP/Data;  ln -s $dirpath Elsepa)
    echo New link
    ls -l $COSMOSTOP/Data/Elsepa
    echo established.
fi
exit

    
	
   
    
