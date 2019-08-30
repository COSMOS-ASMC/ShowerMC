#!/bin/bash
if [ -z  "$LIBLOFT" ]; then
    cat <<EOF
    echo \$LIBLOFT is not yet given
    exit
EOF
fi

for f in *;do
    if [ -L $f ]; then
	rm $f
    fi	
done

	 
for f in $LIBLOFT/DataFromEpi/*; do
#	 echo $f; 
         g=`basename $f`;
	 #	 echo $g;
	 ln -s $f $g;
done

#  next is for making link to use old  $COSMOSTOP/Data
#  But At present not employed.  So use as the old style:  
#for f in  $LIBLOFT/DataFromCos/*; do
#    g=`basename $f`;
#    if [ -e $g ]; then
#	echo $g already exist in $LIBLOFT/DataFromEpi/.
#	echo So link name, ${g}-Epi is used.
#	ln -s $LIBLOFT/DataFromCos/$f  ${g}-Cos
#    else
#	ln -s $LIBLOFT/DataFromCos/$f  $g
#    fi
#done
