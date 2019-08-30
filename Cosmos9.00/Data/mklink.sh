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

	 
for f in $LIBLOFT/DataFromCos/*; do
#	 echo $f; 
         g=`basename $f`;
	 #	 echo $g;
	 ln -s $f $g;
done

#   next is to use $EPICSTOP/Data as $COSMOSTOP/Data
#  but not employed. 
#for f in  $LIBLOFT/DataFromEpi/*; do
#    g=`basename $f`;
#    if [ -e $g ]; then
#	echo $g already exist in $LIBLOFT/DataFromCos/.
#	echo So link name, ${g}-Epi is used.
#	ln -s $LIBLOFT/DataFromEpi/$f  ${g}-Epi
#    else
#	ln -s $LIBLOFT/DataFromEpi/$f  $g
#    fi
#done
