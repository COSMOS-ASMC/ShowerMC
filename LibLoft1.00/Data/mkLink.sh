for f in *; do
    if [ -L $f ]; then
	rm $f
    fi
done

for f in ../DataFromCos/*;do
    g=`basename $f`;
    ln -s ../DataFromCos/$g $g;
done

for f in ../DataFromEpi/*;do
    g=`basename $f`;
    if [ -e $g ]; then
	echo $g aleady exist so we add -epi
	ln -s ../DataFromEpi/$g  ${g}-epi
    else
	ln -s ../DataFromEpi/$f $g;
    fi
done

    
