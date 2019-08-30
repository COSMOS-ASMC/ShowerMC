#!/bin/csh -f
#  remove last ,\ in the gnuplot script 'var.gp' made by
#  splithisto.awk 
#
foreach f(*)
    if( -d $f ) then
      ( cd $f; $0 )
    endif
    if( $f == "var.gp") then
	mv var.gp temp.gp
        sed '$,$s/,\\//' temp.gp >  var.gp
	rm -f temp.gp
    endif
end
