#!/bin/csh -f 
if ( $#argv < 4 ) then
    echo "./gatherAllInTmp.csh hostfile fromdir todir {fileextension}"
    echo "all files or files with specified file extension in fromdir"
    echo "will be copied (not moved) to todir."
    echo "fileextension is such as .hyb. If omitted, all files are target"
    exit
endif
set fold="xx"
@ n = 0
set numb=`awk '$1 != "#" && NF > 0 {print $1}' $1`
set host=`awk '$1 != "#" && NF > 0 {print $2}' $1`
foreach f($numb)
@ n++
    echo "$f, $host[$n] is being inspected"
    if ( $host[$n] == `hostname -s`) then
	echo "this is the current host; skip inspection"
    else
	if ( "x$4" == "x"  ) then
	    scp ${host[$n]}:$2/"*" $3/;
	else
	    scp ${host[$n]}:$2/"*"$4 $3/;
	endif
    endif
end
