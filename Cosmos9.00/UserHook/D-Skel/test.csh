#!/bin/csh -f 
set numb="1 2 3 # 4" 
foreach f($numb)
    echo $f
    if($f == "#" ) continue
    echo $f
end
