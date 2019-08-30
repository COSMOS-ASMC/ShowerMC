#!/bin/csh -f
set fold = "xxxxxx"
foreach f(`awk '{print $2}' Hosts`)
    echo $f
    if( $fold != $f ) then
      set fold = $f
      ssh $f "ps aux |grep $USER"
    endif
end
