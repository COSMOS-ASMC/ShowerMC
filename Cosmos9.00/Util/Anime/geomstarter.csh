#!/bin/csh 
    if( $#argv <= 1 ) then
       echo "Usage: ./geomstarter.csh dir_of_skelfiles loader  [options]"
       echo "options: e.g, -wpos 720,480"
       exit
    endif
    set loader = $2;
    cd $1; 
    shift; shift;

    echo  "(normalization world none)" > starter;
    echo  "(bbox-draw world no)"  >> starter;
    head  $loader | awk 'NR>2 && $0 ~ "load" {print;exit}' >> starter;
    awk 'NR>100 && $0 ~ "load" {print;exit}' $loader  >> starter;
    awk 'NR>300 && $0 ~ "load" {print;exit}'   $loader >> starter;
    awk 'NR>700 && $0 ~ "load" {print;exit}'   $loader >> starter;
    awk 'NR>1000 && $0 ~ "load" {print;exit}'   $loader >> starter;
    tail -3  $loader | awk '$0 ~ "load" {print;exit}' >> starter;
    echo  "(look-recenter)" >> starter;
    echo  "(transform world world world rotate -1.5708 0 0)" >> starter;
    echo  "(transform world world  world  rotate 0 0  -1.5708  )" >> starter;
    echo  "(look)" >> starter;
    echo "(sleep-for 5)" >> starter
#    echo  "(load ground)" >> starter
#    echo  "(sleep-for 2)" >> starter
#    echo  "(transform world world world  scale 1.8 1.8 1.8  2)" >> starter;
#    echo  "(sleep-for 2)" >> starter
#    echo  "(transform world world world  translate 0  0  -1000 2)" >> starter;
#    echo  "(sleep-for 2)" >> starter
#    echo  "(transform world world world  translate  -2000 0 0  2 )" >> starter;
#    echo  "(sleep-for 2)" >> starter
#    echo  "(transform world world focus  rotate  0.15  0   0 2 )" >> starter;
#    echo  "(sleep-for 2)" >> starter
#    echo  "(transform world world world  rotate 0  0  -6.28  10)" >> starter;
     echo "(sleep-for 20)" >> starter
     echo "(delete world)" >> starter;
     echo "(load $loader)" >> starter;
#    $GEOMVIEW  -c starter -wpos 720,486  $argv;
    $GEOMVIEW  -c starter   $argv;


