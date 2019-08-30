#!/bin/csh 
    if( $#argv <= 1 ) then
       echo "Usage: ./geomstarter.csh dir_of_skelfiles loader [options]"
       echo "options: e.g -wpos 720,480"
       exit
    endif
    set loader = $2;
    cd $1; 
    shift; shift;


    echo  "(normalization world none)" > starter;
    echo  "(bbox-draw world no)"  >> starter;
    echo  "(look-recenter)" >> starter;
    echo  "(transform world world focus rotate -1.5708 0 0 5)" >> starter
    echo  "(sleep-for 6)" >> starter;
    echo  "(look)" >> starter;
    echo  "(sleep-for 2)" >> starter;
    echo  "(load $loader)" >> starter;

    $GEOMVIEW   $EPICSTOP/Util/Work/*.list -c starter   $argv;
