#!/bin/bash
if [ $# -le 1 ]; then
    cat <<EOF
       Usage: ./geomstarter.bash dir_of_skelfiles loader  [options]
       options: e.g, -wpos 720,480
EOF
 exit
fi
loader=$2
pwd=`pwd`
cd $1
shift; shift;

    echo  "(normalization world none)" > starter;
    echo  "(bbox-draw world no)"  >> starter;
    head  $loader | awk 'NR>2 && $0 ~ "load" {print;exit}' >> starter;
    awk 'NR>100 && $0 ~ "load" {print;exit}' $loader  >> starter;
    awk 'NR>300 && $0 ~ "load" {print;exit}'   $loader >> starter;
    awk 'NR>700 && $0 ~ "load" {print;exit}'   $loader >> starter;
    awk 'NR>1000 && $0 ~ "load" {print;exit}'   $loader >> starter;
    awk 'NR>1200 && $0 ~ "load" {print;exit}'   $loader >> starter;
    echo  "(load Array/loadarray)" >> starter
#    echo  "(transform world world world rotate -1.0 0 0)" >> starter;    
    echo  "(look-recenter)" >> starter;
    echo "(sleep-for 10)" >> starter
    echo  "(transform world world world rotate -1.57079632679 0 0)" >> starter;
    echo "(sleep-for 10)" >> starter
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
    echo "(backcolor c0 0 0 0)" >>  starter;
    echo "(sleep-for 40)" >> starter
    echo "(delete world)" >> starter;
    echo "(backcolor c0 0 0 0)" >>  starter;
    echo  "(load Array/loadarray)" >> starter    
    echo "(load $loader)" >> starter;
#    $GEOMVIEW  -c starter -wpos 720,486  $argv;
    $GEOMVIEW  -c starter   $argv;
cd $pwd


