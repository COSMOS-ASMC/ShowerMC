#!/bin/bash
if [ $# -lt 2 ]; then
    cat << "EOF"
    Usage: ./ForManyMedia.sh dir list
  where 
     dir-- is a directory where the sampling table files are to
            be stored.
     list-- is a list of  media names in $EPICSTOP/Data/BaseM/
      such as Pb PWO ...
    E.G  ./ForManyMedia.sh /tmp  W PWO BGO Air

    To make all media in  $EPICSTOP/Data/BaseM/ target,
      ./ForManyMedia.sh /tmp  `(cd $EPICSTOP/Data/BaseM/; ls)`
     would work
EOF
exit
fi
dir=$1
if [ ! -d $dir ]; then
    mkdir -p $dir
fi
shift
for f in $*; do
    echo $f
    ./mkSampMCStab0.sh  $f  $dir 
done

