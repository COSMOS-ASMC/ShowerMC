#!/bin/csh 
#
#  make Glauber preinitialization data for files specified
#  by a regular expression. This uses iniGlauber.
#  
if( $#argv < 3 ) then
    cat <<EOF
    Usage: iniGlauberAll.csh outputdir FeOrPb regular-expression 
           This can be invoked at any place.
    where
   outputdir: directory to store the .GLB and .inp files.
      FeOrPb: = 1  .GLB is for projectile upto Fe
              = 2  .GLB is for projectile upto to Pb 
                   (take much longer time than 1)
 regular-exprssion: such as \$EPICSTOP/Data/BaseM/[A-C]* to 
   specify files starting with A,B or C  in \$EPICSTOP/Data/BaseM/.

 So typical usage is:

    iniGlauberAll.csh /tmp/ 2  \$EPICSTOP/Data/BaseM/[A-B]*  > /tmp/out 2>/tmp/err &

    iniGlauberAll.csh /tmp/ 1  \$EPICSTOP/Data/BaseM/{Al,Argas} >/tmp/out 2>/tmp/err &
    etc
EOF
exit 
endif
(cd $COSMOSTOP/Util/DPM; make clean; make -f mkglauber.mk)

set dir = $1
if( ! -d $dir ) then
    echo $dir non existent
    exit 1
endif
set FeOrPb = $2
echo '# of arg' $#argv
shift 
shift

foreach f($argv)
   echo $f
   $EPICSTOP/Scrpt/iniGlauber $f $dir $FeOrPb no
end

