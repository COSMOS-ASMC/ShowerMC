#!/bin/bash
# current atmos model number
modelnum=`(cd $COSMOSTOP/Atmosphere; ls -l AtmosModel | awk 'BEGIN{FS="/"}; {print $2}'| awk 'BEGIN{FS="\:"}; {print $1}')`
exit $modelnum

