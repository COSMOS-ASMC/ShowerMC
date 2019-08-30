#!/bin/bash
if [ $# -ne 1 ]; then
    cat <<EOF
   useage: ./AutoCreateTab.sh  dir
         dir:   ouput directory such as /tmp/$USER/Media
EOF
exit
fi

echo "Be sure to have fixed 'pairvmin' in ../../Mu/MkTab/AutoCreatemuTab"
echo "It may be around 1.0e-4 to 1.0e-3"
echo "Current one is "
awk '$0 ~ "pairvmin" {print $0; exit}' ../../Mu/MkTab/AutoCreatemuTab
echo "you specified output dir: $1" 
echo "May I go; enter yes if so."
read yes
if [ -z "${yes}" ]; then
    exit 1
fi
if [ ${yes} != "yes" ]; then
    exit 1
fi
nohup ./slave.sh  $1  &
