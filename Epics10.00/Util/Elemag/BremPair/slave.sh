#!/bin/bash
# slave of AutoCreateTab.sh
for f in `ls ../../../Data/BaseM/`; do                              
    echo $f;                                                         
     ./AutoCreateTab ../../../Data/BaseM/$f $1
done    
