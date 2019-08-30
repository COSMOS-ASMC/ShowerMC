#!/bin/sh
#  manipulate histograms
source ../setupenv.sh
DIR=$ARENA
export DIR
source ${COSMOSTOP}/UserHook/${DIR}/Smash/setupenv.sh $0
source ${COSMOSTOP}/UserHook/${DIR}/FleshHist/setupenv.sh $0


echo " execid =  $EXECID"
echo " fleshdir= $FLESHDIR"


source ./setupManipHistEnv.sh
source $COSMOSTOP/Scrpt/setarch.sh

here=`pwd`
cd ../../Hist
for job in $JOBTYPE
do
    case $job in
	1) 
	     make -f bin2bin.mk
	     echo ".hist is being converted to .chist"
	    ./bin2bin$ARCH
	     echo "------end of bin to bin-----"
	    ;;
	2) 
	     make -f bin2ascii.mk
	     echo "bin is being converted to ascii"
	    ./bin2ascii$ARCH > $AHISTFILE
	    echo "--------end of bin to ascii----"
	    ;;
	    
	3|5|7)
	     make -f bin2ascii.mk
	    export -n HYBFILE0
	    echo ".chist is being converted to .achist"
	    echo "HYBFILE0 is not needed so it is unexported"
	    ./bin2ascii$ARCH > $AHISTFILE
	    echo "----end of .chist to .achist conversion----"
	    ;;
	4|6)
	    echo "making files for plotting"
	    splitHisto.sh $PLOTDIR $AHISTFILE
	    echo "====end of plotting file creation===="
	    ;;
    esac
done
if [ -d  ../${DIR}/Assemble/${EXECID}-Plot ]; then
    rm -fr ../${DIR}/Assemble/${EXECID}-Plot
fi
mv ${EXECID}*  ../${DIR}/Assemble/


