#!/bin/bash
ShowIntModel () {
    qgs=`(cd  $COSMOSTOP/Import; ls -l QGS) | awk '{print  $NF}' | awk -F / '{print $NF}'`
##    (cd  $COSMOSTOP/Import; ls -l QGS) | awk '{printf("%s %s ", $(NF-2), $(NF-1) )}'
    echo  If IntModel in \$PARAM contains \"qgsjet2\", it  means $qgs

    epos=`(cd  $COSMOSTOP/Import; ls -l EPOS) | awk '{print  $NF}' | awk -F / '{print $NF}'`
#    (cd  $COSMOSTOP/Import; ls -l EPOS) | awk '{printf("%s %s ", $(NF-2), $(NF-1) )}'
    echo If IntModel in \$PARAM contains \"epos\", it  means $epos
  }



if [ "$COSMOSTOP" == "" ]; then
    echo "Env. variable COSMOSTOP is empty"
    exit
fi

cat <<EOF
  ***************************************************************
  Current model selection condition is as follows:
EOF
ShowIntModel

cat <<EOF
If want to use a different model, select a number
     1) qgsjetII-04 or qgsjetII-03
     2) epos-lhc-v3700 or epos199 or epos-lhc-v3400
 other) exit

  If you need to use  both 1 and 2, use intModel.sh 2 times.
EOF
read ans
if [ -z $ans ]; then
    exit
fi

if [ $ans == 1 ]; then
    $COSMOSTOP/Scrpt/fixQGSII.sh
elif [ $ans == 2 ]; then
    $COSMOSTOP/Scrpt/fixEPOS.sh
else
    exit
fi
cat <<EOF
   *********** NOTE **********
   Even if compilation has been done and successfull,
   the procedure so far only updated the library.
   Your application must also be updated by, say, make clean; make
   in your application folder.
EOF








