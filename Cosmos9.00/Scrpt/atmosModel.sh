#!/bin/bash
GetCurrentAtmosModel() {
    model=`(cd $COSMOSTOP/Atmosphere; ls -l AtmosModel | awk 'BEGIN{FS="/"}; {print $NF}')` 
    return
}
ShowAvailableModels() { 
    (cd $COSMOSTOP/Atmosphere/Hidden/; ls | sort)
    return
}

    
GetModelDir(){
    read num
    nummax=`ls $COSMOSTOP/Atmosphere/Hidden/ | wc -l`
    if [ -z $num ]; then
	echo no number
	exit
    fi
    if [ $num -le 0 ] || [ $num -gt $nummax ]; then
	echo invalid number: $num
	exit
    fi
    modeldir=`ls $COSMOSTOP/Atmosphere/Hidden/ | awk 'NR==nth {print;exit}' nth=$num`
    return
}



#----------------- main --------------

if [ "$COSMOSTOP" == "" ]; then
    echo "Env. variable COSMOSTOP is empty"
    exit
fi


GetCurrentAtmosModel
echo " "
if [ -z "$model" ]; then
    echo atmospheric model NOT YET fixed
else
    echo Current atmosphere model is    
fi    

echo "   " "$model"
echo which is one of the following  available models: 
ShowAvailableModels 
echo " "
cat <<EOF
    What do you want to do? 
    Enter a number selecting from the following list
    1) Change the atmospheric model
    2) Show the usage of an atmosphere model 
other) exit
EOF
read num
if [ -z $num ]; then
    exit
fi


if [ $num -eq 1 ]; then
    # change model
    echo To change the model, enter a model number selecting from the available list.
    GetModelDir;
    (cd $COSMOSTOP/Atmosphere/;   rm -f AtmosModel; ln -s Hidden/"$modeldir" AtmosModel)
#    ls $COSMOSTOP/Atmosphere/AtmosModel
    #  remake   library; first remove related component in the library
    $COSMOSTOP/Scrpt/ar-d-atmos.sh
    #  then  make
    (cd $COSMOSTOP/Atmosphere/AtmosModel/; make)
    cat <<EOF
 
   NOTE: You have to re-compile your application; only the library was updated so far.
EOF
elif [ $num -eq 2 ]; then
    #  Show usage of an atmosphere model 
    echo To show the usage of a model, enter a model number from the available list
    GetModelDir
    echo "=============================="
    cat $COSMOSTOP/Atmosphere/Hidden/"$modeldir"/usage
else 
  echo your input $num invalid
fi
exit


