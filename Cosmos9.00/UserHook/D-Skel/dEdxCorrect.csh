#!/bin/csh
if( $#argv != 1  ) then
    echo "Usage: dEdxCorrect.csh fdddatDir"
    exit
endif

foreach f($1/*.nrfai)
  set nrfaifilepath=$f
end
echo ".nrfai file path is $nrfaifilepath"
set nrfaifile=$nrfaifilepath:t
echo ".nrfai file is $nrfaifile"
foreach f($1/*.hyb)
  set hybfilepath=$f
end
echo ".hyb file path is $hybfilepath"
set hybfile=$hybfilepath:t
echo ".hyb file  is $hybfile"

    
set lymax=`tail $nrfaifilepath | awk '$1=="dE/dx" {print $3}'`
echo "FDD layer max is  $lymax"

    @ ly=3
    while ( $ly <= $lymax )
	set ElWeb=`awk  -f getsumdEdxOneLayer.awk layer=$ly $nrfaifilepath`
	if( $ElWeb == 0 ) break
	set Elhyb=`awk '$2==ly {print $12}' ly=$ly $hybfilepath`
	set age=`awk '$2==ly {print $5}' ly=$ly $hybfilepath`
#	echo "ElWeb=" $El  
#        echo "Elhyb=" $Elhyb
	if( x$ElWeb != "x" )  then
	    echo $ElWeb $Elhyb | awk '{print age, $1/$2, ly}' age=$age ly=$ly
	endif
	@ ly++
    end
set yesno="n"
while (x$yesno != "xy" ) 
  echo "Enter enhance factor for correction or put 0 to stop here "
  set enhance=$<
  if( $enhance == 0 ) then
	echo  "good bye"
	exit
  endif
  echo "your input is  $enhance ; is it ok ?. Enter y if ok"
  set yesno=$<
end
if ( $enhance == 0 )  then
    echo "good bye"
    exit
endif
echo "Now correction job starts"

if( -f $1/${nrfaifile}.wrong ) then
    echo "old file .nrfai.wrong already exist"
    echo "Therefore I don't copy current .nrfai to .nrfai.wrong"
    echo "assuming .nrfai.wrong is original file"
    echo "Is it OK?  Enter y if yes"
    set yesno=$<
    if( x$yesno != "xy" ) then
      exit
    endif
else
    echo "Current .nrfai is saved as .nrfai.wrong"
    cp $nrfaifilepath $1/${nrfaifile}.wrong
endif

awk -f correctNrfai.awk enhance=$enhance $1/${nrfaifile}.wrong > $1/$nrfaifile


    


