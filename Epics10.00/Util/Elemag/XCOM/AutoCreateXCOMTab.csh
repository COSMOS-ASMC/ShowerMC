#!/bin/csh -f
#
#      to creat cross-section table for 
#    1) photo-electric effect
#    2) coherent scattering
#    3) incoherent scattering (Compton)
#    4) Pair creattion (by nucl. and atomic elec)
#    5) Absoption coeff.
#  at energies 1 keV to 100 GeV
if( $#argv != 2 ) then
	echo "Usage: AutoCreateXCOMTab.csh  basic_media_table_path  output_table_dir"
	exit 1
endif
if( ! -f $1 )  then
	echo $1 not exisst
	exit 1
endif
if( ! -d $2 ) then
	echo directory $2 not exists
	exit 1
endif
set dir_of_media=$1:h
if( x$dir_of_media == x$2 ) then
	echo "The table must be saved in a different"
        echo "directory than the one in which  basic media exists"
	exit 1
endif
set matter=$1:t
set lastname=$1:t.xcom
if( -f $2/$lastname )  then
    echo $2/$lastname already exists
    echo "May I delete it ?: if yes, enter y"
#   set yesno=$< 
    set yesno=y
    if( x$yesno == "xy" ) then
       rm -f $2/$lastname
    else
       exit
    endif
endif
echo "*********************************************************"
echo "*"                                                   
echo "*  Creating a x-section table for low energy photons  "
echo "*  by Berger & Hubbell's XCOM. "
echo "*  The table is typically used at 1keV< E <<1 MeV,  "
echo "*  althoug we create it for 1 keV to 100 GeV "
echo "*  1) coherent scatt."
echo "*  2) incohrent scatt. (Compton)"
echo "*  3) photo-electric effect "
echo "*  4) Pair creation (by nucl. and atomic elec)"
echo "*  5) total absorption coeff.(with/without) coh.scat" 
echo "* "
echo "*********************************************************"
echo 
echo "1)  make the tab for low energy photons "
echo "2)  skip making tab(later you can make it, if needed)"
echo "Select the number "
#     set maketab=$<
     set maketab=1
if ( $maketab != 1 ) then
  echo "you can make the tab later, if you feel it necessary."
  echo "Use CreateXCOMTab.csh command (in Util/XCOM)"
  echo "The tab name will be H2O.xcom etc where H2O is the"
  echo "basic media name"
  exit
endif
make clean
make

echo
#     put remark at the top; which cross-sections, date, by whom
    set  version=`epicsvn`
    echo "# **********XCOM Tab by Berger & Hubell**************** " >> $2/$lastname
    echo "#     This table  has been made by using $1 " >> $2/$lastname
    echo "#     `date` by $USER using "                   >> $2/$lastname
    echo "#      version $version  "                   >> $2/$lastname
    echo "# "                                   >> $2/$lastname
    echo "#    E        Coh.        InCoh.    P.E      N.Pair    E.Pair    Attn Coeff."  >> $2/$lastname
    echo "#   GeV      1/(g/cm2)     //        //        //        //     with coh without coh" >> $2/$lastname
#    next -------- is mandatory
    echo "#-------------------------------------------"  >> $2/$lastname
source $EPICSTOP/Scrpt/setarch
if( -f epXCOM$ARCH ) then
    setenv MATTER $matter
    setenv BASEPATH  $1
    ./epXCOM$ARCH >> $2/$lastname
else	
    echo "executable epXCOM$ARCH not yet made"
    exit 1
endif
awk -f extract.awk $2/$lastname > $2/tempfile
mv $2/tempfile $2/$lastname

