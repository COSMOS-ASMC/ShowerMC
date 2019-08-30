#!/bin/csh -f
#  this is to change "DisPar" in the files listed in DisParaFiles.
#  You may use this when DisPara is copied to other directory under UserHook.
while (1)
  echo "Enter the dir name of this ancestor"
  set oldname=$<
  echo "It is " $oldname ": enter y if ok"
  set yesno=$<
  if( x$yesno == "xy" ) break
end


foreach g(`cat DisParaFiles`)
#    echo $g
#            g  FleshBasic/setupnv.sh
    set f = $g:t
#            f  setupenv.sh   
    set dir = $g:h
#           dir FleshBasic 
    cp $g  $f.Old
    set x = `pwd`
    set y = $x:t
#    echo $y
    sed s/$oldname/$y/ $g > $f.temp
    mv $f.temp $dir/$f
    set ext = $f:e
    if ( $ext == "sh" ) then
      chmod +x $dir/$f
    endif
end
echo 'You may deletete *.Old files in this directory, if the command succeeded'
echo "You have to use $y  in stead of DisPara in SkeleFlesh "
