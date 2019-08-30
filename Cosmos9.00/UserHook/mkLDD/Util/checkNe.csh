#!/bin/csh
@ n=0
echo " This is check routine for E0=10^20 eV"-
echo " For other E0, you may check next:"
echo " max Ne number is far from 6x10^10 x (E0/10^20 eV)"
echo " for  many showers, the results are suspcted invalid"

foreach f(*.hyb)
    echo $f
    awk -f $COSMOSTOP/UserHook/mkLDD/Util/checkNe.awk $f
    @ n++
    if( $n == 5  ) then
       set dummy=$<
       @ n=0
    endif
end
