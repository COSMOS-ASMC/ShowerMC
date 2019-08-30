#!/bin/csh
#   1:  start file number
#   2:  end
#   3:  overlap frame
#   4:  snapshot file number offset 
if( $#argv  != 4) then
    echo  'front.csh startfile# endfile# overlapframe# offset_for_snapshot#'
    echo  'Say,  front.csh 1 850 2 1000'
    exit
endif
@ n = $1;
@ i = 0;
while ( $i < $3)
      @ fn = $i + $n; 
      awk 'END {printf("(load ts%5.5d.skel)\n", n)}' n=$fn /dev/null
      @ i++;
end
#   1   2  3

@ i = 1;
while ( $i <= $2 )  
   @ n = $i + $3;
   @ m = $n + $4;
   awk 'END {printf("(delete ts%5.5d.skel)\n", n)}' n=$i /dev/null;
   awk 'END {printf("(load ts%5.5d.skel)\n", n)}' n=$n /dev/null;
   awk 'END {printf("(snapshot focus  ts%5.5d.ppm ppmscreen)\n",m)}' m=$m /dev/null;
#   awk 'END {print "(sleep-for 0.01)"}' /dev/null;
   @ i++;
   @ n++; 
end
#      2  3  4


    

