#!/bin/csh
#   1:  start file number
#   2:  end
#   3:  frames to be shown togther for front display
#   4:  snapshot file number offset 
#   5:  track or front only 
if( $#argv  != 5) then
    echo  'timedev.csh startfile# endfile# overlapframe# offset_for_snapshot# track(0)_or_front(1)_disply'
    echo  'Say,  timedev.csh 1 850 2 1000 0'
    exit
endif
@ n = $1;
@ i = 0;
while ( $i < $3)
      @ fn = $i + $n; 
      awk 'END {printf("(load ts%5.5d.skel)\n", n)}' n=$fn /dev/null
      if($5 == 0) then
        @ i = $3; 
      else	
        @ i++;
      endif
end
#   1   2  3

@ i = $1;
while ( $i <= $2 )  
   @ n = $i + $3;
   @ m = $n + $4;
   if($5 == 1)  awk 'END {printf("(delete ts%5.5d.skel)\n", n)}' n=$i /dev/null;
   awk 'END {printf("(load ts%5.5d.skel)\n", n)}' n=$n /dev/null;
   awk 'END {printf("(snapshot focus  ts%5.5d.ppm ppmscreen)\n",m)}' m=$m /dev/null;
#   awk 'END {print "(sleep-for 0.01)"}' /dev/null;
   @ i++;
   @ n++; 
end
#      2  3  4


    

