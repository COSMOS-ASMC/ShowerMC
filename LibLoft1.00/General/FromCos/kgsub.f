      subroutine kgsub(r, t, si, so)
      implicit none
!        replace string r in si by t and put into so
!        old kgsub has bug if t is "" or " " ...
!    moved from Epics/prog/KKlib  to Cosmos/KKlib/ 2014/Sep/27
      character*(*)  r  ! input.
      character*(*)  t  ! input.
      character*(*)  si ! input.
      character*(*)  so ! output.

      integer klena, ir, it, isi, iso, i, j

      ir = len(trim(r))
      it = len(trim(t))
      if(it == 0 ) it = len(t)

      isi =len(trim(si))
      iso = len(so)
      
      i = 1
      j = 0
      so = ' '
      do while ( i .le. isi )
         if(i+ir-1 .le. isi) then
            if( si(i:i+ir-1) .eq. r(1:ir)) then
               j = j + 1
               if(it > 0 )  then
                  so(j:j+it-1) = t(1:it)
                  j = j + it -1
               endif
               i = i + ir
            else
               j = j + 1
               if( it> 0 ) then
                  so(j:j+it-1) = si(i:i)
               else
                  so(j+it-1:J) = si(i:i)
               endif
               i = i + 1
            endif
         else
            j = j + 1 
            if(it > 0 ) then
               so(j:j+it-1) = si(i:i)
            else
               so(j+it-1:j) = si(i:i)
            endif
            i = i +1
         endif
         if(j .gt. iso) then
            write(0,*)
     *      ' output string length is too short: kgsub'
            stop 9999
         endif
      enddo
      end
