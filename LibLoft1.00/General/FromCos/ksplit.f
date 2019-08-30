      subroutine ksplit(a,   m, n, b,  nr)
      implicit none
      character*(*)  a   ! input. e.g   abc xyz xxkxkxkk
      integer m          ! input.  max characters containable in b(i)
      integer n          ! input.  array size of b    
      character*(*) b(n) ! output.  field  value of i=1,2,...nr
      integer nr         ! output.  nr field is obtained

      integer klena, tl, bp, i, nc
      character*1 prev

      tl = klena(a)

      nr = 0
      prev = ' '
      do i = 1, tl
         if( a(i:i) .eq. ' ' ) then
            prev = ' '
         else
            if(prev .eq. ' ') then
               if(nr .ge. n) return  !!!!
               nr = nr + 1
               nc = 0
               prev = 'x'
            endif
            nc = nc + 1
            if(nc .eq. 1)  b(nr) = ' '
            if(nc .le. m) then
               b(nr)(nc:nc) = a(i:i)
            endif
         endif
      enddo
      end
