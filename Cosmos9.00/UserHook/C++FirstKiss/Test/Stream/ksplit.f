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
      write(*, *) ' m,n=', m, n
      write(*, *) ' tl=', tl
      write(*,*) 'b1', b(1)(1:m)
      write(*,*) 'b2', b(2)(1:m)
      write(*,*) 'b3', b(3)(1:m)
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
            if(nc .le. m) then
               b(nr)(nc:nc) = a(i:i)
               write(*,*) ' b=', b(nr)(nc:nc)
            endif
         endif
      enddo
      end
