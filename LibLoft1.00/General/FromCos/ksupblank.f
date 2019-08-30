!      character*120 text
!      text = ' ***  ls 000   0 12   23 s () (   )    |'
!      call ksupblank(text, nc)
!      write(*,*) nc
!      write(*,*) text(1:nc)
!      end


      subroutine ksupblank(text, nc)
      implicit none
      character*(*) text  ! in/out. string
      integer nc          ! output.  resultant string length
!
!         supress two or more blanks into one blank
!      
      integer i, klena, tl 

      i = 1
      do while (i .ne. 0)
         tl = klena(text)
         i = index(text(1:tl), "  ")
         if(i .gt. 0) then
            text(i+1:tl) = text(i+2:tl)
         endif
      enddo
      nc = tl
      end
