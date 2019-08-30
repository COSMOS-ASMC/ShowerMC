!     test kmkDataStm1
!
!     implicit none
!     integer m, j
!     parameter (m=120)
!     real*8 v(m)
!     
!        do j = 1, m
!           v(j) = j
!       enddo

!     call kmkDataStm1(v, m, 'vname', 'f9.4', 9)
!     end


      subroutine kmkDataStm1(v, m,  name, fmt, w)
      implicit none
!        make data statement for a given array data values
!        so that one can use it in fortran90  program
!         (wiht f90 extension).
!
      integer,intent(in):: m    !  input.  see below
      real(8),intent(in):: v(m)  !  input. given data
      character*(*) name    ! input. character str. for array name
      character*(*) fmt    ! input.  basic e.g 'f8.3'
      integer,intent(in):: w            ! input.  width of fmt  e.g 8  for the above case
!           output is on standard output.
!
      integer:: nInRow, nInBlock, colwidth, contLines
      integer::  k1, k2, j, l, cc, k
      parameter (colwidth = 100, contLines = 19 - 2)
      character*40  form
      character*(contLines) contSymb/'123456789abcdefgh'/
      character*3 cm /"','"/
      character*150 cache
      integer klena

      nInRow = colwidth/(w+1)             ! # of data in one row 
      nInBlock = contLines * nInRow       ! # of data in one DATA stm.

      write(form, '( "(", a, ",", a,  ")")'  )  fmt, cm

         do j = 1, m, nInBlock
            k1 = j
            k2 = min(k1 + nInBlock -1, m)
            
            write(*, '(a,a,a,i3,a,i3,a)') &
                 ' data ( ',name, '(i), i=', k1, ', ',k2, ')/ &'

            cc = 0
            do k = k1, k2, nInRow
               cc = cc + 1
               if( k + nInRow -1 .lt. k2) then
                  write(*, form)  (v(l), l = k, min(k + nInRow-1, k2) ), ", &"
               else
                  write(cache, form)  &
                       (v(l), l = k, min(k + nInRow-1, k2) )
                  l = klena(cache)
                  cache(l:l)=' '        ! delete last ,
                  write(*,'(a,"/")') cache
               endif
            enddo
         enddo
       end subroutine kmkDataStm1



               

