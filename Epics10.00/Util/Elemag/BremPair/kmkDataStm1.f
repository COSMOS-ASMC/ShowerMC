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
!        so that one can use it in fortran program.
!
      integer m    !  input.  see below
      real*8 v(m)  !  input. given data
      character*(*) name    ! input. character str. for array name
      character*(*) fmt    ! input.  basic e.g 'f8.3'
      integer w            ! input.  width of fmt  e.g 8  for the above case
!           output is on standard output.
!
      integer nInRow, nInBlock, colwidth, contLines
      integer  k1, k2, j, l, cc, k
      parameter (colwidth = 64, contLines = 19 - 2)
      character*40  form
      character*(contLines) contSymb/'123456789abcdefgh'/
      character*3 cm /"','"/
      character*72 cache
      integer klena

      nInRow = colwidth/(w+1)             ! # of data in one row 
      nInBlock = contLines * nInRow       ! # of data in one DATA stm.

      write(form, '( "(5x,a1,1x,", i2, "(", a, ",", a,
     *                    "))")') nInRow, fmt, cm

         do j = 1, m, nInBlock
            k1 = j
            k2 = min(k1 + nInBlock -1, m)
            
            write(*, *) '          data ( ',
     *      name, '(i), i=', k1, ', ',k2, ')/'
            cc = 0

            do k = k1, k2, nInRow
               cc = cc + 1
               if( k + nInRow -1 .lt. k2) then
                  write(*, form) contSymb(cc:cc),
     *           (v(l), l = k, min(k + nInRow-1, k2) )
               else
                  write(cache, form) contSymb(cc:cc),
     *           (v(l), l = k, min(k + nInRow-1, k2) )
                  l = klena(cache)
                  cache(l:l)=' '        ! delete last ,
                  write(*,'(a)') cache
               endif
            enddo
            write(*, '(a)') '     * /   '

         enddo

      end


               

