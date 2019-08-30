!     test kmkDataStm2a
!
!     implicit none
!     integer m, n, i, j
!     parameter (m=120, n=2)
!     real*8 v(m, n)
!     
!     do i = 1, n
!        do j = 1, m
!           v(j, i) =(i-1) *m + j
!        enddo
!        write(*, *) (v(j, i),j =1, m)
!     enddo
!     call kmkDataStm2a(v, m, n, 'vname', 'g15.3', 15)
!     end


      subroutine kmkDataStm2a(v, m, n, name, fmt, w)
      implicit none
!        make data statement for a given array data values
!        so that one can use it in fortran program.
!
      integer m, n    !  input.  see below
      real*8 v(m, n)  !  input. given data
      character*(*) name    ! input. character str. for array name
      character*(*) fmt    ! input.  basic e.g 'f8.3'
      integer w            ! input.  width of fmt  e.g 8  for the above case
!           output is on standard output.
!
      integer nInRow, nInBlock, colwidth, contLines
      integer  k1, k2, i, j, l, cc, k
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
!cc      write(*, *) form
      do i = 1, n

         do j = 1, m, nInBlock
            k1 = j
            k2 = min(k1 + nInBlock -1, m)
            
            write(*, *) '          data ( ',
     *      name, '(i, ',i,'), i=', k1, ', ',k2, ')/'
            cc = 0

            do k = k1, k2, nInRow
               cc = cc + 1
               if( k + nInRow -1 .lt. k2) then
                  write(*, form) contSymb(cc:cc),
     *           (v(l, i), l = k, min(k + nInRow-1, k2) )
               else
                  write(cache, form) contSymb(cc:cc),
     *           (v(l, i), l = k, min(k + nInRow-1, k2) )
                  l = klena(cache)
                  cache(l:l)=' '        ! delete last ,
                  write(*,'(a)') cache
               endif
            enddo
            write(*, '(a)') '     * /   '

         enddo

      enddo
      end


               

