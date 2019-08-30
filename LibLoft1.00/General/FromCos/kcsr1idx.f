!  test kcsr1idx (comb-sort for 1 dimension real double prec. array)
!      
!      implicit none
!      integer nrec, l
!      parameter (nrec =10000)
!      real*8  a(nrec), u
!      integer idx(nrec)
!         do l = 1, nrec
!              call rndc(u)
!              a(l)= 2*u -1
!         enddo
!         write(*,*) '============ original'
!          do l=1, 10
!               write(*,*) a(l)
!          enddo
! 
!          call kcsr1idx(a,  nrec, idx,  'd')
!          write(*,*) '============ descending  first 10'
!          do l=1, 10
!               write(*, *) a(idx(l))
!          enddo 
!          write(*,*) '   last 10'
!          do  l=nrec-10, nrec
!               write(*, *) a(idx(l))
!          enddo
!          
!          call kcsr1idx(a, nrec, idx, 'a')
!          write(*,*) '============ ascending  first 10'
!          do l=1, 10
!               write(*, *) a(idx(l))
!          enddo 
!          write(*,*) '   last 10'
!          do  l=nrec-10, nrec
!               write(*, *) a(idx(l))
!          enddo
!      end
 

!--------------------------------------------------------------

      subroutine kcsr1idx(a,  nrec, idx,  ad)
      implicit none

      integer nrec
      real*8  a(nrec)    ! input. unchagned
      integer idx(nrec)  ! outut. a(k) is to be moved idx(k)-th  position.
                         ! ,that is, a(idx(k)) is the k-th largest or
                         !  smallest among  a.
      character*(*) ad

      integer j, k, gap, imax, i
      real*8  x, sf/1.30/
      logical exch, more

      gap=nrec
      more=.true.
      do i = 1, nrec
         idx(i) = i
      enddo

      do while( more )
         gap=float(gap)/sf
         if(gap .le. 0) then
              gap=1
         elseif(gap .eq. 9 .or. gap .eq. 10) then
              gap=11
         endif
         imax=nrec - gap
         exch = .false.
         do j=1, imax
              k=j+gap
              if(ad .eq. 'a') then
                   if(a(idx(j))  .gt. a(idx(k))) then
                                x =idx(j)
                                idx(j) = idx(k)
                                idx(k) = x
                                exch = .true.
                   endif
              else
                   if(a(idx(j)) .lt. a(idx(k))) then
                            x = idx(j)
                            idx(j) = idx(k)
                            idx(k) = x
                            exch = .true.
                   endif
              endif     
         enddo
         more=exch .or. gap .ne. 1
      enddo
      end
