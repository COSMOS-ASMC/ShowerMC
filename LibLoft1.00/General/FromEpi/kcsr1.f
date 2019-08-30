!  test kcsr1 (comb-sort for 1 dimension real double prec. array)
!      
!     implicit none
!     integer nrec, l
!     parameter (nrec =10000)
!     real*8  a(nrec), u
!         do l = 1, nrec
!              call rndc(u)
!              a(l)= 2*u -1
!         enddo
!         write(*,*) '============ original'
!         do l=1, 10
!              write(*,*) a(l)
!         enddo
!
!         call kcsr1(a,  nrec,  'd')
!         write(*,*) '============ descending  first 10'
!         do l=1, 10
!              write(*, *) a(l)
!         enddo 
!         write(*,*) '   last 10'
!         do  l=nrec-10, nrec
!              write(*, *) a(l)
!         enddo
!         
!         call kcsr1(a, nrec, 'a')
!         write(*,*) '============ ascending  first 10'
!         do l=1, 10
!              write(*, *) a(l)
!         enddo 
!         write(*,*) '   last 10'
!         do  l=nrec-10, nrec
!              write(*, *) a(l)
!         enddo
!     end
!
!--------------------------------------------------------------
!        Sort a 1 dimensional  real*8 array
!        see for details kcsr2.f
!
!        call kcsr1(a,  nrec, d)
!        
!        a:  real*8 array to be sorted.
!     nrec:  # of records in a.  
!       ad: 'a' or 'd' to specify the ascending or descending sort.
!
      subroutine kcsr1(a,  nrec, ad)
      implicit none

      integer nrec
      real*8  a(nrec)
      character*(*) ad

      integer j, k, gap, imax
      real*8  x, sf/1.30/
      logical exch, more

      gap=nrec
      more=.true.
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
                   if(a(j)  .gt. a( k)) then
                                x = a(j)
                                a(j) = a(k)
                                a(k) = x
                                exch = .true.
                   endif
              else
                   if(a(j) .lt. a(k)) then
                            x = a( j)
                            a(j) = a(k)
                            a( k) = x
                            exch = .true.
                   endif
              endif     
         enddo
         more=exch .or. gap .ne. 1
      enddo
      end
