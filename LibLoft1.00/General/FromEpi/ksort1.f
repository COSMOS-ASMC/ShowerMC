!  test ksort1 (comb-sort for arbitray data structure with 1 key)
!
!     implicit none
!     integer nrec, l, intv, i, k
!     external keycmp
!     common /cmp/ a
!     parameter (intv=4, nrec =10000)
!     real a(intv, nrec), rand
!
!     do  i=1, intv
!          do  l = 1, nrec
!               a(i, l)=rand(0)
!          enddo
!     enddo
!     do  l=1, 10
!          write(*,*) (a(k,l), k=1, intv)
!     enddo
!     write(*,*) ' --------------------'
!         ****************************************
!     call ksort1(0, nrec, 3, 'a', keycmp)
!         ****************************************
!     do l=1, 10
!              write(*, *) (a(k,l), k=1, intv)
!     enddo 
!     write(*,*) ' ------------------'
!     do l=nrec-10, nrec
!              write(*, *) (a(k,l),k=1, intv)
!     enddo 
!     end
! ***************************************************
!     integer function keycmp(i, j, keyp)
!     implicit none    
!     common /cmp/ a
!     integer intv, nrec
!     parameter (intv=4, nrec =10000)
!     real a(intv, nrec), x
!     integer i, j,  keyp
!     integer k
!
!     if(a(keyp, i) .gt. a(keyp, j)) then
!             do k=1, intv
!                x=a(k,i)
!                a(k,i)=a(k,j)
!                a(k,j)=x
!              enddo
!              keycmp=1
!     elseif(a(keyp, i) .eq. a(keyp,j)) then
!              keycmp=0
!     else
!              keycmp=-1
!     endif
!     end  
!*********************************************************************************
! Sort data with an arbitrary structure using 1 key
!  Usage:
!         call ksort1(noff, nrec, keyp, ad, judge)
!
!     noff:  integer input. noff+1 is the first record position where sort
!            is to be started.
!            In a normal application, give 0. (fisrt record in the data)
!     nrec:  integer input. # of records  in the data to be sorted.
!            Data structure may be arbirary.  nrec records from noff+1-th record
!            will be sorted. 
!     keyp:  keyp-th field in a record is the key to be used for sorting.
!            (keyp should be positive integer ; 1,...)
!       ad: 'a' or 'd' to specify the ascending or descending sort.
!    judge:  integer function name which has the following calling sequence.
!            intval=judge(i, j, keyp)
!            This function should compare the i-th and j-th record of the data
!            (for keyp-th field), and 
!               if the i-th record is > the j-the  record then
!                    1) exchange the i-th and j-th record,
!                    2) give 1 to the function value
!               elseif both are equal
!                    1) give 0 to the function value
!               else (i.e., i-th < j-th) 
!                    1) give -1 to the function value
!
!            Example of how to construct judge is shown in the above
!            test program. (The name 'keycmp' is used; it must be declared
!            external)
!
!************************** note ************************************************
!    If the data is 1 dimensional real or integer array,
!       use kcsr1 or kcsi1.  They are much faster and simpler in usage.
!    If the data is multidimensional real or integer array,
!       the user may use kcsr2 or kcsi2.  They are much simpler in usage.
!       However, the speed is almost comparable with this one. (Sometimes,
!       this is faster).
!********************************************************************************
!  Method:  comb-sort (see Nikkei-byte; Nov. 1991)
      subroutine ksort1(noff, nrec, keyp,  ad, judge)
      implicit none

      integer nrec, keyp, judge, noff
      character*(*) ad

      integer i, j, k, gap, imax
      real  sf/1.30/
      logical exch, more

      gap=nrec
!             
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
                   i=judge(j+noff, k+noff, keyp)
              else
                   i=judge(k+noff, j+noff, keyp)
              endif     
              if(i .eq. 1) exch=.true.
         enddo
         more=exch .or. gap .ne. 1
      enddo
      end










