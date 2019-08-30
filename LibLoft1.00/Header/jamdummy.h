         subroutine mydummy
#include "jam1.inc"
#include "jam2.inc"
        common/pjpars/mstp(200),parp(200),msti(200),pari(200)
        parameter(itblsz=100)
        parameter(isgmax=29)
        common /sigma1/sigfit(isgmax,itblsz)
        parameter (itblsz2=100)
        parameter (isgmax2=26)
        common /stbl1/sfits1(isgmax2,itblsz2)
        parameter (itblsz3=100)
        parameter (isgmax3=65)
        common /tbls2/sfits2(isgmax3,itblsz3)
        common/hiparnt/hipr1(100),ihpr2(50),hint1(100),ihnt2(50)
        common/bveg1/xl(10),xu(10),acc,ndim,ncall,itmx,nprn,iprn,ierr
!        common/vnibins/ihist(4),indx(1000),bin(20000)
!	write(0,*) ihist
!	write(0,*) '--------------------'
         write(0,*) '-------------------xu'
         write(0,*) xu
         write(0,*) '-------------------ihpr2'
        write(0,*) ihpr2
         write(0,*) '-------------------hipr1'
        write(0,*) hipr1
         write(0,*) '-----------------sfits2'
        write(0,*) (sfits2(1,i),i=1,itblsz3)
         write(0,*) '-----------------sigfit'
         write(0,*) (sigfit(1,i),i=1,itblsz)
         write(0,*) '-----------------parp'
         write(0,*) parp
         write(0,*) '-----------------mstc'
         write(0,*) mstc
!////////////
!	  call  jamdummy2('jamdummy')
!////////////
         end

!////////////////////////
           subroutine jamdummy2(id)
           implicit none
              character*(*) id
 	integer i

         integer n, npad, k
         real*8  p, v
         common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
         integer j
              write(*,*) '************* k @', trim(id)
         do i= 1, 1000
            do j = 1, 5
             if(k(i,j) /= 0) then
                write(*,*) i, j, k(i,j)
             endif
            enddo
         enddo
        end
