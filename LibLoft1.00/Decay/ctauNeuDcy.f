
!     basically   inclusive treatment for neutau
!      B.R       %
!     tau---> mu + neumu + neutau;   17.4    17.4
!              e  + neue  + neutau;  17.8    35.2
!              pic + neutau          11      46.2
!             pi0+pic + neutau       25.5    71.7
!             2pi0 + pic + neutau     9.5    81.2
!     H +  pic + neutau       rest   18.8     100    H ~3pi.
!     tested by Test/testTauDcy.f
!     (adjust br below: also  Makefile there)
      subroutine ctauNeuDcy(pj, a, np)
      implicit none

#include  "Zptcl.h"
#include  "Zmass.h"
#include  "Zcode.h"
      
      integer np     ! output.  no. of produced particles
      type(ptcl)::pj  ! input. tau
      type(ptcl):: a(*)  ! output.  produced ptcls

      integer,intent(out):: brinfo   ! see last part
! 
      integer,parameter:: nbr=6
      real(8),save:: br(nbr)=(/17.4,17.8,11.0,25.5, 9.5,18.8/)
!        use next to check each branch      
!      real(8),save:: br(nbr)=(/0., 0., 0., 0., 0., 100./)
      real(8),save:: cbr(nbr)


      real(8)::u, w
!         polarization is negleced.  
      integer  i, charge, subcode,  icon
      
      logical,save:: first=.true.
      integer,parameter:: jw=0
      integer,save:: branch
      
      real*8 pab
      if( first ) then
         !  make integrated BR.
         cbr(1) = br(1)
    !         cbr(2:nbr)= cbr(1:nbr-1)+ br(2:nbr)
         do i = 2, nbr
            cbr(i)= cbr(i-1)+ br(i)
         enddo
         !  normalize
         cbr(:) = cbr(:)/100.
         cbr(nbr) = 1.0
         first = .false.
      endif
      

!         make neutau
      subcode = -pj%charge
      call cmkptc(kneutau, subcode, 0, a(1))

      charge = pj%charge
      call rndc(u)
      if( u < cbr(1) ) then
         !  tau->mu + neumu + neutau;   17.4 %
         subcode = pj%charge
         call cmkptc(kneumu, subcode, 0, a(2))
         call cmkptc(kmuon, 0, charge,  a(3))
         np = 3
         branch = 1
      elseif( u < cbr(2) ) then
        !       e  + neue  + neutau;  17.8    35.2
         subcode = -pj%charge
         call cmkptc(kneue, subcode, 0, a(2))
         call cmkptc(kelec, 0, charge,  a(3))
         np = 3
         branch = 2
      elseif( u < cbr(3) ) then
!     pic + neutau          11   %
         call cmkptc(kpion, 0, charge,  a(2))
         np = 2
         branch=3
      elseif( u < cbr(4) ) then
!     pi0+pic + neutau       25.5    71.7!
         call cmkptc(kpion, 0, charge,  a(2))
         call cmkptc(kpion, 0, 0,  a(3))
         np = 3
         branch = 4
      elseif( u < cbr(5) ) then
!     2pi0 + pic + neutau     9.5    81.2
         call cmkptc(kpion, 0, charge,  a(2))
         call cmkptc(kpion, 0, 0,  a(3))
         call cmkptc(kpion, 0, 0,  a(4))
         np = 4
         branch = 5
      else
!  assume heavy 1 ptcl with mass 4 pi ; two body
! decay, take only neutau;  ==> modified 
!     a(2)%mass = 4.0*maspic

         call cmkptc(kpion, 0, charge,  a(2))
         call cmkptc(kpion, 0, 0,  a(3))
         call rndc(u)
         if( u < 0.5) then
!               3pi0 + pic + neutau            
            call cmkptc(kpion, 0, 0,  a(4))
            call cmkptc(kpion, 0, 0,  a(5))
         else
!                  pi0 + pic + pi+ + pi- + neutau  
            call cmkptc(kpion, 0, 1,  a(4))
            call cmkptc(kpion, 0, -1,  a(5))
         endif
         np = 5
         branch = 6
      endif

      if (np >= 3 ) then
         call cnbdcy(np, pj%mass, a, jw,  w, icon)
         if( icon /= 0 ) then
            write(0,*) ' icon =',icon, ' from cnbdcy'
            stop
         endif
         !               boost to lab.
         do i = 1, np
            call cibstPol(i, pj, a(i), a(i) )
         enddo
      else
         ! np = 2
         call c2bdcy(pj, a(1), a(2))
! alredy boosted;  
         np = 2 
      endif
      return
!
      entry ctauNeuDcyBr(brinfo)
      brinfo = branch
      end
