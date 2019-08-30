#include "ZsaveStruc.h"
!c           to test clorep.
!      include 'cgetRotMat4.f'
!      include 'cmkptc.f'
!      include 'cpm2e.f'
!      include 'clorez.f'
!c       ----------------------------
!      implicit none
!      include '../Zptcl.h'
!      include '../Zcode.h'
!      type(ptcl):: p, q, r
!      type(fmom):: gb, gbn
!      real*8  g, gba, bx, by, bz
!      integer i, j
!      g=1.d0
!      do  j=1, 45
!         gba=g*sqrt(1.d0-1.d0/g/g)
!         bx=-sqrt(2.d0)/2.d0
!         by=sqrt(2.d0)/5.d0
!         bz=-sqrt(1.d0 - bx**2 - by**2)
!         gb.p(1)=bx*gba
!         gb.p(2)=by*gba
!         gb.p(3)=bz*gba
!         gb.p(4)=g
!         do i=1, 3
!            gbn.p(i)=-gb.p(i)
!         enddo
!         gbn.t=g     
!         p.fm.p(1)=10.d0
!         p.fm.p(2)=10.d0
!         p.fm.p(3)=10000.d0
!         call cmkptc(knuc, 0, 1, p)
!         call cpm2e(p, p)
!         do i=1, 1
!            call clorep(i, gb, p, q)
!         enddo
!         write(*,*) ' after clorep q=',q.fm.p
!         call clorep(1, gbn, q, r)
!         write(*,*) ' --------g=', g
!         write(*,*)  ( (p.fm.p(i)-r.fm.p(i))/p.fm.p(i), i=1, 4)
!         g = g * 10.d0**.25
!      enddo   
!      end
!       **************************************************************
!       *
!       *    clorep: general Lorentz transformation
!       *           (vector defining axes are paralell in both systems).
!       *
!       **************************************************************
!
! /usage/   call clorep(j, gb, q,  p)
!
!          suppose a system (K') moving with a velocity beta
!          (3-d vector) and gamma factor g relative to another
!          systeme (K).  q is  4 momenta given in
!          the frame K' (of which x,y,z axies are parallel to
!                       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!          those of the system K).  This routine gives
!          ~~~~~~~~~~~~~~~~~~~~~
!          4 momenta, p, seen from the K-system.
!
!          4 moemnta are assumed to be in the order of (px, py,
!          pz, e)
!
!     j: input. integer*4   j=1--->gb are new
!                           j^=1-->gb are the same as previous
!                                  values.
! gb: /fmom/  input.  (g*beta, g).
!  q: /ptcl/  input.   4 momenta and mass of a ptcl
!  p: /ptcl/  output.  transformed 4 one.
!                        (p may be the same one as q)
!
!        If we apply the formula directly for gb*q < 0 case at
!        very high  energy,
!        subtraction mig.p(4) result in complete loss of accuracy.
!        To get rid of this, q is converted to a system whose z-axis
!        coinsides with beta and lorentz transformation is applied
!        for the z-direction.  After that, the vector is ratoted
!        so that z axis be parallel to that of the K-system.
!
!        If beta is completely parallel to the z axis, use
!        clorez(;faster).
!
!   Accuracy is better than 7 dig.p(4)s in normal applicaltion
!      (g upto 10e12).
!
       subroutine clorep(j, gb, q, p)
         implicit none
!----         include '../Zptcl.h'
#include  "Zptcl.h"

         type(fmom):: gb
         type(ptcl):: q, p
!
         real*8 rm(4, 4), rmy(4, 4), rmyi(4, 4), rmz(4, 4),
     *          rmzi(4, 4),  rmi(4, 4), gmin/1.e4/, g
         real*8 fai1, fai2,  tmp,  gbq, a
         type(fmom):: agb
         type(ptcl):: qt

#ifdef   USESAVE
         save agb
#endif

         integer jsv/0/, i, j
         save  rm, rmi, jsv
!
         p = q
         gbq = 0.d0
         do   i=1, 3
            gbq=gbq + gb%p(i)*q%fm%p(i)
         enddo
!
         g = gb%p(4) 
         a=1.d0/(1.d0+g)
         if(gbq .ge. 0.d0 .or. g .lt. gmin) then
!         if(gbq .ge. 0.d0 ) then
              do   i=1, 3
                p%fm%p(i) = q%fm%p(i) + 
     *                      gb%p(i)*(q%fm%p(4) + a*gbq)
              enddo
             p%fm%p(4) = g*q%fm%p(4) + gbq
!               j=1, but matrix is not computed
             if(j .eq.1) jsv=0
         else
!              rotate the axes by atan(beta(y)/beta(x)) around z,
!              then rotate the axes by  atan(beta/beta(z))
!              around y, then the orignal z axis coincide with
!              beta.  apply lorentz trans. there and re-rotate
!                 matrix for z-axis
             if(j .eq. 1 .or. jsv .eq. 0) then
                 if(gb%p(2) .eq. 0. .and. gb%p(1) .eq. 0.) then
                     fai1=0.
                 else
                     fai1= atan2(gb%p(2), gb%p(1))
                 endif
                 call cgetRotMat4(3, fai1, rmz)
!                     matrix for y-axis
                 tmp=gb%p(1)**2 + gb%p(2)**2
                 agb%p(3)= sqrt(tmp + gb%p(3)**2)
                 agb%p(4)=g
                 if(tmp .eq. 0. and. gb%p(3) .eq. 0.) then
                    fai2 = 0.
                 else
                    fai2= atan2(sqrt(tmp), gb%p(3))
                 endif
                 call cgetRotMat4(2, fai2, rmy)
!                     combined rotaion matrix
                 call cmultRotMat4(rmy, rmz, rm)
             endif
!                 do combined rotaion
             qt = q
             call capplyRot4(rm, q%fm, qt%fm)
             qt%fm%p(4)=q%fm%p(4)
!           ////////////
!             call ctestOnShell('q before rot', q)
!             call ctestOnShell('qt after rot', qt)
!           ////////////////                      

!                   lorentz trans. along beta
             call clorez(agb, qt,  qt)
!           /////////
!             call ctestOnShell('after lorez', qt)
!           /////////////
!                   re-rotate; get inverse rotation matrix
             if(j .eq. 1 .or. jsv .eq. 0) then
                 call cinvRotMat4(rmz, rmzi)
                 call cinvRotMat4(rmy, rmyi)
                 call cmultRotMat4(rmzi, rmyi, rmi)
             endif
             call capplyRot4(rmi, qt%fm, p%fm)
             p%fm%p(4) = qt%fm%p(4)
!           ////////////
!             call ctestOnShell('after rot', p)
!           ////////////////                      
             jsv=1
          endif
      end
