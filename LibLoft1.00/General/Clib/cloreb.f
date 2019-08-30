#include "ZsaveStruc.h"
!c              to test cloreb
!        include 'clorez.f'
!        include 'cgetRotMat4.f'
!        include 'clorep.f'
!  -----------------------------
!      implicit none
!      include '../Zptcl.h'
!      type(fmom):: p, q, r, s, gb, gbn
!      real*8  rm(4,4), rmz(4,4), rmy(4,4), g, bx, by, bz
!      real*8  fai1, fai2, tmp, pabs, gba
!      real*8  am/.938/ 
!      integer i, j
! 
!      g=1.1
!      do j=1, 45
!         gba=g*sqrt(1.d0-1.d0/g/g)
!         bx=-sqrt(2.d0)/2.
!         by=sqrt(2.d0)/5.
!         bz=sqrt(1.d0 - bx**2 - by**2)
!         gb.p(1) = bx*gba
!         gb.p(2) = by*gba
!         gb.p(3) = - bz*gba
!         gb.p(4) = g
!         gbn.x = - gb.p(1)
!         gbn.y = - gb.p(2)
!         gbn.z = - gb.p(3)
!         gbn.t= gb.p(4)
!c                 matrix for z-axis
!         fai1=atan2(gb.p(2), gb.p(1))
!         call cgetRotMat4(3, fai1, rmz)
!         tmp=gb.p(1)**2 + gb.p(2)**2
!c               matrix for y-axis
!         fai2=atan2(sqrt(tmp), gb.p(3))
!         call cgetRotMat4(2, fai2, rmy)
!c                 combined rotaion matrix
!         call cmultRotMat4(rmy, rmz, rm)
!         p.px=1.
!         p.py=10.
!         p.pz=10000.
!         call cpxyzp(p, pabs)
!         p.e=sqrt( pabs**2 + am**2)
!         call cloreb(1, gb, p, q)
!c
!         call clorep(1, gbn, q, r)
!c                 do rotaion
!         call capplyRot4(rm, r.p, s.p)
!         s.e = r.e
!         write(*, *) " ------ gamma=", g
!         write(*,*) ( (p.p(i)-s.p(i))/p.p(i), i=1, 4)
!c        write(*,*)   (p.p(i),s.p(i), i=1, 4)
!         g = g * 10.d0**.25
!      enddo   
!      end
!       **************************************************************
!       *
!       *    cloreb: General Lorentz transformation.
!       *            The z axis is on the beta direction.
!       *
!       **************************************************************
!
! /usage/   call cloreb(i, gb, q,  p)
!
!          suppose a system (K') moving with a velocity beta
!          (3-d vector) and gamma factor g relative to another
!          system (K).  q is  4 momentum given in
!          K' (of which the z axis is parallel to beta.)
!                       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!          This routine converts q into p seen from K.
!
!          4 moemnta are assumed to be in the order of (px, py,
!          pz, e)
!
!     i: input. integer*4.  if 1, gb is assumed to be new
!                           else they are the same as previous call.
!    gb: input. type fmom (g*beta, g)
!     q: input. type ptcl. 4 momentum and mass of a particle.
!     p: output. type ptcl transformed one ( p may be q)
!
!
      subroutine cloreb(i, gb, q, p)
         implicit none
#include  "Zptcl.h"
         type(fmom):: gb
         type(ptcl)::  q, p
         integer i

         type(ptcl):: qq
         type(fmom):: a, b, beta
         real*8 betanorm, anorm

#ifdef USESAVE
         save a, b, beta, betanorm
#endif
  
         if(i .eq. 1) then
#ifdef NEXT486
            betanorm = gb%x**2 + gb%y**2 + gb%z**2
            if(betanorm .ne. 0.) then
               betanorm = sqrt(betanorm)
               beta%x = gb%x/betanorm
               beta%y = gb%y/betanorm
               beta%z = gb%z/betanorm
               a%x= beta%y
               a%y = -beta%x
               a%z= 0
               anorm = a%x**2 + a%y**2
               if(anorm .eq. 0.) then
                  a%x = 1.
                  a%y = 0.
                  a%z = 0.
               else
                  anorm = sqrt(anorm)
                  a%x = a%x/anorm
                  a%y = a%y/anorm
               endif
!
               b%x = beta%y*a%z - beta%z * a%y
               b%y = beta%z*a%x - beta%x* a%z
               b%z = beta%x*a%y - beta%y* a%x
            endif
         endif
         if(betanorm .eq. 0.) then
            p = q
            if(gb%z .lt. 0.) then
               p%fm%z = - p%fm%z
               p%fm%y = - p%fm%y
               p%fm%x = - p%fm%x
            endif
         else
!            assume that  qx is  a-direction and qy is b direction

               qq = q
            qq%fm%x = q%fm%x*a%x +  q%fm%y*b%x +
     *         q%fm%z*beta%x
            qq%fm%y = q%fm%x*a%y +  q%fm%y*b%y +
     *         q%fm%z*beta%y
            qq%fm%z = q%fm%x*a%z +  q%fm%y*b%z +
     *         q%fm%z*beta%z
            qq%fm%e = q%fm%e
!              Lorentz boost by gb
            call clorep(i, gb, qq, p)
         endif
#else
            betanorm = gb%p(1)**2 + gb%p(2)**2 + gb%p(3)**2 
            if(betanorm .ne. 0.) then
               betanorm = sqrt(betanorm)
               beta%p(1) = gb%p(1)/betanorm
               beta%p(2) = gb%p(2)/betanorm
               beta%p(3) = gb%p(3)/betanorm
               a%p(1)= beta%p(2)
               a%p(2) = -beta%p(1)
               a%p(3)= 0
               anorm = a%p(1)**2 + a%p(2)**2
               if(anorm .eq. 0.) then
                  a%p(1) = 1.
                  a%p(2) = 0.
                  a%p(3) = 0.
               else
                  anorm = sqrt(anorm)
                  a%p(1) = a%p(1)/anorm
                  a%p(2) = a%p(2)/anorm
               endif
!
               b%p(1) = beta%p(2)*a%p(3) - beta%p(3) * a%p(2)
               b%p(2) = beta%p(3)*a%p(1) - beta%p(1)* a%p(3)
               b%p(3) = beta%p(1)*a%p(2) - beta%p(2)* a%p(1)
            endif
         endif
         if(betanorm .eq. 0.) then
            p = q
            if(gb%p(3) .lt. 0.) then
               p%fm%p(3) = - p%fm%p(3)
               p%fm%p(2) = - p%fm%p(2)
               p%fm%p(1) = - p%fm%p(1)
            endif
         else
!            assume that  qx is  a-direction and qy is b direction
!
            qq = q
            qq%fm%p(1) = q%fm%p(1)*a%p(1) +  q%fm%p(2)*b%p(1) +
     *           q%fm%p(3)*beta%p(1)
            qq%fm%p(2) = q%fm%p(1)*a%p(2) +  q%fm%p(2)*b%p(2) + 
     *           q%fm%p(3)*beta%p(2)
            qq%fm%p(3) = q%fm%p(1)*a%p(3) +  q%fm%p(2)*b%p(3) +
     *           q%fm%p(3)*beta%p(3)
            qq%fm%p(4) = q%fm%p(4)
!              Lorentz boost by gb
            call clorep(i, gb, qq, p)
         endif
#endif
      end





