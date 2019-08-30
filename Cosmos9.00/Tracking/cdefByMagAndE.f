      module modBEfield
      real(8),save:: Bfld(3)
      real(8),save:: Efld(3)
      real(8),save:: eval=0.29979   ! e:  Ze*E => in GeV/c/s
                        !  Ze vB  ==> in GeV/c/s
      real(8),save:: ptclmass     ! mc^2 in GeV
      real(8),save:: ptclZ        ! charge
      end       module modBEfield
#include "Zcondc.h"
      subroutine cgetEfield(aTrack)
      use modBEfield
      use modEfield
      implicit none
#include "Zglobalc.h"
#include "Ztrack.h"
! #include "Zmagfield.h"      
#include "Ztrackv.h"
#include "Zmass.h"

      type(track)::aTrack ! input. current  track
      if(HowEfield == 1 ) then
         call cEfield0(aTrack, Efld)
      else
#if defined (MYEFIELD) 
!        not simple Efield; user must define next
         call cmyEfield(aTrack, Efld)
#else
         write(0,*) "HowEfield=",HowEfield," invalid"
         write(0,*) "You must define MYEFIELD  in Zcondc%h"
         write(0,*)
     *    "and prepare cmyEfield subroutine in your apps.",
     *    "Interface is cmyEfield(aTrack,Efout); see",
     *    "manual or cEfield0 in Cosmos/Module"
         stop
#endif
      endif
      end subroutine  cgetEfield

      subroutine cdefByMagAndE(aTrack,  leng,  dispmr, dispmd,
     *  newmom)
      use modBEfield
      use modEfield
      implicit none
#include "Zglobalc.h"
#include "Ztrack.h"
! #include "Zmagfield.h"
#include "Ztrackv.h"
#include "Zmass.h"

      type(track)::aTrack ! input. current  track
      real(8),intent(in):: leng !  moved length in m
      type(coord)::dispmr !  output displacement 
      type(coord)::dispmd !  output new direction cos
      real(8),intent(out):: newmom(3) ! chnaged momentum

      integer,parameter:: nsim=6
      integer,parameter::nstep=1  ! if 2, accuracy may increse
      real(8):: t(0:nstep), u(6,0:nstep)
!xxxxxxxxxx next equivalence does not work (ifort)
!      real(8):: r(3,0:nstep), p(3,0:nstep)
!      equivalence( r(1,0), u(1,0) )
!      equivalence( p(1,0), u(4,0) )
      real(8):: dt  ! time step in s. 
      real(8),external:: fBEdeflection ! function to give f of u'=f
         ! where  u=(r,p), u'=(v, Ze(E+VxB)
      real(8):: p2, E, gamma, bt
      integer::n

!      here we assume B is given by Mag (const).
      Bfld(:) =(/Mag%x, Mag%y, Mag%z/)
!      Bfld(:) =(/0.3d-4, 0.d0, 0.d0/)  ! test

      call cgetEfield(aTrack)  ! ans is in Efld

      u(1:3, 0) = 0.   ! r=0 we need only dr 
!      u(1:3, 0) =  aTrack.pos.xyz.r(1:3)

      u(4:6, 0) = aTrack%p%fm%p(1:3)  ! px,py,pz in GeV/c
      ptclZ = aTrack%p%charge
      ptclmass =aTrack%p%mass
      E = aTrack%p%fm%p(4)
      gamma = E/ptclmass
      bt = sqrt( (gamma-1)*(gamma+1))/gamma
      dt = leng/c/bt
      forall(n=0:nstep) t(n) = n*dt
      call kSRunge_Kutta4(t, u,  nstep, nsim, dt,
     *     fBEdeflection)
      dispmr%r(:) = u(1:3,nstep) - u(1:3,0)
      p2 =dot_product( u(4:6,nstep), u(4:6,nstep) )
      if( p2 > 0.) then
         dispmd%r(:) = u(4:6,nstep)/ sqrt( p2) 
      else
         dispmd%r(:) = (/0.,0.,1./)
      endif
      newmom(:) = u(4:6,nstep)
      end    subroutine cdefByMagAndE
      
      function fBEdeflection(t, u, nv)  result(ans)
!     charged ptcl motion in E and B
!     dp/dt =Ze(E+ VxB)=ZeE+Ze(VyBz-VzBy, VzBx-BxBz, VxBy-VyBx)
!     u=(x,y,z, px,py,pz) u'=(v, force)
!     v = pc/E
      use modBEfield
      implicit none
#include "Zglobalc.h"
      integer,intent(in):: nv   ! # of variable; simultaneous eq
      real(8),intent(in):: t  ! given time
      real(8),intent(in):: u(nv)  ! value of u for du/dt=f
      real(8)::ans(nv)  !  f
!               f(1:3) = dr/dt = beta c = pc/E
!               f(4;6) = dp/dt = Ze(E+vxB) 
      real(8)::p2, E
      real(8)::vp(3), temp(3)
      vp(:) = u(4:6)  ! p
      p2= dot_product( vp, vp)
      E = sqrt(ptclmass**2 + p2)
!      gamma = E/m
!      beta = sqrt( (gamma-1)*(gamma+1))/gamma
      ans(1:3) = vp(:)*c/E  ! pc/E = c*beta => v
      call kvec_prod3(ans(1:3), Bfld, temp(1:3))  ! vxB
      ans(4:6) = ptclZ*eval*(temp(1:3) + Efld(1:3))

      end function fBEdeflection







  


