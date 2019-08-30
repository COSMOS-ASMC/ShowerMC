      module BPLPM
      implicit none
      character(8),save::name=' '
      private
      real(8),save:: s1, alogs1, sconst, x0g, conv2mb, Z,
     *                Eesave, Egsave
      real(4),save:: s
      real,parameter:: pi12=37.699112
      real,parameter:: pi6=18.849556
      real,parameter:: er= 1.e-3
      real,parameter:: eps= 1.e-4
      public sconst, s
      public smigb, smigp, epBremH, epPairH, epBPLPMconst
      public epBremSH, epPairSH, smigb0
      contains

!     ****************************************************************
!      cross-sections is given in mb/target.
!            That is, compound or molecule is regarded as
!        one atom with some effective (A,Z).
!

!     interface to epPairH and epBremH
!     which are in epBPfuncH.f   for the LPM region
!    *************************
      function epPairSH(media, Eg, x) result(ans)
!    *************************
      implicit none
#include "Zmedia.h"
       type(epmedia)::  media  ! input target media
      real(8),intent(in):: Eg  ! Eg  GeV
      real(8),intent(in)::  x    !  Ee/Eg
      real(8):: ans   !     ds/dx in mb

      call epBPLPMconst(media)  
      ans = epPairH(media, Eg, x) ! s is obtained within epPairH
                             ! since, this is used only at high E.
      end     function epPairSH

!    *************************
      function epBremSH(media, Ee, x) result(ans)
!    *************************
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
       type(epmedia)::  media  ! input target media
      real(8),intent(in):: Ee  ! Ee GeV
      real(8),intent(in):: x     !  Eg/Ee
      real(8):: ans       ! ds/dx mb


!      write(0,*) ' epBPLPMconst'
      call epBPLPMconst(media)
!      Eesave = Ee
      call smigb0(media, Ee, x, s) ! s is  put in BPLPM
!      write(0,*) ' epBremH'
      ans = epBremH(media, Ee, x)
!      write(0,*) ' exiting BremSH'
      end     function epBremSH

      subroutine smigb0(media, Ee, x, sout)
      implicit none
#include "Zmedia.h"
       type(epmedia)::  media
      real(8),intent(in):: Ee  ! electron energy GeV
      real(8),intent(in):: x   ! Eg/Ee
      real(4),intent(out):: sout !Migdal  s value obtained. s>1==>LPM works
 
      !      
      real(4):: E, v
                  !v, E must be single
      v = x    
      if(v >= 1.0) then
!         for x= 0.999999.. v becomes 1
         sout = 10.
      else
         E = Ee

         sout = smigb(v, E, 1.0, er)
      endif
      end subroutine smigb0

!     ****************************************************************
!     *                                                              *
!     * smigb:  get root s, from recursive relation                  *
!     *                                                              *
!     ****************************************************************
!
!
      function smigb(v, E, sin, er) result(ans)
      implicit none
!
      real,intent(in):: v, E, sin, er
      real:: ans

      real:: ss, s
!      real::  sbrem2   ! internal
!
      s = sin
      if( v >= 1.0 ) then
         ss = 10.
      else
         do
            ss=sqrt(sbrem2(v,E,s))
            if(abs((s-ss)/ss) < er) exit
            s=ss
         enddo
      endif
      ans = ss
      end function smigb
!
!     ***********
      function smigp(v, E, sin, er) result(ans)
!     ***********
      implicit none
!
      real,intent(in):: v, E,  sin, er
      real:: ans

      real:: ss, s

!      real:: spair2   ! internal

      s = sin

      do 
         ss = sqrt( spair2(v, E, s)  )
         if(abs((s-ss)/ss) < er) exit
         s=ss
      enddo
      ans = ss
      end function smigp
!     ****************************************************************
!     *                                                              *
!     * sbrem2:  auxliary function for brem with landau effect       *
!     * spair2:  //                    pair                          *
!     *                                                              *
!     ****************************************************************
!
!
      function sbrem2(v,e,s) result(ans)
      implicit none
!
      real,intent(in):: v, E,  s
      real:: ans

      real:: tmp
!      real::  gzai  ! internal

      tmp=sconst*v
      ans = tmp/(1.-v)/E/gzai(s)
      end function sbrem2

      function  spair2(v,E,s) result(ans)
      implicit none
      real,intent(in):: v, E,  s
      real:: ans

      real:: tmp
!      real:: gzai  ! internal

      tmp=sconst/v
      ans = tmp/(1.-v)/E/gzai(s)
      end function spair2
!     ****************************************************************
!     *                                                              *
!     * gzai:  gzai function which appear in ladanu effect           *
!     *                                                              *
!     ****************************************************************
!
!
      function gzai(s) result(ans)
      implicit none
      real,intent(in):: s
      real:: ans

      if(s > 1.) then
         ans = 1.
      elseif(s > s1) then
         ans =  log(s)/alogs1+1.
      else
         ans = 2.
      endif
      
      end function gzai
!     ****************************************************************
!     *                                                              *
!     * gmigdl:  g(s) function which appear in landau effect         *
!     * psimig:  pis(s) //                                           *
!     *                                                              *
!     ****************************************************************
!
!             .... psiim is needed.....
!
      function gmigdl(s,eps) result(ans)
      implicit none
      real,intent(in):: s, eps
      real::ans

      real  psiim  ! exernal
!
      ans = (pi12*s-48.*s*s*psiim(s+0.5,s,0,eps))*s
      end function gmigdl
!
!     ************
      function psimig(s,eps) result(ans)
!     ************
      implicit none
      real,intent(in):: s, eps
      real:: ans
      real psiim   ! external
!

      ans = ( (psiim(s,s,1,eps)*s*24.-pi6)*s+6.) *s
      end function psimig
!     ****************************************
!     *                                                             
!     * epBPLPMconst:  LPM const is calculated

!      For LPM, we use effective Z / media but not Z /atom or Z/molecule
! 


!
      subroutine epBPLPMconst(media)
      implicit  none
#include "Zmedia.h"      
#include "Zmass.h"      
!
       type(epmedia):: media  ! input
 
!      s1=    ( z**0.3333333/ 183 )**2
      if(name /= media%name) then
         name = media%name
         s1 = media%s1
         alogs1 =  media%logs1
         conv2mb = 1.d0/media%mbtoPX0 
!
!        
         Z = media%Z
         x0g = media%X0g
!         const in eq.60 of migdal's paper. phys. rev. vol 103 1956
!         energy is in gev
         sconst=( 1.37e3 ) **2 * media%X0 * masele
      endif
      end subroutine epBPLPMconst


      function epBremH(media, Ee, vv) result(ans)
      implicit none
#include "Zglobalc.h"
#include "ZbasicCnst.h"
#include "Zmass.h"      
#include "Zmedia.h"      
       type(epmedia)::  media
      real(8),intent(in)::  Ee  ! Ee GeV
      real(8),intent(in)::   vv   ! Eg/Ee

      real(8):: ans    ! ds/dx in mb
      integer i


      real(8):: Tomb   ! to mb conversion
      parameter (Tomb = 1.d27/N0)
!       The original LPM cross section is given in unit of prob/r.l,
!       With  N being the number density(/g),
!             L the length and S the cross-section
!         Note description below dose not apply to LPM case
!       NLS=1.    
!       N =  N0/A (/g), 1/L (/r.l) = 1/( L*X0g) (/(g/cm^2))
!       so that
!       S = prob/r.l * (1/X0g) * A/N0  (X0g is the radiaiton length
!       in unit of (g/cm^2).  Since X0g is propotional to A we may
!       express X0g= A*x0g where x0g is the r.l for A=1 with the same
!       Z for X0g. Then
!       S= prob/r.l / (x0g * N0) (cm^2).  To get this in mb
!       10^27 must be multiplied.
!
       real:: v
!      real:: smigb, gzai, gmigdl, psimig  ! internal

      real(8)::f,  epCompScrBr, epBremS

      
      v=vv


      if(vv .ge. .99999999d0) then
         ans = (v*v+2.*(1.+(1.-v)*(1.-v)))/v/3.
!         ans = epCompScrBr(Zin, vv)     ! or
      elseif(vv .eq. 0.) then  
         ans = 0.
      else

 
!         s = smigb(v, E, 1.0, er)  
!          use s in module
!///////////
!         write(0,*) 'in bremH E=',Ee, ' s=',s, ' v=',vv
!/////////////
         if( s > 1.0) then
!          ans =   epCompScrBr(Zin, vv)   ! or
!            f=epBremS(vv)  ! alredy in mb
!                            or
            ans = (v*v+2.*(1.+(1.-v)*(1.-v)))/v/3.  *  conv2mb  
         else
!            write(0,*)  ' smigb 2:  v E er =', v, E, er
!            s = smigb(v, E, s, er)
            ans = 
     *       gzai(s)/v*(v*v*gmigdl(s,eps)+2.*(1.+(1.-v)*(1.-v))*
     *       psimig(s, eps))/3.

            ans =  ans * conv2mb
!
!           note that as v-->0, gzai(s) becomes 2 and
!           epBremH---> 2/v *( v*v*12pi*s**2 + 2*(1+(1-v)**2 )* 6 s) )/3
!           where s---> sqrt( sconst*v/2/e/(1-v))
!           so that epBremH--->8*sqrt(2*sconst/v/e)
!      
!           at s =1 ,  gzai=1, gmigdl=1, psimig=1
!           so that epBremH becomes the same as 
!               (v*v+2.*(1.+(1.-v)*(1.-v)))/v/3. 
!           we normalize this to be Tsai's CS or PS function
!           value at v=vc which gives s=1 at  given Ee 
!
         endif
      endif
      end function epBremH

!     ***********
      function epPairH(media, Eg, vv) result(ans)
      implicit none
#include "Zglobalc.h"
#include "ZbasicCnst.h"
#include "Zmass.h"
#include "Zmedia.h"
       type(epmedia)::  media
      real(8),intent(in):: Eg  ! input  Eg  GeV
      real(8),intent(in)::  vv   ! input  Ee/Eg, Ee is the higher energy  of pair.
      real(8):: ans

      real(8)::  epCompScrPrs
      real E, v, s
!      real  gzai, gmigdl, psimig, smigp  ! internal
      real*8 cmTomb
      parameter (cmTomb = 1.d27/N0)


      E = Eg
      v=vv
      if(vv .ge. 0.9999d0 ) then
         ans = (1.+2.*(v*v+(1.-v)*(1.-v)))/3.* conv2mb
      elseif(vv .eq. 0.) then
         ans = 0.
      else
!         s = spair2(v,E,s)
         s = 2.0
         s = smigp(v,E,s,er)    ! this is always used here
        if(s .gt. 1.) then
           ans = (1.+2.*(v*v+(1.-v)*(1.-v)))/3.* conv2mb
!             above one and next one give the same result
!        within a very very small diff.  But if you use
!        Media.Z or Media.Zeff with epCompScrPr, the 
!        result will be very much different for non sigle
!        element media such as PWO or BGO etc. So we
!        use above simple formula (for W, Pb epCompScrPr
!        can be used).   
!         ans = epCompScrPrs(media, vv)
        else
!            s = smigp(v,E,s,er)
            ans = 
     *       gzai(s)/3.*(gmigdl(s,eps)+2.*(v*v+(1.-v)*(1.-v))*
     *       psimig(s, eps))
            ans = ans * conv2mb
        endif
      endif

      end function epPairH
      end module BPLPM
