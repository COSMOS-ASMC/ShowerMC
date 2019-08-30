!     **********************************
      function epBrgene(media, force, Ee, x) result(ans)
      use BPLPM  ! to use LPM const and s, Xc
      implicit none
!
!     Generic bremsstrahlung function for entire energy region.
!     The same note as for epBrgenex
!
!
#include "Zglobalc.h"
#include "ZbasicCnst.h"
#include "Zmedia.h"
#include "Zmass.h"
#include "ZepTrackp.h"
       type(epmedia)::  media  ! input. media 
      character(8),intent(in):: force
      real*8  Ee  ! input. electron energy in GeV.
      real*8  x   ! input. x = Eg/Ee.
      real(8):: ans


      real*8  epBremS,  epBrSfs
      real*8  epCompScrBrs
      external epBremS, epBrSfs, epCompScrBrs

      real(8):: Eeme 

      integer i
      real*8 f, Xc, xn
      real,parameter::er=1.e-4
      logical LPMworks
      real(8):: fXc, fNLPMx,  fLPMx, fLPMXc

      Eeme = Ee/masele
      LPMworks = .false.
         
      if( LPMeffect .and. Ee > Flpm* media%cnst%BremEeminLPM ) then
         call smigb0(media, Ee, x, s)  ! s is in BPLPM
         if( s < 1.) then
!               LPM may  work at small x  

!            Xc = 1./(1.+sconst/Ee) !Xc=x for s=1 so shuould be v<Xc 
!                       above one is correct according to Migdal's formula
!             but there appear some bump at x's < Xc, at low energies,
!             which seems to be  unnatural. Also at high energies, 
!             continuation is somewhat awkward.
!             So we set it little bit smaller value. Next factor 2 is
!             employe.     Factor 1.5,1.75 is still  not enough
            Xc = 1./(1.+ 2.0*sconst/Ee) !Xc=x for s=1 so shuould be v<Xc 
            LPMworks = x < Xc
         endif
      endif
      if( LPMworks ) then
         fLPMx = epBremSH(media, Ee, x) * media%cnst%NormSH
         fLPMXc = epBremSH(media, Ee, Xc) * media%cnst%NormSH
      endif

      if( (force == "?" .and. Ee <= media%cnst%BrEemaxS2 ) .or.
     *  force == "seltzer" ) then
!                   Seltzer region 
        if( LPMworks ) then
!          !  normalize at Xc to Seltzer
           fXc = epBrSfs(media, Eeme, Xc)* media%cnst%NormS
           f = fLPMx * fXc/fLPMXc
!           f = fLPMx/fNLPMx *
!     *        epBrSfs(media, Eeme, x)* media.cnst.NormS
        else
           if( Ee > EpartialSC .and. force == "?" ) then
!                 special request
              f = epBremS(media, Eeme, x)*media%cnst%NormPS
           else
              f = epBrSfs(media, Eeme, x)*media%cnst%NormS
           endif
        endif
      elseif( ( force == "?" .and. Ee <= media%cnst%BrScrE) .or.
     *   force(1:2) == "ps" ) then
!                    Z=1:  150 GeV ~ Z=90: 30 GeV
!          Partial screening
         if( LPMworks ) then
            fXc = epBremS(media, Eeme, Xc)*media%cnst%NormPS
            f = fLPMx * fXc/fLPMXc
!            f = epBremS(media,Eeme,x)*fLPMx/fNLPMx
!     *           * media.cnst.NormPS
         else
            f = epBremS(media, Eeme, x)*media%cnst%NormPS
         endif
      elseif( (force=="?" .and.  Ee >  media%cnst%CompScrE)  .or.
     *     force == "TsaiCS" ) then  ! note CompScrE= BrScrE
         if( LPMworks ) then
            fXc =  epCompScrBrs(media,  Xc)*media%cnst%NormCS
            f = fLPMx * fXc/fLPMXc
         else            
            f =  epCompScrBrs(media, x)*media%cnst%NormCS
         endif
      elseif( force == "CS+LPM") then
!            not used now.  included above
         if( LPMworks ) then
            f =  fLPMx 
         else
            f =  epCompScrBrs(media, x)*media%cnst%NormCS
         endif
      else
         write(0,*) ' force=',force,' invalid'
         stop
      endif
      ans = f
      end  function epBrgene
!     ************************************
      subroutine epBrgeneTX(xmin, xmax, tx)
      implicit none
#include "Zglobalc.h"
#include "ZbasicCnst.h"
#include "Zmedia.h"
#include "Zmass.h"
#include "ZBPgene.h"

      real*8  xmin  ! input.  Eg/Ee min.
      real*8  xmax  ! input. Eg/Ee max.
      real*8  tx    ! output.  Integral of Brems function from
                    !         x= xmin to xmax
      external epBrgenex
      real*8   epBrgenex, v1, v2, ans1
      real*8  Ee


      Ee = Eeme *masele
!      tx = 0.
      v2 = xmax
      v1=max( v2*0.9d0,  xmin)
      if(v1 >= v2 ) then
         tx = 0.
      else
         call k16pGaussLeg(epBrgenex, v1, v2, 16,  tx)
      endif
      v2 = v1
      do while (v1 .ne. xmin)
!         v1=max( v2/10.d0,  xmin)
         v1 = max( v2/4.d0,  xmin)
         if(v1 >= v2 ) then
            tx = 0.
         else
            call k16pGaussLeg(epBrgenex, v1, v2, 16,  ans1)
            tx = tx + ans1
            v2=v1
         endif
      enddo
!         vc = media.cnst.BremEgmin/(Ee-masele)
!         tx = 4. * ar02 * ( ( log(1.d0/vc) -
!     *                 (1.d0-vc)) * media.cScrMain +
!     +                 (1.d0-vc)* (1.d0+vc)/2 * media.cScrC1 )
      end subroutine epBrgeneTX
!     **************
      real*8 function epBrgeneSolv(v)
!     **************
      implicit none
!
!          used to solve total-cross-section * u = integral of
!          brem function from min to v.
!
!
!
      common/upsic/upsi,vmax
      real*8 upsi, vmax
      real*8 v

      real*8 ans
      call epBrgeneTX(v, vmax, ans)
!     epBrgeneSolv =1.0d0- ans/upsi
      epBrgeneSolv = log(ans/upsi)  ! ifans 0 ? seems ok
      end function epBrgeneSolv
