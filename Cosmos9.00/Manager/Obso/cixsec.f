      subroutine cixsec
!       setup Air media for cross-section calculations
      use modXsecMedia
      use modBremPairAng
      implicit none
#include "Zglobalc.h"
#include "Ztrackp.h"
#include "Ztrack.h"
#include "Ztrackv.h"
#include "Zevhnv.h"
#include "Zcode.h"

!     c
      
!      real*8 shp, xsecmin, xsecstep
!      parameter ( xsecmin = 10.d0, xsecstep = 10 )
      integer i, j, mbindex, jcon
      integer,parameter::noOfElem = 3
      real(8),save::elemA(noOfElem)=(/14.0d0, 16.0d0, 40.0d0/)
      real(8),save::elemZ(noOfElem)=(/7.0d0, 8.0d0, 18.0d0/)
!      real(8),save::No(noOfElem)=(/1.56d0, 0.42d0, 0.01d0/)
!        if next is not normalized, correct xsec will be obtained.
!        if normalized, same as above. correct xs may be obtained
!        by multiling ~2.55 in that case.  mfp is correct in
!        either case.
      real(8),save::No(noOfElem)=(/2.0d0, 0.5385d0, 0.0128d0/)

!      elemS(noOfElem, nxsec), sigma(nxsec)
!      real*8  cumsigma(noOfElem, nxsec)

      real*8  u, tgA, tgZ

!      TargetMassN = 14.55   !  sum( Ai*Noi)  / sum (Noi )
!      TargetAtomicN = 7.266  !  sum( Zi*Noi) / sum(Noi)
!      TargetZ2 = 53.54  (cf. <Z>^2 =  52.8)
      call cAllocXsMediaArray(1)
      call cAllocXsElemArray(1,noOfElem)
      media(1)%elem(:)%OrigNo = No(:)
      media(1)%sumNo = sum( No(:) )
        ! normalized
      No(:) = No(:) /media(1)%sumNo

      TargetMassN =sum( No(:)*elemA(:) )
      TargetAtomicN = sum( No(:)*elemZ(:) )

      media(1)%name = "Air"
      media(1)%noOfElem = noOfElem
      media(1)%elem(:)%A = elemA(:)

      media(1)%elem(:)%Z = elemZ(:)
      media(1)%elem(:)%No = No(:)
!      mfp = A2deninv*TargetMassN/xs ! in kg/m2
      media(1)%mbtoPkgrm = 1.d0/(A2deninv*TargetMassN)
!      then, mfp   = 1/(mbtoPkgrm*xs)  ! in kg/m2
      TargetZ2 = sum( No(:)*elemZ(:)**2 )
!      TargetZ2_3rd = TargetAtomicN**0.66666
      TargetZ2_3rd = sum( No(:)*elemZ(:)**(2.d0/3.d0) )
!      TargetZ1_3rd = TargetAtomicN**0.33333
      TargetZ1_3rd = sum( No(:)*elemZ(:)**(1.d0/3.d0) )
!//////////
!     <A>=   14.5527638190955     
!     <Z>=   7.26633165829146     
!     <Z^2>=   53.5477386934673     
!     <Z>^2=   52.7995757682887     
!     <Z>^(2/3)=   3.75149534225373     
!     <Z>^(1/3)=   1.93687773033140     
!     <Z^(2/3)>=   3.74733626838206     
!     <Z^(1/3)>=   1.93486427511770     
!      write(0,*) "<A>=", TargetMassN
!      write(0,*) "<Z>=", TargetAtomicN 
!      write(0,*) "<Z^2>=", TargetZ2
!      write(0,*) "<Z>^2=", TargetAtomicN**2
!      write(0,*) "<Z^(2/3)>=", TargetZ2_3rd
!      write(0,*) "<Z^(1/3)>=", TargetZ1_3rd
!//////////// 

      mediumNo=1  ! for Air, fixed
      TargetNucleonNo = 14   !  these 2  will be updated 
      TargetProtonNo =7      !  for each collision (also TargetXs)
          !  Target... were used before xsmedia was introduced
          ! in Cosmos.  We can live without Target.. but
          ! are kept for easy modifcation.
!      do j = 1, noOfElem
!         shp = xsecmin          ! (say 10mb)
!         do i = 1, nxsec
!            call cxp2xAXsec(elemA(j), shp,elemS(j,i))
!            shp = shp + xsecstep
!         enddo
!      enddo


!      do i = 1, nxsec
!         sigma(i) = 0.
!         do j = 1, noOfElem
!            sigma(i) = sigma(i) + No(j)* elemS(j,i)
!         enddo
!      enddo

!      do i = 1, nxsec
!         do j = 1, noOfElem
!            cumsigma(j, i) =
!     *           No(j)*elemS(j,i)/sigma(i)
!         enddo
!      enddo
!
!      do i = 1, nxsec
!         do j = 2, noOfElem
!            cumsigma(j, i) = 
!     *         cumsigma(j, i) + cumsigma(j-1, i)
!         enddo
!c          for safety
!         cumsigma(noOfElem, i) = 1.0
!      enddo
!  &&&&&&
!      do i = 1, nxsec
!         write(*,*) ' sigma=', sigma(i)
!         write(*,*) (cumsigma(j, i), j=1,noOfElem)
!      enddo
!  &&&&&&


!!!!!!!!!!!!!!!
!
!      exerpt from epGetEffZA;  only for rigorous Moliere scattering
!         
!      sum = 0.
!      media.A = 0.
!      media.Z = 0.
!      media.Z2 = 0.
!      media.ZZ1 = 0.
!      media.MoliereExpb = 0.
!      media.Z1_3rd =0.
!      media.Z2_3rd =0.
!      do i = 1, nc
!         media.A = media.A + media.No(i)* media.elem(i).A
!      media(1)%A = sum( media(1)%elem(:)%No* media(1)%elem(:)%A )
!      TargetMassN =sum( No(:)*elemA(:) )
      media(1)%A =   TargetMassN 

!         media.Z = media.Z + media.No(i)* media.elem(i).Z
!      media(1)%Z = sum( media(1)%elem(:)%No* media(1)%elem(:)%Z )
!      TargetAtomicN = sum(No(:)*elemZ(:))
      media(1)%Z =   TargetAtomicN

!         media.Z2 = media.Z2 + media.No(i)* media.elem(i).Z**2

      media(1)%Z2 = sum( No(:)* elemZ(:)**2 )
!         media.ZZ1 = media.ZZ1 + 
!     *     media.No(i)* media.elem(i).Z*(media.elem(i).Z + 1.)
!      media(1)%ZZ1 = sum(
!     *   media(1)%elem(:)%No*
!     *   media(1)%elem(:)%Z*(media(1)%elem(:)%Z + 1.)
!     *     )

      media(1)%ZZ1 = sum( No(:)* elemZ(:)*(elemZ(:) + 1) )


!            for beta=1 and z=1; Moliere's  exp(b) = MoliereExb*t
!         media.MoliereExpb = media.MoliereExpb +
!     *       media.No(i)*media.elem(i).Z**0.333333 *
!     *        (media.elem(i).Z+1.d0)/
!     *              (1.0d0 + 3.327d0*(media.elem(i).Z/137.d0)**2)

      media(1)%MoliereExpb = sum(
     *      media(1)%elem(:)%No*media(1)%elem(:)%Z**0.333333*
     *      (media(1)%elem(:)%Z+1.d0)/
     *        (1.0d0 + 3.327d0*(media(1)%elem(:)%Z/137.d0)**2)
     *    )

!         media.Z1_3rd =  media.Z1_3rd 
!     *          +  media.No(i)* media.elem(i).Z**(1./3.)
!         media.Z2_3rd =  media.Z2_3rd 
!     *          +  media.No(i)* media.elem(i).Z**(2./3.)
!      enddo
!                only for e+/e- 
      media(1)%MoliereForXc2 =0.6011 * media(1)%ZZ1/media(1)%A

      media(1)%MoliereExpb = 6702. * media(1)%MoliereExpb/
     *     media(1)%A
!     for    brems/pair angle  v7.646
      call cBremPairAngInit(1, media(1)%A, media(1)%Z)
      end
