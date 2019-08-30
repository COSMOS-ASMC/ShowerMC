!
!           initialize xsection calculation.
!  many intermediate values such as npercm3 etc are
!  dep. on rhoc but final quantities are indep. of rhoc
!
      subroutine epixsec(media)
      implicit none
#include "Zmedia.h"
!
       type(epmedia):: media   ! in/out.  
!

!      
      real*8 shp
      integer i, j

!           the output media.elem(j).sigma(i) will not
!       be used anywhere !!!
      do j = 1, media%noOfElem
         shp = xsecmin          ! (say 10mb)
         do i = 1, nxsec
            call cxp2xAXsec(media%elem(j)%A, shp,
     *                      media%elem(j)%sigma(i))
            shp = shp + xsecstep
         enddo
      enddo


      end
!!!!!!!!   use cGetXsec
      subroutine epgetxs(model, pj, media, xs, mfp)
      
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
#include "Zmedia.h"
!  #include "ZepTrackp.h"
#include "Zevhnp.h"
!   
      character*16 model    ! input. Active interaction model.
            ! for x-section calclation
       type(ptcl)::  pj      ! input.  a particle (hadronic)
       type(epmedia)::  media  ! input
      real*8 xs             ! output. obtained cross-section for the
                            !  media. in mb
                    !    if xs==smallxs, no collision
                    !       xs==largexs, instant collision
      real*8 mfp    ! output. obtained  mean free path in kg/m2


      real*8  mbTocm2
      parameter (mbTocm2=1.d-27)

      if( model == "jam" ) then
         call epxsJam(pj, media)
      elseif( model == "phits" ) then
         call epxsPhits(pj, media)
      elseif( model == "dpmjet3" ) then
         call epxsDpmjet3(pj, media)
      elseif( model == "qgsjet2" ) then
         call epxsQgsjet2(pj, media)
      elseif( model == "epos" ) then
         call epxsEPOS(pj, media)
      elseif( model == "sibyll" ) then
         call epxsSibyll(pj, media)
      elseif( model .eq. 'gheisha' ) then
         call epxsGheisha(pj, media)
      else
         call epxsOther( pj, media)
      endif
      if(media%sumns /= smallxs .and. media%sumns /= largexs) then
!           prop to rhoc so next two are indep. of rhoc
         xs = media%sumns/media%nd ! mb/cm3 /(1/cm3)
         if( xs <= 0. ) then
            xs = smallxs
            mfp = largexs
         else
            mfp = 10.*media%rho/(media%sumns*mbTocm2) ! in kg/m2. Not in g/cm2
         endif
                                          ! to be consistent with other  part.
!/////////////
!         if( pj.code == knuc) then
!            write(*,*)' E=',pj.fm.p(4),  ' charge=',pj.charge,
!     *      ' xs=',xs, ' mb',  ' mfp=', mfp*0.1
!         endif
!////////////
         
      elseif( media%sumns == smallxs ) then
         xs = smallxs
         mfp = largexs
      else
         xs = largexs
         mfp = smallxs
      endif
      end

      subroutine epxsJam(pj,  media)
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
#include "Zmedia.h"
!#include "ZepTrackp.h"
#include "Ztrackp.h"
#include "Zevhnp.h"
       type(ptcl)::  pj  ! input ptcl
       type(epmedia)::  media    ! input/output

      integer i
      real(8):: sumns,  xs
!//////////
!      real(8):: temp
!      integer,save:: ncall = 0
!//////////
      sumns=0. 
      do i =1, media%noOfElem
         if(pj%code >= kpion .and. pj%code <= knuc ) then
            if( JamXs == 1 ) then
               call ctotx(pj, media%elem(i)%A,  xs)
            elseif( JamXs == 0 ) then
               call cinelx(pj, media%elem(i)%A, media%elem(i)%Z, xs)
            else
               write(0,*)
     *           ' JamXs=',JamXs, ' not usable in epixsec%f'
               stop
            endif
         else
            call cinelx(pj, media%elem(i)%A,  media%elem(i)%Z, xs)
         endif
         if( xs == smallxs .or. xs == largexs ) then
            sumns = xs
            exit
         else
            media%nsigma(i) = xs*media%npercm3(i)
            sumns = sumns + media%nsigma(i)
         endif
      enddo
      media%sumns = sumns
      end

      subroutine epxsPhits(pj,  media)
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
#include "Zmedia.h"
#include "Zevhnp.h"
       type(ptcl)::  pj  ! input ptcl
       type(epmedia)::  media    ! input/output

      integer i, icon
      real(8):: sumns,  xs, u, elaxs
      integer::ka, subc, ia, iz

      sumns=0. 
      ka = pj%code
      subc = pj%subcode
      do i =1, media%noOfElem
         ia = media%elem(i)%A + 0.5
         iz = media%elem(i)%Z
         if( ( ka == knuc .and. subc /= antip )
     *        .or. ka == kgnuc) then
!c            call rndc(u)
!c            if(u <  media.elem(i).A - ia) then
!c               ia = ia +1
!c            endif
!            if(u < media.elem(i).Z- iz ) then
!               iz = iz + 1
!            endif
            call cphitsXs(pj, ia, iz, elaxs,xs, icon)
!c             no need to add. xs is already total 2010.Nov.16
!cc            xs = xs + elaxs  ! phits elaxs for heavy is 0
         else
!            if( ka >= kpion .and. ka <= knuc ) then
!               if( (pj.fm.p(4)-pj.mass)  < 10.d0 )  then ! include elastic
!                  call ctotx(pj, media.elem(i).A,  xs)
!               else
                  call cinelx(pj, media%elem(i)%A, media%elem(i)%Z, xs)

!               endif
!            else
!               call cinelx(pj, media.elem(i).A, media.elem(i).Z, xs)
!            endif
         endif
         if( xs == smallxs .or. xs == largexs ) then
            sumns = xs
            exit
         else
            media%nsigma(i) = xs*media%npercm3(i)
            sumns = sumns + media%nsigma(i)
         endif
      enddo
      media%sumns = sumns

      end

      subroutine epxsDpmjet3(pj, media)
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
#include "Zmedia.h"
!  #include "ZepTrackp.h"
#include "Zevhnp.h"
       type(ptcl)::  pj  ! input ptcl
       type(epmedia)::  media    ! input/output

      integer i
      real(8):: sumns,  xs

      sumns=0. 
      do i =1, media%noOfElem
         if( pj%code >= kpion .and. pj%code <= knuc ) then
            if( pj%fm%p(4) .lt.  4.1d0 ) then  ! Et is better than Ek
               call ctotx(pj, media%elem(i)%A,  xs)
            else
               call cinelx(pj, media%elem(i)%A,  media%elem(i)%Z, xs)
            endif
         else
            call cinelx(pj, media%elem(i)%A, media%elem(i)%Z, xs)
         endif
         if( xs == smallxs .or. xs == largexs ) then
            sumns = xs
            exit
         else
            media%nsigma(i) = xs*media%npercm3(i)
            sumns = sumns + media%nsigma(i)
         endif
      enddo
      media%sumns = sumns
      end


      subroutine epxsQgsjet2(pj, media)
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
#include "Zmedia.h"
!  #include "ZepTrackp.h"
#include "Zevhnp.h"
       type(ptcl)::  pj  ! input ptcl
       type(epmedia)::  media    ! input/output

      integer i, ia
      real(8):: sumns,  xs, tga, u, tgz


      sumns=0. 
      do i =1, media%noOfElem
         tga = media%elem(i)%A
         tgz = media%elem(i)%Z
         ia =tga 
         call rndc(u)
         if(u .lt.  tga - ia ) then
            ia = min(ia + 1, 209)
         endif
         if( (pj%code >= kpion .and. pj%code <= knuc) .or. 
     *        pj%code == kgnuc ) then
            call cxsecQGS(pj, ia, xs)               
         else
            call cinelx(pj, tga, tgz, xs)
         endif

         if( xs == smallxs .or. xs == largexs ) then
            sumns = xs
            exit
         else
            media%nsigma(i) = xs*media%npercm3(i)
            sumns = sumns + media%nsigma(i)
         endif
      enddo
      media%sumns = sumns
      end

      subroutine epxsEPOS(pj, media)
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
#include "Zmedia.h"
!  #include "ZepTrackp.h"
#include "Zevhnp.h"
       type(ptcl)::  pj  ! input ptcl
       type(epmedia)::  media    ! input/output

      integer i, ia, iz
      real(8):: sumns,  xs, tga, u, tgz

       type(ptcl)::  tg

      sumns=0. 
      do i =1, media%noOfElem
         tga = media%elem(i)%A
         tgz = media%elem(i)%Z
         ia =tga 
         iz =tgz 
         call rndc(u)
         if(u .lt.  tga - ia ) then
            ia = min(ia + 1, 209)
         endif

         if(ia > 1 ) then
            call cmkptc(kgnuc, ia, iz, tg)
         else
            call cmkptc(knuc, ia, iz, tg)
         endif
         tg%fm%p(1:3) = 0.
         tg%fm%p(4) = tg%mass
         call ceposIniOneEvent(pj, tg, xs)
         if( xs == smallxs .or. xs == largexs ) then
            sumns = xs
            exit
         else
            media%nsigma(i) = xs*media%npercm3(i)
            sumns = sumns + media%nsigma(i)
         endif
      enddo
      media%sumns = sumns
      end      subroutine epxsEPOS

      subroutine epxsSibyll(pj, media)
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
#include "Zmedia.h"
!  #include "ZepTrackp.h"
#include "Zevhnp.h"
       type(ptcl)::  pj  ! input ptcl
       type(epmedia)::  media    ! input/output

      integer i, ia, iz
      real(8):: sumns,  xs, tga, u, tgz

       type(ptcl)::  tg

      sumns=0. 
      do i =1, media%noOfElem
         tga = media%elem(i)%A
         tgz = media%elem(i)%Z
         ia =tga 
         iz =tgz 
         call rndc(u)
         if(u .lt.  tga - ia ) then
            ia = min(ia + 1, 209)
         endif

         if(ia > 1 ) then
            call cmkptc(kgnuc, ia, iz, tg)
         else
            call cmkptc(knuc, ia, iz, tg)
         endif
         tg%fm%p(1:3) = 0.
         tg%fm%p(4) = tg%mass
         if( media%name == "Air" ) then
            tg%subcode = 0   ! sibyll can do for almos Air target only. 
         endif
         call csibyllXs(pj, tg, xs)
         if( xs == smallxs .or. xs == largexs ) then
            sumns = xs
            exit
         else
            media%nsigma(i) = xs*media%npercm3(i)
            sumns = sumns + media%nsigma(i)
         endif
      enddo
      media%sumns = sumns
      end

      subroutine epxsGheisha(pj, media)
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
#include "Zmedia.h"
!  #include "ZepTrackp.h"
#include "Zevhnp.h"
       type(ptcl)::  pj  ! input ptcl
       type(epmedia)::  media    ! input/output

      integer i, ia
      real(8):: sumns,  xs, tga, tgz


      sumns=0. 
      do i =1, media%noOfElem
         tga = media%elem(i)%A
         tgz = media%elem(i)%Z
         if(pj%code  >= kpion .and. pj%code <= knuc ) then
            call cxsecGheisha(pj, tga,  tgz, xs)
         else
            call cinelx(pj, tga, tgz, xs)
         endif
         if( xs == smallxs .or. xs == largexs ) then 
            sumns = xs
            exit
         else
            media%nsigma(i) = xs*media%npercm3(i)
            sumns = sumns + media%nsigma(i)
         endif
      enddo
      media%sumns = sumns
      end

      subroutine epxsOther(pj, media)
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
#include "Zmedia.h"
! #include "ZepTrackp.h"
#include "Zevhnp.h"
          
       type(ptcl)::  pj      ! input.  a particle (hadronic)
       type(epmedia)::  media  ! input
 
      integer i
      real(8):: sumns, xs

      sumns = 0.
        do i = 1, media%noOfElem
         call cinelx(pj, media%elem(i)%A, media%elem(i)%Z, xs)
!//////////////
!         write(0,*) 'in othert;  xs=',xs, i,  media.noOfElem
!////////////
         if( xs == smallxs .or. xs == largexs ) then
            sumns = xs
            exit
         else
            media%nsigma(i) = xs*media%npercm3(i)
            sumns = sumns + media%nsigma(i)
         endif
      enddo
      media%sumns = sumns
      end


      subroutine epfixTarget(model, media)
      implicit none
#include "Zmedia.h"
#include "Zevhnp.h"
!         
!
!       Fix the target element, 
!  
      character*16 model    ! input. Active interaction model.
       type(epmedia):: media   ! input/output.  colElem, colA, colZ
      real*8 xs   ! output. x-section on selected target element (mb)
!      integer ia, iz        ! output target mass number and charge; now  
!      xs                    they are put colA and colZ colXs
      real*8 u, csigma
      integer  j, ia


!      if( media.noOfElem .eq. 1 ) then
      if(  media%sumns == smallxs .or.  media%sumns == largexs ) then
         j = 1
         media%nsigma(j) = media%sumns
      elseif( media%noOfElem .eq. 1 ) then
         j = 1
      else
         call rndc(u)
         u = u* media%sumns
         csigma = 0.
         do j = 1, media%noOfElem
            csigma = csigma + media%nsigma(j)
            if(u .le. csigma)  goto 10
         enddo
         write(0,*) 'medianame=', media%name, ' model=', model 
         write(0,*) 'media%sumns=',media%sumns
         write(0,*) 'media%noOfElem=', media%noOfElem
         write(0,*) 'media%nsigma=',  media%nsigma(1:j)
         write(0,*) ' u=',u, ' csigma=',csigma, ' j=',j
         call cerrorMsg('should not come here; epfixTarget',0)
      endif
 10   continue
      media%colElem = j
!          int value is taken.
!c      if(model .eq. "dpmjet3" ) then
         ia = media%elem(j)%A  + 0.5
!c      else
!c         call rndc(u)
!c         ia = media.elem(j).A
!c         if(u .lt. media.elem(j).A -ia ) then
!c            ia = ia + 1
!c         endif
!c      endif
      media%colZ = media%elem(j)%Z
      media%colA = ia
      media%colXs = media%nsigma(j)/media%npercm3(j)
      end
      subroutine epfixTarget2(model, media)
      implicit none
#include "Zmedia.h"
#include "Zevhnp.h"
      character*16 model    ! input. Active interaction model.
!          This may be used for small basic cross
!     section case so that the xA cross section is
!     propotinal to A. (for gamma A/ mu A)
!      for gammma, >=v9.16, epgetPhotoPxs is used and
!      this is no more used. (epfixTarget is used).
!         

       type(epmedia):: media   ! input/output  colElem, colA, colZ are output
      integer ia
      real*8 u, sum
      integer j
      if(  media%sumns == smallxs .or.  media%sumns == largexs ) then
         j = 1
         media%nsigma(j) = media%sumns
      elseif(media%noOfElem .eq. 1) then
         j = 1
      else
         call rndc(u)
         sum = 0.
         do j = 1, media%noOfElem
            sum = sum + media%w(j)
            if(u .le. sum) goto 10
         enddo
         call cerrorMsg('should not come here; epfixTarget2',0)
      endif
 10   continue
      media%colElem = j
!          int value is taken.
!      if(model .eq. "dpmjet3" ) then
         ia = media%elem(j)%A  + 0.5
!      else
!         call rndc(u)
!         ia = media.elem(j).A 
!         if(u .lt. media.elem(j).A -ia ) then
!            ia = ia + 1
!         endif
!      endif
      media%colZ = media%elem(j)%Z
      media%colA = ia
      media%colXs = media%nsigma(j)/media%npercm3(j)
      end


      subroutine epgetCaprate( media)
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
#include "Zmedia.h"
#include "Zevhnp.h"
          
       type(epmedia)::  media  ! input
 
      integer i
      real(8):: sumns, capr

      sumns = 0.
      do i = 1, media%noOfElem
         call cmucap( int(media%elem(i)%A), int(media%elem(i)%Z), 
     *            capr)
         media%nsigma(i) = capr*media%npercm3(i)
         sumns = sumns + media%nsigma(i)
      enddo
      media%sumns = sumns
      end

      subroutine epgetPhotoPxs(model, pj, media, xs, mfp)
!        epgetxs for photon hadron production
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
#include "Zmedia.h"
!  #include "ZepTrackp.h"
#include "Zevhnp.h"
!   
      character*16 model    ! input. Active interaction model.
            ! for x-section calclation .  at present not used
       type(ptcl)::  pj      ! input.  a photon
       type(epmedia)::  media  ! input
      real(8),intent(out):: xs   !  obtained cross-section for the
                            !  media. in mb
                    !    if xs==smallxs, no collision
                    !       xs==largexs, instant collision
      real(8),intent(out):: mfp  !  obtained  mean free path in kg/m2

      real*8  mbTocm2
      parameter (mbTocm2=1.d-27)

      integer i
      real(8):: sumns

      sumns = 0.
      do i = 1, media%noOfElem
         call cgpXsec(media%elem(i)%A, pj%fm%p(4), xs)
         if( xs == smallxs .or. xs == largexs ) then
            sumns = xs
            exit
         else
            media%nsigma(i) = xs*media%npercm3(i)
            sumns = sumns + media%nsigma(i)
         endif
      enddo
      media%sumns = sumns

      if(media%sumns /= smallxs .and. media%sumns /= largexs) then
         xs = media%sumns/media%nd ! mb/cm3 /(1/cm3)                 
         if( xs <= 0. ) then
            xs = smallxs
            mfp = largexs
         else
            mfp = 10.*media%rho/(media%sumns*mbTocm2) ! in kg/m2. 
         endif
      elseif( media%sumns == smallxs ) then
         xs = smallxs
         mfp = largexs
      else
         xs = largexs
         mfp = smallxs
      endif
      end
