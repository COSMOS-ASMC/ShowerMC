      subroutine cGetXsec2(model, pj, media, xs, mfp)
       ! normally this is not to be used. use cGetXsec
!      use modXsecMedia, xmedia=>media
      implicit none
#include "Zcode.h"
#include "Zptcl.h"
#include "Zmedia.h"      
#include "Zevhnp.h"
      character(*),intent(in):: model
      type(ptcl)::pj           ! input.   projectile 
!     type(xsmedia),intent(inout):: media ! media
      type(epmedia),intent(inout):: media ! media 
      real(8),intent(out):: xs  ! collision xsection in mb
                    !    if xs==smallxs, no collision
                    !       xs==largexs, instant collision
      real(8),intent(out):: mfp ! mean free path in kg/m2

      call cGetXsec(model, pj, media, xs, mfp)
      if(xs /= largexs .and. xs /= smallxs) then
         xs = xs * media%sumNo
      endif
      end subroutine cGetXsec2

      subroutine cGetXsec(modelin, pj, media, xs, mfp)
! @JAXA, xs=0 leads to 0-divide error and after 10 repeats
! stop is made.  At other systems, 0 divide
! in this case has no bad effect. 
! So change is made as follows.
!    xs == smallxs is made to xs <= smallxs.
!
!         get hadronic interaction cross-section for
!    given interaction model,  projectile and medium.
!    The media information must have been set via
!    cXsecMedia.f90
!             xmedia=>media is to avoid name
!       collision of media  in modXsecMedia and
!       media argument in the subroutine def.
!      use modXsecMedia, xmedia=>media
      implicit none
#include "Zcode.h"
#include "Zptcl.h"
#include "Zmedia.h"      
#include "Zevhnp.h"
      character(*),intent(in):: modelin
      type(ptcl)::pj           ! input.   projectile 
!     type(xsmedia),intent(inout):: media ! media
      type(epmedia),intent(inout):: media ! media 
      real(8),intent(out):: xs  ! collision xsection in mb
                    !    if xs==smallxs, no collision
                    !       xs==largexs, instant collision
      real(8),intent(out):: mfp ! mean free path in kg/m2

      character(8):: model
      
      if( pj%code == knuc ) then
         if( pj%subcode == antip .and. pj%fm%p(4) - pj%mass <= 0.) then
            xs = largexs
            media%xs = xs
            mfp = 0. 
            return  ! *******************
         endif
      elseif( pj%code == kneue .or. pj%code == kneumu .or.
     *         pj%code == kneutau ) then
         xs = smallxs
         media%xs = xs
         mfp  = largexs
         return                 ! *******************
      endif

      model = modelin      
      if( model == "special" ) then
!     model might be given in next call on return by the user;
!     (pj and target are not suited for  giving xs so that
!     some another relevant model  might have been given there for XS calc.
         call cxsSpecial(pj, media, model)
      endif
       
      if( model /= "special" ) then
        select case(model)
          case( "phits" ) 
            call cxsPhits(pj, media)
         case( "dpmjet3" )
            call cxsDpmjet3(pj, media)
         case( "jam" ) 
            call cxsJam(pj, media)
         case( "qgsjet2" )
            call cxsQgsjet2(pj, media)
         case( "epos" )
            call cxsEPOS(pj, media)
         case( "sibyll" )
            call cxsSibyll(pj, media)
         case( "gheisha" ) 
            call cxsGheisha(pj, media)
         case("incdpm3")   
            call cxsIncdpm3(pj, media)
         case default
            call cxsOther( pj, media)
         end select
      endif
      if(media%xs < smallxs) then
         media%xs = smallxs
      endif
      if( media%xs /= smallxs .and. media%xs /= largexs) then
         mfp = 1.0d0/( media%mbtoPkgrm * media%xs)
      elseif( media%xs <= smallxs ) then
         mfp = largexs
      else
         mfp = 0.
      endif
      xs = media%xs

      end subroutine cGetXsec


      subroutine cxsJam(pj,  media)
!      use modXsecMedia, xmedia=>media
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
#include "Zmedia.h"      
#include "Ztrackp.h"
#include "Zevhnp.h"
      type(ptcl)::pj  ! input ptcl
!     type(xsmedia),intent(inout):: media    ! input/output
      type(epmedia),intent(inout):: media ! input/output


      integer i, iA
      real(8):: sumns,  xs, iAR
      sumns=0.

      do i = 1, media%noOfElem
         iA = media%elem(i)%A + 0.5d0
         iAR =  iA
         if(pj%code >= kpion .and. pj%code <= knuc ) then
            if( JamXs == 1 ) then
               call ctotx(pj, iAR,  xs)
            elseif( JamXs == 0 ) then
               call cinelx(pj, iAR, media%elem(i)%Z, xs)
            else
               write(0,*)
     *          ' JamXs=',JamXs, ' not usable in cxsJam'
               stop
            endif
         else
            call cinelx(pj, iAR,  media%elem(i)%Z, xs)
         endif
         if( xs <  smallxs ) then
            xs = smallxs
         endif
         if( xs == smallxs .or. xs == largexs ) then
            sumns = xs
            exit
         else
            media%nsigma(i) = xs*media%No(i)
            sumns = sumns + media%nsigma(i)
         endif
      enddo
      media%xs = sumns
      end subroutine cxsJam

      subroutine cxsPhits(pj,  media)
!      use modXsecMedia, xmedia=>media
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
#include "Zmedia.h"      
#include "Zevhnp.h"
      type(ptcl)::pj  ! input ptcl
!     type(xsmedia),intent(inout):: media    ! input/output
      type(epmedia),intent(inout):: media    ! input/output      


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
            call cphitsXs(pj, ia, iz, elaxs,xs, icon)
!c             no need to add. xs is already total 2010.Nov.16
!cc            xs = xs + elaxs  ! phits elaxs for heavy is 0
         else
!            if( ka >= kpion .and. ka <= knuc ) then
!               if( (pj.fm.p(4)-pj.mass)  < 10.d0 )  then ! include elastic
!                  call ctotx(pj, media%elem(i)%A,  xs)
!               else
            call cinelx(pj, media%elem(i)%A, media%elem(i)%Z, xs)

!               endif
!            else
!               call cinelx(pj, media%elem(i)%A, media%elem(i)%Z, xs)
!            endif
         endif
         if( xs <= smallxs .or. xs == largexs ) then
            sumns = xs
            exit
         else
            media%nsigma(i) = xs*media%No(i)
            sumns = sumns + media%nsigma(i)
         endif
      enddo
      media%xs = sumns
      end subroutine cxsPhits

      subroutine cxsDpmjet3(pj, media)
!      use modXsecMedia, xmedia=>media
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
#include "Zmedia.h"      
#include "Zevhnp.h"
      type(ptcl)::pj  ! input ptcl
!     type(xsmedia),intent(inout):: media    ! input/output
      type(epmedia),intent(inout):: media ! input/output

      integer i,  iA
      real(8):: sumns,  xs, iAR

      sumns=0. 
      do i =1, media%noOfElem
         iA = media%elem(i)%A + 0.5d0
         iAR = iA
         if( pj%code >= kpion .and. pj%code <= knuc ) then
            if( pj%fm%p(4) .lt.  4.1d0 ) then ! Et is better than Ek
!              call ctotx(pj, media%elem(i)%A,  xs)  !2017/sep/20 
               call ctotx(pj, iAR,  xs)
!//////////////////
!               if( pj.code == 6 .and. pj.charge == -1  .and.
!     *              (pj.fm.p(4)- pj.mass) < 1d-3 ) then
!                  write(0,*) '********* ', pj.fm.p(4)-pj.mass,
!     *             xs
!                  write(0,*) ' largexs=',largexs
!               endif
!///////////////////
            else
!     call cinelx(pj, media%elem(i)%A,  media%elem(i)%Z, xs)
               call cinelx(pj, iAR,  media%elem(i)%Z, xs)
            endif
         else
            call cinelx(pj, iAR, media%elem(i)%Z, xs)
!            call cinelx(pj, media%elem(i)%A, media%elem(i)%Z, xs)
         endif
         if( xs <= smallxs .or. xs == largexs ) then
            sumns = xs
            exit
         else
            media%nsigma(i) = xs*media%No(i)
            sumns = sumns + media%nsigma(i)
         endif
      enddo
      media%xs = sumns
      end subroutine cxsDpmjet3


      subroutine cxsQgsjet2(pj, media)
!      use modXsecMedia, xmedia=>media
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
#include "Zmedia.h"
#include "Zevhnp.h"
      type(ptcl)::pj  ! input ptcl
!     type(xsmedia),intent(inout):: media    ! input/output
      type(epmedia),intent(inout):: media ! input/output      


      integer i, ia
      real(8):: sumns,  xs, tga, u, tgz


      sumns=0. 
      do i =1, media%noOfElem
         tga = media%elem(i)%A
         tgz = media%elem(i)%Z
!         ia =tga 
!         call rndc(u)
!         if(u .lt.  tga - ia ) then
!            ia = min(ia + 1, 209)
!         endif
!     above  method is not so good. should  statistically be ok
!     For H2,  A  is 1.008 so some times ia becomes 2. If one makes
!     a table, 2 is always used. and dangerous.  To be accurate
!     media must be mixed A=1 andn A=2...
!     Now we treat such case roughly.  by round off (四捨五入)
         ia = tga + 0.5d0
         ia = min(ia, 209)
         if( (pj%code >= kpion .and. pj%code <= knuc) .or. 
     *        pj%code == kgnuc ) then
            call cxsecQGS(pj, ia, xs)               
         else
            call cinelx(pj, tga, tgz, xs)
         endif

         if( xs <= smallxs .or. xs == largexs ) then
            sumns = xs
            exit
         else
            media%nsigma(i) = xs*media%No(i)
            sumns = sumns + media%nsigma(i)
         endif
      enddo
      media%xs = sumns
      end subroutine cxsQgsjet2

      subroutine cxsEPOS(pj, media)
!      use modXsecMedia, xmedia=>media
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
#include "Zmedia.h"      
#include "Zevhnp.h"
      type(ptcl)::pj  ! input ptcl
!     type(xsmedia),intent(inout):: media    ! input/output
      type(epmedia),intent(inout):: media ! input/output

      integer i, ia, iz
      real(8):: sumns,  xs, tga, u, tgz

      type(ptcl)::tg

      sumns=0. 
      do i =1, media%noOfElem
         tga = media%elem(i)%A
         tgz = media%elem(i)%Z
         ia =tga + 0.5d0
         iz =tgz
!          see qgs part         
!         ia = tga
!         call rndc(u)
!         if(u .lt.  tga - ia ) then
!            ia = min(ia + 1, 209)
!         endif

         if(ia > 1 ) then
            call cmkptc(kgnuc, ia, iz, tg)
         else
            call cmkptc(knuc, ia, iz, tg)
         endif
         tg%fm%p(1:3) = 0.
         tg%fm%p(4) = tg%mass
         call ceposIniOneEvent(pj, tg, xs)
         if( xs <= smallxs .or. xs == largexs ) then
            sumns = xs
            exit
         else
            media%nsigma(i) = xs*media%No(i)
            sumns = sumns + media%nsigma(i)
         endif
      enddo
      media%xs = sumns
      end      subroutine cxsEPOS

      subroutine cxsSibyll(pj, media)
!      use modXsecMedia, xmedia=>media
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
#include "Zmedia.h"
#include "Zevhnp.h"
      type(ptcl)::pj  ! input ptcl
!     type(xsmedia),intent(inout):: media    ! input/output
      type(epmedia),intent(inout):: media ! input/output      

      integer i, ia, iz
      real(8):: sumns,  xs, tga, u, tgz

      type(ptcl)::tg

      sumns=0. 
      do i =1, media%noOfElem
         tga = media%elem(i)%A
         tgz = media%elem(i)%Z
         iz =tgz
         ia = tga + 0.5d0
!     see QGS part
!                 ia =tga 
!         call rndc(u)
!         if(u .lt.  tga - ia ) then
!            ia = min(ia + 1, 209)
!         endif
         ia = min(ia, 209)
         if(ia > 1 ) then
            call cmkptc(kgnuc, ia, iz, tg)
         else
            call cmkptc(knuc, ia, iz, tg)
         endif
         tg%fm%p(1:3) = 0.
         tg%fm%p(4) = tg%mass
         if( media%name == "Air" )then
            if( i == 1) then 
               tg%subcode = 0   !sibyll can do for almost Air target only.
               call csibyllXs(pj, tg, xs)
            else
               xs =0.
            endif
         else
            call csibyllXs(pj, tg, xs)
            if( xs <= smallxs .or. xs == largexs ) then
               sumns = xs
               exit
            endif
         endif
         if( media%name ==  "Air") then
!           only 1 virtual element 'air' so weight should be 1
!     In the event generator, should not choose N target
!     but use "air" target.  i>1 xs=0
            media%nsigma(i) = xs
         else
            media%nsigma(i) = xs*media%No(i)
         endif
         sumns = sumns + media%nsigma(i)
      enddo
      media%xs = sumns
      end subroutine cxsSibyll

      subroutine cxsGheisha(pj, media)
!      use modXsecMedia, xmedia=>media
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
#include "Zmedia.h"       
#include "Zevhnp.h"
      type(ptcl)::pj  ! input ptcl
!     type(xsmedia),intent(inout):: media ! input/output
      type(epmedia),intent(inout):: media ! input/output

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
         if( xs <= smallxs .or. xs == largexs ) then 
            sumns = xs
            exit
         else
            media%nsigma(i) = xs*media%No(i)
            sumns = sumns + media%nsigma(i)
         endif
      enddo
      media%xs = sumns
      end subroutine cxsGheisha

      subroutine cxsIncdpm3(pj, media)
!      use modXsecMedia, xmedia=>media
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
#include "Zmedia.h"      
#include "Zevhnp.h"
          
      type(ptcl)::pj      ! input.  a particle (hadronic)
!     type(xsmedia),intent(inout):: media    ! input/output
      type(epmedia),intent(inout):: media    ! input/output      

      real(8):: ek, crossint
      integer:: kinc

      ek = pj%fm%p(4)- pj%mass
      if( ek > 0.2d0 ) then
!            special for inclusive treatment.  target is always air   
!         *********************************                           
         call cccode2hcode(pj, kinc)
         media%xs = crossint(kinc, ek)
      else
         media%xs = smallxs
      endif
      end subroutine cxsIncdpm3

      subroutine cxsOther(pj, media)
!      use modXsecMedia, xmedia=>media
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
#include "Zmedia.h"      
#include "Zevhnp.h"
          
      type(ptcl)::pj      ! input.  a particle (hadronic)
!     type(xsmedia),intent(inout):: media    ! input/output
      type(epmedia),intent(inout):: media ! input/output      

 
      integer i, iA
      real(8):: sumns, xs, iAR

      sumns = 0.
      
      do i = 1, media%noOfElem
         iA =  media%elem(i)%A + 0.5d0
         iAR = iA
         call cinelx(pj, iAR, media%elem(i)%Z, xs)
!//////////////
!         write(0,*) 'in othert;  xs=',xs, i,  media%noOfElem
!////////////
         if( xs <= smallxs .or. xs == largexs ) then
            sumns = xs
            exit
         else
            media%nsigma(i) = xs*media%No(i)
            sumns = sumns + media%nsigma(i)
         endif
      enddo
      media%xs = sumns
      end subroutine cxsOther

      subroutine cfixTarget(media)
!     use modXsecMedia, xmedia=>media
!     use modColMediaInfo
      use modColInfo
      implicit none
#include "Zmedia.h"      
#include "Zevhnp.h"

!
!       Fix the target element, 
!  
!     type(xsmedia),intent(inout):: media  ! input/output
      type(epmedia),intent(inout):: media ! input/output
                        !   xsmedia information.


      real*8 u, csigma
      integer  j, ia


      if(  media%xs == smallxs .or.  media%xs == largexs ) then
         j = 1
         media%nsigma(j) = media%xs
      elseif( media%noOfElem .eq. 1 ) then
         j = 1
      else
         call rndc(u)
         u = u* media%xs
         csigma = 0.
         do j = 1, media%noOfElem
            csigma = csigma + media%nsigma(j)
            if(u <= csigma)  goto 10
         enddo
         write(0,*) 'media name=', media%name
         write(0,*) 'media%xs=',media%xs
         write(0,*) 'media%noOfElem=', media%noOfElem
         write(0,*) 'media%nsigma(1:j)=',
     *          media%nsigma(1:j)
         write(0,*) ' u=',u, ' csigma=',csigma, ' j=',j
         call cerrorMsg('should not come here; cfixTarget',0)
      endif
 10   continue
      colElemNo = j
!          int value is taken.
!c      if(model .eq. "dpmjet3" ) then
      ia = media%elem(j)%A  + 0.5
!c      else
!c         call rndc(u)
!c         ia = media%elem(j).A
!c         if(u .lt. media%elem(j).A -ia ) then
!c            ia = ia + 1
!c         endif
!c      endif
      !  next 3lines are used now in unified version
      media%colZ = media%elem(j)%Z
      media%colA = ia
      media%colXs = media%nsigma(j)/media%No(j)
!     in many cases from now, next will not be used now (see above)
!    ( in chAcol.f cheavyint.f use above)
      TargetProtonNo = media%elem(j)%Z
!      TargetNucleonNo = media%elem(j)%A
      TargetNucleonNo = ia
      TargetXs =  media%nsigma(j) / media%No(j)
!       Nxt needs not be called if !MacIFC  
!      (see chAcol.f and cheavyInt.f)
#if defined (MacIFC)      
!      call  cworkaround(TargetNucleonNo, TargetProtonNo, TargetXs,
!     *   colElemNo)
#endif
      end   subroutine cfixTarget

!!      subroutine cworkaround(A, Z, xs, nelem)  ! not needed now
!!!     With MacIFC,
!!!     TargetNucleonNo, TargetProtonNo, TargetXs
!!!     colElemNo
!!!     cannot be recognized at link time
!!!     from chAcol.f cheavyInt.f so we trnasfer info.
!!!     thru common.  (VERY STRANGE; compiler problem)
!!      implicit none
!!      integer,intent(in):: A  ! TargetNucleonNo
!!      integer,intent(in):: Z  ! TargetProtonNo
!!      real(8),intent(in):: xs ! TargetXs
!!      integer,intent(in):: nelem ! element # of the target elemnt
!!                      ! in the mediumNo
!!!    !#include "Zworkaround.h"
!!
!!      TargetNucleonNo = A
!!      TargetProtonNo  = Z
!!      TargetXs   = xs
!!      colElemNo = nelem
!!      end      subroutine cworkaround

      
      subroutine cfixTargetMuNI(media)
!        fix target for muon nuclear interaction
!      In  the case of muon n.i,  x-section for each
!      element is not computed and elem(:)%nsigma is
!      not ready. So we roughly compute it and fix
!      the  target 
!      use modXsecMedia, xmedia=>media
      implicit none
#include "Zmedia.h"      
#include "Zevhnp.h"

!     type(xsmedia),intent(inout):: media  ! input/output
      type(epmedia),intent(inout):: media ! input/output
                        !   xsmedia information.


      real*8 u, csigma, s0
      integer  j, ia

      s0=media%xs/sum( media%No(:) * media%elem(:)%A )
      media%nsigma(:) =
     *     s0 *media%No(:) * media%elem(:)%A
      call cfixTarget(media)
      end   subroutine cfixTargetMuNI

      subroutine cgetCaprate( media)
!      use modXsecMedia, xmedia=>media
      implicit none
!#include "Zptcl.h"
!#include "Zcode.h"
#include "Zmedia.h"
!#include "Zevhnp.h"
          
!     type(xsmedia),intent(inout):: media    ! input
      type(epmedia),intent(inout):: media    ! input      
 
      integer i
      real(8):: sumns, capr

      sumns = 0.
      do i = 1, media%noOfElem
         call cmucap( int(media%elem(i)%A), int(media%elem(i)%Z), 
     *            capr)
         media%nsigma(i) = capr*media%No(i)
         sumns = sumns + media%nsigma(i)
      enddo
      media%xs = sumns   ! this is not mb x-sec. but is used
                 ! to fix the target (with media%nsigma(i)
      end subroutine cgetCaprate

      subroutine cgetPhotoPxs(model, pj, media, xs, mfp)
!        cgetxs for photo-hadron production
!      use modXsecMedia, xmedia=>media
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
#include "Zmedia.h"
#include "Zevhnp.h"
!   
      character(*):: model    ! input. Active interaction model.
            ! for x-section calclation .  at present not used
      type(ptcl)::pj      ! input.  a photon
!     type(xsmedia),intent(inout):: media  ! input
      type(epmedia),intent(inout):: media  ! input      
      real(8),intent(out):: xs   !  obtained cross-section for the
                            !  media% in mb
                    !    if xs==smallxs, no collision
                    !       xs==largexs, instant collision
      real(8),intent(out):: mfp  !  obtained  mean free path in kg/m2

      integer i
      real(8):: sumns

      sumns = 0.
      do i = 1, media%noOfElem
         call cgpXsec(media%elem(i)%A, pj%fm%p(4), xs)
         if( xs < smallxs ) then
            xs = smallxs
         endif
         if( xs == smallxs .or. xs == largexs ) then
            sumns = xs
            exit
         else
            media%nsigma(i) = xs*media%No(i)
            sumns = sumns + media%nsigma(i)
         endif
      enddo
      media%xs = sumns

      if(media%xs /= smallxs .and. media%xs /= largexs) then
         if( media%xs <= 0. ) then
            xs = smallxs
            mfp = largexs
         else
            xs = media%xs
            mfp =1.0d0/( media%mbtoPkgrm *media%xs)
         endif
      elseif( media%xs == smallxs ) then
         xs = smallxs
         mfp = largexs
      else
         xs = largexs
         mfp = smallxs
      endif
      end      subroutine cgetPhotoPxs
