!     generate Cerenkov light from a given charged track
!         
      subroutine epLightCerenkov( path )
      use modepLight
      use modepLightPty
      implicit none
#include "Zcode.h"
#include "ZepTrackp.h"
#include "ZepTrackv.h"
#include "Zcnfig.h"
#include "ZsepManager.h"


       type(epTrack)::  path  ! charged ptcl segment for Cerekov


      real(8):: wl, wl0    ! wave length 
!      real(8), pointer:: wl1, wl2
      real(8):: wl1, wl2
      real(8)::  E1, E2, E
      real(8):: photons ! <no. of photons from scinti>
      integer:: nphotons ! with poisson fluctuatin
      integer:: npreal !  really sampled #
      integer::Nsmp  ! max # of really sampled
       type(epTrack)::  cerenlight
      integer:: i
      real(8):: u, betax
      real(8):: refracn, maxf
      real(8):: n1, n2, ans, cost, sint, cs, sn

       type(epDirec)::  w
      integer:: icon 
      real(8):: weight

      call epLightCerenTh( path, wl1, wl2, E1, E2, n1, n2, betax,icon)
!///////
!      write(0,*) ' aft th;  wl1, wl2, E1, E2, n1, n2, betax,icon '
!      write(0,'(1p,7g12.3,0p,i3)')  
!     *    wl1, wl2, E1, E2, n1, n2, betax,icon 
!/////////
      if( icon == 0 ) then
         if( n1 == n2 ) then
            !  assume n is independent of wl ; so integral (1/n^2) from
            !  E2 to E1 is (E1-E2)x 1/n^2
            photons =
     *        370.*path%p%charge**2 *(E1-E2)*(1.-1/(n1*betax)**2)
     *      * path%wl   ! wl is path lenght in this case
!///////////
!            write(0,*)  ' photons=',photons, ' path.wl=', path.wl
!            write(0,*) 'E1-E2=', E1-E2, ' charge=',path.p.charge
!/////////
         else   
            ! we must integrate  1/n^2 dE from E2 ~ E1.  (E1 is fixed
            ! already in epLightPreCeren.)
!///////
!            call Lcompchk(' cc ',  cLcompNo)
!            call Lcompchk(' ccp ',  cPtyNo)
!////////
            call csampAFinteg( comInfo( cPtyNo )%invN2intAF, E2, ans)
            photons = 370.*path%p%charge**2 *((E1-E2) -
     *            ans/betax**2) * path%wl  
         endif
         if( photons < usePoisson ) then
            call kpoisn(photons, nphotons)
         else
            nphotons = photons + 0.5
         endif

            
         cerenlight = path   ! inherit some from charged ptcl track
          !   make light:       Ceren's subcode is 2
         call cmkptc(klight, kceren, 0, cerenlight%p)
         Nsmp = comInfo( cPtyNo )%NpSample 
         if(nphotons > Nsmp) then
!               no starter no stopper will be needed
!            cerenlight.p.subcode= cerenlight.p.subcode + 100
            weight = 1. + real(nphotons -Nsmp)/Nsmp
!////////////
!            write(0,*) ' cw=', weight
!//////////
!            call eppush(cerenlight)   ! only act as a stopper 
            npreal = Nsmp
         else
            npreal = nphotons
            weight = 1.
         endif
         sumni = sumni + npreal
         sumniwi = sumniwi + npreal*weight

         cerenlight%wgt = weight
         cerenlight%p%subcode=kceren

         do i = 1, npreal
              !    light generated uniformly along path
            call rndc(u)
            cerenlight%pos%x = path%pos%x + u*path%wl*path%w%x
            cerenlight%pos%y = path%pos%y + u*path%wl*path%w%y
            cerenlight%pos%z = path%pos%z + u*path%wl*path%w%z
              ! sample energy from 1/(1- 1/n^betax^2) dE in E2~E1
            if(n1 == n2 ) then
               call rndc(u)
               E = (E1-E2)*u + E2
               call epLightE2wl(E, 1.d0,  wl0, wl)
               refracn = n1
            else
               maxf = 1./(1.0 - 1.0/(n1*betax)**2 )
                   ! use rejection method
               do while (.true.)
                  call rndc(u)
                  E = (E1-E2)*u + E2
                  call  epLightE2wl(E, 1.d0,  wl0, wl)
!///////////////
!                  call Lcompchk(' dd ', cLcompNo)
!                  call Lcompchk(' ddp ', cPtyNo)
!/////////
                  call  epLightwl2N(wl, 
     *              Lcomp(cLcompNo), comInfo(cPtyNo), 
     *              Media(MediaNo)%gasF, refracn)
                  call rndc(u)
                  if(u < 
     *             (1./(1.0 - 1.0/(refracn*betax)**2))/ maxf ) then
                     exit
                  endif
               enddo
            endif
               !  sample angle
            cost = 1./refracn/betax
            if( cost > 1.0 ) then
               write(0,*) 'warning:  Cerenkov angle cos =',cost, '>1'
               write(0,*) ' n=',refracn, ' beta=',betax
               cost = 1.
            endif
               
            sint = sqrt(1.0 - cost**2)
            call kcossn(cs,sn)
            w%x = cs*sint
            w%y = sn*sint
            w%z = cost
            call eptransVect(path%w, w, cerenlight%w)
            cerenlight%wl = wl
            cerenLight%p%fm%p(4) = E
            call epe2p(cerenlight) ! energy is eV but ok since m=0 
!            if(i == npreal) then
!               if( nphotons > Nsmp ) then
!                  cerenlight.p.subcode =  cerenlight.p.subcode +10  ! starter flag
!               endif
!            endif
            call eppush(cerenlight)
         enddo
      endif
      end
!   
!       see if Cerenkov light can be generated
      subroutine epLightCerenTh( path, wl1, wl2, E1, E2, n1, n2, betax,
     *                           icon)
      use modepLight
      use modepLightPty
      implicit none
#include "Zcode.h"
#include "ZepTrackp.h"
#include "ZepTrackv.h"
#include "Zcnfig.h"
#include "ZsepManager.h"


       type(epTrack)::  path  ! input. charged ptcl segment for Cerekov
                             ! note: path.wl should have path length in cm


!      real(8), pointer:: wl1, wl2
      real(8),intent(out):: wl1, wl2  ! Cerenkov w.l  region (wl1< wl2)( nm)
      real(8),intent(out):: E1, E2    ! corresponding energy in eV. (E1>E2)
      real(8),intent(out):: n1, n2  !  corresponding refracton index.
      real(8),intent(out):: betax   ! beta of the ptcl
      integer,intent(out):: icon   ! 0--> Cerenkov possilbe. 1 --> no Ceren

      !----------------------
      real(8):: gamma

      real(8):: wl, wl0    ! wave length 

      integer:: gasF
!////////
      if( cLcompNo == 0 ) then
         write(0,*) ' in cerenth; cLcompNo=0 ; cn=',
     *      path%cn,' code=', path%p%code, ' pos=',
     *      path%pos%x,  path%pos%y,  path%pos%z,
     *      path%w%x,  path%w%y,  path%w%z
      endif
      
!      call Lcompchk(' ee ',  cLcompNo )
!      call Lcompchk(' eep ', cPtyNo)
!/////////

      wl1 = comInfo(cPtyNo)%minmaxWL(1)
      wl2 = comInfo(cPtyNo)%minmaxWL(2)
      gamma = path%p%fm%p(4)/path%p%mass
      betax = sqrt(1.-1./gamma**2)
      call epLightgasF( path%cn, gasF)
      call epLightwl2N(wl1,
     *    Lcomp(cLcompNo), comInfo(cPtyNo), gasF, n1)
      if( n1*betax <= 1.)  then   ! < 1. no Cerenkov
         icon = 1
      else
         icon = 0
         E1 = comInfo(cPtyNo)%ElightH
         call epLightwl2N(wl2, 
     *    Lcomp(cLcompNo), comInfo(cPtyNo),gasF, n2)
         if( n2*betax <= 1.) then 
               ! close to the threshold
            n2 = 1./betax
            call epLightN2wl(n2,
     *       Lcomp(cLcompNo), comInfo(cPtyNo), gasF, wl)
            !     should be  wl1<  wl < wl2 
            if( wl1 > wl .or. wl > wl2 ) then
               ! stange
               write(0,*) '*** warning in epLightCerenkov:'
               write(0,*) ' wl1=', wl1, ' wl2=',wl2
               write(0,*) ' neglect  Cerekov'
               icon = 1
               return !**************
            endif
            wl2 = wl
            call epLightwl2E(wl2, 1.d0, wl0, E2) 
         else
            E2=comInfo(cPtyNo)%ElightL
         endif
      endif

      end


