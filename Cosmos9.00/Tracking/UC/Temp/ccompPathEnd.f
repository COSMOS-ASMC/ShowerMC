!!!  note   #ifdef XXXXXX  before cputDeflection  for debug
!     ccompPathEnd
!       compute path end information. including multiple scattering and mag.
!       deflection, Electric field effect.
!
      subroutine ccompPathEnd
      use modSetIntInf
      implicit none

#include  "Ztrack.h"
#include  "Ztrackv.h"
#include  "Ztrackp.h"
#include  "ZmediaLoft.h"      
!
      real(8):: Rx(3), dir(3), len
      integer:: tobecut
      integer:: nodeX
!
!         get endpoint coord. assuming no scat, no mag. effect. at first
!     MovedTrack has complete information at path end.
!!!!!!!!!!
      write(0,*) ' in comppathend       IntInfArray(ProcessNo)%length '
      write(0,*)   IntInfArray(ProcessNo)%length
!!!!!!!!!!!!!!!      
       MovedTrack = TrackBefMove   !   copy TrackBefMove into MovedTrack. 
       if(IntInfArray(ProcessNo)%length .gt. 0.0 ) then
          call cmoveStreight(IntInfArray(ProcessNo)%length, 
     *         TrackBefMove%vec%w)
!!!!!!!!!
          write(0,*) ' back from movestreight'
!!!!!!!!!!          
!         If a chargde ptcl, compute energy loss by dE/dx  and  reset energy
!         in MovedTrack. If it is too large, path is truncated. Also compute
!         deflection.

          if(TrackBefMove%p%charge .ne. 0) then
!            if(.not. (mod(HowGeomag, 10) .eq. 1 .and.     this existed from uv4.92 to 5.13
!     *                Zfirst .eq. 0.)) then
             call cputEnergyLoss
!!!!!!!!!
             write(0,*) ' back from put Eloss; MoveStat', MoveStat
!!!!!!!!!!          
!               endif
!                if "if" is omitted, segmentation violation
!                may take place occasionally
             if(MoveStat .ne. Dead) then
!!!!!!!!
                write(0,*) ' enter cputDeflection'
!!!!!!!!!!             
                call cputDeflection
             endif
          endif
          if( NoOfMedia > 1 ) then
!              see if goes into diff. medium. before MovedTrack
             call cmanageMediaChange(tobecut, nodeX, Rx, dir,len)
             if( tobecut == 0  ) then
!              goes into diff. media.  Rx is closer than MovedTrack.
!     amove TrackBefMov by len length with dir direction
!     dir may be diff. from direction in TrackBefMove
                call ccutAndAdjust(dir, len)
             endif
          endif      
       else
          EnergyLoss = 0.
       endif
       end
!     --------------------------------------------------------
       subroutine cmoveStreight(leng, dir)
!        move a track by leng m
!         after this is called,  MovedTrack has inf for moved pos.
!         MovedTrack.pos.where may not be correct. It must be
!         corrected by seeing if the track passes acrosss a
!         observation layer.
       implicit none

#include  "Ztrack.h"
#include  "Ztrackp.h"
#include  "Ztrackv.h"
      


       real*8  leng, cnewcos
       type(coord)::dir        ! input.
!

!            new coord in xyz system.

       MovedTrack%pos%xyz%r(:) = TrackBefMove%pos%xyz%r(:)
!     *                      + TrackBefMove.vec.w.r(1) * leng
     *                      + dir%r(:) * leng

       call csetPos(MovedTrack%pos) ! set height, radial length, thickness
!     get new cos


!        cgetZenith  or next one ; which is better ?       
       MovedTrack%vec%coszenith = cnewcos(MovedTrack%pos%radiallen,
     *  TrackBefMove%vec%coszenith, leng)
       if(TimeStructure ) then
          call cgetBeta(MovedTrack%p, Beta)  
!             note righthand is not MovedTrack.t
          if(Beta .gt. 0.) then
             MovedTrack%t = TrackBefMove%t + leng/Beta
          endif
       endif
       end
!      -----------------------------------
       subroutine cputEnergyLoss
!     observation layer.
       use modEMcontrol
       use modSetIntInf
       implicit none

#include  "Zcode.h"
#include  "Ztrack.h"
#include  "Ztrackv.h"
#include  "Ztrackp.h"
#include  "ZmediaLoft.h" 
!     
       real*8  rho, cvh2den, dedt, dedtF,  thick, leng, cupsilon, Eloss
       save dedt, dedtF
       real*8  csyncTELoss
       integer jcut
       real*8  dedx, dedxF  ! output.  dedt, dedtF is put


       real*8  dedxmu, dedxmuF

!
!          first consider the synchrotron loss.
!
       EnergyLoss = 0.
       if( TrackBefMove%p%fm%p(4) .gt. MagBremEmin) then
          if( MagBrem .eq. 1 .and. TrackBefMove%p%code .eq. kelec) then
                            !  Cosmos upsilon is used. Mag is Zmagfield.h
	    Upsilon = cupsilon(TrackBefMove%p, Mag)
	    if(Upsilon .gt. UpsilonMin) then
!                 compute energy loss due to sychrotron rad.
                 EnergyLoss = csyncTELoss(Upsilon) * 
     *           IntInfArray(ProcessNo)%length
                 if(Reverse .eq. 0) then
                    MovedTrack%p%fm%p(4) = 
     *                MovedTrack%p%fm%p(4) -  EnergyLoss
                 elseif(Reverse .eq. 2) then
                    MovedTrack%p%fm%p(4) = 
     *                MovedTrack%p%fm%p(4) +  EnergyLoss
                 endif
!                 don't  worry about death
	    endif
         endif
       endif 
!            
!           dE/dX
!
       rho = cvh2den(TrackBefMove%pos%height)

       Media(MediaNo)%rhoc = rho/Media(MediaNo)%rho
!     rhoc is used to get rho inside epdeedxe for the dendity
!     effect in dE/dx
       if( TrackBefMove%p%code == kelec) then
          call epdedxe(Media(MediaNo), TrackBefMove%p, dedt, dedtF)
       else
          call epdedxNone(Media(MediaNo), TrackBefMove%p, dedt, dedtF)
       endif
                ! in GeV/(g/cm2). we need GeV/(kg/m2) = GeV/(1000g/10000cm2)
!        =>10 GeV/(g/cm2).  i.e,  1 GeV/(g/cm2) = GeV/(kg/m2) / 10.
       dedt= dedt/10.
       dedtF = dedtF/10.
       
!      call cdedxInAir(TrackBefMove%p, rho, dedt, dedtF )  ! dedt; GeV/(kg/m2)
       if(TrackBefMove%p%code .eq. kmuon ) then

!     dE/dx due to muon pair, brem, nuc.i
!     next is not important if we treat pair, brems, ni
!       as stacasting descrete process.          
          call epmudEdx(Media(MediaNo)%mu%MuNI,
     *       Media(MediaNo)%mu%MuBr,  Media(MediaNo)%mu%MuPr,
     *         Media(MediaNo), TrackBefMove%p%fm%p(4), dedxmu)
          
!          call cmudEdx(MuNI, MuBr, MuPr, TrackBefMove%p%fm%p(4),
!     *         dedxmu, dedxmuF)   ! dedxmu in GeV/(g/cm2)
          dedxmu = dedxmu /10.  !  GeV/(kg/m2)
!     dedxmuF = dedxmuF /10.  !  GeV/(kg/m2)
!     MuNI, etc: 0. the interaction is neglected
!          others.  v=Eemited/Emu.
!     dedexmu is calculated as           
!     1.     Emu  vds/dv  is integrated from vc to vmax 
!     2.                      //        from vc to vmin.  v>vmin is
!            treated as stacatics process but product is tracking is neglected
!     3      same as 2 and tracking of product is included.
!           
          dedt = dedt + dedxmu
!
!         dedtF = dedtF + dedxmuF      dedxmu  becomes  dedxmuF value 
!         only if  MuNI, MuBR, MuPr are all 1. but such setting is
!     only for test puarpose and not  good sice no fluctuation is included.
!     In versions before 9, corresponding values for MuNI etc= 1 are alwaays
!     computed and given in dedxmF, together with dedxmu
       endif

       ELoss = dedt * IntInfArray(ProcessNo)%thickness
       EnergyLoss = EnergyLoss + Eloss
       if(Reverse .eq. 0) then
           MovedTrack%p%fm%p(4) = MovedTrack%p%fm%p(4) -  ELoss          
!           see if <=mass
           if(MovedTrack%p%fm%p(4) .lt. MovedTrack%p%mass) then
              EnergyLoss =( TrackBefMove%p%fm%p(4) - MovedTrack%p%mass )
              if(dedt .gt. 0.) then
                 thick =EnergyLoss / dedt
              else
                 thick = 0.
              endif
              MovedTrack%p%fm%p(4) = MovedTrack%p%mass

!              call cthick2len(TrackBefMove.pos.height, 
!     *          TrackBefMove.vec.coszenith, thick, leng, thick, jcut)
!                   Mar.18
              call cthick2len(TrackBefMove, thick, leng, thick, jcut)
!               reset position information in MovedTrack
              call cmoveStreight(leng, TrackBefMove%vec%w)
              MoveStat = Truncated
           endif
!                reset 3 momenta px, py, pz 
!                 (but assume direction is unchanged)
           call ce2p(MovedTrack)
       elseif(Reverse .eq. 2) then
          MovedTrack%p%fm%p(4) = MovedTrack%p%fm%p(4) +  EnergyLoss
!                reset 3 momenta px, py, pz 
!                 (but assume direction is unchanged)
          call ce2p(MovedTrack)
!          MoveStat = Truncated
       else
!          MoveStat = Truncated
       endif 
       return
!      **************
       entry cqElossRate(dedx, dedxF)
!      **************
!         inquire the dedt; which may be used when particle cross 
!     an obs. level and recompute the energy
!     New Cosmos version 10. dedxF should not be used.
!     For muon case,  it does not include dedexF of muons.
!     For electrons,  dedxF is only a measure of energy loss including
!     high energy delta rays without fluctuation.
       dedx = dedt 
       dedxF = dedtF
       end
!     --------------------------------------------------------------
!    !   #ifdef XXXXXX
       subroutine cputDeflection
       implicit none

#include  "Ztrack.h"
#include  "Ztrackp.h"
#include  "Ztrackv.h"
#include  "Zelemagp.h"

!
       real*8 dt, dispx, dispy
       logical nodef
       if(Reverse .eq. 0) then
          nodef = OneDim .ne. 0  .or.
     *    (mod(HowGeomag, 10) .eq. 1 .and. Zfirst%pos%depth .eq. 0.)
       else
          nodef = .false.
       endif
       if(.not. nodef) then
!
!         compute   magnetic deflection first. independently of scattering.
!     system is xyz;  if Efield exists, together with it
!!!!!!!!!!!!!!
          write(0,*) ' enter cmagDef'
!!!!!!!!!!!!!!!!!!!          
          call cmagDef
       endif
!!!!!!!!!!!
       write(0,*)  ' basck from cmagDef'
!!!!!!!!!!!       
       if(.not. nodef .and. Reverse .eq. 0) then
!     this is for multiple scattering
!!!!!!!!!!!!
          write(0,*) ' enter celecDef'
!!!!!!!!!!          
!!!!!!!!!!!!   ***** test commentout       call celecDef(dispx, dispy)  ! effect alrady put in
! MovedTrack. dispx, y are dummy
!!!!!!!!!!
          write(0,*) ' back from celecDef'
!!!!!!!!!!!          
       endif
       if(.not. nodef) then
          call csetPos(MovedTrack%pos)
          call cgetZenith(MovedTrack, MovedTrack%vec%coszenith)
!          reset momentum to be compatible with direction cos.
          call cresetMom(MovedTrack)
!
          if(TimeStructure .and. Reverse .eq. 0) then  ! only for mul. scat
!                compute excess path length to be added to streight path
             call cexcessLen(dispx, dispy, dt)
             if(Beta .gt. 0.) then
                MovedTrack%t = MovedTrack%t + dt/Beta
             endif
          endif
       endif
       end
!      *************************
#include "ZsubstRec.h"
       subroutine cmagDef
       use modEfield
       use modSetIntInf
       implicit none
#include  "Zobs.h"
#include  "Zobsp.h"
#include  "Ztrack.h"
#include  "Ztrackp.h"
#include  "Ztrackv.h"
#include  "Zelemagp.h"

!

!              by  Geomag  (dr and ddirec)
       type(coord)::dispmr, dispmd
       type(magfield)::tempmag

       real*8 temp1, temp2, leng
!cc    *   ,norm
       integer i

       integer icon
       type(coord):: middle
       logical usemiddle, high, usecmiddle
       real*8 cheight
       data cheight/30.d3/  
       real*8 newmom(3), dE, E1, E2, p1sq, p2sq, newE, norm
       save cheight
       real(8):: dEout, lengout
       real(8):: dEsave, lengsave
       leng = IntInfArray(ProcessNo)%length
!!!!!!!!!!!!!!!
       write(0,*) ' HowEfield ',HowEfield
       write(0,*) ' leng =',leng
!!!!!!!!!!       
       
       if( HowEfield >= 1 ) then
          call cdefByMagAndE(TrackBefMove,  leng,  dispmr, dispmd,
     *     newmom)
          p1sq =dot_product(TrackBefMove%p%fm%p(1:3),
     *           TrackBefMove%p%fm%p(1:3))
          E1 = sqrt( p1sq + TrackBefMove%p%mass**2)
          p2sq =dot_product(newmom(1:3),newmom(1:3))
          E2=sqrt(p2sq +  TrackBefMove%p%mass**2)
          dE=E2-E1  !  may be + or -  depending on Ef and charge
          ! dE/dx loss has been put in MovedTrack. adjust it again
          
          newE = MovedTrack%p%fm%p(4) + dE
          if( newE <  TrackBefMove%p%mass ) then
             newE =  TrackBefMove%p%mass
             newmom(:) = 0
             dE = MovedTrack%p%fm%p(4) -MovedTrack%p%mass 
          else
             ! keeping the direction,adjust new momentum to be
             ! consistent with newE
             if( p2sq <= 0. ) then
                write(0,*) ' p2sq = ',p2sq
                write(0,*)
     *          ' TrackBefMove%p=',TrackBefMove%p%fm%p(1:4)
                write(0,*) ' dE, newE=',dE, newE
                write(0,*) ' leng =', leng
                write(0,*) 'dispmr=',dispmr%r(:)
                write(0,*) 'dispmd=',dispmd%r(:)
                write(0,*)
     *            'code, chg=',TrackBefMove%p%code,
     *            TrackBefMove%p%charge

                newE = TrackBefMove%p%mass
                newmom(:) = 0
                dispmd%r(:) = (/0.,0.,1./)
             else
                norm = sqrt( (newE**2 - TrackBefMove%p%mass**2)/p2sq)
                newmom(:) = newmom(:) *norm
             endif
          endif

          dEsave =dE
          lengsave = leng

          
          MovedTrack%p%fm%p(4) = newE
!          note. dispmr is r(new)-r(old) vector and different
!          from other routines below. 
!           dispmd is new dir and set at 100
          MovedTrack%pos%xyz%r(:) =  TrackBefMove%pos%xyz%r(:) +
     *         dispmr%r(:) 
          MovedTrack%p%fm%p(1:3) = newmom(:)
          goto 100
       endif
!
!            UseRungeKutta         Height       
!                0                  any         Mag and cmagneticDef

!                1                 >cheight     middle Mag and cmagneticDef
!                1                 <            Mag and cmagneticDef

!                2                 >            middle Mag and cmagneticDef 
!                                               or cbDefByRK2        
!                2                 <            Mag and cmagneticDef 

!                3                 >            middle Mag and cmagneticDef 
!                                               or cbDefByRK        
!                3                 <            Mag and cmagneticDef 

!                4                 >            use Mag at curved middle point
!                                               and estimate end point by 
!                                               cmagneticDef or cbDefByRK2
!                4                 <            Mag and cmagneticDef   

!                5                 >            use Mag at curved middle point
!                                               and estimate end point by 
!                                               cmagneticDef or cbDefByRK
!                5                 <            Mag and cmagneticDef   

!                6                 >            cbDefByRK2
!                6                 <            Mag and cmagneticDef 

!                7                 >            cbDefByRK
!                7                 <            Mag and cmagneticDef 
!
!                8      at any height           cbDefUser ; interface is
!                                               same as cbDefByRK
!!!!!!!!!!!!!!!!!
      write(0,*) ' UseRungeKutta =', UseRungeKutta
!!!!       
       if(UseRungeKutta .eq. 8 ) then
          call cbDefUser(TrackBefMove,  leng,  dispmr, dispmd,
     *    MovedTrack%p%fm%p(1:3))
#if defined SUBSTREC             
             MovedTrack%pos%xyz = dispmr ! this is not a dispalcement vector
#else 
             call csubstcoord(dispmr, MovedTrack%pos%xyz)
#endif
          goto 100
       endif
       high = TrackBefMove%pos%height .gt. cheight
       usemiddle = ( UseRungeKutta .ge. 1 .and.
     *               UseRungeKutta .le. 3 .and.
     *               high )
       usecmiddle = ( UseRungeKutta .ge. 4 .and.
     *               UseRungeKutta .le. 5 .and.
     *               high )

!!!!!!!!!
       write(0,*)  ' high =', high
       write(0,*) ' usemddle, usecmidl=',   usemiddle,  usecmiddle
!!!!!!!!!!!       
       if(usemiddle) then
          do i = 1, 3
             middle%r(i) = TrackBefMove%pos%xyz%r(i)+ 
     *          leng/2 * TrackBefMove%vec%w%r(i)
          enddo
          middle%sys='xyz'

          call cgeomag(YearOfGeomag,  middle,
     *                tempmag, icon)
          call ctransMagTo('xyz', middle,
     *        tempmag, tempmag)
       elseif(usecmiddle) then
!            get curved middle point
          call cmagneticDef(TrackBefMove, Mag, leng/2.0d0,
     *     dispmr, dispmd)  
          do i = 1, 3
             middle%r(i) = TrackBefMove%pos%xyz%r(i) + dispmr%r(i)
          enddo
          middle%sys ='xyz'
          call cgeomag(YearOfGeomag,  middle,
     *                tempmag, icon)
          call ctransMagTo('xyz', middle,
     *        tempmag, tempmag)
       endif

       if(usemiddle .and. UseRungeKutta .eq. 1) then
          Mag = tempmag
       elseif( usemiddle .and.  UseRungeKutta .eq. 2 ) then
          temp1 = tempmag%x**2 + tempmag%y**2 + tempmag%z**2
          temp2 = Mag%x**2 + Mag%y**2 + Mag%z**2
!           if  dB/B~ dB^2/B^2/2 > 0.4 %, use RungeKutta.
          if( abs(temp1-temp2)/temp1/2.0 .gt. 4.0d-3) then
             call cbDefByRK2(TrackBefMove,  leng, dispmr, dispmd,
     *         MovedTrack%p%fm%p(1:3))
#if defined SUBSTREC             
             MovedTrack%pos%xyz = dispmr ! this is not a dispalcement vector
#else 
             call csubstcoord(dispmr, MovedTrack%pos%xyz)
#endif
             goto 100
          else
             Mag = tempmag
          endif
       elseif(usemiddle .and. UseRungeKutta .eq. 3) then
          temp1 = tempmag%x**2 + tempmag%y**2 + tempmag%z**2
          temp2 = Mag%x**2 + Mag%y**2 + Mag%z**2
!           if  dB/B~ dB^2/B^2/2 > 0.45 %, use RungeKutta.
!                  next 4.5d-3 is very sensitive to cpu time
!                 if it was 4.0d-3, cpu time becomes twice.
          if( abs(temp1-temp2)/temp1/2.0 .gt. 4.5d-3) then
             call cbDefByRK(TrackBefMove,  leng, dispmr, dispmd,
     *         MovedTrack%p%fm%p(1:3))
#if defined SUBSTREC             
             MovedTrack%pos%xyz = dispmr ! this is not a dispalcement vector
#else 
             call csubstcoord(dispmr, MovedTrack%pos%xyz)
#endif
             goto 100
          else
             Mag = tempmag
          endif
       elseif(usecmiddle .and. UseRungeKutta .eq. 4) then
          temp1 = tempmag%x**2 + tempmag%y**2 + tempmag%z**2
          temp2 = Mag%x**2 + Mag%y**2 + Mag%z**2
!           if  dB/B~ dB^2/B^2/2 > 0.1 %, use RungeKutta.
          if( abs(temp1-temp2)/temp1/2.0 .gt. 4.0d-3) then
             call cbDefByRK2(TrackBefMove,  leng, dispmr, dispmd,
     *         MovedTrack%p%fm%p(1:3))
#if defined SUBSTREC             
             MovedTrack%pos%xyz = dispmr ! this is not a dispalcement vector
#else 
             call csubstcoord(dispmr, MovedTrack%pos%xyz)
#endif
             goto 100
          else
             Mag = tempmag
          endif
       elseif(usecmiddle .and. UseRungeKutta .eq. 5) then
          temp1 = tempmag%x**2 + tempmag%y**2 + tempmag%z**2
          temp2 = Mag%x**2 + Mag%y**2 + Mag%z**2
!           if  dB/B~ dB^2/B^2/2 > 0.45 %, use RungeKutta.
          if( abs(temp1-temp2)/temp1/2.0 .gt. 4.5d-3) then
             call cbDefByRK(TrackBefMove,  leng, dispmr, dispmd,
     *         MovedTrack%p%fm%p(1:3))
#if defined SUBSTREC             
             MovedTrack%pos%xyz = dispmr ! this is not a dispalcement vector
#else 
             call csubstcoord(dispmr, MovedTrack%pos%xyz)
#endif
             goto 100
          else
             Mag = tempmag
          endif
       endif
!!!!!!   default comes here
       if(UseRungeKutta .le. 5 .or. .not. high ) then
!!!!!!!!!!
          write(0,*) ' enter  cmagneticDef'
!!!!!!!!!!!!!          
          call cmagneticDef(TrackBefMove, Mag, leng,
     *     dispmr, dispmd)
!!!!!!!!!!
          write(0,*) ' back  cmagneticDef'
          write(0,*) 'dispmr, dispmd ',  dispmr, dispmd
!!!!!!!!!!!!!          

          do i = 1,  3
             MovedTrack%pos%xyz%r(i) = TrackBefMove%pos%xyz%r(i) +
     *         dispmr%r(i)
          enddo
       elseif(UseRungeKutta .eq. 6) then
          call cbDefByRK2(TrackBefMove,  leng, dispmr, dispmd,
     *         MovedTrack%p%fm%p(1:3))
#if defined SUBSTREC
          MovedTrack%pos%xyz = dispmr  ! this is not a dispalcement vector
#else 
          call csubstcoord(dispmr, MovedTrack%pos%xyz)
#endif
                                      !  but the vector itself.
       else
          call cbDefByRK(TrackBefMove,  leng, dispmr, dispmd,
     *         MovedTrack%p%fm%p(1:3))
#if defined SUBSTREC
          MovedTrack%pos%xyz = dispmr  ! this is not a dispalcement vector
                                      !  but the vector itself.
#else 
          call csubstcoord(dispmr, MovedTrack%pos%xyz)
#endif
       endif
 100   continue
#if defined SUBSTREC
       MovedTrack%vec%w = dispmd
#else
       call csubstcoord(dispmd, MovedTrack%vec%w)
#endif
!!!!!!!!!!!!!!!
       write(0,*) ' before entry'
!!!!!!!!!!!!!!
       return
!!!!!!!!!!!!next entry why ??       
       entry cqdEbyEfield(dEout, lengout)
       dEout = dEsave
       lengout = lengsave
       end
!        for IBM AIX
      subroutine csubstcoord(c1, c2)
      implicit none
#include "Zcoord.h"
      type(coord)::c1, c2   
      c2 = c1
      end

!      ==============================================================
      subroutine celecDef(dx, dy)
      use modSetIntInf
      implicit none

#include  "Ztrack.h"
#include  "Ztrackp.h"
#include  "Ztrackv.h"
#include  "Zelemagp.h"
#include  "ZmediaLoft.h"
      real*8 dx, dy             ! output. scatttering displacement
!     
      type(coord)::dircos
!              by  Multiple Scattering
      type(coord)::dsa   ! dire ccos of scattering angle
      type(coord)::w
      real*8  sint, cs, sn, tmp, avx, avy, disp, dt, dl
      real*8 r, g1, g2, gf1, gf2, beta2, tetarms
      real*8 theta
      
!     call cmulScat(theta)
      call epmulScat(TrackBefMove%p, MovedTrack%p,
     *     IntInfArray(ProcessNo)%length*100.d0,    ! cm
     *     IntInfArray(ProcessNo)%thickness*0.1, ! g/cm2
     *     Media(MediaNo), theta)
      
      if(theta .lt. 0.01d0) then
!                 cos
!     dsa%z = 1.-theta**2/2
         dsa%r(3) = 1.-theta**2/2
         sint = theta
      else
!     dsa%z = cos(theta)
         dsa%r(3) = cos(theta)
         sint = sin(theta)
      endif
!        azimuthal angle
      call kcossn(cs, sn)
!     dsa%x = sint * cs
!     dsa%y = sint * sn
!     dsa%z = cos(theta)
      dsa%r(1) = sint * cs
      dsa%r(2) = sint * sn
      dsa%r(3) = cos(theta)    
      dt = IntInfArray(ProcessNo)%thickness/Media(MediaNo)%X0 ! r%l
      dl = IntInfArray(ProcessNo)%length !  m

!      if(Moliere .ge. 0) then
!     if( ALateCor ) then       ! >v7.645
!         v>=7.655      
      if( (Moliere == 0 .and. AlateCor>=1) .or. 
     *     AlateCor == 2) then
!     ALateCor D=1.  if Gaussian (Moliere=0), correlation is taken
!     into account.
!     If 0, no correlation is considered.
!        2, correlation based on Gaussian formula is forced even if
!           Moliere is not 0.
         
!         sample displacement correlated to theta
!                 this is the same as P.D.B though look like
!                diff.
        tmp = dl/2.d0
        avx = tmp * dsa%r(1)
        avy = tmp * dsa%r(2)
!                dispersion
        gf1 = TrackBefMove%p%fm%p(4)/MovedTrack%p%mass
        gf2 = MovedTrack%p%fm%p(4)/MovedTrack%p%mass
        beta2 = 1.d0 - 1.d0/gf1/gf2
        if(beta2 .le. 0.) then
           disp = 0.d0
        else
           if(dt .gt. 1.d-3) then
              tetarms = Es/TrackBefMove%p%fm%p(4)*
     *            abs(MovedTrack%p%charge)*
     *            sqrt(dt)*(1.0 + 0.038*log(dt))
           else
              tetarms = Es/TrackBefMove%p%fm%p(4)*
     *                  abs(MovedTrack%p%charge)*
     *              sqrt(dt)
           endif
           disp=tetarms/sqrt(6.d0*beta2)*dl/2.d0
!               sample 2 independent gaussian variables
!             with mean 0 and var 1
        endif
        call kgauss2(0.d0, 1.0d0, g1, g2)
        dx = g1 * disp + avx
        dy = g2 * disp + avy
!                  displacement
        r=sqrt(dx*dx+dy*dy)     ! in m
!              direction cos of vector r in original sys.
        if(r .ne. 0.) then
!           w%x = dx/r
!           w%y = dy/r
!           w%z = 0.
           w%r(1) = dx/r
           w%r(2) = dy/r
           w%r(3) = 0.d0
           
!                 transform wx,wy,wz to original sys.
!                    TrackBefMove is better
           call ctransVectZ(TrackBefMove%vec%w, w, w)
!               r is already in m.
!              add scattering effect.
!              r*w is displacement by scattering
!           MovedTrack%pos%xyz%x = r*w%x + MovedTrack%pos%xyz%x
!           MovedTrack%pos%xyz%y = r*w%y + MovedTrack%pos%xyz%y
!           MovedTrack%pos%xyz%z = r*w%z + MovedTrack%pos%xyz%z
           MovedTrack%pos%xyz%r(1:3) =
     *         r*w%r(1:3)+ MovedTrack%pos%xyz%r(1:3)
        endif
      else
         dx = 0.
         dy = 0.
      endif
!        convert scattering angle at end of path to
!        original system . MovedTrack is better since
!        mag. def is contained there already.
      call ctransVectZ(MovedTrack%vec%w, dsa,
     *      MovedTrack%vec%w)

      end
!!!   #endif

      
      subroutine cmanageMediaChange(tobecut, nodeX, Rx, dir,len)
!
!     This is called when TrackBefMove is moved to MovedTrack
!     where some interaction is going to occure or the path is
!     too long and cut there.  This sub. is to check if MovedTrack
!     already enters into a new medium part or not. If it enters into the
!     new medium, the path should be cut there.  
!    
      use modAtmosDef
      implicit  none
#include "Ztrack.h"
#include "Ztrackv.h"      
!  #include "Zatmos.h"

      integer,intent(out):: tobecut ! == 0   ==> to be cut at Rx, nodeX.
!     /= 0   ==> not cross with diff. medium, so dont worry; use
!     current MovedTrack.  nodex, Rx, len are undefined.
      integer,intent(out):: nodeX ! the cut position is inside the entrance
!     of nodeX-th layer.  (if down going, at the top, if going up
!     at the bottom.
      real(8),intent(out):: Rx(3) !  x-sing point. should be in the
!     nodeX node.
      real(8),intent(out):: dir(3) ! diection cose of TracBefMove to
            !       moved tracke
      real(8),intent(out):: len ! TraceBefMove-to-Rx length in m.
      

      integer:: node, mno, nodeB, nodeU
      integer:: icon
      real(8):: vh
      real(8):: Hx

      tobecut = 1  ! most probably not enter into new media
!        get height of TrackBefMove
      vh = TrackBefMove%pos%height
!        get its node
      call kdwhereis(vh, atmos%nodes, atmos%z, 1,  node)
!     see if there is  diff medium from the one  at node.
!     Get nearest one below: nodeB and upper: nodeU
!      if 0, non existent      
      call cnearestDiffMedia(node, nodeB, nodeU)
      if( TrackBefMove%vec%coszenith > 0.) then
!             down going
         if( nodeB > 0 ) then
!            node and nodeB have diff. Media.  if nodeB is closer to the
!            Moved pos, cut the  path  at the entrance of nodeB. so see
!           where is xing point with nodeB
            Hx= atmos%z(nodeB) + Eradius  ! radius of candidate node
            call cgetNodeXP(TrackBefMove%pos%xyz%r(:),
     *           TrackBefMove%pos%radiallen, MovedTrack%pos%xyz%r(:),
     *           MovedTrack%pos%radiallen, Hx, Rx(:), len, tobecut)
            if( tobecut == 0 ) then
               nodeX = nodeB
            endif
         elseif( nodeU > 0 ) then
!          rare case but possible. down doing but might cross with
!          an upper node.  so see the top of one node below the candidate
            Hx= atmos%z(nodeU+1) + Eradius ! 
            call cgetNodeXP(TrackBefMove%pos%xyz%r(:),
     *           TrackBefMove%pos%radiallen, MovedTrack%pos%xyz%r(:),
     *           MovedTrack%pos%radiallen, Hx, Rx(:), len, tobecut)
            if( tobecut == 0 ) then
               nodeX = nodeU
            endif
         endif
      else  ! upgoint or horizontal
         if( nodeU > 0 ) then
!            seems ptcl goes upper node
            Hx= atmos%z(nodeU+1) + Eradius ! radius of bottom of the
!                            candidate node
            call cgetNodeXP(TrackBefMove%pos%xyz%r(:),
     *           TrackBefMove%pos%radiallen, MovedTrack%pos%xyz%r(:),
     *           MovedTrack%pos%radiallen, Hx, Rx(:), len, tobecut)
            !  Rx and len are adjusted to be in nodeU 
            if( tobecut == 0 ) then
               nodeX = nodeU
            endif
         elseif( nodeB > 0 ) then
            ! no possibilty but check
            Hx= atmos%z(nodeB) + Eradius  ! radius of candidate node
            call cgetNodeXP(TrackBefMove%pos%xyz%r(:),
     *           TrackBefMove%pos%radiallen, MovedTrack%pos%xyz%r(:),
     *           MovedTrack%pos%radiallen, Hx, Rx(:), len, tobecut)
            if( tobecut == 0 ) then
               nodeX = nodeB
               write(0,*) ' strange diff. medium for upgoing ptcl'
               write(0,*) ' it is below current node '
               write(0,*) ' TrackBef ', TrackBefMove%pos%xyz%r(:)
               write(0,*) ' Movedf ',MovedTrack%pos%xyz%r(:)
               write(0,*) ' len, Rx ', len, Rx(:)
               stop
            endif
         endif
      endif
      if( tobecut == 0 ) then
         if( len > 1.d-2 ) then
            dir(:) = ( Rx(:)-TrackBefMove%pos%xyz%r(:)) /len
         else
            dir(:)= TrackBefMove%vec%w%r(:)
         endif
      endif
      end  subroutine cmanageMediaChange

!!      program main  ! test cgetNodeXP(R1, H1, R2, H2, Hx,  Rx, len, icon)
!!      implicit none
!!      real(8):: R1(3),  R2(3), Rx(3), len
!!      integer:: icon
!!      real(8),parameter:: pi=asin(1.d0)*2
!!      real(8),parameter::
!!     *       H1 = 6800.d3, teta1=35.0d0*pi/180.d0, fai1=50.0d0*pi/180.d0
!!      real(8),parameter::
!!!     *       H2 = 6750.d3, teta2=35.0d0*pi/180.d0, fai2=50.0d0*pi/180.d0
!!     *       H2 = 6750.d3, teta2=37.0d0*pi/180.d0, fai2=53.0d0*pi/180.d0
!!      real(8)::Hx = 6760.0d3
!!
!!      call cgetNodeXPsetEps(.true.,1.d-7)
!!      
!!      R1(1)=H1 * sin(teta1) *cos(fai1)
!!      R1(2)=H1 * sin(teta1) *sin(fai1)
!!      R1(3)=H1 * cos(teta1)
!!      
!!      R2(1)=H2 * sin(teta2) *cos(fai2)
!!      R2(2)=H2 * sin(teta2) *sin(fai2)
!!      R2(3) = H2* cos(teta2)
!!
!!      call  cgetNodeXP(R1, H1, R2, H2, Hx,  Rx, len, icon)
!!      write(0,*) R1(:)
!!      write(0,*) R2(:)     
!!
!!      write(0,*) ' icon=', icon
!!      write(0,*) ' len, Rx=', len, Rx(:)
!!      write(0,*) ' |Rx|=',sqrt( dot_product(Rx(:), Rx(:)))
!!      write(0,*) 'H2=',H2
!!      end
      
      subroutine cgetNodeXP(R1, H1, R2, H2, Hx,  Rx, len, icon)
      implicit none
!     Earth x-y-z system
      real(8),intent(in):: R1(3) ! postion vector  (abs(R1) = R0 + h1= H1)
                                 !  TrackBefMove%pos%xyz%r(3)
      real(8),intent(in):: H1   !  R0 + h1 
      real(8),intent(in):: R2(3) ! MovedTrack%vecxyz%r(3)
      real(8),intent(in):: H2   !  R0 + h2
      real(8),intent(in):: Hx  !  R0+hx of the candidate node boundary
!          the node is the nearest node of which medium is diff from current medium
!          (if exist)
      real(8),intent(out):: Rx(3) !  vector of crossing point at (modified) Hx.
      real(8),intent(out):: len !  length from R1 to Rx.  0< len <= 1
!       is first obtaiend and to force the  xp is in  side the new media
!       acutual len could be > 1.0 (little bit)
      integer,intent(out):: icon !  if Rx is obtained 0. else 1.

!     ===================
!     Rx(:) = R1(:) + alfa(R2(:)-R1(:)) is  obtained with some redundant factor so
!     that it is surely  in the new medium but very close to the boundary. For this we use
!     two parameters modifiable by call cgetNodeXPsetEps(rel, eps)
!     First we get alfa and Rx as accurately as  possible. Even
!     there exists a crossing point, the position may or may not be in
!     the new medium region. For our purse, it's better it is inside the
!     new medium region but very close to the bundary.
!   
      logical,intent(in):: rel  !if T, original alfa is modified to be
!          alfa=min((1+eps)*alfa,1) and len(gth) from R1 to Rx is bit
!          enlarged. len= |alfa(R2-R1)|
!     If F, len(gth) obtained with original alfa is modified to be len=len+eps
  
      real(8),intent(in):: eps !  see above
!     Defaults are rel=F, eps=5.0d-2 (5cm).
!      
      logical,save:: relSave=.false. ! epsSave is  absolute value
      real(8),save:: epsSave=5.0d-2   !in m  5cm. 

!       
!
      real(8)::  R1sq, R2sq, R1R2, R(3), dR(3)
      real(8):: alfa, a, b, c, D, sd
!     get alfa for Rx(:) - R1(:) = alfa*(R2(:)-R1(:)) to get
!         Rx(:).       
!     where   Rx^2 = (R0+hx)^2= Hx^2
!     Rx(:)= R1(:) + alfa*(R2(:)-R1(:)); define  R(:)=R2(:)-R1(:)
!     Rx^2  = R1^2 + 2 alfa * R1(:).R(:) + alfa^2 R^2
!        a=R^2  b= 2R1(:).R(:)  c= R1^2-Hx^2

      
      R1sq = dot_product( R1(:), R1(:) )
      R2sq = dot_product( R2(:), R2(:) )
      R1R2 = dot_product( R1(:), R2(:) )
      R(:) = R2(:) - R1(:)
      a  = dot_product( R(:), R(:) )
      write(0,*) ' a, R1sq=', a, R1sq
      write(0,*) ' R2sq, H2sq=',R2sq, H2**2
      b =  R1R2 - H1**2
             
      c =  (H1 -Hx) *( H1 + Hx)
!!!!
!      write(0,*) 'a,b,c=',a,b,c
!!!!!!      
!            a*x^2 + 2bx + c =0;  
      D = b**2 - a*c
!!!!!!
!      write(0,*) 'D=',D
!!!!
      if (D > 0 ) then
!     skim  is neglected  ; take nearer one (positive L)
!    only 0< alfa <1 is obtained         
         sd = sqrt(D)
!!!!!!
!         write(0,*) ' -b =', -b, ' sd=',  sd, ' -b-sd=',-b-sd
!!!!!!!!         
         if( -b > sd) then
            alfa = (-b - sd)/a
!!!!!!
!            write(0,*) ' alfa=', alfa
!!!!!!!!!            
            if( alfa < 1.0 ) then

               if( relSave ) then
                  dR(:) = min(alfa*(1.d0+epsSave), 1.d0)*R(:)
                  len = sqrt( dot_product(dR(:), dR(:)) )
               else
                  dR(:) =  alfa*R(:)
                  len = sqrt( dot_product(dR(:), dR(:)) )
                  dR(:) = alfa* R(:)*(len+epsSave)/len
                  len = len + epsSave
               endif
               Rx(:) = R1(:) + dR(:)
               icon = 0
            else
               icon = 1
            endif
         else
            icon = 1
         endif
      else
         icon = 1
      endif
      return

      entry cgetNodeXPsetEps(rel, eps)
      epsSave = eps
      relSave = rel
      end   subroutine cgetNodeXP
      
      subroutine ccutAndAdjust(dir, len)
      use modSetIntInf
      implicit none
#include "Ztrack.h"      
#include "Ztrackv.h"
      real(8),intent(in):: dir(3) ! direction cos in E-xyz
      real(8),intent(in):: len  !  length in m


      real(8):: leninkgm2, dedx, dedxF,  lenm
      real(8),external:: clen2thick
!        TrackBefMove direction could be diff. from dir so
!        update it
      TrackBefMove%vec%w%r(:) = dir(:) 
      call cgetZenith(TrackBefMove, TrackBefMove%vec%coszenith)
      
      call cmoveStreight(len,  TrackBefMove%vec%w)
      MoveStat = Truncated
       
      leninkgm2 = clen2thick(TrackBefMove%pos%height,
     *     TrackBefMove%vec%coszenith, len)
       ! energy loss rate has been computed in cputEloss; so get it
!      if( HowEfield >= 1) then
!         call cqdEbyEfield(ELoss, lenm)
!         EnergyLoss = ELoss*len/lenm   ! maybe >0 or <0
!      else
!         EnergyLoss = 0
!      endif
!     call cqElossRate(dedx, dedxF)
!       adjust EnergyLoss prop. to the length.         
!     EnergyLoss= dedx * leninkgm2  + EnergyLoss
      EnergyLoss = EnergyLoss * len/IntInfArray(ProcessNo)%length
!        update length      
      IntInfArray(ProcessNo)%length = len
      IntInfArray(ProcessNo)%thickness = leninkgm2
      
      MovedTrack%p%fm%p(4) = TrackBefMove%p%fm%p(4) -  EnergyLoss
!     see if <=mass;  should not happen because
!     path is shorter than to MovedTrack
      if(MovedTrack%p%fm%p(4) .lt. MovedTrack%p%mass) then
         call checkstat('from ccutAndAdjust')
         EnergyLoss =( TrackBefMove%p%fm%p(4) - MovedTrack%p%mass )
         MovedTrack%p%fm%p(4) = MovedTrack%p%mass
         MovedTrack%p%fm%p(:)=(/0.,0.,0.,0./)
         stop
      else
         call ce2p(MovedTrack)
      endif
      end
