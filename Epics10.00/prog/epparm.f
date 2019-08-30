      subroutine epprmr(dsn)
      use modV1ry
      use SampMolReducedA
      use modMCSparam
      implicit none
#include  "ZepTrackp.h"
#include  "ZepManager.h"
!
!             read epics file
!
            character*(*) dsn
            integer io,  icon, klena, i
            character*190 msg
            logical epgetParmN
            character*24 vname
            character*100 vvalue
            real*8 EpsLengX
            integer mediadirc
            
            mediadirc = 0 

            call copenf(iowk, dsn, icon)
            if(icon .ne. 0) then
               call cerrorMsg(dsn, 0)
               call cerrorMsg('above file not exits', 1)
            endif

!            read baisc parameters; find separater
            call afsep(iowk)
            do while( epgetParmN(iowk, vname, vvalue ) )
               if(vname .eq. 'AngleB') then
                  call arprml(vvalue, AngleB)
               elseif( vname .eq. 'DtMax' ) then
                  call arprmr(vvalue, DtMax)
               elseif( vname .eq. 'Eabsorb') then
                  call arprmi(vvalue,  Eabsorb)
               elseif( vname .eq. 'Eanihi' ) then
                  call arprmr(vvalue,  Eanihi)
               elseif( vname .eq. 'Ecut') then
                  call arprmr(vvalue, Ecut)
               elseif( vname .eq. 'EdepDedx' ) then
                  call arprml(vvalue, EdepdEdx)
               elseif( vname .eq. 'ElowerBndPair') then
                  call arprmr( vvalue, ElowerBndPair)
               elseif( vname .eq. 'EminElec' ) then
                  call arprmr(vvalue, EminElec)
               elseif( vname .eq. 'EminGamma') then
                  call arprmr(vvalue, EminGamma)
               elseif( vname .eq. 'AutoEmin') then
                  call arprmi(vvalue, AutoEmin)
               elseif( vname .eq. 'ModifyFile') then
                  call arprmc(vvalue, ModifyFile)
               elseif( vname .eq. 'EminH') then
                  call arprmr(vvalue, EminH)
               elseif( vname .eq. 'Es' ) then
                  call arprmr(vvalue, Escat)
               elseif( vname == 'Zcorrec') then
                  call arprml(vvalue, Zcorrec)
               elseif( vname .eq. 'EupperBndCS') then
                  call arprmr(vvalue, EupperBndCS)
               elseif( vname .eq. 'EpartialSC') then
                  call arprmr(vvalue, EpartialSC)
                  if( EpartialSC < 0.1d0 .or.
     *                EpartialSC > 10.d0 ) then
                     write(0,*) 'EpartialSC=',EpartialSC
                     write(0,*) ' invalid in epicsfiel'
                     write(0,*) ' must be 0.1~10 '
                     stop
                  endif        
               elseif( vname == 'TargetElecBrems' ) then
                  call arprmi(vvalue, TargetElecBrems)
                  if( TargetElecBrems < 0 .or. TargetElecBrems >7 ) then
                     write(0,*) 'TargetElecBrems=',TargetElecBrems,
     *                 ' invalid'
                     stop
                  endif
               elseif( vname .eq. 'Flpm') then
                  call arprmr(vvalue, Flpm)
                  if( Flpm < 1.0 ) then
                     write(0,*) 'Flpm=',Flpm, ' invalid'
                     write(0,*) 'must be >= 1 (<500)'
                     stop
                  elseif( Flpm > 500.) then
                     write(0,*) 'Flpm=',Flpm, ' too large'
                     write(0,*) ' normal range is 1~100'
                     write(0,*) '======================'
                     write(0,*) '======================'
                     write(0,*) 
     *               ' ***** Experimental purpose ?****'
                     write(0,*) '======================'
                     write(0,*) '======================'
                  endif
               elseif( vname .eq. 'HowNormBrems') then
                  call arprmi(vvalue, HowNormBrems)
               elseif( vname .eq. 'HowQuench') then
                  call arprmi(vvalue, HowQuench)
               elseif( vname .eq. 'IncGp') then
                  write(0,*) 'Sorry: IncGp is managed by param file'
                  write(0,*) 
     *       'If needed, use HowPhotoP(D=1)in &HPARAM part of param'
                  stop
               elseif( vname .eq. 'IoTrace') then
                  call arprmi(vvalue, IoTrace)
               elseif( vname .eq. 'Knckon' ) then
                  call arprml(vvalue, Knckon)
               elseif( vname .eq. 'Molier' .or.
     *                 vname == 'Moliere'  ) then  ! v9.170
                  call arprmc(vvalue, msg)
                  if( msg == 't' .or. msg == 'T') then
                     Moliere = 1
                  elseif( msg == 'f' .or. msg == 'F' ) then
                     Moliere = 0
                  else
                     call arprmi(vvalue, Moliere)
                  endif
               elseif( vname == 'MCSmodel' ) then
                  call arprmc(vvalue, MCSmodel)
               elseif( vname == 'MCSdir' ) then
                  call arprmc(vvalue, MCSdir)
               elseif( vname == 'MCSparam' ) then
                  call arprmc(vvalue, MCSparam)
!             elseif( vname == 'MCSxRange' ) then
!                  if(MCSxCond < MCSxyzCondmx) then
!                     MCSxCond = MCSxCond+1
!                     call arprmm(vvalue, MCSxRange(MCSxCond))
!                  else
!                     write(0,*) ' too many MCSxRange'
!                     stop
!                  endif
                  
!                  if(MCSyCond < MCSxyzCondmx) then
!                     MCSyCond = MCSyCond+1
!                     call arprmm(vvalue, MCSyRange(MCSyCond))
!                  else
!                     write(0,*) ' too many MCSyRange'
!                     stop
!                  endif
             elseif( vname == 'MCSzRange' ) then
                  if(MCSzCond < MCSxyzCondmx) then
                     MCSzCond = MCSzCond+1
                     call arprmm(vvalue, MCSzRange(MCSzCond))
                  else
                     write(0,*) ' too many MCSzRange'
                     stop
                  endif
                  
               elseif( vname == 'MCSnumRange' ) then
                  if(MCSnumCond < MCSnumCondmx) then
                     MCSnumCond = MCSnumCond+1
                     call arprmm(vvalue, MCSnumRange(MCSnumCond))
                  else
                     write(0,*) ' too many MCSnumRange'
                     stop
                  endif
               elseif( vname == 'MCSandor' ) then
                  call arprmc(vvalue, MCSandor )
                  if( MCSandor /= 'and' .and. MCSandor /= 'or' ) then
                     write(0,*) ' MCSandor value =',MCSandor
                     write(0,*) ' invalid. '
                     stop
                  endif
               elseif( vname == 'MCSrevert' ) then
                  call arprml(vvalue, MCSrevert)
               elseif( vname == 'MCSdebug' ) then
                  call arprml(vvalue, MCSdebug)
               elseif( vname == 'sbMin' ) then ! v9.201;don't use
                  write(0,*)
     *            'NOTE****: sbMin in epicsfile obsolete; not used'
!                  call arprmr(vvalue, sbMin)
               elseif( vname .eq. 'Photo') then
                  call arprml(vvalue, Photo)
               elseif( vname .eq. 'RecoilKeMin') then
                  call arprmr(vvalue, RecoilKEmin)
               elseif( vname .eq. 'Tcoef' ) then
                  call arprmr(vvalue, Tcoef)
               elseif( vname .eq. 'Tmin') then
                  call arprmr(vvalue, Tmin)
               elseif( vname .eq. 'Trace') then
                  if(.not. Trace) then
                     call arprml(vvalue, Trace)
                  endif
               elseif( vname .eq. 'TimeStruc') then
                  call arprml(vvalue, TimeStruc)
               elseif( vname .eq. 'MagField' ) then
                  call arprmi(vvalue, MagField)
               elseif( vname .eq. 'Bxu' ) then
                  call arprmr(vvalue, Bxu)
               elseif( vname .eq. 'Byu' ) then
                  call arprmr(vvalue, Byu)
               elseif( vname .eq. 'Bzu' ) then
                  call arprmr(vvalue, Bzu)
               elseif( vname .eq. 'ElecField') then
                  call arprmi(vvalue, ElecField)
               elseif( vname .eq. 'Exu') then
                  call arprmr(vvalue, Exu)
               elseif( vname .eq. 'Eyu') then
                  call arprmr(vvalue, Eyu)
               elseif( vname .eq. 'Ezu') then
                  call arprmr(vvalue, Ezu)
               elseif( vname .eq. 'FreeC' ) then
                  call arprml(vvalue, FreeC)
               elseif( vname == "FollowV1ry") then
                  call arprml(vvalue, FollowV1ry)
               elseif( vname .eq. 'EpsLeng') then
                  call arprmr(vvalue, EpsLengX)
                  if(EpsLengX < 0. ) then
                     EpsLeng = - EpsLengX
                  endif
               elseif( vname .eq. 'ALateCor') then
                  call arprmi(vvalue,ALateCor)
               elseif( vname .eq. 'Sync' ) then
                  call arprmi(vvalue, Sync)
               elseif( vname .eq. 'SyncLoop') then
                  call arprmr(vvalue, SyncLoop)
               elseif( vname .eq. 'MagPair' ) then
                  call arprmi(vvalue, MagPair)
               elseif( vname .eq. 'MuNI') then
                  call arprmi(vvalue, MuNI)
               elseif( vname .eq. 'MuBr') then
                  call arprmi(vvalue, MuBr)
               elseif( vname .eq. 'MuPr') then
                  call arprmi(vvalue, MuPr)
               elseif( vname .eq. 'KEmin' ) then
                  call arprmr(vvalue, KEmin)
               elseif( vname .eq. 'MsgLevel') then
                  call arprmi(vvalue, MsgLevel)
               elseif( vname .eq. 'MediaDir' ) then
                  mediadirc = mediadirc + 1
                  if(mediadirc .gt. MaxMediaDir) then
                     write(0,*) ' too many MediaDir'
                     stop 9789
                  endif
                  call arprmc(vvalue, MediaDir(mediadirc))
               elseif( vname .eq. 'TraceDir') then
                  call  arprmc(vvalue, TraceDir)
               elseif( vname .eq. 'Excom1' ) then
                  call arprmr(vvalue, Excom1) 
               elseif( vname .eq. 'Excom2' ) then
                  call arprmr(vvalue, Excom2) 
               elseif( vname .eq. 'LPMeffect' ) then
                  call arprml(vvalue, LPMeffect)                 
               elseif( vname == "StoppingPw" ) then
                  call arprmi(vvalue, StoppingPw)
               elseif( vname == "SrimEmax" ) then
                  call arprmr(vvalue, SrimEmax)
               elseif( vname == "PhitsXs" ) then
                  write(0,*) 'PhitsXs is not used now'
                  write(0,*) 'Please drop it from epicsfile'
                  stop
               elseif( vname == "JamXs" ) then
                  write(0,*) 'Sorry: JamXs is managed by param file'
                  write(0,*) 
     *       'If needed, put it(D=0) in &HPARAM part of param'
                  stop
               else
                  write(0,*) ' epicsfile parameter: ', vname,
     *                 ' is undefined '
                  if( vname .eq. "epHooks" ) then
                     write(0,*)
     *                "Now, this must be placed in sepicsfile"
                  endif
                  stop 0000
               endif
            enddo
            close(iowk)
            write(msg,*) ' epics parameters have been read from ',
     *      dsn(1:klena(dsn))
            call cerrorMsg(msg, 1)
            return
!      *************
       entry epprmw(io)
            write(io,*)'----------------------'
            call awprml(io,'AngleB', AngleB)
            call awprmr(io,'DtMax', DtMax)
            call awprmi(io,'Eabsorb', Eabsorb)
            call awprmr(io,'Eanihi', Eanihi)
            call awprmr(io,'Ecut ', Ecut)
            call awprml(io,'EdepDedx', EdepdEdx)
            call awprmr(io,'ElowerBndPair ', ElowerBndPair)
            call awprmr(io,'EminElec', EminElec)
            call awprmr(io,'EminGamma', EminGamma)
            call awprmi(io,'AutoEmin', AutoEmin)
            call awprmc(io, 'ModifyFile', ModifyFile)
            call awprmr(io,'EminH', EminH)
            call awprmr(io,'Es   ', Escat)
            call awprml(io,'Zcorrec', Zcorrec)
            call awprmr(io,'EupperBndCS ', EupperBndCS)
            call awprmr(io,'EpartialSC ', EpartialSC)
            call awprmi(io,'TargetElecBrems', TargetElecBrems)
            call awprmr(io,'Flpm ', Flpm)
            call awprmi(io,'HowNormBrems', HowNormBrems)
            call awprmi(io,'HowQuench', HowQuench)
!            call awprmi(io,'IncGp', IncGp)
            call awprmi(io,'IoTrace', IoTrace)
            call awprml(io,'Knckon', Knckon)
            call awprmi(io,'Molier', Moliere)
            call awprml(io,'Photo', Photo)
            call awprmr(io,'RecoilKeMin', RecoilKEmin)
            call awprmr(io,'Tcoef' , Tcoef)
            call awprmr(io,'Tmin' , Tmin)
            call awprml(io,'Trace' , Trace)
            call awprml(io,'TimeStruc' , TimeStruc)
            call awprml(io,'FreeC' , FreeC)
            call awprmr(io, 'EpsLeng', EpsLeng)
            call awprmi(io, 'ALateCor', ALateCor)             
            call awprmi(io, 'Sync',     Sync)
            call awprmr(io, 'SyncLoop',     SyncLoop)
            call awprmi(io, 'MagPair',     MagPair)
            call awprmi(io, 'MuNI', MuNI)
            call awprmi(io, 'MuBr', MuBr)
            call awprmi(io, 'MuPr', MuPr)
            call awprmr(io, 'KEmin', KEmin)
            call awprmi(io, 'MsgLevel', MsgLevel)
            call awprmr(io, 'Excom1', Excom1)
            call awprmr(io, 'Excom2', Excom2)
            call awprml(io, 'LPMeffect', LPMeffect)
            call awprmi(io, 'StoppingPw', StoppingPw)
!            call awprmi(io, 'PhitsXs', PhitsXs)
!            call awprmi(io, 'JamXx', JamXs)
            do i = 1, MaxMediaDir
               call awprmc(io, 'MediaDir', MediaDir(i))
            enddo
            call awprmc(io, 'TraceDir', TraceDir)

      end
!      ************
       subroutine  epcmp1
!      ************
       use epModify
       use SampMolReducedA
       implicit none
#include "ZepTrackp.h"
#include "ZepTrackv.h"
#include "ZepManager.h"
#include "Zmass.h"
!             compute something using basic parameters
!              v 9.201
       if( sbMin == 0. .and. Moliere == 1 ) then
          sbMin=3.395d0
       elseif(sbMin == 0. .and. Moliere == 2 ) then
          sbMin=1.5d0
       endif

       if(ALateCor .eq. 2 .or.
     *        ( ALateCor .eq. 1 .and.  Moliere == 0)) then
!          angular latelaral correlation will be taken at
!          multiple scattering. use standard trucation
          Tcoefx = Tcoef
          Tminx = Tmin
       else
!            no correlation will be considered so make
!          truncation length shorter
          Tcoefx = Tcoef/2
          Tminx = Tmin/2
       endif
       if(KEmin .eq. 0.) then
          KEmin = max(EminElec - masele, 100.d-6)  ! 100 keV
       endif
       if( EminH == 0. ) then  ! for neutron
          EminH = 20.d-3
       else
!             use EminH as it is.  Mar. 18, 2014 
!          EminH = min( EminH, 20.d-3)  
       endif
       if( RecoilKEmin .eq. 0.) then
          RecoilKEmin =EminElec - masele
       endif
!           save;  to be used in epModify;  in GeV
       EminGsave = EminGamma
       EminEsave = EminElec
       RecoilEsave = RecoilKEmin
       KEminsave = KEmin 
       Enminsave = EminH 
       ElecMass = masele  


       if( ModifyFile /= " " ) then  ! special modification requested
          call epfixModifier(ModifyFile)
       endif
!        get info from  Cosmos
       call epInfoPhotoP( IncGp )
! 
!       call epInfoPhitsJamXs( PhitsXs, JamXs )
       end
       subroutine epInfoPhotoP(gpinfo)
   !    if this is placed in epcosIntF.f, epgen must be linked
   !    which in tern requires ephook. So cannot be used
   !   from other main.
   !          copy the info for how to treat photo-had.prog
   !          from cosmos: HowPhotoP
      implicit none
#include "Ztrackp.h"
      integer,intent(out):: gpinfo
      gpinfo = HowPhotoP
      end
