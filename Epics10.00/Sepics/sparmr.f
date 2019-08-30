!     *****************************************************************
!     *                                                               *
!     *  sparmr: read       sepics parameter
!     *  sparmw  write                                                *
!     *****************************************************************
!
!
       subroutine sparmr(dsn)
       use modShiftInciPos, only: ShiftInciPos
       implicit none
#include  "ZsepManager.h"
#include  "ZepPos.h"
#include  "ZepDirec.h"
#include  "Zsparm.h"
#include  "ZepTrackp.h"
#include  "ZepManager.h"

           character*(*) dsn
           integer io, icon, i
           character*100 msg
           character*24 vname
           character*100 vvalue
           logical epgetParmN
           integer hookcc, hookic, hookrc
           hookcc = 0
           hookic = 0
           hookrc = 0
                                                      
!

           call copenf(iowk, dsn, icon)
           if(icon .ne. 0) then
              call cerrorMsg(' file not exist', 1)
              call cerrorMsg(dsn, 0)
           endif
!               skip separator
           call afsep(iowk)
!               read each item
           do while ( epgetParmN(iowk, vname, vvalue) )
              if( vname .eq. 'BaseTime' ) then
                 call arprmr(vvalue, BaseTime)
              elseif( vname .eq. 'CosNormal') then
                 call arprmm(vvalue, CosNormal)
              elseif( vname .eq. 'DCInpX') then
                 call arprmr(vvalue, DCInp%x)
              elseif( vname .eq. 'DCInpY') then
                 call arprmr(vvalue, DCInp%y)
              elseif( vname .eq. 'DCInpZ') then
                 call arprmr(vvalue, DCInp%z)
              elseif( vname .eq. 'DeadLine') then
                 call arprmc(vvalue, DeadLine)
              elseif( vname .eq. 'E0Max') then
                 call arprmr(vvalue, E0Max)
              elseif( vname .eq. 'Hwhm' ) then
                 call arprmr(vvalue, Hwhm)
              elseif( vname .eq. 'InputA') then
                 call arprmc(vvalue, InputA)
              elseif( vname .eq. 'InputP') then
                 call arprmc(vvalue, InputP)
              elseif( vname .eq. 'Ir1(1)') then
                 call arprmi(vvalue, Ir1(1))
              elseif( vname .eq. 'Ir1(2)') then
                 call arprmi(vvalue, Ir1(2))
              elseif( vname .eq. 'JobTime') then
                 call arprmi(vvalue, JobTime  )
              elseif( vname .eq. 'Nevent') then
                 call arprmi(vvalue, Nevent)
              elseif( vname .eq. 'PrimaryFile') then
                 call arprmc(vvalue, PrimaryFile)
              elseif( vname .eq. 'ProfR') then
                 call arprmr(vvalue, ProfR)
              elseif( vname == 'ShiftInciPos') then
                 call arprmc(vvalue, ShiftInciPos)
              elseif( vname .eq. 'Xinp') then
                 call arprmr(vvalue, PosInp%x)
              elseif( vname .eq. 'Yinp') then
                 call arprmr(vvalue, PosInp%y)
              elseif( vname .eq. 'Zinp') then
                 call arprmr(vvalue, PosInp%z)
              elseif( vname .eq. 'LogIr' ) then
                 call arprml(vvalue,LogIr) 
              elseif( vname .eq. 'Xrange' ) then
                 call arprmm(vvalue, Xrange)
              elseif( vname .eq. 'Yrange' ) then
                 call arprmm(vvalue, Yrange)
              elseif( vname .eq. 'Zrange' ) then
                 call arprmm(vvalue, Zrange)
              elseif( vname .eq. 'epHooks') then
                 call arprmc(vvalue,  msg)
                 read(msg, *)  epHooks
              elseif( vname .eq. 'epHookc') then
                 hookcc = hookcc + 1
                 if(hookcc .gt. epHooks(1)) then
                    write(0,*)
     *             ' warning: number of userdefined char const>',
     *             ' epHooks(1)=',epHooks(1)
                 endif
                 call arprmc(vvalue,  epHookc(hookcc))

              elseif( vname .eq. 'epHooki') then
                 hookic = hookic + 1
                 if(hookic .gt. epHooks(2)) then
                    write(0,*)
     *               ' warning: number of userdefined int  const>',
     *               ' epHooks(2)=',epHooks(2)
                 endif
                 call arprmi(vvalue,  epHooki(hookic))

              elseif( vname .eq. 'epHookr') then
                 hookrc = hookrc + 1
                 if(hookrc .gt. epHooks(3)) then
                    write(0,*)
     *              ' warning: number of userdefined double  const>',
     *              ' epHooks(3)=',epHooks(3)
                 endif
                 call arprmr(vvalue,  epHookr(hookrc))
              elseif( vname .eq. 'Trace' ) then
                 call arprml(vvalue,Trace) 
              elseif( vname .eq. 'TraceDir' ) then
                 call arprmc(vvalue, TraceDir)
              elseif(vname .eq. 'TraceErg') then
                 call arprmra(vvalue, TraceErg, 6)
              elseif(vname .eq. 'Light' ) then
                 call arprmi(vvalue, Light)
              elseif(vname .eq. 'LightDir' ) then
                 call arprmc(vvalue, LightDir)
              elseif(vname .eq. 'StackDiskFile') then
                 call arprmc(vvalue, StackDiskFile)
              elseif(vname .eq. 'StackDiskFileNo' ) then
                 call arprmi(vvalue, StackDiskFileNo)
              elseif(vname .eq. 'OutPrimaryFile') then
                 call arprmc(vvalue, OutPrimaryFile)
              elseif(vname .eq. 'OutPrimaryFileNo' ) then
                 call arprmi(vvalue, OutPrimaryFileNo)
              else
                 write(0,*) ' sepicsfile parameter: ', vname,
     *           ' is undefined '
                 stop 0000
              endif
           enddo
           close(iowk)
           if(hookcc .le. EPMAX_UHOOKC ) then
              epHooks(1) = hookcc
           else
              write(0,*) 'too many user defined char const',
     *         ' in sepicsfile'
              write(0,*) 'Increase EPMAX_UHOOKC'
              stop  1111
           endif
           if(hookic .le. EPMAX_UHOOKI ) then
              epHooks(2) = hookic
           else
              write(0,*) 'too many user defined int const',
     *         ' in sepicsfile'
              write(0,*) 'Increase EPMAX_UHOOKI'
              stop  2222
           endif
           if(hookrc .le. EPMAX_UHOOKR ) then
              epHooks(3) = hookrc
           else
              write(0,*) 'too many user defined double const',
     *         ' in sepicsfile'
              write(0,*) 'Increase EPMAX_UHOOKR'
              stop  3333
           endif
           
           write(msg,*) ' sepics parameters have been read from '
           call cerrorMsg(msg, 1)
           call cerrorMsg(dsn, 1)
           return
       entry sparmw(io)
           write(io,*)'----------------------'
           call awprmr(io,'BaseTime', BaseTime)
           call awprmm(io,'CosNormal', CosNormal)
           call awprmr(io,'DCInpX', DCInp%x)
           call awprmr(io,'DCInpY', DCInp%y)
           call awprmr(io,'DCInpZ', DCInp%z)
           call awprmc(io,'DeadLine', DeadLine)
           call awprmr(io,'E0Max', E0Max)
           call awprmr(io,'Hwhm', Hwhm)
           call awprmc(io,'InputA', InputA)
           call awprmc(io,'InputP', InputP)
           call awprmi(io,'Ir1(1)', Ir1(1))
           call awprmi(io,'Ir1(2)', Ir1(2))
           call awprmi(io,'JobTime  ', JobTime  )
           call awprmi(io,'Nevent',Nevent)
           call awprmc(io, 'PrimaryFile', PrimaryFile)
           call awprmr(io,'ProfR', ProfR)
           call awprmr(io,'ShiftInciPos', ShiftInciPos)
           call awprmr(io,'Xinp', PosInp%x)
           call awprmr(io,'Yinp', PosInp%y)
           call awprmr(io,'Zinp', PosInp%z)
           call awprml(io,'LogIr', LogIr)
           call awprmm(io,'Xrange', Xrange)
           call awprmm(io,'Yrange', Yrange)
           call awprmm(io,'Zrange', Zrange)
           call awprml(io,'Trace', Trace)
           call awprmc(io,'TraceDir', TraceDir)
           call awprmra(io, 'TraceErg', TraceErg, 6)
           call awprmi(io, 'epHooks', epHooks)
           do i =1, epHooks(1)
              call awprmc(io, 'epHookc', epHookc(i))
           enddo
           do i =1, epHooks(2)
              call awprmi(io, 'epHooki', epHooki(i))
           enddo
           do i =1, epHooks(3)
              call awprmr(io, 'epHookr', epHookr(i))
           enddo
           call awprmi(io, 'Light', Light)
           call awprmc(io, 'LightDir', LightDir)
           call awprmi(io, 'StackDiskFileNo', StackDiskFileNo)
           call awprmc(io, 'StackDiskFile', StackDiskFile)
           call awprmi(io, 'OutPrimaryFileNo', OutPrimaryFileNo)
           call awprmc(io, 'OutPrimaryFile', OutPrimaryFile)
       end
