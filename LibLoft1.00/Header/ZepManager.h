#include "ZepMaxdef.h"
       integer Nevrun   ! event counter generated in this run
       logical DoUpdateNo ! if t, update Nevrun in epe1ev
       integer iowk     ! logical  device number used for input disk file
       integer IoTrace  ! logical device number of trace info. 
                        ! If >0, only charged ptcl tracks are recorded
           	        ! If <0, abs is used. neutral ptcl tracks also recorded
			      
       character*128 TraceDir   ! The director where you take trace info
       character*128 ModifyFile ! some component must have non default
	 ! quenching, Emin  etc.  This file contains such info.
!->ZmediaLoft.h  integer MaxMediaDir  ! max number of directories where
                            ! media file is searched for.
       integer MsgLevel     ! 0,1,or 2.  0-->only fatal error is shown
                            ! 1--> some usefule/warning level info.
                            ! 2--> more details.
       integer epHooks(3)   ! number of c, i, r variables user wants to 
                            ! define and input from epicsfile(epHookc/i/r)

!->ZmediaLoft.h    parameter(MaxMediaDir = MAX_MEDIA_DIR)
       parameter( iowk = 11 )
!->ZmediaLoft.h     character*128 MediaDir(MaxMediaDir)
       character*100 epHookc(EPMAX_UHOOKC)
       integer epHooki(EPMAX_UHOOKI)
       real*8 epHookr(EPMAX_UHOOKR)
       common /epManag/ epHookr, epHooki, epHooks,
     *  Nevrun, MsgLevel, IoTrace, DoUpdateNo

!->ZmediaLoft.h     common /epManagC/ MediaDir
        common /epManagC/  TraceDir, ModifyFile, epHookc
