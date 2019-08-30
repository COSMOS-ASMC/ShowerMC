#include "Zmaxdef.h"
!   Parameters   needed  for the Launcher.
!
!	(->	------------------------------------

         integer ErrorOut    !2 Error output logical  dev number.
         character*128  PrimaryFile  !1  Primary Spectrum data file (full or relative path)
         character*128  CutOffFile   !1  Geomagnetic cut-off file
         character*128  ContFile     !1  Job continuation information file  (full or relative path).
                                   !   default is "ContInfo".  This will be created when job
                                   !   is finished normally.
         character*128  GeomagFile   !2  IGRF or WMM file path which contains geomagnetic field expansion 
                                   !   coefficients.  Their format is the same one as given in their web 
                                   !   page.  If ' ' (default), Cosmos/Data/Geomag/igrf is used.
         character*128  SkeletonFile !1   Skeleton information file (full or relative path). created if Job =
                                   !    'skeleton'. Default is 'skeletonParam'.  This is the Namelist data
                                   !     referred by Cosmos automatically  if Job='flesh' is specified. For
                                   !     Job='flesh', you have to modify some part of  this file.
        character*128  DpmFile      !2  control card to specify the dpmjet execution conditions. If ' ',
                                    !   Cosmos/Data/DPM/atmos.inp is assumed.  
        character*10  Job          !1  What kind of job you are going to do.\newline
                                   !   =' ' (default).  nothing special.\newline
                                   !   ='skeleton'.  Makes skeleton. \newline
                                   !   ='flesh'. Flesh skeleton events.  See manual.\newline
                                   !   ='newskel'   \newline
                                   !   ='newflesh'  see manual. \newline
        character*128  SeedFile     !1   File to  contain the initial random numbers for those events to 
                                   !    which you want to flesh. You can create the file by calling
                                   !     cwriteSeed in a user hook routine (say, in chookEnEvent) at 
                                   !     skeleton making time. Default is 'Seed'.  For a normal run with
                                   !      Job=' ', if SeedFile is not ' ',  two integer initial random numbers
                                   !      and the event number are  automatically output on the speicfied disk file.
        integer       SeedFileDev  !2   logical device number of SeedFile.
        logical       Cont         !1  If T, continuation from a previous job is assumed. Contfile content is used.
        integer       InitRN       !1  Initial random number seed. 2 integers. If InitRN(1) $<$ 0, file dev  \# 14
                                   !    is  assumed to have  pairs of IR in each row, and they are read to
                                   !    initialize each event.  This feature is ignored when Job = 'flesh' or 
                                   !    'newflesh'. The \# 14 file should be opened by the user routine
                                   !    (chookBgRun). This is almost debug purpose.\newline
                                   !   If InitRn(2)$<$0, timer, hostname and process number are used for the 
                                   !    initialization.
         integer       EventNo      !2  cumulative event number counter.(excluding discarded ones due to cutoff).
         integer       EventsInTheRun !2  Counter for event number in the run. Internal use.
                                     !          (excluding discarded ones due to cutoff).
         integer       DestEventNo    !1 2 integers: Final event no. to be generated and events to be generated
                                     !  in the current run.  If negative, their absolute is used and counting 
                                     !  includes discarded ones due to rigidity cutoff.
                                     !  If DestEventNo(2)=0, DestEventNo(1) is used. If it is negative, only
                                     !  DestEventNo(2) is checked to see events in the current run. For the
                                     !  flux calculation, negative ones are better.
         logical       Hidden         !1  Make T, if hidden parameters are to be written.
         integer	      TempDev	   !2  Logical Dev. number for temporary disk use.
         integer       PrevEventNo  !2  The event number already finished.  System use for Cont job.
                                   !        (excluding discarded ones due to cutoff).
         character*8   DeadLine     !1  The dead line before which the job should terminate.
                                   !   Should be given like '10.11.15' which means the nearest 10th, 11 O'clock,
                                   !   15 min.  Not used if Within has non zero value.  
        integer       Within       !1  The job should end within this minutes from now.  Default is 99999.
                                   !   If 0 is given,  DeadLine is used.
        real*8        BaseTime     !1  Rough cpu time needed for completing one event (say, for protons, or
                                   !   gamma rays) with energy BaseErg.  The cpu time estimation is based on 
                                   !   A * ( E1ry par nucleon )**BasePower / BaseErg * BaseTime, where A is mass number
                                   !  (for nucleus; otherwise 1).
        real*8        BaseErg      !2  See BaseTime.  The default is  1000 (GeV).
        real*8        BasePower    !2  See BaseTime.   Default is 1.0
        character*100 UserHookc    !2  array size is MAX\_USERHOOKC(=5). Usage is left for the user. To get the i-th
                                   !   component, the use may 'call  cqUHookc(i, cv)' in the userHook routine, 
                                   !   where cv is a character variable to receive the data.
        real*8        UserHookr    !2  array size is MAX\_USERHOOKR(=10). Usage is left for the user. To get the i-th
                                   !   component, the use may 'call cqUHookr(i, rv)' in the userHook routine,
                                   !   where rv is a real*8 variable to receive the data.
        integer       UserHooki    !2  array size is MAX\_USERHOOKI(=10). Usage is left for the user.  To get the i-th
                                   !   component, the use may 'call ccqUHooki(i, iv)' in the userHook rouitne,
                                   !   where iv is an integer varialbe to receive the data.
       character*128 AtmosFile    !2  path to the atmospheric data as in 'Cosmos/Data/Atmos/stdatmos1.d' 
                                   ! Normally  this may be blank. Then, standard atmosphere
                                   ! is employed.  If another data with the same format
                                   ! atmosphere model is available, it can be specified.
                                   ! To change the atmosphere model, use shell command "atmosModel.sh"
                                   !   --> moved to moAtmosDef
        integer:: AtmosModel       !2  default is 1. which means "1: Standard"
                                   !  This one is not intended to change the atmospheric model, but
                                   !  to check the consistency with the current model in the library
                                   !  and the user's inent.
	                           !  
        integer:: NRL_period       !2 see Cosmos/Import/NRL/Util or manual for NRL atmsophere
        character*32  AtEnv        !2  If this is non blank, an environmental variable with that name is
                                   !   assumed to exist and Cosmos tries to get the value of that env variable.
                                   !   If the value is obtained, the \verb/@/ in \verb/@_/ or \verb/@./
                                   !   expressing a part of a file name is replaced by that value. 
                                   !   (default is blank and in that case the \verb/@/ is replaced by
                                   !    the host name where  the job runs.)

        character*32 SharpEnv      !2  If this is non blank, an environmental variable with that name is
                                   !   assumed to exist and Cosmos tries to get the value of that env variable.
                                   !   If the value is obtained, the \verb/#/ in \verb/#_/ or \verb/#./ 
                                   !   expressing a  part of a file name is replaced by that value. 
                                   !   (default is blank and in that case the \verb/#/ is replaced by
                                   !    the process number of the run).

        character*32 PercentEnv    !2  If this is non blank, an environmental variable with that name is
                                   !   assumed to exist and Cosmos tries to get the value of that env variable.
                                   !   If the value is obtained, the \verb/%/ in \verb/%_/ or \verb/%./ 
                                   !   expressing a  part of a file name is replaced by that value. 
                                   !   (default is blank and in that case the \verb/%/ is replaced by
                                   !    the USER name).
        character(128):: ObjFile ! default is " ". If want to define non
!     earth (say, Sun) , give a  file name in which namelist, ObjParam, is
!    defined to give  the radius, etc  of the sun.

!	<-)	-------------------------------------
        common /Zmanagerpc/ &
       BaseTime,  BaseErg, BasePower, Within, UserHookr(MAX_USERHOOKR), &
       ErrorOut, Cont, InitRN(2), UserHooki(MAX_USERHOOKI),  &
       EventsInTheRun, DestEventNo(2), NRL_period(4),  &
       Hidden, TempDev, AtmosModel, &
       PrevEventNo, SeedFileDev, EventNo


       common /Zmanagerpc2/  &
      UserHookc(MAX_USERHOOKC), PrimaryFile, ObjFile, &
      CutOffFile,  Job, ContFile, AtmosFile, GeomagFile, &
      SkeletonFile, SeedFile, DpmFile, DeadLine, SharpEnv, &
      PercentEnv, AtEnv
 
