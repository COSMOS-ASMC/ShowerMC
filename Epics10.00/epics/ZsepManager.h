        integer Ioprim       ! primary container.  code subcode ...etc
        integer IoScratch    
        integer Inpxyz, Inpdir, Inperg, Inpsubcode, Inptime, 
     *  Inp1ry, Inpmul, Inpdisk, Inpuser, Inpcn, Inpwgt, Inppol,
     *  Inpwl, Inpmass, InpLight,  Out1ry
	integer IoCont, Ir1,  Nevent
	real*8  JobTime, BaseTime
        character*8 pdatE0, ptimE0, Ddate, DeadLine
        character*128 LightDir
        character*256 EpicsFile, SepicsFile, ConfigFile
        character*256  PrimaryFile, StackDiskFile
        character*256 OutPrimaryFile
        character*200  buf
        integer StackDiskFileNo, OutPrimaryFileNo, OutPrimEff
!          OutPrimEff assigned to OutPrimaryFileNo or IoScratch
 

       parameter(  Ioprim=15, IoScratch=16 )
!             sepics common
         common /ZsepManegc/ EpicsFile, SepicsFile, ConfigFile,
     *  PrimaryFile, LightDir,  StackDiskFile,
     *  OutPrimaryFile, buf,
     *  ptimE0,  pdatE0, Ddate, DeadLine

         common /ZsepManeg/ BaseTime,
     *   IoCont, JobTime,
     *   Ir1(2), Nevent, 
     *   Inpxyz, Inpdir, Inperg, Inpsubcode, Inptime, Inp1ry,
     *   Inpmul, Inpdisk, Inpuser, Out1ry,
     *   Inpcn, Inpwgt, Inppol, Inpwl, Inpmass, InpLight,
     *   StackDiskFileNo,  OutPrimaryFileNo, OutPrimEff

!
! IoCont: fortran dev # for cont dataset.
!  BaseTime: base time for completing 1 gev event.
!        to is used to judge time needed for an event
! pdatE0: to keep first exEcution date
! ptimE0: //                      time
!        integer StackDiskFileNo
              ! file # for direct access disk file for 
	      ! ptcl stacking.  If not given, 13 will be used
	      ! = DdiskNoDefault in epLightMaxDef.f90
              ! actual variable is DdiskNo in modepStack
!    character*128 StackDiskFile ! file name path of the disk  for above.
			      
