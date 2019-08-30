#if defined (KEKA) || defined (KEKB)
#define DOMPI
#endif
!     ****************************************************************
!     * cmanager: Managing Cosmos Simulation
!     ****************************************************************
      subroutine cmanager
!      use modAtmosDef
      implicit none
#include "Zevhnp.h"
#include "Zmanagerp.h"      
      integer:: icon
#if defined (DOMPI)
#include "mpif.h"
#include "Ztrack.h"
#include "Zmpi.h"


      integer err, intdata
      character*120 parampath  ! absolute param path given by job script 
                      ! via envrionment


      
      call mpi_init(err)
      call mpi_comm_size(MPI_COMM_WORLD, mpisize, err)
      call mpi_comm_rank(MPI_COMM_WORLD, mpirank, err)
!        this is to treat exising name list as it is.
!         no delim="apstrophe" is needed
#if defined (KEKB) || defined (KEKA)      
      call setrteopts("namelist=old")
#endif
!!!      call kgetenv2("PARAMPATH", parampath) ! this is OK ??
      err = kgetenv2("PARAMPATH", parampath)
!           to avoid simultaneous access to HD, we
!           read namelist from rank0, rank1,... 
      if(mpirank .eq. 0) then
         write(0,*) ' mpisize=',mpisize
         open(TempDev, file=parampath, action="read")
         call creadParam(TempDev)
         close(TempDev)
         if( mpisize .gt. 1) then
            call MPI_SEND(mpirank, 1, MPI_INTEGER, 1, 1,
     *       MPI_COMM_WORLD, err)
         endif
      else
         call MPI_RECV(intdata, 1, MPI_INTEGER, mpirank-1, 1,
     *      MPI_COMM_WORLD, mpistat, err)
         open(TempDev, file=parampath, action="read")
         call creadParam(TempDev)
         close(TempDev)
!              mpisize=5 and mpirank=4  ==> no more
!                            mpirank=3  ==> send
         if( mpirank .lt. mpisize-1 ) then
            call MPI_SEND(mpirank, 1, MPI_INTEGER, mpirank+1, 1,
     *       MPI_COMM_WORLD, err)
         endif
      endif
!      write(0,*)' rank=', mpirank, ' Job =', Job
!//////////////
#else
#if defined (MacIFC)
      write(0,*) ' '
#endif
      call creadParam(5)        ! read execution conditions from stdin
#endif
#ifdef NEXT486
      if( index(IntModel, 'fritiof1.6') .gt. 0 ) then
         call cerrorMsg('fritiof1.6 cannot be used'//
     *  ' for NEXT-Absoft Fortran', 0)
      endif
#endif
      
      if( ObjFile /= " " ) then
!      non Earth; e.g, sun 
         call copenNLf(TempDev, ObjFile, icon)
         if(icon == 0) then
            call creadObjParam(TempDev)
            close(TempDev)
         else
            write(0,*) ' "ObjFile=', trim(ObjFile), ' cannot be opened'
            stop
         endif
      endif
      
      call cbeginRun            ! initialize the simulation
      call ceventLoop           ! begin simulation and enter  event Loop
      call cendRun              ! close the simulation
      end
