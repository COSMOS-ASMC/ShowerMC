         integer mpisize, mpirank, mpierr
         integer mpistat(MPI_STATUS_SIZE) 
         integer maxmpisize
         parameter ( maxmpisize = 2048 )
         common /commmpi/  mpisize, mpirank, mpierr, mpistat
