!            we note the last index may be actual max used
!
      do k = 1, ansites
         do j = 1, 4
            do l = 1,nfai
               do i = 1, nrbin
                  nrfaiRecA(i, l, j, k)=0.
               enddo
            enddo
         enddo
      enddo
      do k = 1, MaxNoOfSites
         do j = 1, 4
            do l = 1,nfai
               do i = 1, nrbin
                  nrfaiAllA(i, l, j, k)=0.
               enddo
            enddo
         enddo
      enddo
      do k = 1,  MaxNoOfsites
         do l = 1,nfai
            do i = 1, nrbin
               dErfaiA(i, l, k) = 0.
            enddo
         enddo
      enddo
!      if(mpirank .eq. 0 ) then
      call mpi_reduce(nrfaiRec, nrfairecA, ansites*4*nfai*nrbin, 
     *   MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, icon)
!      call mpi_reduce(MPI_IN_PLACE, nrfairec, ansites*4*nfai*nrbin, 
!     *   MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, icon)

      call mpi_reduce(nrfaiAll, nrfaiAllA, MaxNoOfSites*4*nfai*nrbin, 
     *   MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, icon)

!         call mpi_reduce(MPI_IN_PLACE, nrfaiAll,
!     *         MaxNoOfSites*4*nfai*nrbin, 
!     *   MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, icon)

      call mpi_reduce(dErfai, dErfaiA, MaxNoOfSites*nfai*nrbin, 
     *   MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, icon)
!         call mpi_reduce(MPI_IN_PLACE, dErfai, MaxNoOfSites*nfai*nrbin, 
!     *   MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, icon)
!      else

!         call mpi_reduce(nrfaiAll,  nrfaiAll, MaxNoOfSites*4*nfai*nrbin, 
!     *   MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, icon)

!         call mpi_reduce(dErfai, dErfai, MaxNoOfSites*nfai*nrbin, 
!     *   MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, icon)

!      endif
