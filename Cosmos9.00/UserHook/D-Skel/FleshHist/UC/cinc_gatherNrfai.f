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
      call mpi_reduce(nrfaiRec, nrfairecA, ansites*4*nfai*nrbin, 
     *   MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, icon)

      call mpi_reduce(nrfaiAll, nrfaiAllA, MaxNoOfSites*4*nfai*nrbin, 
     *   MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, icon)

      call mpi_reduce(dErfai, dErfaiA, MaxNoOfSites*nfai*nrbin, 
     *   MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, icon)
!          gatherred rec data is broadcasted; not needed; we work at rank 0 
!      call mpi_bcast(nrfairecA, ansites*4*nfai*nrbin, MPI_REAL,
!     *      0, MPI_COMM_WORLD, icon)
