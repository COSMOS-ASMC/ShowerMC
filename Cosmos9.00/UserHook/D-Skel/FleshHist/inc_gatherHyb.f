
      do i = 1, NoOfAsSites
         call mpi_gather(ASObsSites(i).age, 1, MPI_REAL8,
     *       ager, 1, MPI_REAL8, 0, MPI_COMM_WORLD, icon)
         call mpi_gather(ASObsSites(i).esize, 1, MPI_REAL8,
     *       esizer, 1, MPI_REAL8, 0, MPI_COMM_WORLD, icon)
         if(mpirank .eq. 0) then
            wage = 0.
            wsize = 0.
            do j = 1, mpisize
               wage = wage + ager(j)*esizer(j)
               wsize = wsize + esizer(j)
            enddo
            if(wsize .gt. 0.) then
               ASObsSites(i).age =wage/wsize
               ASObsSites(i).esize = wsize
            else
               ASObsSites(i).age =0.
               ASObsSites(i).esize =0.
            endif
         endif
      enddo
!       MPI_IN_PLACE dose not work at KEKB
      if( mpirank .eq. 0 ) then
         do i = 1, NoOfSites
            NgA(i) =0.
            NeA(i) =0.
            NmuA(i) = 0.
            NhadA(i) = 0.
            SumElossA(i) = 0.
         enddo
      endif

!      if(mpirank .eq. 0) then
        call mpi_reduce(Ng, NgA, NoOfSites,
     *   MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, icon)
!         call mpi_reduce(MPI_IN_PLACE, Ng, NoOfSites,
!     *   MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, icon)

        call mpi_reduce(Ne, NeA, NoOfSites,
     *   MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, icon)
!         call mpi_reduce(MPI_IN_PLACE, Ne, NoOfSites,
!     *   MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, icon)

       call mpi_reduce(Nmu, NmuA, NoOfSites,
     *   MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, icon)
!         call mpi_reduce(MPI_IN_PLACE, Nmu, NoOfSites,
!     *   MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, icon)

       call mpi_reduce(Nhad, NhadA, NoOfSites,
     *   MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, icon)
!         call mpi_reduce(MPI_IN_PLACE, Nhad, NoOfSites,
!     *   MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, icon)

       call mpi_reduce(SumEloss, SumElossA, NoOfSites,
     *   MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, icon)
!         call mpi_reduce(MPI_IN_PLACE, SumEloss, NoOfSites,
!     *   MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, icon)

!     else

!         call mpi_reduce(Ng, Ng, NoOfSites,
!     *   MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, icon)

!         call mpi_reduce(Ne, Ne, NoOfSites,
!     *   MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, icon)

!         call mpi_reduce(Nmu, Nmu, NoOfSites,
!     *   MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, icon)

!         call mpi_reduce(Nhad, Nhad, NoOfSites,
!     *   MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, icon)

!         call mpi_reduce(SumEloss, SumEloss, NoOfSites,
!     *   MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, icon)

!      endif
