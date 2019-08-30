!        here we write limit; at reduce process, non enhanced value
!        is needed (at MPI env, limit is real limit now)
#if defined (DOMPI)
      msg =dir(1:lengdir)//"/"//execid(1:lengid)//".nrfai"               
#define  nrfaiAllXXX nrfaiAllA
#define  dErfaiXXX   dErfaiA
#else
      msg =dir(1:lengdir)//"/"//execid(1:lengid)//
     *     "-@."//numb(1:lengn)//".nrfai"

#define  nrfaiAllXXX nrfaiAll
#define  dErfaiXXX   dErfai
#endif

      call copenfw2(fnonrfai, msg, 1, icon)
      if(icon .gt. 1) then
         write(0,*) ' icon=', icon 
         call cerrorMsg(msg, 1)
         call cerrorMsg('could not be opened', 0)  
      endif            
      write(fnonrfai, 
     *  '(i2,1pE11.3, 0p,i3, f8.4, 1pE11.3, 0p, 4i4,4f10.0)')
     *   EventNo, E0, NN, cosz, KEminObs, nrbin, nfai, ansites,
     *   NoOfSites, limit   ! this is not exist in the older version 
!              KEminObs above was limit (non array). 
!
      do i = 1, ansites
         l = indivdep(i)
         do j = 1, 4
            do k = 1, nfai
               write(fnonrfai, '("rec",f7.1, 4i4)' )
     *              DepthList(l)*0.1, l, i, j, k
               write(fnonrfai, '(1p10E11.3)')
     *             ( nrfaiRec(m,k,j,i), m=1,nrbin )
            enddo
         enddo
      enddo

      do i = 1, NoOfSites
         do j = 1, 4
            do k = 1, nfai
               write(fnonrfai, '("all",f7.1, 4i4)' )
     *         DepthList(i)*0.1, i, i, j, k
               write(fnonrfai, '(1p10E11.3)')
     *             ( nrfaiAllXXX(m,k,j,i), m=1,nrbin )
            enddo
         enddo
      enddo
!           dE/dx lateral
      do i = 1, NoOfSites
         do k = 1, nfai
            write(fnonrfai, '("dE/dx",f7.1, 3i4)' )
     *         DepthList(i)*0.1, i, i, k
            write(fnonrfai, '(1p10E11.3)')
     *             ( dErfaiXXX(m,k,i), m=1,nrbin )
!                 next is not yet assumed in MPI
            if(SeeLowdE) write(fnonrfai, '(1p10E11.3)') 
     *             ( dErfai2(m,k,i), m=1,nrbin )
         enddo
      enddo

      write(fnonrfai, *)
