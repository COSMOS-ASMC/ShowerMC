      implicit none
      real(8)::p1(4), p2(4)
      integer irej, i, j
      real(8):: eelec, prom, eprot, egam, sigmax, fac

      integer NMXHEP
      PARAMETER (NMXHEP=4000)

      INTEGER NEVHEP,NHEP,ISTHEP,IDHEP,JMOHEP,JDAHEP
      DOUBLE PRECISION PHEP,VHEP
      COMMON /POEVT1/ NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &                JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),
     &                VHEP(4,NMXHEP)
!  extension to standard particle data interface (PHOJET specific)              
      INTEGER IMPART,IPHIST,ICOLOR
      COMMON /POEVT2/ IMPART(NMXHEP),IPHIST(2,NMXHEP),ICOLOR(2,NMXHEP)



!      call pho_inp(-1, irej)
      call pho_init(-1, irej)

      prom =0.938
!      eprot=820.
      eprot=prom
      p1(1)=0
      p1(2) = 0
      p1(3) = sqrt(eprot**2 - prom**2)
      p1(4) = eprot


      eelec =100
      egam = 0.45*eelec
      p2(1) = 0
      p2(2) = 0
      p2(3)=-egam
      p2(4)=egam

      call pho_setpar(1, 2212, 0, 0.d0)
      call pho_setpar(2, 22, 0, 0.d0)      

      call pho_event(-1, p1, p2,sigmax, irej)

      write(0,*) ' sigmax=',sigmax, ' ireg=',irej


      do i = 1, 10

         call pho_event(3, p1, p2, fac, irej)
         write(0,*) ' fac=',fac, ' ireg=',irej

        write(0,*) ' nhep=',NHEP
        do  j = 1, NHEP
           if( ISTHEP(j) == 1 ) then
              write(0,'(2i5, 1p, 5g13.3)')
     *             ISTHEP(j), IDHEP(j), PHEP(:, j)
           endif
        enddo
      enddo
      
      call pho_event(-2, p1, p2, sigmax,irej)
      end
