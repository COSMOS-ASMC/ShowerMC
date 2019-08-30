c....Example for pi0 and eta decay.

      include 'jam1.inc'
      include 'jam2.inc'
c...PYTHIA common block
      common/pjjets/njet,npad,kjet(1000,5),pjet(1000,5),vjet(1000,5)
      character frame*8,proj*8,targ*8,cwin*15


c...Initialize jam.
      nevent=10
      frame='nn'         ! comp. frame
      proj='p'           ! projectile
      targ='12c'         ! target
      cwin='12gevc'      ! incident energy
      call jaminit(nevent,0d0,-1d0,100d0,1,frame,proj,targ,cwin)

c...Loop for Simulation
      do iev=1,nevent

c....Simulate one JAM event.
          call jamevt(iev)

c......save mdcy for JAM simulation.
          mdcsave1=mdcy(jamcomp(111),1)
          mdcsave2=mdcy(jamcomp(221),1)
          mdcy(jamcomp(111),1)=1          ! make pi0 decay
          mdcy(jamcomp(221),1)=1          ! make eta decay

c...Loop for particle from jam.
        do ip=1,nv
          if(k(1,ip).gt.10) goto 1000     ! skip dead particle.
          if(k(2,ip).eq.111.or.k(2,ip).eq.211) then
            njet=1
            kjet(1,1)=1
            kjet(1,2)=k(2,ip)
            kjet(1,3)=0
            kjet(1,4)=0
            kjet(1,5)=0
            pjet(1,1)=p(1,ip)
            pjet(1,2)=p(2,ip)
            pjet(1,3)=p(3,ip)
            pjet(1,4)=p(4,ip)
            pjet(1,5)=p(5,ip)

c..........pi0/eta decay.
            call pjdecy(1,icon)
            call pjlist(1)         ! list particles.

c...........do something for analysis


          endif

 1000    end do

         mdcy(jamcomp(111),1)=mdcsave1
         mdcy(jamcomp(221),1)=mdcsave2 


      end do         

      end
