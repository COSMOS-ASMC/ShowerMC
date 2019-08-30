         integer maxErgNodeForMuPair
         parameter( maxErgNodeForMuPair= MAXERGNODE_MUPAIR)
       type mubpn 
       sequence
!                  Z dependent  const for muon pair creation
            real*8 pa  
            real*8 ra
            real*8 pb
            real*8 qb
!                 used for brems sampling
            real*8 Ak
            real*8 Akm
            real*8 Akm2
            real*8 logf0
!                 usef for muon nuc. inte, if Kobayakawa
!                 model is used
            real*8 Shadow
            real*8 PointLike
            integer MuNI     ! see MuNI in ZepTrackp.h ==>modMuNucontrol
            integer MuBr     ! These are for individual media. 
            integer MuPr     ! First MuNI in ZepTrackp.h is copied
                             ! to mu.MuNI. If table is not found for
                             ! the media, mu.MuNI will be reset.
                             ! These will be set when media file is read

              !  sampling of emitted energy of  pair electrons  by muon 
              ! is not simple due to complex cross-section form
              ! at low energies (<1 TeV)  so  we give exact ds/dx
              !  (x=Epair/Emu) at several energy nodes. 
              ! For those energies, we use csampAF.  For energies
              ! in between, we use two eneriges at the nodes 
              ! with a weight of log-distance in energy.
              ! E1 <=  E  <= E2;  
              ! E1 is used with the weight of  log(E2/E)/log(E2/E1).
              ! if sample x is < possible xmin, retry is made
!         this is not used now since all direct pair spectrum has the
!         same form independent of Z. (cross-section itself is diff.)
!            so epmuPrsmp has another common PairId which is common
!            to all media.
!            integer muPairId(maxErgNodeForMuPair)  ! csampAF's id
!
            real*8  muPairErg(maxErgNodeForMuPair) ! energy node
            real*8  logmuPairErg(maxErgNodeForMuPair) ! log energy node
!                   this is also not used as explained a few lines above
!            integer NoOfErgNodeForMuPair         ! actual # of  id's 
!
                 ! cleared when muon Tab is read.  fixed when
                  ! muon pair creation is requested to be
                  ! maxErgNodeForMuPair. Erg and logErg are fixed
                  ! in epRdmuTab.
       end type mubpn

