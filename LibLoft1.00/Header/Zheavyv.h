!
        integer SkipPtclGen

         common /Zheavyvc2/ SkipPtclGen

!               in heavy particle collision, particle generation by 
!           interacting nucleons are skipped if this is 1.
!           if this is 0, generation is tried normallay .  This is set to 1
!           if Generate has 'qas': quick as generation. so that interacting
!           nucleons are replaced by component showers and AS generation
!           is managed by those  componen showers. No further tracking
!           of interacting nucleons are done. 
!
