!        Zmanager.h
         
         integer SeedSave(2) ! to store initial seed of random number 
                            ! for each event.  got at cevenLoop.f
        logical RefreshIR   ! becomes t if InitRN(1) < 0 && Job != 'flesh'
                            ! each event IR is read from #14
                            ! the file must be opened by the user routine
        character*80 PrefixConf ! Epics config file directory
        character*80 TopDir ! $LIBLOFT
        character*6  CosOrEpi ! cosmos or epics; from of cintModels is set.
        integer  TopDirLeng ! length of TopDir 
        integer  PrefixLeng ! lenght of PrefixConf
         common /Zmanager/ SeedSave, RefreshIR, TopDirLeng, PrefixLeng
        common /Zmanagerc/ TopDir, PrefixConf, CosOrEpi



