            real*8  dEbydEdx(0:MaxNoOfSites+1)   ! real energy loss by dE/dx
            real*8  dEbyDeath(0:MaxNoOfSites+1)  ! Energy contained in all dead ptcls
            real*8  dEbyDeathG(0:MaxNoOfSites+1) ! in gamma
            real*8  dEbyDeathE(0:MaxNoOfSites+1) ! in electron
            real*8  dEbyDeathMuPiK(0:MaxNoOfSites+1)  ! in mu, pi, K
            real*8  dEbyDeathNeu(0:MaxNoOfSites+1)    ! in neutrino
            real*8  dEbyDeathP(0:MaxNoOfSites+1)    ! in proton
            real*8  dEbyDeathNut(0:MaxNoOfSites+1)    ! in Nutron 
            real*8  dEbyDeathO(0:MaxNoOfSites+1)      ! in other ptcls
            real*8  Espace(7),  Ecrash(7)
            real*8  MaxEbreak(2), MaxRelEbreak(2)
            real*8  SumEdiff, SumAbsEdiff 
            type(track):: inci
            type(coord):: angle
            common /Zabsob/  inci, angle,
     *      dEbydEdx, dEbyDeath, dEbyDeathG, dEbyDeathE, dEbyDeathMuPiK,
     *      dEbyDeathNeu, dEbyDeathP, dEbyDeathO, deByDeathNut, Espace,
     *      Ecrash, MaxEbreak, MaxRelEbreak, SumEdiff, SumAbsEdiff 
