!            klast; max ptcl codes in the system.  (except krare)
      integer  kphoton, kelec, kmuon, kpion,  kkaon, knuc,
     1   kneue, kneumu,  knnb, kddb, kdmes, krho,
     2   komega, kphi, keta, kgnuc, kalfa, klibe, kcno, khvy, kvhvy,
     3   kiron, khvymax, klast, klambda, ksigma, kgzai, kbomega,
     4   ktriton, klambdac, krare, klight, kEdepo, kchgPath,
     5   kdeuteron,  kds, kXic, komeC0,  ktau, kneutau, ketap,
     6   kDelta, kXic0,
     7   kseethru
!            subcode
        integer
     1  regptcl, antip, k0s,  k0l, kneutron,
     2  kneutronb, kd0, kd0b, kdirectg, kcasg,
     3  kscinti, kceren, ksync, kChaX
  !
        parameter(
     1  kphoton=1, kelec=2, kmuon=3, kpion=4,  kkaon=5,
     2  knuc=6,
     3  kneue=7, kneumu=8, kgnuc=9, kalfa=10, klibe=11, kcno=12, 
     4  khvy=13, kvhvy=14, kiron=15, kdmes=16, 
!          next line added Nov. 17,'95. '
     5  ktriton=17, klambda=18, ksigma=19, kgzai=20, klambdac=21,
     6  kbomega=22, knnb=23, kddb=24, krho=25,komega=26, kphi=27,
     7  keta=28,
     9  khvymax=kiron, kdeuteron=29,
     a  kds=30, kXic=31, komeC0=32,  ktau=33, kneutau=34,
     b  ketap=35, kDelta=36, kXic0=37,
     c  krare = 0, 
     d  klight=-1, kEdepo=-2, KchgPath=-3, kseethru=-4,
     e  klast=37+4 ) ! 

!       kindmx=kbomega  not used now

        parameter(
     1  regptcl=-1, antip=1,
     2  kdirectg=2, kcasg=3, 
     3  k0s = 4,  k0l= 5, kscinti=1, kceren=2, ksync=3, kChaX=4,
     4  kneutron=regptcl,
     5  kneutronb = antip, kd0 =-8, kd0b =-kd0)
!             for heavy next are not used. ( to give
!             isotope, iso Z,A	may be used 				     
         integer maxHeavyMassN, maxHeavyCharge, maxHeavyG
         parameter (maxHeavyMassN = 56, maxHeavyCharge = 26, 
     1             maxHeavyG = 7)
!       kphoton: gamma ray 
!        kelec: electron, positiron
!        kmuon: muon
!        kpion: pion
!        kkaon: kaon
!        knuc: neucleon
!        kneue: electron neutrino
!       kneumu: muon neutrino
!        kgnuc: general nucleus(A>=2.)
!        kalfa: alpha  (heliunm)
!        klibe: Li, Be, B
!         kcno: C, N, O 
!         khvy: heavy such as, Na/Mg/Si
!        kvhvy: very heavy such as S/Cl/Ar
!        kiron: iron group
!        regptcl: particle index
!        antip: anti-particle index
!        klight: light normally 100 nm~1000 nm
!             subcode: kscinit scintillation light
!                      kceren  Cerekov light
!                      ksycn   synchrotron light
!        kEdepo: energy deposit in a small cell from whcih
!                scintillation lightis produced.
!        kchgPath: charged particle path form which Cerenkov
!               light is generated.
!        krare:  used to set very rare particle code
!                which might come from imported soft.
!                They are neglected in Cosmos. 
