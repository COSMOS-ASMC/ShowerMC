         subroutine ccos2gheCode(aPtcl, ghecode)
!        convert cosmos particle code  to Gheisha code
!        neutrinos are assumed to be not comming here, since
!        no interactions are assumed.  (Gheisha dose not
!        recognize the flavor).
!
         implicit none
#include "Zcode.h"
#include "Zptcl.h"
         type (ptcl):: aPtcl  ! input. a cosmos particle 
        integer ghecode  ! output.  Gheisha particle code.
!	
         integer code, subcode, charge, kdeut

         integer cos2ghe(klast, -1:1), anticos2ghe(klast, -1:1)

         character*100 msg
         parameter (kdeut =9)  ! legacy definition.

!       data cos2ghe( cosmos ,charge: )/ # /  ! subcode Gheisha
       data cos2ghe( kphoton ,0 )/ 1 /  ! 0 GAMMA
       data cos2ghe( kelec ,1 )/ 2 /  ! antip POSITRON
       data cos2ghe( kelec ,-1 )/ 3 /  ! regptcl ELECTRON
!       data cos2ghe( kx ,0 )/ 4 /  ! kx NEUTRINO
       data cos2ghe( kmuon ,1 )/ 5 /  ! 0 MUON+
       data cos2ghe( kmuon ,-1 )/ 6 /  ! 0 MUON-
       data cos2ghe( kpion ,0 )/ 7 /  ! 0 PION0
       data cos2ghe( kpion ,1 )/ 8 /  ! 0 PION+
       data cos2ghe( kpion ,-1 )/ 9 /  ! 0 PION-
!            need programing
       data cos2ghe( kkaon ,0 )/ 10 /  ! k0l KAON0LONG
       data cos2ghe( kkaon ,1 )/ 11 /  ! 0 KAON+
       data cos2ghe( kkaon ,-1 )/ 12 /  ! 0 KAON-
!            need programing
       data cos2ghe( knuc ,0 )/ 13 /  ! regptcl NEUTRON
       data cos2ghe( knuc ,1 )/ 14 /  ! regptcl PROTON
       data cos2ghe( knuc ,-1 )/ 15 /  ! antip ANTIPROTON

!            need programing
!       data cos2ghe( kkaon ,0 )/ 16 /  ! k0s KAON0SHORT
       data cos2ghe( keta ,0 )/ 17 /  ! 0 ETA
!            need programing

       data cos2ghe( klambda ,0 )/ 18 /  ! regptcl LAMBDA
!            need programing 
       data cos2ghe( ksigma ,1 )/ 19 /  ! regptcl SIGMA+
!            need programing 

       data cos2ghe( ksigma ,0 )/ 20 /  ! regptcl SIGMA0
!            need programing 
       data cos2ghe( ksigma ,-1 )/ 21 /  ! regptcl SIGMA-
!            need programing 
       data cos2ghe( kgzai ,0 )/ 22 /  ! regptcl XI0
!            need programing 

       data cos2ghe( kgzai ,-1 )/ 23 /  ! regptcl XI-
       data cos2ghe( kbomega, -1 )/ 24 /  ! kx OMEGA-
!            need programing 
       data anticos2ghe( knuc ,0 )/ 25 /  ! antip ANTINEUTRON
!            need programing 
       data anticos2ghe( klambda ,0 )/ 26 /  ! antip ANTILAMBDA
!            need programing 
       data anticos2ghe( ksigma ,-1 )/ 27 /  ! antip ANTISIGMA-
!            need programing 
       data anticos2ghe( ksigma ,0 )/ 28 /  ! antip ANTISIGMA0
!            need programing 
       data anticos2ghe( ksigma ,1 )/ 29 /  ! antip ANTISIGMA+
!            need programing 
       data anticos2ghe( kgzai ,0 )/ 30 /  ! antip ANTIXI0
!            need programing 
       data anticos2ghe( kgzai ,1 )/ 31 /  ! antip ANTIXI+
       data cos2ghe( kbomega, 1)/ 32 /  ! kx ANTIOMEGA+
!       data cos2ghe( kx ,1 )/ 33 /  ! kx TAU+
!       data cos2ghe( kx ,-1 )/ 34 /  ! kx TAU-
       data cos2ghe( kdmes ,1 )/ 35 /  ! 0 D+
       data cos2ghe( kdmes ,-1 )/ 36 /  ! 0 D-
!            need programing c        

       data cos2ghe( kdmes ,0 )/ 37 /  ! regptcl D0
!            need programing 
       data anticos2ghe( kdmes ,0 )/ 38 /  ! antip ANTID0
!       data cos2ghe( kx ,1 )/ 39 /  ! kx DS+
!       data cos2ghe( kx ,-1 )/ 40 /  ! kx DS-
!            need programing 
       data cos2ghe( klambdac ,1 )/ 41 /  ! regptcl LAMBDAC+
!       data cos2ghe( kx ,1 )/ 42 /  ! kx W+
!       data cos2ghe( kx ,-1 )/ 43 /  ! kx W-
!       data cos2ghe( kx ,0 )/ 44 /  ! kx Z0
        data cos2ghe( kdeut, 1 )/ 45 /  ! 0 DEUTERON
       data cos2ghe( ktriton ,1 )/ 46 /  ! 0 TRITON
!       data cos2ghe( kalfa , 1 )/ 47 /  ! 0 ALPHA
!       data cos2ghe( kx ,0 )/ 48 /  ! kx GEANTINO

       save cos2ghe, anticos2ghe

         code = aPtcl%code
         subcode = aPtcl%subcode
         charge = aPtcl%charge
!       	knuc, kkaon, klambda, ksigma, kgzai, kdmes
        if(code .eq. kpion) then
           ghecode = cos2ghe(kpion, charge)
        elseif(code .eq. kkaon) then
           if(charge .ne. 0) then
              ghecode = cos2ghe(kkaon, charge)
           elseif(subcode .eq. k0s ) then
              ghecode =16
           else
              ghecode = 10
           endif
         elseif(code .eq. knuc) then
           if(charge .ne. 0) then
              ghecode = cos2ghe(knuc, charge)
            elseif(subcode .eq. antip) then
               ghecode = anticos2ghe(knuc, 0)
           else	
               ghecode = cos2ghe(knuc, 0)
           endif	
        elseif(code .eq. ksigma ) then
           if(subcode .eq. antip) then
              ghecode = anticos2ghe(ksigma, charge)
           else
              ghecode = cos2ghe(ksigma, charge)
           endif
         elseif( code .eq. kgnuc .and. subcode .eq. 2 ) then
!	elseif(code .eq. kdeut) then
            ghecode = cos2ghe(kdeut, 1)
         elseif(code .eq. ktriton) then
            ghecode = cos2ghe(ktriton, 1)
         elseif(code .eq. kalfa) then
            ghecode = 47
        elseif(code .eq. klambda) then
           if(subcode .eq. antip) then
              ghecode = anticos2ghe(klambda, 0)
           else
              ghecode = cos2ghe(klambda, 0)
           endif
        elseif(code .eq. kgzai) then
           if(subcode .eq. antip) then
              ghecode =anticos2ghe(kgzai, charge)
           else
              ghecode =cos2ghe(kgzai, charge)
           endif
        elseif(code .eq. kdmes) then
           if(charge .ne. 0) then
              ghecode =cos2ghe(kdmes, charge)
           elseif(subcode .eq. antip) then
              ghecode = anticos2ghe(kdmes, 0)
           else
              ghecode = cos2ghe(kdmes, 0)
           endif
        elseif(code .eq. keta) then
            ghecode = cos2ghe(keta, 0)
         elseif(code .eq. kbomega) then
            ghecode = cos2ghe(kbomega, charge)
         elseif(code .eq. klambdac) then
            ghecode = cos2ghe(klambdac, charge)
        else
           write(msg, *)
     *     ' cosmos code=', code, 'charge=',charge,
     *     'subcode=',subcode,' cannot be put to Gheisha code'
           call cerrorMsg(msg, 0)
        endif
        end
