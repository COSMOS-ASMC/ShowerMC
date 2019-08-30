/*
c             ptcl kind code; kindmx is the no. of observable ptcls
c             klast; max ptcl code in the system.
c
*/
const int kphoton = 1;
const int kelec = 2;
const int kmuon = 3;
const int kpion = 4;
const int kkaon = 5;
const int knuc = 6;
const int kneue = 7;
const int kneumu = 8;



const int kdmes = 16;


 



const int kgnuc = 9;
const int kalfa = 10;
const int klibe = 11;
const int kcno = 12;
const int khvy = 13;
const int kvhvy = 14;
const int kiron = 15;
const int khvymax = kiron;


const int klambda = 18;
const int ksigma = 19;
const int kgzai = 20;
const int kbomega =  22;
const int ktriton = 17;
const int klambdac =21;
const int krare = 0;
const int kindmx = kbomega;
const int knnb = kindmx+1;
const int kddb = knnb+1;
const int krho = kddb+1;
const int komega =krho+1;
const int kphi = komega+1;
const int keta = kphi+1;
const int klast = keta;
//            subcode
const int regptcl = -1; 
const int antip = 1;
const int k0s = 4;    
const int k0l = 5;
const int kneutron = regptcl;
const int kneutronb = antip;
const int kd0 = -8;
const int kd0b = -kd0;
const int kdirectg = 2;
const int kcasg = 3;
const int maxheavymassn = 56;
const int maxheavycharge = 26;
const int maxheavyg = 7;
/*
 c       kphoton: gamma ray 
 c        kelec: electron, positiron
 c        kmuon: muon
 c        kpion: pion
 c        kkaon: kaon
 c        knuc: neucleon
 c        kneue: electron neutrino
 c       kneumu: muon neutrino
 c        kgnuc: general nucleus(A>=2.)
 c        kalfa: alpha  (heliunm)
 c        klibe: Li, Be, B
 c         kcno: C, N, O 
 c         khvy: heavy such as, Na/Mg/Si
 c        kvhvy: very heavy such as S/Cl/Ar
 c        kiron: iron group
 c        regptcl: particle index
 c        antip: anti-particle index
 c        krare:  used to set very rare particle code
 c                which might come from imported soft.
 c                such as tau. They are neglected in
 c                Cosmos.
*/
