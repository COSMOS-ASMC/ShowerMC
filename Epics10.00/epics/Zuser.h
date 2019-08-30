!             for user interface
!     containes variables for a current particle
!
         common /Zuser/ xbndry, ybndry,  zbndry, tm,
     1   x, y, z, e, wx, wy, wz, de,  am, SumDe,
     2   x1st, y1st, z1st, pap, cn, ic, k
         common /Zuserc/ proc1
!
         real*8 tm, x, y, z, e, wx, wy, wz, de,
     1          x1st, y1st, z1st, am, SumDe, xbndry,
     2          ybndry, zbndry

	integer pap

         integer  cn, ic, k
         character*4 proc1

!
!   xbndry, ybndry:  bondary x,y at component boundary.
!        tm: sum of track length(cm)/beta
!   x, y, z: particle position
!         e: total energy in gev
!  wx,wy,wz: direction cosines
!        de: energy loss in gev
!        am: mass in gev
!        cn: component # where the particle resides
!        ic: charge
!         k: particle kind =  kgamma--> gamma
!                             kelec-->  electron
!                             kmuon --> muon
!                             kpion --> pion
!                             kkon -->  kon
!                             knuc -->  nucleon
!                             kneue ---> electron neutrino
!                             kneumu---> muon neutrino
!                              these are defined in Zepdef
!   proc1 :  if '  ' no collision nor decay
!            otherwise, 'comp' --> compton
!                       'col'  --> collision
!                       etc
!  x1st, y1st, z1st:  first interaction possition.
!               have meaning if proc1 !='   '
!
