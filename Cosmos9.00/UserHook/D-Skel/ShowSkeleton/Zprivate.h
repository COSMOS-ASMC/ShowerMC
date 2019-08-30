#include "Zprivate0.h"
       integer w2hl(nsites), w2il(nsites)
       integer histdep(nsites)  !  histdep(j) = k >0 ==> at depth k, histograminng 
                                !  is tried.
       integer indivdep(nsites)  ! indivdep(j) = k>0 ==> at depth k, 
                                       !  individual particle info. is output. 
       integer ansites, hnsites
       integer fno, fnoB, fnonrfai
       parameter ( fno=2, fnoB=3, fnonrfai=4 ) 

       logical  recxy  !   T --> record (x,y) or r
!            These must be double.  because addition is done 1 by 1
!            and reach 10^7 or more.(at least  Ng, Ne, at E0
!            >10^16eV.)
       real*8  SumEloss(maxNoOfASSites),
     *   Ng(maxNoOfASSites), Ne(maxNoOfASSites),  Nmu(maxNoOfASSites),
     *   Nhad(maxNoOfASSites)
       real*8 dfai
       parameter ( dfai = 360.d0/nfai )
       real recprob(nrbin,  4, MaxNoOfASsites),
     *      nptcls(nrbin,  4,  MaxNoOfASsites),
     *      nrfaiRec(nrbin, nfai, 4, MaxNoOfASsites),
     *      nrfaiAll(nrbin, nfai, 4, MaxNoOfASsites),
     *      dErfai(nrbin, nfai, MaxNoOfASsites)
       logical tklat      ! take lateral or not
       logical tkelosslat ! take energy loss lateral or not

       logical tkarspec  ! azimuthal angle dependent lateral
       logical tkrespec

       logical  tkrzspec

       logical tkzfspec
  
       logical tkrfspec
       
       logical tkefspec

       logical tkrtspec

       logical tkretspec

       logical tkrezspec

       logical tkrzfspec

       logical tkrefspec


       real*8  rbin(nrbin+1)  
       common /Zprivatec/   Ng, Ne, Nmu, Nhad, SumEloss,
     *    rbin, recprob, nptcls, nrfaiRec, nrfaiAll,
     *    dErfai,
     *    histdep, indivdep,  ansites, hnsites,
     *    recxy, tkelosslat,  tklat, tkarspec,
     *    tkrespec, tkrzspec, tkzfspec, tkrfspec, tkefspec,
     *    tkrtspec,
     *    tkretspec, tkrezspec, tkrzfspec, tkrefspec


  
       
