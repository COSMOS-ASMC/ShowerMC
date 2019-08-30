#include "Zprivate0.h"
       integer w2hl(nsites), w2il(nsites)
       integer histdep(nsites)  !  histdep(j) = k >0 ==> at depth k, histograminng 
                                !  is tried.
       integer indivdep(nsites)  ! indivdep(j) = k>0 ==> at depth k, 
                                       !  individual particle info. is output. 
       integer ansites, hnsites
       integer fno, fnoB, fnonrfai
       parameter ( fno=2, fnoB=3, fnonrfai=4 ) 



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
     *      dErfai(nrbin, nfai, MaxNoOfASsites),
     *      dErfai2(nrbin, nfai, MaxNoOfASsites)


       real*8  CosRot, SinRot 	 
       real*8  rbin(nrbin+1)  
         common /Zprivatec/  SumEloss, CosRot, SinRot,
     *    rbin, recprob, nptcls, nrfaiRec, nrfaiAll,
     *    dErfai, dErfai2,
     *    histdep, indivdep,  ansites, hnsites,
     *    Ng, Ne, Nmu, Nhad

