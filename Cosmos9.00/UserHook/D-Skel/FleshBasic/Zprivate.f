       integer nsites ! # of sites to be histogramed.
       real bin, rmin  !  rbin in log10.  rmin. for lateral in Moliere u.
       parameter (bin=0.1, rmin=0.01)
       parameter (nsites=35)
       integer nrbin
       parameter ( nrbin = 42 )

       integer w2hl(nsites), w2il(nsites)
       integer histdep(nsites)  !  histdep(j) = k >0 ==> at depth k, histograminng 
                                !  is tried.
       integer indivdep(nsites)  ! indivdep(j) = k>0 ==> at depth k, 
                                       !  individual particle info. is output. 
       integer ansites, hnsites
       save ansites
       integer fno, fnoB
       parameter ( fno=2, fnoB=3)

       logical  recxy  !   T --> record (x,y) or r
!            These must be double.  because addition is done 1 by 1
!            and reach 10^7 or more.(at least  Ng, Ne, at E0
!            >10^16eV.)
       real*8  SumEloss(maxNoOfASSites),
     *   Ng(maxNoOfASSites), Ne(maxNoOfASSites),  Nmu(maxNoOfASSites),
     *   Nhad(maxNoOfASSites)

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
     *    rbin, 
     *    histdep, indivdep,  ansites, hnsites,
     *    recxy, tkelosslat,  tklat, tkarspec,
     *    tkrespec, tkrzspec, tkzfspec, tkrfspec, tkefspec,
     *    tkrtspec,
     *    tkretspec, tkrezspec, tkrzfspec, tkrefspec


  
       
