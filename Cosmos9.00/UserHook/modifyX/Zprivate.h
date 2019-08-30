       integer nsites ! # of sites to be histogramed.
       real bin, rmin  !  rbin in log10.  rmin. for lateral in Moliere u.
       integer binw
       parameter (bin=0.1, rmin=0.01)
       parameter (nsites=3)  ! max number of histogram layers
       
       integer histdep(nsites)  !  histdep(j) = k >0 ==> at depth k, histograminng 
                                !  is tried.

!            These must be double.  because addition is done 1 by 1
!            and reach 10^7 or more.(at least  Ng, Ne, at E0
!            >10^16eV.)
       real*8  SumEloss(maxNoOfSites),
     *  Ng(maxNoOfSites), Ne(maxNoOfSites),
     *  Nmu(maxNoOfSites), Nhad(maxNoOfSites)



       logical tkarspec
       logical tkrtspec
       logical tkweb
       character*192 basefilename, basefilename2, filename

       integer fnoT, fnoL,  fnoB, fnoN
       integer nrbin, nfai
       parameter ( nrbin = 42, nfai=12, 
     *   fnoT=42, fnoL=43, fnoN=44, fnoB=45) 
       real*8  rbin(nrbin), 
     * webmin(nrbin,  nfai, MaxNoOfSites)
       real dErfai(nrbin, nfai, MaxNoOfSites)
       real dECent( MaxNoOfSites )
       real*8 dfai
       parameter ( dfai = 360.d0/nfai )


       common /Zprivatec/   Ng, Ne, Nmu, Nhad, SumEloss,
     *   rbin,   dErfai,
     *   histdep,  tkrtspec, tkarspec, tkweb

       common/ Zprivatec2/ basefilename

