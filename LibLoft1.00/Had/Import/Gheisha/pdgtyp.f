
       function pdgtyp(ight)
c function to go from gheish/geant types to pdg number
       integer pdgtyp 
      
      integer types(49)
      data types/
c  gamma, e+, e-, nu(mu), mu+, mu-
     *      22, -11, 11, 14, -13, 13, 

c   pi-0, pi+, pi-, kl, k+, k-
     *      111, 211, -211, 130, 321, -321,

c   n, p, p(bar), ks, eta, lamda
     *      2112, 2212, -2212, 310, 0, 0, 
c   sigma+, sigma0, sigma-, xi0, xi-, 
     *       0, 0, 0, 0, 0,

c   omega-, n(bar), lambda(bar), 
     *       0, -2112, 0,

c   sigma(bar)-, sigma(bar)0, sigma(bar)+, xi(bar)0, xi(bar)+,
     *       0, 0, 0, 0, 0,

c   omega(bar)+, tau+, tau-, 
     *       0,  0, 0, 

c   d+, d-, d0, d0(bar), ds+, ds-, lambdac
     *       421, -421, 411, -411, 0, 0, 0,

c   w+, w-, z0, 
     *         0, 0, 0,
 
c   deutron, triton, alpha, 
     *         0, 0, 0,

c   geantino, psuedo
     *         0, 0/
       pdgtyp = 0
       if(ight.gt.0 .and. ight.le.48) pdgtyp = types(ight)
       return
       end
