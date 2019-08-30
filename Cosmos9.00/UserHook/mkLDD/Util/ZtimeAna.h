#include "ZprivateSub.h"
       integer nfrac
       parameter( nfrac=11)
       real*8 frac(nfrac), tf(nfrac)
       data frac/0.05d0, 0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0, 0.6d0,
     *             0.7d0, 0.8d0, 0.9d0, 0.95d0/ 
       real*8 tfary0(nrbin, nfrac)
       real*8 tfary(nrbin,  nfrac)
