!            subset of Zprivate.f
       real bin, rmin  !  rbin in log10.  rmin. for lateral in Moliere u.
       integer smooth
       parameter (bin=0.1, rmin=0.01)
       integer nrbin, nfai
       parameter ( nrbin = 42, nfai=12)
       real*8  rbin(nrbin)
       common /ZprivateSub/  rbin, smooth
