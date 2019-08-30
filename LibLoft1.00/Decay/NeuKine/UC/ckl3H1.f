      real*8 function ckl3H1(gzai, el, ml, mpi, mk)
      implicit none
!          K-> mu + numu + pi
      real**  gzai  ! input  f-/f+
      real*8  el  ! input. Energy of muon or electron
      real*8  ml ! input  lepton mass (mu or electron)
      real*8  mpi  ! input pion mass (pi0 or pic)
      real*8  mk   ! input  kaon mass (K0 or kc)
!   function value : output.  eq.(4.3b)/f+**2

      real*8 delta, dm2mke, ml2, mkel

      delta = mk**2 + ml**2

      dm2mke = delta - 2*mk*el
      ml2= ml*ml
      mkel = mk*el
      ckl3H1 = (mkel + ml2/4 * (mkel-ml2) +
     *    gzai**2*ml2/4 *(mkel-ml2)/dm2mke +
     *    2*gzai*ml2/4 *(2*mk**2 + ml2-3*mkel)/
     *    dm2mke  ) * (dm2mke -mpi**2)**2/dm2mke
      end
      real*8 function ckl3lepSpec

