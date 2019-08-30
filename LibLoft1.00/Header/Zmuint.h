         real*8  Z          !  Z
         real*8  A          !  A
         real*8  A3         !  A**(1/3)
         real*8  Z2         !  Z^2
         real*8  Zp2        !  Z(Z+1.2)
         real*8  D          !  alpha*(2r0Me/Mmu)**2
         real*8  LogZ, Z3   !  ln(Z), Z^(1/3)
         real*8  Shadow     !  shadowing factor A**(-0.1)
         real*8  PointLike  !  fraction of point like interaction 
         real*8  Ak         !  log factor 189 (or 191 by Rozental)
         real*8  Akm        !  Ak * MubyMe/Z3
         real*8  Akm2       !  Ak * sqrt(e)
         real*8  Gzai       !  gzai factor. longitudinally / transvers
                            !         polarized meson absorption XS.
                            !  
!


         common /zmucom/
     *   Z, A, Z2, Zp2, D, LogZ, Z3, Shadow,
     *   PointLike, Ak, Akm, Akm2, A3, Gzai
