!              constants thru Cosmos
            real*8 pi,          !  3.14..
     *             sqrtpi,      !  sqrt(pi)
     *             Torad,       !  if multiplied to deg  --> radian
     *             Todeg,       !  1/Torad
     *             c,           !  light velocity m/sec
     *             Infty,       !  infinty
     *             Togpcm2,     !  kg/m2 *Togpcm2 --> g/cm2
     *             Tokgpm2,     !  g/cm2 *Tokgpm2    --> kg/m2
     *             Tom,         !  cm *Tom  --> m
     *             Tocm,        !  m*Tocm  --> cm
     *             Tokgpm3,     !  g/cm3 * Tokgpm3 -->  kg/m3
     *             Togpcm3,     !  kg/m3 * Togpcm3 --> g/cm3
     *	           Tonsec,      !  sec *Tonsec --> nsec 
     *             Bcr,         !  Tesla. m^2c^3/eh = 4.414x10^13G=4. x10^9
     *             SyncConvR    !  alpha(mc^2/Lc). GeV/m. conv. rate of synch.r
       real*8 Avogn,            ! Avogadro #. /mol
     *        A2deninv          ! mfp * n* xs = 1;  mfp in kg/m2
                                ! xs in mb.  1/n = A2deninv*A  


        parameter(pi = 3.141592653589793238d0, 
     *    sqrtpi = 1.772453850905516d0, Torad = pi/180.d0, 
     *    Todeg = 180.d0/pi, c=2.998d8, Infty=1.d50, Tom = 1.d-2,
     *    Tocm = 1.d0/Tom,
     *    Togpcm2 = 0.1d0, Tokgpm2 = 1.d0/Togpcm2,
     *	  Tokgpm3 = Tokgpm2/Tom, Togpcm3 = 1.d0/Tokgpm3,
     *    Tonsec = 1.d9, Bcr = 4.414d9, SyncConvR=9.657d6)
       parameter (Avogn=6.0221415d23, A2deninv = 1.d28/Avogn)




  
