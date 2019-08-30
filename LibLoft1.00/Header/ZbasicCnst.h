!             Basic constant
          real*8 N0  !  Avogadro number
          real*8 r0  !  classical elecron radius. m
          real*8 alpha  ! fine structure const.
          real*8 m2Tomb ! m2 to mb conversion
          real*8 cm2Tomb ! cm2 to mb conversion
          real*8 ar02   ! alpha * r0**2 in mb
          real*8 pir02  ! pi x r0**2 in mb 
       parameter(
     * N0 = 6.0221367d23,
     * r0 = 2.81794092d-15,
     * alpha = 1./137.0359895d0,
     * m2Tomb = 1./1.0d-31,
     * ar02 = alpha * r0**2 * m2Tomb,
     * pir02 = pi * r0**2 * m2Tomb,
     * cm2Tomb = m2Tomb*1.d-4
     * )


 
