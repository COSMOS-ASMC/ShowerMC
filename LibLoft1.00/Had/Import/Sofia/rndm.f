      DOUBLE PRECISION FUNCTION RNDM(IDUMMY)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       IMPLICIT INTEGER (I-N)
c       original one is relaced by rndc
       real(8):: u
       call rndc(u)
       RNDM = u
      END
