      DOUBLE PRECISION FUNCTION GASDEV(Idum)
C***********************************************************************
C     Gaussian deviation
C***********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON /RNDMGAS/ ISET
      SAVE
      DATA ISET/0/      
      gasdev=idum
      IF (ISET.EQ.0) THEN
1       V1=2.D0*S_RNDM(0)-1.D0
        V2=2.D0*S_RNDM(1)-1.D0
        R=V1**2+V2**2
        IF(R.GE.1.D0)GO TO 1
        FAC=SQRT(-2.D0*LOG(R)/R)
        GSET=V1*FAC
        GASDEV=V2*FAC
        ISET=1
      ELSE
        GASDEV=GSET
        ISET=0
      ENDIF
      RETURN
      END
