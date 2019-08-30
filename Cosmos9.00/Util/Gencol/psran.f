!          defined in QGSJETII-03 but not in II-04
       real*8 function psran(seed)
       implicit none
       real*8 seed
       real*8 u
       call  rndc(u)
       psran = u
       end
!=======================================================================
!            GAMFUN is defined in qgsjetII-03 but not in II-04.
      DOUBLE PRECISION FUNCTION GAMFUN(Y)
! Gamma function : See Abramowitz, page 257, form. 6.4.40
!-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION
     +     Y,R,S,T,AFSPL,X,
     +     COEF(10),PI,ZEROD,HALFD,ONED,TWOD,TEND
      SAVE
!
      DATA COEF/8.3333333333333334D-02,-2.7777777777777778D-03,
     .          7.9365079365079365D-04,-5.9523809523809524D-04,
     .          8.4175084175084175D-04,-1.9175269175269175D-03,
     .          6.4102564102564103D-03,-2.9550653594771242D-02,
     .          0.1796443723688306    ,-0.6962161084529506    /
      DATA PI/  3.141592653589793D0/
      DATA ZEROD/0.D0/,HALFD/0.5D0/,ONED/1.D0/,TWOD/2.D0/,TEND/10.D0/
!
      X=Y
      AFSPL=ONED
      N=INT(TEND-Y)
      DO 10 I=0,N
        AFSPL=AFSPL*X
        X=X+ONED
10    CONTINUE
      R=(X-HALFD)* LOG(X)-X+HALFD* LOG(TWOD*PI)
      S=X
      T=ZEROD
      DO 20 I=1,10
        T=T+COEF(I)/S
        S=S*X**2
20    CONTINUE
      GAMFUN = EXP(R+T)/AFSPL
      END
