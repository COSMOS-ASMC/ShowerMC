      subroutine csibSetStblPtcl
      implicit none
!     !       specify short life patciles which should not be
!     decayed inside sibyll.
!     as IDB(x)   = - abs(IDB(x))
!
      INTEGER          IDB,KDEC,LBARP
      REAL             CBR
      COMMON /S_CSYDEC/IDB(49), CBR(102), KDEC(612), LBARP(49)
 


!
!
!      
!  --------------------------------------------------------
!  Particle    SIB PID      SIB2PDG      SIB2PDG^-1    MASS
!  --------------------------------------------------------
!    gam          1            22           0         0.000
!    e+           2           -11           0         0.001
!    e-           3            11           0         0.001
!    mu+          4           -13           0         0.106
!    mu-          5            13           0         0.106
!    pi0          6           111           0         0.135
!    pi+          7           211           0         0.140
!    pi-          8          -211           0         0.140
!    k+           9           321           0         0.494
!    k-          10          -321           0         0.494
!    k0l         11           130           0         0.498
!    k0s         12           310           0         0.498
!!!      IDB(4:12) = -IDB(4:12)   ! in 2.3c this is set
!                            better to take abs and negate
!  
!    p           13          2212           0         0.938
!    n           14          2112           0         0.940
!    nue         15            12           0         0.000
!    nueb        16           -12           0         0.000
!    num         17            14           0         0.000
!    numb        18           -14           0         0.000
!    pbar        19         -2212           0         0.938
!    nbar        20         -2112           0         0.940
!    k0          21           311           0         0.498
!    k0b         22          -311           0         0.498
!    eta         23           221           0         0.548
!         make eta stable
      IDB(23)  = -abs(IDB(23))
!    etap        24           331           0         0.958
!    rho+        25           213           0         0.767
!    rho-        26          -213           0         0.767
!    rho0        27           113           0         0.768
!    k*+         28           323           0         0.892
!    k*-         29          -323           0         0.892
!    k*0         30           313           0         0.896
!    k*0b        31          -313           0         0.896
!    omeg        32           223           0         0.783
!    phi         33           333           0         1.019
      !!     force decay etap etc
      IDB(24:33)  = abs(IDB(24:33))
!    SIG+        34          3222           0         1.189
!    SIG0        35          3212           0         1.193
!    SIG-        36          3112           0         1.197
!    XI0         37          3322           0         1.315
!    XI-         38          3312           0         1.322
!    LAM         39          3122           0         1.116
!       make stable  SIGMA, XI, LAMBDA
!      IDB(34:36) = -IDB(34:36)  ! this is  in 2.3c
!      IDB(39) = -IDB(39)    !          //
      
!    DELT++      40          2224           0         1.231
!    DELT+       41          2214           0         1.235
!    DELT0       42          2114           0         1.234
!    DELT-       43          1114           0         1.233
!    SIG*+       44          3224           0         1.383
!    SIG*0       45          3214           0         1.384
!    SIG*-       46          3114           0         1.387
!    XI*0        47          3324           0         1.532
!    XI*-        48          3314           0         1.535
!    OME-        49          3334           0         1.672
!           force decay      
      IDB(40:49)= abs(IDB(40:49))
      end
