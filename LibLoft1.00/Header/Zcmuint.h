!         muon interaction related variables
!

!          
      real*8 Zeff, Zeff3  ! Zeff**(1/3)
      real*8 muPrVmin, muPrdETX, muPrdE, muPrEmin,
     *       muPrEmax, muPrdU,  muPrEmax1
      integer    muPrUsize, muPrEsize,  muPrTXT
!
      real*8 muBrVmin, muBrdETX, muBrdE, muBrEmin, muBrEmax,
     *       muBrdU, muBrEmax1
      integer muBrUsize, muBrEsize, muBrTXT
!
      real*8 muNVmin, muNdETX, muNdE, muNEmin,
     *  muNEmax, muNdU, muNEmax1
      integer muNUsize, muNEsize, muNTXT
!
      real*8 muNpwtx, muNpwdEdx0, muNpwdEdxt

      real*8 MuPrTX(38), MuPrdEdx0(38), MuPrdEdxt(38)
      real*8 MuBrTX(36), MuBrdEdx0(36), MuBrdEdxt(36)
      real*8 MuNTX(34), MuNdEdx0(34), MuNdEdxt(34)

      real*8 MuNTbl(101, 17)

      real*8  mupa, mura, mupb, muqb, muAk, muAkm, muAkm2,
     *        muPointLike, muShadow, mulogf0


      real*8  muNLEmin, muPrLEmin, muBrLEmin
      common /muintc/ MuNTbl,  MuPrTX, MuPrdEdx0, MuPrdEdxt,
     *    MuBrTX, MuBrdEdx0, MuBrdEdxt, MuNTX, MuNdEdx0, MuNdEdxt,
     *    muNpwtx, muNpwdEdx0, muNpwdEdxt, muPrVmin, muPrdETX, 
     *    muPrdE, muPrEmin, muPrEmax, muPrdU,  muPrEmax1,
     *    muBrVmin, muBrdETX, muBrdE, muBrEmin, muBrEmax,
     *    muBrdU, muBrEmax1, muNVmin, muNdETX, muNdE, muNEmin,
     *    muNEmax, muNdU, muNEmax1,  muNLEmin, muPrLEmin, muBrLEmin,
     *    mupa, mura, mupb, muqb, muAk, muAkm, muAkm2,
     *    muPointLike, muShadow, mulogf0, Zeff, Zeff3,
!
     *    muPrUsize, muPrEsize,  muPrTXT, muBrUsize, muBrEsize,
     *    muBrTXT, muNUsize, muNEsize, muNTXT
