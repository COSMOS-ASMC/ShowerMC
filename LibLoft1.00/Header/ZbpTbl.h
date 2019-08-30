      integer mxBrTXL, mxBrTblLA, mxBrTblLB
      integer mxPrTXL, mxPrTblLA, mxPrTblLB

      integer mxPrTXH, mxPrTblH
      integer mxBrTXS, mxBrTblSA,  mxBrTblSB

!           for muon interaction
      integer mxMuNTX, mxMuBrTX, mxMuPrTX
      integer mxMuNTbl, mxMuBrTbl, mxMuPrTbl
      integer mxMuNUsize, mxMuNEsize, mxMuBrUsize, 
     *        mxMuBrEsize, mxMuPrUsize,mxMuPrEsize


      parameter ( mxBrTXL = 60, mxBrTblLA=mxBrTXL*51,  ! v8.0
     *            mxBrTblLB = mxBrTblLA)
      parameter ( mxPrTXL = 60, mxPrTblLA=mxPrTXL*51,  ! v8.0
     *            mxPrTblLB = mxPrTblLA)
!      parameter ( mxBrTXH = 60, mxBrTblHA=mxBrTXH*51,   ! v8.0
!     *            mxBrTblHB = mxBrTblHA)
	integer,parameter:: mxBrTXH = 81            ! v9.15
	integer,parameter:: mxBrTblHA=mxBrTXH*51        ! //
        integer,parameter:: mxBrTblHB = mxBrTblHA       ! //

      parameter ( mxPrTXH = 60, mxPrTblH=mxPrTXH*51)    ! v8.0

      parameter ( mxBrTXS = 85, mxBrTblSA = mxBrTXS*51,  ! v8.0
     *            mxBrTblSB = mxBrTblSA)
!         next is when Seltzer table is used upto 10GeV
      integer,parameter:: mxBrTXS2=41, mxBrTblSA2=mxBrTXS2*51, 
     *   mxBrTblSB2=mxBrTblSA2  !!! v9.135


  
!////////////
!              for muon
      parameter ( mxMuNTX =51, mxMuBrTX=51, mxMuPrTX= 51 )
      parameter ( mxMuNUsize = 101, mxMuNEsize = 51, 
     *            mxMuBrEsize = 51, mxMuPrEsize= 51 )
      parameter ( mxMuBrUsize = 1, mxMuPrUsize = 1 )  ! not used
      parameter ( mxMuNTbl = mxMuNUsize* mxMuNEsize)
      parameter ( mxMuBrTbl = mxMuBrUsize* mxMuBrEsize)
      parameter ( mxMuPrTbl = mxMuPrUsize* mxMuPrEsize)

       type bpTbl     
       sequence
!          real*8 BrdEdx0(mxBrTXL)        ! for x*dE/dx(v<100 keV)
	real*8 BrTXL(mxBrTXL,2)
          real*8 BrSTLA(mxBrTblLA, 1)
          real*8 BrSTLB(mxBrTblLB, 1)

	real*8 BrTXH(mxBrTXH,2) 
          real*8 BrSTHA(mxBrTblHA, 1)
          real*8 BrSTHB(mxBrTblHB, 1)


	real*8 PrTXL(mxPrTXL)
          real*8 PrSTLA(mxPrTblLA, 1)
          real*8 PrSTLB(mxPrTblLB, 1)

	real*8 PrTXH(mxPrTXH)
          real*8 PrSTH(mxPrTblH, 1)     
        
	real*8 BrTXS(mxBrTXS,2)
          real*8 BrSTSA(mxBrTblSA, 1)
          real*8 BrSTSB(mxBrTblSB, 1)
!              v9.135
	real*8 BrTXS2(mxBrTXS2,2)
          real*8 BrSTSA2(mxBrTblSA2, 1)
          real*8 BrSTSB2(mxBrTblSB2, 1)
!              muon
          real*8 MuNTX(mxMuNTX)
          real*8 MuNdEdx0(mxMuNTX)
          real*8 MuNdEdxt(mxMuNTX)


          real*8 MuBrTX(mxMuBrTX)
          real*8 MuBrdEdx0(mxMuBrTX)
          real*8 MuBrdEdxt(mxMuBrTX)

          real*8 MuPrTX(mxMuPrTX)
          real*8 MuPrdEdx0(mxMuPrTX)
          real*8 MuPrdEdxt(mxMuPrTX)

!          real*8 MuNTbl(mxMuNUsize, mxMuNEsize)
          real*8 MuNTbl(mxMuNTbl, 1)
!          real*8 MuBrTbl(mxMuBrUsize, mxMuBrEsize)
          real*8 MuBrTbl(mxMuBrTbl, 1)
!          real*8 MuPrTbl(mxMuPrUsize, mxMuPrEsize)
          real*8 MuPrTbl(mxMuPrTbl, 1)
          
       end type bpTbl

  
