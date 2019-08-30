c       by this you can learn how to call cross-section subroutines
c       in phits.  the corss-section calculation part is
c       extracted and saved as phitsXsec.f in Com/ directoryr
c  
      subroutine getflt(fpl)
*                                                                      *
*       flight path length sampling                                    *
*       last modified by K.Niita on 10/05/2000                         *
*       for proton, total cross section is the max value.              *
*                                                                      *
*     output:                                                          *
*                                                                      *
*       fpl     : sampling flight path length (1/cm)                   *
*                                                                      *
*       totttl  : total flight length (1/cm)                           *
*       tothyd  : Hydrogen flight length (1/cm)                        *
*       totdcy  : decay fligt length (1/cm)                            *
*                                                                      *
*       siggn(ksign+lem) : nonelastic flight length for (lem,mat)      *
*       sigge(ksige+lem) : elastic flight length for (lem,mat)         *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param.inc'
      include 'mmbank.inc'

      parameter ( rpmass = 938.27, rnmass = 939.58 )
      parameter ( nemd = 20 )

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp

      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)
      common /tlgeom/ iblz1,iblz2
      common /mathzn/ mathz, mathn, jcoll, kcoll
      common /qparm/  ielas,icasc,iqstep,lvlopt,igamma
      common /xtotal/ totttl,tothyd,totdcy
      common /eparm/  esmax, esmin, emin(20)
      common /ndemax/ dnmax(20)
      common /rcmini/ rcmin(20)
      common /emode/ ge1, ge2, iemode, nevhin, nlowlv, nefiss, ntwidt

      common /celdn/  denr(kvlmax), denm(kvlmax), denc(kvlmax)
      common /regcm/  icmg(kvlmax)

      common /kmat1g/ kmat(kvlmax)
      common /xgeosm/ ksig(kvlmax)
      common /xinels/ ksige, ksign

      dimension dnel(1), denh(1), zz(1), a(1), den(1)
      equivalence ( das, dnel, denh, zz, a, den )
      dimension sigg(1), siggm(1)
      equivalence ( das, sigg, siggm )
      dimension siggn(1), sigge(1)
      equivalence ( das, siggn, sigge )

*-----------------------------------------------------------------------

      common /clustt/ nclsts, iclusts(nnn)
      common /clustw/ jclusts(0:7,nnn),  qclusts(0:12,nnn)

      dimension crsen(kvlmax), crstl(kvlmax)

*-----------------------------------------------------------------------

      data rcc / 2.997925d+10 /
      data pi / 3.1415926535897932d0 /

*-----------------------------------------------------------------------
*        initial values and constants
*-----------------------------------------------------------------------
               fpl    = 1.0d+10

               totttl = 0.d0
               tothyd = 0.d0
               totdcy = 0.d0

               ein    = e(ibke+no)
               kf     = nkf(ibknkf+no)

               tme    = abs(t(ibkt+no)) / 10.0d0
               icl    = idgr(iblz1)

               ieletr = 0
c/////////////
      write(0,*) 'getflt; ein=', ein, ' kf=',kf, ' ityp=',ityp
c///////////// 

*-----------------------------------------------------------------------
*        forced collision or reaction cutoff option
*-----------------------------------------------------------------------

            if( nfcs(ibknfc+no) .eq. 1 ) return

            if( ityp .lt. 15 ) then

               emint = rcmin(ityp)

            else

               emint = rcmin(ityp) * dble( mtyp )

            end if
c//////
            write(0,*) ' getflt; emint=',emint
c//////////////
            if( ein .lt. emint ) return

*-----------------------------------------------------------------------
*        totdcy : decay flight legth
*-----------------------------------------------------------------------

            totdcy = dcyfl(kf,ein,rtyp)

            if( totdcy .gt. 1.0d+10 ) then

               fpl = 1.0d-10
               return

            end if

*-----------------------------------------------------------------------
*     for no void
*-----------------------------------------------------------------------

      if( mat .gt. 0 ) then

               lem   = nint( dnel(kmat(mat)+1) )
               hydro = denh(kmat(mat)+2)

*-----------------------------------------------------------------------
*        non nucleus
*-----------------------------------------------------------------------

         if( ityp .lt. 15 ) then

*-----------------------------------------------------------------------
*        neutron for nuclear data
*-----------------------------------------------------------------------

            if( ityp .eq. 2 .and. ein .lt. dnmax(2) ) then

                  if( iemode .eq. 0 ) then

                     jcoll = 6

                  else

                     jcoll = 10

                  end if

                     mk = mat
                     rh = denm(mat)

               call acetot(0,sigt,sigaa,icl,ein,tme,mk)

                     totttl = totttl + sigt * rh

*-----------------------------------------------------------------------
*        photon for nuclear data
*-----------------------------------------------------------------------

            else if( ityp .eq. 14 .and. ein .lt. dnmax(14) ) then

                     jcoll = 7

                     mk = mat
                     rh = denm(mat)

               call photot(0,sigt,sigaa,ein,mk)

                     totttl = totttl + sigt * rh

*-----------------------------------------------------------------------
*        electron for nuclear data
*-----------------------------------------------------------------------

            else if( ( ityp .eq. 12 .or. ityp .eq. 13 ) .and.
     &                 ein .lt. dnmax(12) ) then

                     jcoll = 8

                     mk = mat
                     rh = denm(mat)

               call eletot(pmff,icl,ein,mk,rh)

                     ieletr = 1

                     totttl = totttl + 1.0 / pmff

*-----------------------------------------------------------------------
*        for high energy part
*-----------------------------------------------------------------------

            else

*-----------------------------------------------------------------------
*              Hydrogen cross section by JAM
*-----------------------------------------------------------------------

                     sig   = 0.0
                     sigel = 0.0

               if( hydro .gt. 0.0d0 ) then

                     kf1 = kf
                     kf2 = 2212

                     call sigjam(kf1,kf2,ein,sig,sigel,signo)

               end if

                     tothyd = sig * hydro

*-----------------------------------------------------------------------
*              proton cross sections by sigrc or nuclear data
*-----------------------------------------------------------------------

               if( ityp .eq. 1 ) then

                  if( nfcs(ibknfc+no) .ne. 2 ) then
*                       --- proton pseudo for all ---

                     do i = 1, lem

                        siggn(ksign+i) = 0.0
                        sigge(ksige+i) = 0.0

                        totttl = totttl + siggm(ksig(mat)+(i-1)*2+2)

                     end do

                  else if( ein .ge. dnmax(1) ) then
*                       --- proton real for high energy ---

                     do i = 1, lem

                        itz = nint( zz(kmat(mat)+(i-1)*3+3) )
                        ita = nint(  a(kmat(mat)+(i-1)*3+4) )

                        call sigrc(ityp,ein,ita,itz,sigt,signe,sigel)

                        if( ielas .le. 1 ) sigel = 0.0

                        signe = signe * den(kmat(mat)+(i-1)*3+5)
                        sigel = sigel * den(kmat(mat)+(i-1)*3+5)
                        sigt  = signe + sigel

                        siggn(ksign+i) = signe
                        sigge(ksige+i) = sigel

                        totttl = totttl + sigt

                     end do

                  else if( ein .lt. dnmax(1) ) then
*                       --- proton real for nuclear data ---

                        jcoll = 9

                        mk = mat
                        rh = denm(mat)

                        call sig_tot(ityp,sigt,sigaa,ein,mk)

                        totttl = totttl + sigt * rh

                  end if

*-----------------------------------------------------------------------
*              total, elastic and nonelastic cross sections by sigrc
*-----------------------------------------------------------------------

               else if( ityp .eq. 2 ) then
*                       --- neutron ---

                     do i = 1, lem

                        itz = nint( zz(kmat(mat)+(i-1)*3+3) )
                        ita = nint(  a(kmat(mat)+(i-1)*3+4) )

                        call sigrc(ityp,ein,ita,itz,sigt,signe,sigel)

                        if( ielas .eq. 0 ) sigel = 0.0

                        signe = signe * den(kmat(mat)+(i-1)*3+5)
                        sigel = sigel * den(kmat(mat)+(i-1)*3+5)
                        sigt  = signe + sigel

                        siggn(ksign+i) = signe
                        sigge(ksige+i) = sigel

                        totttl = totttl + sigt

                     end do

               else if( ityp .eq. 6 .or. ityp .eq. 7 ) then
*                       --- muon ---

                     do i = 1, lem

                        siggn(ksign+i) = 0.0
                        sigge(ksige+i) = 0.0

                     end do

               else if( ityp .eq. 12 .or. ityp .eq. 13 ) then
*                       --- electron or positron ---

                     do i = 1, lem

                        siggn(ksign+i) = 0.0
                        sigge(ksige+i) = 0.0

                     end do

               else if( ityp .eq. 14 ) then
*                       --- photon ---

                     do i = 1, lem

                        siggn(ksign+i) = 0.0
                        sigge(ksige+i) = 0.0

                     end do

               else if( ityp .eq. 11 .and. abs(kf) .lt .100 ) then
*                       --- lepton ---

                     do i = 1, lem

                        siggn(ksign+i) = 0.0
                        sigge(ksige+i) = 0.0

                     end do

               else
*                       --- other hadron ---

                  do i = 1, lem

                        siggn(ksign+i) = sigg(ksig(mat)+(i-1)*2+1)
                        sigge(ksige+i) = 0.0

                        totttl = totttl + siggn(ksign+i)

                  end do

               end if

            end if

*-----------------------------------------------------------------------
*        transport particle is nucleus
*-----------------------------------------------------------------------

         else if( ityp .ge. 15 .and. ityp .le. 19 ) then

                     jta = mtyp
                     jtz = jtyp

*-----------------------------------------------------------------------
*           --- nucleus pseudo ---
*-----------------------------------------------------------------------

            if( nfcs(ibknfc+no) .ne. 2 ) then

*-----------------------------------------------------------------------
*              nucleus - proton
*-----------------------------------------------------------------------

                     signe = 0.0
                     sigel = 0.0

               if( hydro .gt. 0.0d0 ) then

                     itypn = 1
                     einn  = rpmass / rtyp * ein

                  if( einn .gt. esmin ) then

                        edel = log( einn / esmin ) / ( nemd - 1 )

                     do k = 1, nemd

                        ei = exp( ( k - 1 ) * edel ) * esmin

                        call sigrc(itypn,ei,jta,jtz,ssigt,
     &                             ssigne,ssigel)

                        if( ssigne .gt. signe ) signe = ssigne

                     end do

                  end if

               end if

                     tothyd = signe * hydro

*-----------------------------------------------------------------------
*              nucleus - nucleus
*-----------------------------------------------------------------------

                     ap = dble( jta )
                     zp = dble( jtz )

               if( ein / ap .gt. esmin ) then

                     edel = log( ein / esmin / ap ) / ( nemd - 1 )

                  do i = 1, lem

                     zt = zz(kmat(mat)+(i-1)*3+3)
                     at =  a(kmat(mat)+(i-1)*3+4)

                     signe = 0.0
                     sigel = 0.0

                     do k = 1, nemd

                        ei = exp( ( k - 1 ) * edel ) * esmin * ap

                        call sighi(ap,zp,ei,at,zt,ssigne,ssigel,bmax)

                        if( ssigne .gt. signe ) signe = ssigne
                        if( ssigel .gt. sigel ) sigel = ssigel

                     end do

                     if( ielas .le. 2 ) sigel = 0.0

                     signe = signe * den(kmat(mat)+(i-1)*3+5)
                     sigel = sigel * den(kmat(mat)+(i-1)*3+5)
                     sigt  = signe + sigel

                     siggn(ksign+i) = signe
                     sigge(ksige+i) = sigel

                     totttl = totttl + sigt

                  end do

               end if

*-----------------------------------------------------------------------
*           --- nucleus real ---
*-----------------------------------------------------------------------

            else

*-----------------------------------------------------------------------
*              nucleus - proton
*-----------------------------------------------------------------------

                     signe = 0.0
                     sigel = 0.0

               if( hydro .gt. 0.0d0 ) then

                     itypn = 1
                     einn  = rpmass / rtyp * ein

                     call sigrc(itypn,einn,jta,jtz,sigt,signe,sigel)

               end if

                     tothyd = signe * hydro

*-----------------------------------------------------------------------
*              nucleus - nucleus
*-----------------------------------------------------------------------

                     ap = dble( jta )
                     zp = dble( jtz )

                  do i = 1, lem

                     zt = zz(kmat(mat)+(i-1)*3+3)
                     at =  a(kmat(mat)+(i-1)*3+4)

                     call sighi(ap,zp,ein,at,zt,signe,sigel,bmax)

                     if( ielas .le. 2 ) sigel = 0.0

                     signe = signe * den(kmat(mat)+(i-1)*3+5)
                     sigel = sigel * den(kmat(mat)+(i-1)*3+5)
                     sigt  = signe + sigel

                     siggn(ksign+i) = signe
                     sigge(ksige+i) = sigel

                     totttl = totttl + sigt

                  end do

            end if

*-----------------------------------------------------------------------

         end if

      end if

*-----------------------------------------------------------------------
*     flight path length by exp. random number
*     distance mesh for electrons
*-----------------------------------------------------------------------

      if( ieletr .eq. 0 ) then

         if( totttl + tothyd + totdcy .le. 0.0d0 ) return

         if( nfcs(ibknfc+no) .eq. 0 ) then

            fpl = exprnf(dummy) / ( totttl + tothyd + totdcy )

         else

            fpl = xfcs(ibkxfc+no)

         end if

      else

            fpl = pmff

      end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      function dcyfl(kf,ein,rtyp)
*                                                                      *
*       decay flight length                                            *
*       last modified by K.Niita on 2002/10/15                         *
*                                                                      *
*     output:                                                          *
*                                                                      *
*       dcyfl : decay flight length (cm)                               *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      data rcc / 2.997925d+10 /

      common /dcyd/ rtdcy(50), kfdcy(50), ndcy

*-----------------------------------------------------------------------

               dcyfl = 0.0d0
               tlife = -1.0

*-----------------------------------------------------------------------

         do k = 1, ndcy

            if( kf .eq. kfdcy(k) ) then

               tlife = rtdcy(k)
               goto 100

            end if

         end do


  100    continue

*-----------------------------------------------------------------------
*        decay width
*-----------------------------------------------------------------------

         if( tlife .gt. 0.0 ) then

               dcyfl = 1.0 / tlife / rcc
     &               / dsqrt((( ein / rtyp ) + 2.0 )
     &               * ( ein / rtyp ) )

         end if


*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      function idcyc(kf)
*                                                                      *
*       check decay particle                                           *
*       idcyc = 1 : decay particle                                     *
*       idcyc = 0 : no decay particle                                  *
*       last modified by K.Niita on 2002/10/15                         *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      common /dcyd/ rtdcy(50), kfdcy(50), ndcy

*-----------------------------------------------------------------------

         idcyc = 0

*-----------------------------------------------------------------------

         do k = 1, ndcy

            if( kf .eq. kfdcy(k) ) then

               idcyc = 1
               return

            end if

         end do

      return
      end


************************************************************************
*                                                                      *
      function pseudpr(ityp)
*                                                                      *
*       for proton, check pseudo collision                             *
*       last modified by K.Niita on 01/06/2001                         *
*                                                                      *
*     output:                                                          *
*                                                                      *
*       pseudpr  = .true.  : pseudo collision                          *
*                  .false. : real collision                            *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      logical pseudpr

      common /xtotal/ totttl,tothyd,totdcy

*-----------------------------------------------------------------------

            pseudpr = .false.

*-----------------------------------------------------------------------
*        proton
*-----------------------------------------------------------------------

         if( ityp .eq. 1 ) then

            tottts = totttl
            tothys = tothyd

            call getfltp

            pseud = ( totttl + tothyd + totdcy )
     &            / ( tottts + tothys + totdcy )

            if( unirn(dummy) .gt. pseud ) pseudpr = .true.

*-----------------------------------------------------------------------
*        nucleus
*-----------------------------------------------------------------------

         else if( ityp .ge. 15 ) then

            tottts = totttl
            tothys = tothyd

            call getfltn

            pseud = ( totttl + tothyd + totdcy )
     &            / ( tottts + tothys + totdcy )

            if( unirn(dummy) .gt. pseud ) pseudpr = .true.

         end if

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine getfltp
*                                                                      *
*       for proton, get total flight path length at ec(no)             *
*       last modified by K.Niita on 11/04/2000                         *
*                                                                      *
*     output:                                                          *
*                                                                      *
*       totttl  : total flight length (1/cm)                           *
*                                                                      *
*       siggn(ksign+lem) : nonelastic flight length for (lem,mat)      *
*       sigge(ksige+lem) : elastic flight length for (lem,mat)         *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'mmbank.inc'

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp
      common /tlgeom/ iblz1,iblz2
      common /mathzn/ mathz, mathn, jcoll, kcoll
      common /qparm/  ielas,icasc,iqstep,lvlopt,igamma
      common /xtotal/ totttl,tothyd,totdcy
      common /kmat1g/ kmat(kvlmax)
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /xinels/ ksige, ksign
      common /ndemax/ dnmax(20)

      common /celdn/  denr(kvlmax), denm(kvlmax), denc(kvlmax)
      common /regcm/  icmg(kvlmax)

      dimension dnel(1), denh(1), zz(1), a(1), den(1)
      equivalence ( das, dnel, denh, zz, a, den )
      dimension siggn(1), sigge(1)
      equivalence ( das, siggn, sigge )

*-----------------------------------------------------------------------
*        only for proton
*-----------------------------------------------------------------------

         if( ityp .ne. 1 ) return

*-----------------------------------------------------------------------
*        initial values and constants
*-----------------------------------------------------------------------

               totttl = 0.0d0
               tothyd = 0.0d0

               if( mat .le. 0 ) return

*-----------------------------------------------------------------------

               ein = ec(ibkec+no)

*-----------------------------------------------------------------------
*        cross section by sigrc
*-----------------------------------------------------------------------

         if( ein .ge. dnmax(1) ) then

*-----------------------------------------------------------------------
*           Hydrogen cross section by JAM
*-----------------------------------------------------------------------

                  sig   = 0.0
                  sigel = 0.0
                  hydro = denh(kmat(mat)+2)

               if( hydro .gt. 0.0d0 ) then

                  kf1 = 2212
                  kf2 = 2212

                  call sigjam(kf1,kf2,ein,sig,sigel,signo)

               end if

                  tothyd = sig * hydro

*-----------------------------------------------------------------------
*           total, elastic and nonelastic cross sections by sigrc
*-----------------------------------------------------------------------

                  lem  = nint( dnel(kmat(mat)+1) )

            do i = 1, lem

                  itz = nint( zz(kmat(mat)+(i-1)*3+3) )
                  ita = nint(  a(kmat(mat)+(i-1)*3+4) )

               call sigrc(ityp,ein,ita,itz,sigt,signe,sigel)

                  if( ielas .le. 1 ) sigel = 0.0

                  signe = signe * den(kmat(mat)+(i-1)*3+5)
                  sigel = sigel * den(kmat(mat)+(i-1)*3+5)
                  sigt  = signe + sigel

                  siggn(ksign+i) = signe
                  sigge(ksige+i) = sigel

                  totttl = totttl + sigt

            end do

                  jcoll = 0

*-----------------------------------------------------------------------
*        cross section by nuclear data
*-----------------------------------------------------------------------

         else if( ein .lt. dnmax(1) ) then

                  icl  = idgr(iblz1)

                  mk = mat
                  rh = denm(mat)

                  call sig_tot(ityp,sigt,sigaa,ein,mk)

                  totttl = totttl + sigt * rh

                  jcoll = 9

         end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine getfltn
*                                                                      *
*       for nucleus, get total flight path length at ec(no)            *
*       last modified by K.Niita on 2002/01/15                         *
*                                                                      *
*     output:                                                          *
*                                                                      *
*       totttl  : total flight length (1/cm)                           *
*                                                                      *
*       siggn(ksign+lem) : nonelastic flight length for (lem,mat)      *
*       sigge(ksige+lem) : elastic flight length for (lem,mat)         *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'mmbank.inc'

      parameter ( rpmass = 938.27, rnmass = 939.58 )

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp

      common /qparm/  ielas,icasc,iqstep,lvlopt,igamma
      common /xtotal/ totttl,tothyd,totdcy

      common /kmat1g/ kmat(kvlmax)
      common /xinels/ ksige, ksign

      dimension dnel(1), denh(1), zz(1), a(1), den(1)
      equivalence ( das, dnel, denh, zz, a, den )
      dimension siggn(1), sigge(1)
      equivalence ( das, siggn, sigge )

*-----------------------------------------------------------------------
*        only for nucleus
*-----------------------------------------------------------------------
c//////////
      write(0,*) ' getfltn , ityp=', ityp
c////////////
         if( ityp .lt. 15 ) return

*-----------------------------------------------------------------------
*        initial values and constants
*-----------------------------------------------------------------------

               totttl = 0.0d0
               tothyd = 0.0d0

               if( mat .le. 0 ) return

               ein    = ec(ibkec+no)
               lem    = nint( dnel(kmat(mat)+1) )
               hydro  = denh(kmat(mat)+2)

*-----------------------------------------------------------------------

               jta = mtyp
               jtz = jtyp
c////////////
               write(0,*) ' ein=',ein,  ' mat =',mat,
     *         ' lem=',lem, 'hydro=',hydro
c//////////
*-----------------------------------------------------------------------
*              nucleus - proton
*-----------------------------------------------------------------------

                     signe = 0.0
                     sigel = 0.0

               if( hydro .gt. 0.0d0 ) then

                     itypn = 1
                     einn  = rpmass / rtyp * ein

                     call sigrc(itypn,einn,jta,jtz,sigt,signe,sigel)

               end if

                     tothyd = signe * hydro

*-----------------------------------------------------------------------
*              nucleus - nucleus
*-----------------------------------------------------------------------

                     ap = dble( jta )
                     zp = dble( jtz )

                  do i = 1, lem

                     zt = zz(kmat(mat)+(i-1)*3+3)
                     at =  a(kmat(mat)+(i-1)*3+4)

                     call sighi(ap,zp,ein,at,zt,signe,sigel,bmax)

                     if( ielas .le. 2 ) sigel = 0.0

                     signe = signe * den(kmat(mat)+(i-1)*3+5)
                     sigel = sigel * den(kmat(mat)+(i-1)*3+5)
                     sigt  = signe + sigel

                     siggn(ksign+i) = signe
                     sigge(ksige+i) = sigel

                     totttl = totttl + sigt

                  end do

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine sigrc(incp,emev,ia,iz,sigt,sigr,sigs)
*                                                                      *
*        calculates total, nonelastic and elastic cross-sections       *
*                                                                      *
*     input:                                                           *
*        incp   : =1, proton, =2, neutron                              *
*        emev   : incident nucleon energy (MeV)                        *
*        ia     : mass number                                          *
*        iz     : charge number                                        *
*                                                                      *
*     output:                                                          *
*       sigt    : total cross-section (b)                              *
*       sigr    : nonelastic cross-section (b)                         *
*       sigs    : sigt-sigr=elastic scattering cross-section (b)       *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)
      dimension       par(4),c(4,3)

      data (c(1,k),k=1,3)/1.353,0.4993,1.076/
      data (c(2,k),k=1,3)/0.6653,0.9900,1.040/
      data (c(3,k),k=1,3)/2.456,2.060,0.9314/
      data (c(4,k),k=1,3)/0.1167,1.623,0.9736/

*-----------------------------------------------------------------------
c////////////
      write(0,*)' sigrc: incp,emev,ia,iz=',  incp,emev,ia,iz
c////////////
         rz = 0.14

         a = float(ia)
         z = float(iz)

         a13 = a**0.3333333
         a23 = a**0.6666667
         alg = log(a)

         ecst = 82.0
         acst = 238.0

*-----------------------------------------------------------------------

         c1 =  0.0825
         c2 = -0.0057
         c3 =  0.14
         c4 = -0.2
         c5 =  2.72
         c6 =  1.62
         c7 = -5.3

*-----------------------------------------------------------------------

         ap = 1.0
         pi = 3.1415927

*-----------------------------------------------------------------------
*     modified by Niita
*-----------------------------------------------------------------------

            cmev = 0.0575 * a + 12.31

            facp = min( 1.0d0, 0.684 + 1.327e-3 * a )
            facc = 1.0 - facp

*-----------------------------------------------------------------------

            g1 = 1.0 - 0.62 * exp( -cmev / 200. )
     &                      * sin( 10.9 / cmev**0.28 )
            g1 = g1 * 1.1
            g2 = 1.0 + 0.016 * sin( 5.3 - 2.63 * log(a) )

            sigpa = 0.045 * a**0.7 * g1 * g2

*-----------------------------------------------------------------------

            eroot = sqrt(cmev)
            rad   = rz * a**0.3333
            wvl   = 0.1 * 1.22 * ( a + ap ) / a * sqrt(14.1) / eroot

            sigcr = pi * ( rad + wvl )**2

*-----------------------------------------------------------------------

            sigpc = facp * sigcr + facc * sigpa

            ftpa = sigpc / sigpa
            ftcr = sigpc / sigcr

            enpa = ftpa * 1.1 - 1.0

*-----------------------------------------------------------------------

         do k=1,4

            par(k) = log(c(k,1))+log(c(k,2))*alg+log(c(k,3))*alg**2
            par(k) = exp(par(k))

         end do

            epk   = par(3) * a13
            fexp1 = par(2) * log(epk/cmev)

            sigtp  = sigpc * ( 1.0 + par(4) )
     &             + par(1) * a13 * exp( -fexp1**2 )

            esub = ecst * a13 / acst**0.3333
            epk2 = epk - esub

         if( epk2 .gt. 0.0 ) then

            fexp2 = par(2) * log( epk2 / cmev )

            sigtp = sigtp + par(1) * a13 * exp( -fexp2**2 )

         end if

*-----------------------------------------------------------------------

            if( eroot .lt. 3.0 ) as = 1.0 - c1 * ( 3.0 - eroot )**2
            if( eroot .ge. 3.0 ) as = 1.0

            p = c2 * eroot + c3
            q = c4 * eroot + c5

            if( eroot .ge. 4.0 ) r = 1.2
            if( eroot .lt. 4.0 ) r = c6 * eroot + c7

            sigtc = 2.0 * sigpc * ( as - p * cos( q * a**0.3333 - r ) )

*-----------------------------------------------------------------------

            facp = min( 1.0d0, 0.578 + 1.77e-3 * a )
            facc = 1.0 - facp

            sigtpc = facp * sigtc + facc * sigtp

            fttp = sigtpc / sigtp
            fttc = sigtpc / sigtc

*-----------------------------------------------------------------------
*     calculation of sigr as
*     J.Letaw etal., Astrophys. J. Supp. series 51,271-276(1983).
*     sigrb = asymptotic high energy reaction cross-section
*-----------------------------------------------------------------------

         sigrb = 0.045 * a**0.7

*-----------------------------------------------------------------------
*     energy dependent factor
*-----------------------------------------------------------------------

         f1 = 1.0 - 0.62 * exp( -emev / 200. )
     &                   * sin( 10.9 / emev**0.28 )

*-----------------------------------------------------------------------
*     Pearlstein added low energy enhancement factor of 10%
*-----------------------------------------------------------------------

         f1 = f1 * ( 1.0 + enpa
     &      * exp( -min( 50.d0, ( emev - cmev ) / 10. )) )


c        f1 = f1 * ( 1.0 + 0.1 * exp( -( emev - 20. ) / 10. ) )

*-----------------------------------------------------------------------
*     mass dependent factor
*-----------------------------------------------------------------------

         f2 = 1.0 + 0.016 * sin( 5.3 - 2.63 * log(a) )

*-----------------------------------------------------------------------

         sigr = sigrb * f1 * f2

*-----------------------------------------------------------------------
*     sigt according to fit by S. Pearlstein, June 86.
*-----------------------------------------------------------------------

      if( emev .ge. cmev ) then

         do 210 k=1,4

            par(k) = log(c(k,1))+log(c(k,2))*alg+log(c(k,3))*alg**2
            par(k) = exp(par(k))

  210    continue

            epk   = par(3) * a13
            fexp1 = par(2) * log(epk/emev)

            sigt  = sigr * ( 1.0 + par(4) )
     &            + par(1) * a13 * exp( -fexp1**2 )

            ecst = 82.0
            acst = 238.
            esub = ecst * a13 / acst**0.3333
            epk2 = epk - esub

         if( epk2 .gt. 0.0 ) then

            fexp2 = par(2) * log( epk2 / emev )

            sigt = sigt + par(1) * a13 * exp( -fexp2**2 )

         end if

            factp = 1.0 + ( fttp - 1.0 )
     &            * exp( - min( 50.d0,( emev - cmev ) / 10.0 ))

            sigt = sigt * factp

*-----------------------------------------------------------------------
*     calculates xsects according to Ramsauer effect
*     ref: Angeli and Csikai, NP A170,577-583(1971)
*     sigt/signe=as-p*cos(q*a**0.3333-r)
*-----------------------------------------------------------------------

      else

            eroot = sqrt(emev)

            if( eroot .lt. 3.0 ) as = 1.0 - c1 * ( 3.0 - eroot )**2
            if( eroot .ge. 3.0 ) as = 1.0

            p = c2 * eroot + c3
            q = c4 * eroot + c5

            if( eroot .ge. 4.0 ) r = 1.2
            if( eroot .lt. 4.0 ) r = c6 * eroot + c7

            rad = rz * a**0.3333
            wvl = 0.1 * 1.22 * ( a + ap ) / a * sqrt(14.1) / eroot

         facr = 1.0 + ( ftcr - 1.0 )
     &        * exp( - min( 50.d0,( cmev - emev ) / 10.0 ))

            sigr = pi * ( rad + wvl )**2 * facr
            sigt = 2.0 * sigr * ( as - p * cos( q * a**0.3333 - r ) )

            factc = 1.0 + ( fttc - 1.0 )
     &            * exp( - min( 50.d0, ( cmev - emev ) / 10.0 ))

            sigt = sigt * factc

      end if

*-----------------------------------------------------------------------
*     coulomb factor
*     ec, coulomb barrier half-height, similar to s. pearlstein,
*     j. nuc. energy 23,87(1975) using optical model inverse p cs's.
*-----------------------------------------------------------------------

      if( incp .eq. 1 .and. emev .lt. 200.0 ) then

            ec = 1.44 * z / ( 8.2 + 0.68 * a13 )

            w1    = 3.816 + 0.1974 * z
            fexc1 = exp( max( -50.d0, ( ec - emev ) / w1 ))
            qfac1 = 1.0 / ( 1.0 + fexc1 )

            w2    = 0.07246 * z + 6.058
            if( z .lt. 10.0 ) w2 = 12.0
            qfac2 =  1.0 - exp( - min( 50.d0, ( emev / w2 )**2 ))

            w3    = 2.0
            fexc3 = exp( max( -50.d0, ( ec - emev ) / w3 ))
            qfac3 = 1.0 / ( 1.0 + fexc3 )

            fcoul = qfac1 * qfac2 * qfac3

            sigt = sigt * fcoul
            sigr = sigr * fcoul

      end if

*-----------------------------------------------------------------------

            sigs = sigt - sigr

*-----------------------------------------------------------------------
c/////////////
            write(0,*) 'sigrc exit: sigt, sigr, sigs=',
     *     sigt, sigr, sigs
c///////////
      return
      end


************************************************************************
*                                                                      *
      subroutine seldsd(incp,ex,ata,atz,csthcm,ifirst)
*                                                                      *
*        determine csthcm according to the elastic angular             *
*        distribution by K. Niita's systematics                        *
*                                                                      *
*     input:                                                           *
*        incp     : 1->proton, 2->neutron                              *
*        ex       : kinetic energy of incident nucleon in lab (MeV)    *
*        ata      : target mass number                                 *
*        atz      : target charge                                      *
*        ifirst   : frag of multiple calculations for the same system  *
*                                                                      *
*     output:                                                          *
*        csthcm   : cosine of scattering angle in cm                   *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      parameter ( pi = 3.1415926535897932 )

*-----------------------------------------------------------------------

      common /paraj/ mstz(200), parz(200)

      dimension       dsiga(1000), dsigp(1000)

      save            dsiga, dsigp

*-----------------------------------------------------------------------
*        mesh for angular distribution
*-----------------------------------------------------------------------

              ianm  = mstz(24)
              andif = 180.0 / ianm

*-----------------------------------------------------------------------
*        ifirst = 1 ; normal,  ne= 1 ; repeatition of the same system
*-----------------------------------------------------------------------

         if( ifirst .ne. 1 ) goto 200

*-----------------------------------------------------------------------

               ia = nint(ata)
               iz = nint(atz)

            call sigrc(incp,ex,ia,iz,sigt,sign,sige)

               sek = 0.0

            do ian = 1, ianm

               aplow = float(ian-1) * andif
               aphig = float(ian-1) * andif + andif
               apoin = ( aplow + aphig ) /2.0
               aprad = apoin * pi / 180.0
               apsin = sin( aprad )
               adrad = ( aphig - aplow ) * pi / 180.0

               sigs  = sign
               icct  = 2
               icm   = 0

               call dsdarc(icct,icm,incp,
     &                     ex,ata,apoin,sigs,dsigs,escat,etarg,angle)

               dsigp(ian) = aplow
               dsiga(ian) = dsigs * 2. * pi * apsin * adrad

               sek = sek + dsiga(ian)

            end do

               fnorm = sek

               sek = 0.0

            do ian = 1, ianm

               sek = sek + dsiga(ian) / fnorm

               dsiga(ian) = sek

            end do

*-----------------------------------------------------------------------
*        random number 0 < ramx < 1
*-----------------------------------------------------------------------

  200 continue

               ram1 = unirn(dummy)
               ram2 = unirn(dummy)

*-----------------------------------------------------------------------

            do ian = 1, ianm

               if( dsiga(ian) .gt. ram1 ) goto 100

            end do

  100       continue

               angcm = ( dsigp(ian) + ram2 * andif ) * pi / 180.0

               csthcm = cos( angcm )

*-----------------------------------------------------------------------
*           high energy > 1000 GeV  => forward only
*-----------------------------------------------------------------------

                  ee1 =    1000.0
                  ee2 = 1000000.0

                  ec1 = log(ee1)
                  ec2 = log(ee2)

               if( ex .gt. ee1 .and. ex .le. ee2 ) then

                  aa  = 1.0 / ( ec1 - ec2 )
                  bb  = - ec2 * aa

                  csthcm = 1.0 - ( 1.0 - csthcm )
     &                   * ( aa * log(ex) + bb )**8

               else if( ex .gt. ee2 ) then

                  csthcm = 1.0

               end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine dsdarc(icc,icm,incp,
     &                  pnrg,atarg,angle,sigs,dsigs,escat,etarg,angcm)
*                                                                      *
*        compute elastic scattering cross-section using                *
*        Fraunhoffer diffraction (b) for scattering formulae,          *
*        see S. Pearlstein, Nuc. Sci. Eng., 49, p.162-171(1972)        *
*                                                                      *
*        Hutcheon PRL 47,315(81) data 200-400 mev p+pb data            *
*                                                                      *
*     input:                                                           *
*        icc      : 0 -> Pearlstein                                    *
*                 : 1 -> Pearlstein + Hutcheon                         *
*                 : 2 -> Niita                                         *
*        icm      : 0 -> cm, 1 -> lab  input and output                *
*        incp     : 1->proton, 2->neutron                              *
*        pnrg     : kinetic energy of incident nucleon in lab (MeV)    *
*        atarg    : target mass number                                 *
*        angle    : scattering angle (deg) cm or lab                   *
*        sigs     : non-elastic scattering cross-section (b)           *
*                                                                      *
*     output:                                                          *
*        dsigs    : differential elastic cross-section (b/sr)          *
*        escat    : lab kinetic energy of scattered nucleon (MeV)      *
*        etarg    : =pnrg-escat, lab k.e. of scattered target (MeV)    *
*        angcm    : cm scattering angle (deg)                          *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------
*     constants
*-----------------------------------------------------------------------

            rz = 0.14
            pi = 3.1415927

*-----------------------------------------------------------------------
*     transformation lab to cm
*-----------------------------------------------------------------------

            call syscrc(icm,incp,pnrg,atarg,angle,vk,angzlm,escat,conv)

            angcm = angzlm
            etarg = pnrg - escat

            a  = atarg
            ia = nint( a )

*-----------------------------------------------------------------------
*        original elastic cross section
*-----------------------------------------------------------------------

         if( icc .ne. 2 ) then

            r = ( rz * a**0.3333
     &          + 0.122 * ( a + 1.0 ) / a ) * 1.d-12

            del = 0.05 * r

            area = pi * r * r * 1.d+24
            vkr2 = 2.0 * vk * r

            call bessel(0,vkr2,besz)
            call bessel(1,vkr2,bes1)

            sig = area * ( 1.0 - besz**2 - bes1**2 )

*-----------------------------------------------------------------------
*        new parametrization by K. Niita
*-----------------------------------------------------------------------

         else

            r = ( sqrt( sigs * 100. ) / pi + 0.0 ) * 1.d-13

            r = r * ( -43.78 / ( a + 136.6 ) + 2.23 )

            r = r * ( 1.0
     &            + ( ( 0.8 + 0.2
     &            *  exp( - min( 50.d0,( a / 100.0 )**4 )) ) - 1.0 )
     &            *  exp( - min( 50.d0,( pnrg / 40.0 )**2 )) )

            del = 0.45 * 1.d-13

            rfac = 0.66 * exp( - min( 50.d0,( 541.4 / pnrg )**2) )
            dfac = 1.80 * exp( - min( 50.d0,( 61.43 / a )**2) )

            r = r * ( 1.0 + rfac )

            del = del * ( 1.0 + dfac * rfac * 1.35 )

         end if

*-----------------------------------------------------------------------
*     differential cross section
*-----------------------------------------------------------------------

            y0 = ( vk * r * r )**2 * 1.d+24

            u = cos( angzlm / 180.0 * pi )
            s = sqrt( ( 1.0 - u ) / 2.0 )

            if( u .eq. 1.0 ) s = 0.0

            xav = 2.0 * vk * r * s

         if( xav .le. 1.d-2 ) then

            xj1x = 0.5
            y1   = xj1x * xj1x
            y2   = 1.0

         else

            xpl = 2.0 * vk * ( r + del ) * s
            xmi = 2.0 * vk * ( r - del ) * s

            call bessel(1,xpl,xj1pl)
            call bessel(1,xmi,xj1mi)

            if( icc .ne. 2 ) then

               fp1 = 1.0
               fp2 = 1.0

            else

               fp1 = 1.4
               fp2 = 2.0 - fp1

            end if

            y1 = 0.5 * (
     &           fp1 * ( xj1pl / xpl )**2 +
     &           fp2 * ( xj1mi / xmi )**2 )

*-----------------------------------------------------------------------
*           following fits Hutcheon PRL 47,315(81)
*           data 200-400 mev p+pb data
*           icc = 0 -> Pearlstein
*                 1 -> Pearlstein + Hutcheon
*                 2 -> Niita
*-----------------------------------------------------------------------

            if( icc .eq. 0 ) then

               fac0 = 1.0

               y2 = fac0

            else if( icc .eq. 1 ) then

               cf1  = 0.2 * ( 208. / a )**0.3333
               fac1 = exp( - min( 50.d0,
     &                            cf1 * xav + cf1 / 40. * xav**2 ) )

               y2 = fac1

            else if( icc .eq. 2 ) then

               cf2 = 0.1 * ( 100. / a )**0.3333
               facc = exp( - min( 50.d0,
     &                            cf2 * xav + cf2 / 40. * xav**2 ) )

               emas = ( -14512.0 / ( a + 103.5 ) + 146.7 )
               fac2 = 1.0 + ( facc - 1.0 )
     &              * ( 1.0 - exp( - min( 50.d0,
     &                                  ( pnrg / emas )**4 ) ) )

               cf3 = 2.5 * ( 80. / a )**0.3333
               angc = 160.0 * 50.0 / ( pnrg + 50.0 ) + 20.0

               angd = - 15.0 / 190.0 * a + 20.79
               ange =   15.0 / 190.0 * a + 49.2

              angc = ( 180.0 - angd ) * ange / ( pnrg + ange ) + angd

               faca = exp( - min( 50.d0, cf3 * ( angzlm / angc )**2 ) )

               fac3 = 1.0 + ( faca - 1.0 )
     &              * ( 1.0 - exp( - min( 50.d0,
     &                                  ( pnrg / 60.0 )**2 ) ) )

               y2 = fac2 * fac3

            end if


         end if

*-----------------------------------------------------------------------

            dsigs = conv * y0 * y1 * y2

            if( abs(dsigs) .le. 1.d-12 ) dsigs = 0.0


*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine syscrc(icm,incp,e,a,angle,vk,angzlm,escat,conv)
*                                                                      *
*        to convert from lab parameters                                *
*        to zero linear momentum (zlm) frame                           *
*        for relativistic formulae,                                    *
*        see J.L. Fowler, J.E. Brolley, Jr.,                           *
*        Rev. Mod. Phys. 28,103-134(1956)                              *
*                                                                      *
*     input:                                                           *
*        icm      : 0 -> cm, 1 -> lab  input and output                *
*        incp     : 1->proton, 2->neutron                              *
*        e        : incident energy in lab (MeV)                       *
*        angle    : scattering angle (deg) cm or lab                   *
*                                                                      *
*     output:                                                          *
*        vk       : wave number (1/cm) in lab                          *
*        angzlm   : scattering angle in cm system (deg)                *
*        escat    : elastic scattering energy in lab (MeV)             *
*        conv     : factor to convert cross-section from c.m.          *
*                   to lab system                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      pi = 3.1415927

*-----------------------------------------------------------------------

            if( incp .eq. 1 ) then

               rmass = 938.3

            else if( incp .eq. 2 ) then

               rmass = 939.58

            end if

            ang = angle * pi / 180.
            rat = e / rmass
            gz  = 1. + rat
            s   = 1. + 2. * a * gz + a**2

*-----------------------------------------------------------------------
*     if nonrelativistic, wave number in lab, vk
*-----------------------------------------------------------------------

            vk = 2.197e+12 * dsqrt(e)

*-----------------------------------------------------------------------
*     corrected for relativistic and zlm motion
*-----------------------------------------------------------------------

            add = 1.
            vk  = vk * a / ( a + add )
     &          * dsqrt( ( 1. + rat / 2. ) / ( 1. + rat )**2 )

*-----------------------------------------------------------------------

         if( icm .eq. 0 ) then

            angzlm = angle

            cosa = cos( angzlm * pi / 180. )
            conv = 1.0

*-----------------------------------------------------------------------

         else

               gcg=(a+gz)/dsqrt(s)
               g1=(1.+a*gz)/dsqrt(s)
               a1=(1.+a*gz)/(a+gz)/a
               c3=gcg**2
               c1=-a1*c3
               c2=(a1**2-1.)*c3
               if(dabs(ang-pi/2.).ge.1.e-2) go to 10
               cosa=c1/c3

               go to 20

   10          continue
               tanx=dsin(ang)/dcos(ang)
               sign=1.
               if(ang.gt.(pi/2.)) sign=-1.
               x2=tanx**2
               cosa=(c1*x2+sign*dsqrt(1.-c2*x2))/(1.+c3*x2)
               if(angle.ge.179.) cosa=-1.
               if(angle.le.0.1) cosa=1.

   20          continue

               angzlm=acos(cosa)

               angzlm=angzlm*180./pi

               if(dabs(ang-pi/2.).ge.0.02) go to 30

               ang1=ang-0.1
               ang2=ang+0.1
               cosx1=dcos(ang1)
               cosx2=dcos(ang2)
               x31=cosx1**3
               x32=cosx2**3
               x21=(dsin(ang1)/dcos(ang1))**2
               x22=(dsin(ang2)/dcos(ang2))**2
               conv=1./2.*((-2.*c1/x31+sign*c2/x31
     &              /dsqrt(1.-c2*x21))/(1.+c3*x21)
     &              +2.*c3/x31*(c1*x21
     &              +sign*dsqrt(1.-c2*x21))/(1.+c3*x21)**2
     &              +(-2.*c1/x32-sign*c2/x32/dsqrt(1.-c2*x22))
     &              /(1.+c3*x22)+2.*c3/x32
     &              *(c1*x22+sign*dsqrt(1.-c2*x22))/(1.+c3*x22)**2)

               go to 40

   30          continue

               cosx=dcos(ang)
               x3=cosx**3

               conv=(-2.*c1/x3+sign*c2/x3/dsqrt(1.-c2*x2))/(1.+c3*x2)+
     &           2.*c3/x3*(c1*x2+sign*dsqrt(1.-c2*x2))/(1.+c3*x2)**2

   40          continue

               if(dabs(conv).le.1.e-6) conv=0.
cKN
               if(conv.le.0.0) conv=1.
cKN
               if(dabs(a-1.).le.0.01.and.ang.ge.(pi/2.-1.e-4)) conv=0.

         end if

*-----------------------------------------------------------------------

         escat = rmass * ( gz - 1. )
     &         / s * ( 1. + a * ( gz - 1. + ( gz + 1. ) * cosa )
     &                    + a**2 )

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine bessel(n,x,xjn)
c ----------------------------------------------------------------------
c     bessel integer functions, see ams-55, handbook of
c     mathematical functions, chapter 9.
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c
      dimension a(7),b(7),c(7),d(7),e(7),g(7)
      data a/1.0,-2.2499997, 1.2656208,-.3163866, .0444479,
     *  -.0039444, .0002100/
      data b/ .79788456,-.00000077,-.00552740,-.00009512, .00137237,
     *  -.00072805, .00014476/
      data c/-.78539816,-.04166397,-.00003954, .00262573,-.00054125,
     *  -.00029333, .00013558/
      data d/ .5,-.56249985, .21093573,-.03954289, .00443319,
     * -.00031761, .00001109/
      data e/ .79788456, .0000156, .01659667, .00017105,-.00249511,
     *  .00113653,-.00020033/
      data g/-2.35619449, .12499612, .00005650,-.00637879, .00074348,
     *  .00079824,-.00029166/
c ----------------------------------------------------------------------
c
      if(n-1) 25,210,25
   25 if(x-3.0) 50,50,100
   50 if(x-0.) 55,55,60
   55 xjzr=1.
      go to 200
   60 xjzr=a(1)
      do 70 m=2,7
      k=2*(m-1)
   70 xjzr=xjzr+a(m)*(x/3.)**k
      go to 200
  100 fzer=b(1)
      do 110 m=2,7
  110 fzer=fzer+b(m)*(3./x)**(m-1)
      thetz=x
      do 120 m=1,7
  120 thetz=thetz+c(m)*(3./x)**(m-1)
      xjzr=fzer*dcos(thetz)/dsqrt(x)
  200 if(n-1) 201,210,210
  201 xjn=xjzr
      return
  210 if(x-3.) 250,250,300
  250 if(x-1.e-6) 255,255,260
  255 xjne=0.
      go to 400
  260 xjne=d(1)
      do 270 m=2,7
      k=2*(m-1)
  270 xjne=xjne+d(m)*(x/3.)**k
      xjne=x*xjne
      go to 400
  300 fzne=e(1)
      do 310 m=2,7
  310 fzne=fzne+e(m)*(3./x)**(m-1)
      thetn=x
      do 320 m=1,7
  320 thetn=thetn+g(m)*(3./x)**(m-1)
      xjne=fzne*dcos(thetn)/dsqrt(x)
  400 if(n-1) 401,401,408
  408 if(x-1.e-6) 409,409,410
  401 xjn=xjne
      return
  409 xjn=0.
      return
  410 do 420 i=2,n
      xn=i-1
      xjn=2.0*xn*xjne/x-xjzr
      xjzr=xjne
      xjne=xjn
  420 continue
      return
      end


************************************************************************
*                                                                      *
      subroutine sighi(ap,zp,ep,at,zt,signe,sigel,bmax)
*                                                                      *
*        calculates reaction cross-sections of nucleus-nucleus         *
*        choose the models                                             *
*        last modified by K.Niita on 2004/03/04                        *
*                                                                      *
*     common:                                                          *
*       icrhi   : choice of cross section formula                      *
*               : 0; Shen's formula                                    *
*               : 1; NASA's formula                                    *
*                                                                      *
*     input:                                                           *
*        ap     : mass number of projectile                            *
*        zp     : charge number of projectile                          *
*        ep     : incident nucleus total energy (MeV)                  *
*        at     : mass number of target                                *
*        zt     : charge number of target                              *
*                                                                      *
*     output:                                                          *
*       signe   : reaction cross-section (b)                           *
*       sigel   : elastic cross-section (b)                            *
*       bmax    : correspond impact parameter (fm)                     *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      common /crshi/  bplus, icrhi, ijudg, imadj, iqmax

      dimension att(10), ztt(10)

*-----------------------------------------------------------------------


         if( icrhi .eq. 0 ) then

            call shen(ap,zp,ep,at,zt,signe,sigel,bmax)

         else

            call nasa(ap,zp,ep,at,zt,signe,sigel,bmax)

         end if

*-----------------------------------------------------------------------
*        for debug to comment out the next return
*-----------------------------------------------------------------------

         return

*-----------------------------------------------------------------------

         ap = 56
         zp = 26

         ap = 20
         zp = 10

         ap = 20
         zp = 10

         at = 208
         zt = 82

         att(1) =  1
         ztt(1) =  1
         att(2) = 12
         ztt(2) =  6
         att(3) = 27
         ztt(3) = 13
         att(4) = 64
         ztt(4) = 29
         att(5) = 119
         ztt(5) =  50
         att(6) = 181
         ztt(6) =  73
         att(7) = 208
         ztt(7) =  82


         write(6,'(''wt:'')')
         write(6,'('' ap ='',i4,'',  zp ='',i3)') nint(ap), nint(zp)
c        write(6,'('' at ='',i4,'',  zt ='',i3)') nint(at), nint(zt)
         write(6,'(''e:'')')
c        write(6,'(''x: Eneergy/nucleon (MeV)'')')
         write(6,'(''x: mass number'')')
c        write(6,'(''y: Cross Section (mb)'')')
         write(6,'(''y: bmax (fm)'')')
c        write(6,'(''p: xlog afac(0.8)'')')
c        write(6,'(''p: ylog xlin afac(0.8)'')')
         write(6,'(''p: ylin xlin afac(0.8)'')')
c        write(6,'(''h:   x   y(SHEN),l0r   y(NASA),l0b'')')
c        write(6,'(''h:   x   y(NASA),l0r   y(bmax),l0b'')')
         write(6,'(''h: x y(SHEN),l0r y(NASA),l0b y(Sato),l0g'',
     &             '' y(bmax),l0'')')

         emin = 100
         emax = 400.0
         nemd = 16
         edel = log( emax / emin ) / ( nemd - 1 )
         edel = ( emax - emin ) / ( nemd - 1 )

         nemd = 1
         nemd = 16

         do i = 1, nemd

c          ep = exp( ( i - 1 ) * edel ) * emin * ap
           ep = ( ( i - 1 ) * edel + emin ) * ap

c          at = att(i)
c          zt = ztt(i)

           at = 12
           zt = 6

c          ep = 1053 * ap
c          ep = 585 * ap
c          ep = 400 * ap

            call shen(ap,zp,ep,at,zt,signe,sigel,bmax)

               sigsh = signe
               bmax1 = bmax

            call nasa(ap,zp,ep,at,zt,signe,sigel,bmax)

               signa = signe
               bmax2 = bmax

c           bmax0 = 1.189 * ( ap**(1./3.) + at**(1./3.) ) - 0.9612
            bmax3 = 1.2 * ( ap**(1./3.) + at**(1./3.) ) - 0.5

c            write(6,'(3e13.5)') ep/ap, sigsh*1000, signa*1000
c            write(6,'(3e13.5)') at, sigsh*1000, signa*1000
c            write(6,'(5e13.5)') at, bmax1, bmax2, bmax0, bmax3
             write(6,'(3e13.5)') ep/ap, sigsh*1000, signa*1000

         end do

        stop 666

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine shen(ap0,zp0,ep0,at0,zt0,signe,sigel,bmax)
*                                                                      *
*        calculates reaction cross-sections of nucleus-nucleus         *
*        this parametrization is taken from                            *
*        Nucl. Phys. A491 (1989) 130 by SHEN Wen-qing, et.al.          *
*                                                                      *
*        coded by H. Iwase, on 2002/01/10                              *
*                                                                      *
*     input:                                                           *
*        ap     : mass number of projectile                            *
*        zp     : charge number of projectile                          *
*        ep     : incident nucleus total energy (MeV)                  *
*        at     : mass number of target                                *
*        zt     : charge number of target                              *
*                                                                      *
*     output:                                                          *
*       signe   : reaction cross-section (b)                           *
*       sigel   : elastic cross-section (b)                            *
*       bmax    : correspond impact parameter (fm)                     *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      parameter ( pi = 3.1415926535897932d0 )

c----------------------------------------------------------------------
c     this routine requires that at > ap
c----------------------------------------------------------------------

      if( ap0 .gt. at0 ) then

         ep = ep0 / ap0 * at0
         ap = at0
         zp = zt0
         at = ap0
         zt = zp0

      else

         ep = ep0
         ap = ap0
         zp = zp0
         at = at0
         zt = zt0

      end if

*-----------------------------------------------------------------------

      signe = 0.0
      sigel = 0.0
      bmax  = 0.0

      ecm = ep * at / ( ap + at )

      call keibeta(ap,zp,at,zt,bbb)

      if( ecm .le. bbb ) return

      call keirrr(ap,zp,at,zt,ep,rrr)

      signe = 10.0 * pi * rrr**2 * ( 1.0 - bbb / ecm )
      signe = signe / 1000.

      bmax = rrr * sqrt( 1.0 - bbb / ecm )

      return
      end


************************************************************************
*                                                                      *
        subroutine keibeta(ap,zp,at,zt,bbb)
*                                                                      *
*        ap,zp : projectile mass and charge                            *
*        at,zt : target mass and charge                                *
*        bbb   : barrier                                               *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      b   = 1.0
      rrp = (1.12*(ap**(1.0/3.0)))-(0.94*(ap**(-1.0/3.0)))
      rrt = (1.12*(at**(1.0/3.0)))-(0.94*(at**(-1.0/3.0)))
      r   = rrp+rrt+3.2

      bbb = ((1.44*zt*zp)/r)-b*((rrt*rrp)/(rrt+rrp))

      return
      end


************************************************************************
*                                                                      *
      subroutine keirrr(ap,zp,at,zt,ep0,rrr)
*                                                                      *
*        ap,zp : projectile mass and charge                            *
*        at,zt : target mass and charge                                *
*        ep0   : projectile total energy                               *
*        rrr   : reaction radius                                       *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      common /cedata/ ee(55),ce(55)

*-----------------------------------------------------------------------

      r0=1.1
      be=0.176

      ep = ep0 / ap

      do i=1,55
         if(ee(i).ge.ep) then
            ced=ce(i)
            go to 10
         end if
      end do
            ced=ce(55)
 10   continue

      za=1.85*(at**(1.0/3.0))*(ap**(1.0/3.0))
      zb=(at**(1.0/3.0)+ap**(1.0/3.0))

      dd=za/zb
      dda=at**(1.0/3.0)+ap**(1.0/3.0)+dd-ced
cKN
      ddb=((at-2.0*zt)*zp)/(at*ap)
cKN
      ecm=ep0*at/(ap+at)

      zc=(ecm**(-1.0/3.0))*(ap**(1.0/3.0))*(at**(1.0/3.0))
      zd=(at**(1.0/3.0)+ap**(1.0/3.0))
      ddc=zc/zd
      rrr=dda*r0+ddb+be*ddc

      return
      end


************************************************************************
*                                                                      *
      block data readce
*                                                                      *
*       transparency coefficients for Shen formula                     *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      common /cedata/ ee(55),ce(55)
      data ee/1.71343E+00, 2.33903E+00, 2.92134E+00, 3.77688E+00,
     +        4.35853E+00, 5.28443E+00, 6.31251E+00, 7.42979E+00,
     +        8.78789E+00, 1.01410E+01, 1.11927E+01, 1.27256E+01,
     +        1.46115E+01, 1.59679E+01, 1.72794E+01, 1.88830E+01,
     +        2.02318E+01, 2.13589E+01, 2.28843E+01, 2.36854E+01,
     +        2.54983E+01, 2.70476E+01, 2.92652E+01, 3.14451E+01,
     +        3.49550E+01, 3.77599E+01, 4.07804E+01, 4.40631E+01,
     +        5.02993E+01, 5.53659E+01, 6.06408E+01, 6.51767E+01,
     +        7.17530E+01, 7.90173E+01, 8.95879E+01, 1.02118E+02,
     +        1.16990E+02, 1.34726E+02, 1.52204E+02, 1.78908E+02,
     +        2.07299E+02, 2.31066E+02, 2.66530E+02, 2.99893E+02,
     +        3.29279E+02, 3.56168E+02, 3.92895E+02, 4.29124E+02,
     +        4.27030E+02, 4.75548E+02, 5.37448E+02, 6.22204E+02,
     +        7.13148E+02, 8.13585E+02, 9.23209E+02/
      data ce/2.77909E-02, 4.44954E-02, 5.85097E-02, 8.40984E-02,
     +        1.01228E-01, 1.24064E-01, 1.58599E-01, 1.87341E-01,
     +        2.24817E-01, 2.50693E-01, 2.94175E-01, 3.28836E-01,
     +        3.78048E-01, 4.30302E-01, 4.68003E-01, 5.26087E-01,
     +        5.81306E-01, 6.24900E-01, 6.83035E-01, 7.38341E-01,
     +        8.48942E-01, 9.42087E-01, 1.02060E+00, 1.15135E+00,
     +        1.26087E+00, 1.33189E+00, 1.41181E+00, 1.47392E+00,
     +        1.52465E+00, 1.59636E+00, 1.67101E+00, 1.71886E+00,
     +        1.78464E+00, 1.83855E+00, 1.90151E+00, 1.94670E+00,
     +        1.98597E+00, 2.01341E+00, 2.02590E+00, 2.01784E+00,
     +        1.99782E+00, 1.96277E+00, 1.92790E+00, 1.90774E+00,
     +        1.88745E+00, 1.87895E+00, 1.86759E+00, 1.86211E+00,
     +        1.86208E+00, 1.86264E+00, 1.86030E+00, 1.87293E+00,
     +        1.89440E+00, 1.90695E+00, 1.93727E+00/
      end


************************************************************************
*                                                                      *
      subroutine nasa(ap0,zp0,ep0,at0,zt0,signe,sigel,bmax)
*                                                                      *
*        calculates reaction cross-sections of nucleus-nucleus         *
*                                                                      *
*        this parametrization is taken from                            *
*        NIM B 117 (1996) 347 by R.K. Tripathi, et.al.,                *
*        NIM B 129 (1997)  11 by R.K. Tripathi, et.al.,,               *
*        and                                                           *
*        NIM B 155 (1999) 349 by R.K. Tripathi, et.al.                 *
*                                                                      *
*        coded by H. Iwase, on 2004/02/06                              *
*                                                                      *
*     input:                                                           *
*        ap     : mass number of projectile                            *
*        zp     : charge number of projectile                          *
*        ep    : incident nucleus total energy (MeV)                  *
*        at     : mass number of target                                *
*        zt     : charge number of target                              *
*                                                                      *
*     output:                                                          *
*       signe   : reaction cross-section (b)                           *
*       sigel   : elastic cross-section (b)                            *
*       bmax    : correspond impact parameter (fm)                     *
*                                                                      *
************************************************************************

      implicit none

      real*8 signe, pi, r0, deltaE, B, Ecm, R, rp, rt, rmsp, rmst,
     $     radius, S, CE, rmsc, rca, E, rap, rat, rac, D, Rc, fact,
     $     G, T1, SL, X1, Xm, sigel, bmax, ep

      real*8 Ap, At, Zp, Zt, Ac
      real*8 Ap0, At0, Zp0, Zt0, ep0
      integer i

      parameter (pi    = 3.14159265358979d0)
      parameter (r0    = 1.1)   ! (fm)
      parameter (Ac    = 12.)

c----------------------------------------------------------------------
c     this routine requires that at > ap
c----------------------------------------------------------------------

      if( ap0 .gt. at0 ) then

         ep = ep0 / ap0 * at0
         ap = at0
         zp = zt0
         at = ap0
         zt = zp0

      else

         ep = ep0
         ap = ap0
         zp = zp0
         at = at0
         zt = zt0

      end if

c----------------------------------------------------------------------

      signe = 0.0
      sigel = 0.0
      bmax  = 0.0

      E    = Ep / Ap
      Ecm  = Ep * At / ( Ap + At )

      rmsp  = radius(Ap,Zp)
      rmst  = radius(At,Zt)

      fact = sqrt(5./3.) ! = 1.29

      rp  = fact * rmsp
      rt  = fact * rmst
      rca = fact * 2.471

      rap = 3. * Ap  / ( 4. * pi * rp**3. )
      rat = 3. * At  / ( 4. * pi * rt**3. )
      rac = 3. * 12. / ( 4. * pi * rca**3. )

c---------------------------------------------------- D selection start

      if(ap.eq.1.and.zp.eq.0)then
c     neutron + X systems
         T1 = 18.
         D  = 1.85 + ( 0.16 ) / ( 1. + exp( (500.-E) / 200. ))

      elseif(ap.eq.1.and.zp.eq.1.and.at.eq.4.and.zt.eq.2)then
c     proton + alpha
         T1 = 40.
         D  = 2.05

      elseif(ap.eq.1.and.zp.eq.1.and.at.eq.3.and.zt.eq.2)then
c     proton + 3He
         T1 = 58.
         D  = 1.70

      elseif(ap.eq.1.and.zp.eq.1.and.at.eq.6.and.zt.eq.3)then
c     proton + 6Li
         T1 = 40.
         D  = 2.05

      elseif(ap.eq.1.and.zp.eq.1.and.at.eq.7.and.zt.eq.3)then
c     proton + 7Li
         T1 = 37.
         D  = 2.15

      elseif((ap.eq.1.and.zp.eq.1).and.At.le.7)then
c     proton + light nucleus
         T1 = 23.
         D  = 1.85 + ( 0.16 ) / ( 1. + exp( (500.-E) / 200. ))

      elseif(ap.eq.1.and.zp.eq.1)then
c     proton + others
         T1 = 40
         D  = 2.05

      elseif(ap.eq.2.and.zp.eq.1.and.at.eq.4.and.zt.eq.2)then
c     deuteron + alpha
         T1 = 23
         D  = 1.65 + ( 0.22 ) / ( 1. + exp( (500.-E) / 200.  ))

      elseif(ap.eq.2.and.zp.eq.1)then ! deuteron + X systems
         T1 = 23.
         D  = 1.65 + ( 0.1 ) / ( 1. + exp( (500.-E) / 200.  ))

      elseif(ap.eq.3.and.zp.eq.2 .or. at.eq.3.and.zt.eq.2)then
         T1 = 40.
         D  = 1.55


      elseif((ap.eq.4.and.zp.eq.2).and.(At.eq.4.and.Zt.eq.2))then
c     alhpa + alpha
         T1 = 40.
         G  = 300.
         D = 2.77 - 8.0E-3 * At + 1.8E-5 * At * At
     $        -0.8 / ( 1. + exp( (250.-E) /G ) )

      elseif((ap.eq.4.and.zp.eq.2).and.(zt.eq.4))then
c     alpha + Be
         T1 = 25.
         G  = 300.
         D = 2.77 - 8.0E-3 * At + 1.8E-5 * At * At
     $        -0.8 / ( 1. + exp( (250.-E) /G ) )

      elseif((ap.eq.4.and.zp.eq.2).and.(zt.eq.7))then
c     alpha + N
         T1 = 40.
         G  = 500.
         D = 2.77 - 8.0E-3 * At + 1.8E-5 * At * At
     $        -0.8 / ( 1. + exp( (250.-E) /G ) )

      elseif((ap.eq.4.and.zp.eq.2).and.(zt.eq.13))then
c     alpha + Al
         T1 = 25.
         G  = 300.
         D = 2.77 - 8.0E-3 * At + 1.8E-5 * At * At
     $        -0.8 / ( 1. + exp( (250.-E) /G ) )

      elseif((ap.eq.4.and.zp.eq.2).and.(zt.eq.26))then
c     alpha + Fe
         T1 = 40.
         G  = 300.
         D = 2.77 - 8.0E-3 * At + 1.8E-5 * At * At
     $        -0.8 / ( 1. + exp( (250.-E) /G ) )

      elseif(ap.eq.4.and.zp.eq.2 .or. at.eq.4.and.zt.eq.2)then
c     alpha + other systems
         T1 = 40.
         G  = 75.
         D = 2.77 - 8.0E-3 * At + 1.8E-5 * At * At
     $        -0.8 / ( 1. + exp( (250.-E) /G ) )

      elseif(zp.eq.3 .or. zt.eq.3)then
c     Li + X
         T1 = 40.
         D = 1.75 * ( rap + rat ) / ( rac + rac ) /3

      else
c     others
         T1 = 40.
         D = 1.75 * ( rap + rat ) / ( rac + rac )

      endif

c---------------------------------------------------- D selection end


c----------------------------------------------------------------------

      CE = D * ( 1 - exp( -E / T1) )  - 0.292 * exp( -E / 792 ) *
     $     cos( 0.229 * E**(0.453) )

      S = ( Ap**(1./3.) * At**(1./3.) ) / ( Ap**(1./3.) + At**(1./3.) )

      deltaE = 1.85 * S + ( 0.16 * S / Ecm **(1./3.) ) - CE + 0.91 * (
     $     At - 2 * Zt ) * Zp / ( At * Ap )


      R = rp + rt + 1.2 * ( Ap**(1./3.) + At**(1./3.) ) / ( Ecm**(1./3.)
     $     )


      B = 1.44 * Zp * Zt / R

c----------------------------------------------------------------------


c--------------------------------------------------- Rc selection start

      if((ap.eq.1.and.zp.eq.1).and.(at.eq.2.and.zt.eq.1))then
         Rc = 13.5              ! p + d

      elseif((ap.eq.1.and.zp.eq.1).and.(at.eq.3.and.zt.eq.2))then
         Rc = 21                ! p + 3He

      elseif((ap.eq.1.and.zp.eq.1).and.(at.eq.4.and.zt.eq.2))then
         Rc = 27                ! p + 4He

      elseif((ap.eq.1.and.zp.eq.1).and.(zt.eq.3))then
         Rc = 2.2               ! p + Li

      elseif((ap.eq.1.and.zp.eq.1).and.(zt.eq.6))then
         Rc = 3.5               ! p + C

      elseif((ap.eq.2.and.zp.eq.1).and.(at.eq.2.and.zt.eq.1))then
         Rc = 13.5              ! d + d

      elseif((ap.eq.2.and.zp.eq.1).and.(at.eq.4.and.zt.eq.2))then
         Rc = 13.5              ! d + 4He

      elseif((ap.eq.2.and.zp.eq.1).and.(zt.eq.6))then
         Rc = 6.0               ! d + C

      elseif((ap.eq.4.and.zp.eq.2).and.(zt.eq.73))then
         Rc = 0.6               ! alpha + Ta

      elseif((ap.eq.4.and.zp.eq.2).and.(zt.eq.79))then
         Rc = 0.6               ! alpha + Au

      else
         Rc = 1.0
      endif

c------------------------------------------------- Rc selection end


c---------------------------------------------- Xm calculation start

      if (ap.eq.1.and.zp.eq.0)then ! neutron + X systems

         SL = 1.2 + 1.6 * ( 1.0 - exp( -E/15. ) )


         if((at.eq.4.and.zt.eq.2) )then
c        n + alpha
            X1 = 5.2

         else
c        others
            X1 = 2.83 - 3.1E-2 * At + 1.7E-4 * At * At

         endif

         Xm = 1. - X1 * exp( -E / ( X1 * SL ) )

      else
         Xm = 1.0
      endif

c---------------------------------------------- Xm calculation end


c------------------------------------------- signe calculation start

      signe = pi * r0 * r0 * ( Ap**(1./3.) + At**(1./3.) + deltaE )**(2.
     $     ) * ( 1 - Rc * B / Ecm ) * Xm

      if(signe.le.0) then
         signe=0
         return
      endif

      signe = signe * 10        ! convert the unit from [fm**2] to [mb]
      signe = signe / 1000      ! convert the unit from [mb] to [b]

      bmax =  r0 * ( Ap**(1./3.) + At**(1./3.) + deltaE )
     $        * sqrt(1. -  Rc * B / Ecm)

c------------------------------------------- signe calculation end

      end


c----------------------------------------------------------------------
      function radius(a,z)
c----------------------------------------------------------------------
c Purpose   : to obtain the "r_rms,i" values
c References: Atomic Data adn Nuclear Data Tables 36, 495-536 (1987)
c             and
c             NIM-B 152(1999)425-431
c             by H.Iwase Fri Feb  6 2004
c----------------------------------------------------------------------

      implicit none

      real*8 rms(300,2)
      real*8 a, z, radius
      integer i
      integer irnm

c----------------------------------------------------------------------
c     these data are refered from
c     Atomic Data adn Nuclear Data Tables 36, 495-536 (1987)
c     data from page 503 to 510 are used
c----------------------------------------------------------------------

       data irnm / 135 /

       data rms(  1,1), rms(  1,2) /     1, 0.3407/
       data rms(  2,1), rms(  2,2) /  1001, 0.8507/
       data rms(  3,1), rms(  3,2) /  1002, 2.1055/
       data rms(  4,1), rms(  4,2) /  1003, 1.7200/
       data rms(  5,1), rms(  5,2) /  2003, 1.8990/
       data rms(  6,1), rms(  6,2) /  2004, 1.6810/
       data rms(  7,1), rms(  7,2) /  3006, 2.5567/
       data rms(  8,1), rms(  8,2) /  3007, 2.4000/
       data rms(  9,1), rms(  9,2) /  4009, 2.5095/
       data rms( 10,1), rms( 10,2) /  5010, 2.4500/
       data rms( 11,1), rms( 11,2) /  5011, 2.3950/
       data rms( 12,1), rms( 12,2) /  6012, 2.4690/
       data rms( 13,1), rms( 13,2) /  6013, 2.4400/
       data rms( 14,1), rms( 14,2) /  6014, 2.5600/
       data rms( 15,1), rms( 15,2) /  7014, 2.5480/
       data rms( 16,1), rms( 16,2) /  7015, 2.6537/
       data rms( 17,1), rms( 17,2) /  8016, 2.7283/
       data rms( 18,1), rms( 18,2) /  8017, 2.6620/
       data rms( 19,1), rms( 19,2) /  8018, 2.7270/
       data rms( 20,1), rms( 20,2) /  9019, 2.9000/
       data rms( 21,1), rms( 21,2) / 10020, 3.0120/
       data rms( 22,1), rms( 22,2) / 10022, 2.9690/
       data rms( 23,1), rms( 23,2) / 11023, 2.9400/
       data rms( 24,1), rms( 24,2) / 12024, 3.0467/
       data rms( 25,1), rms( 25,2) / 12025, 3.0565/
       data rms( 26,1), rms( 26,2) / 12026, 3.0600/
       data rms( 27,1), rms( 27,2) / 13027, 3.0483/
       data rms( 28,1), rms( 28,2) / 14028, 3.1140/
       data rms( 29,1), rms( 29,2) / 14029, 3.1045/
       data rms( 30,1), rms( 30,2) / 14030, 3.1760/
       data rms( 31,1), rms( 31,2) / 15031, 3.1880/
       data rms( 32,1), rms( 32,2) / 16032, 3.2423/
       data rms( 33,1), rms( 33,2) / 16034, 3.2810/
       data rms( 34,1), rms( 34,2) / 17035, 3.3880/
       data rms( 35,1), rms( 35,2) / 16036, 3.2780/
       data rms( 36,1), rms( 36,2) / 18036, 3.3270/
       data rms( 37,1), rms( 37,2) / 17037, 3.3840/
       data rms( 38,1), rms( 38,2) / 19039, 3.4040/
       data rms( 39,1), rms( 39,2) / 18040, 3.4320/
       data rms( 40,1), rms( 40,2) / 20040, 3.4703/
       data rms( 41,1), rms( 41,2) / 20048, 3.4605/
       data rms( 42,1), rms( 42,2) / 22048, 3.6550/
       data rms( 43,1), rms( 43,2) / 22050, 3.5730/
       data rms( 44,1), rms( 44,2) / 24050, 3.6690/
       data rms( 45,1), rms( 45,2) / 23051, 3.5975/
       data rms( 46,1), rms( 46,2) / 24052, 3.6467/
       data rms( 47,1), rms( 47,2) / 24053, 3.7260/
       data rms( 48,1), rms( 48,2) / 24054, 3.7127/
       data rms( 49,1), rms( 49,2) / 26054, 3.6957/
       data rms( 50,1), rms( 50,2) / 25055, 3.6800/
       data rms( 51,1), rms( 51,2) / 26056, 3.7503/
       data rms( 52,1), rms( 52,2) / 26058, 3.7750/
       data rms( 53,1), rms( 53,2) / 28058, 3.7683/
       data rms( 54,1), rms( 54,2) / 27059, 3.8130/
       data rms( 55,1), rms( 55,2) / 28060, 3.7953/
       data rms( 56,1), rms( 56,2) / 28061, 3.8060/
       data rms( 57,1), rms( 57,2) / 28062, 3.8263/
       data rms( 58,1), rms( 58,2) / 29063, 3.9187/
       data rms( 59,1), rms( 59,2) / 28064, 3.8673/
       data rms( 60,1), rms( 60,2) / 30064, 3.9370/
       data rms( 61,1), rms( 61,2) / 29065, 3.9440/
       data rms( 62,1), rms( 62,2) / 30066, 3.9580/
       data rms( 63,1), rms( 63,2) / 30068, 3.9667/
       data rms( 64,1), rms( 64,2) / 30070, 4.0077/
       data rms( 65,1), rms( 65,2) / 32070, 4.0565/
       data rms( 66,1), rms( 66,2) / 32072, 4.0550/
       data rms( 67,1), rms( 67,2) / 32074, 4.0750/
       data rms( 68,1), rms( 68,2) / 32076, 4.0810/
       data rms( 69,1), rms( 69,2) / 38088, 4.2060/
       data rms( 70,1), rms( 70,2) / 39089, 4.2500/
       data rms( 71,1), rms( 71,2) / 40090, 4.2707/
       data rms( 72,1), rms( 72,2) / 40091, 4.3090/
       data rms( 73,1), rms( 73,2) / 40092, 4.2970/
       data rms( 74,1), rms( 74,2) / 42092, 4.3047/
       data rms( 75,1), rms( 75,2) / 41093, 4.3205/
       data rms( 76,1), rms( 76,2) / 40094, 4.3235/
       data rms( 77,1), rms( 77,2) / 42094, 4.3340/
       data rms( 78,1), rms( 78,2) / 40096, 4.3960/
       data rms( 79,1), rms( 79,2) / 42096, 4.3640/
       data rms( 80,1), rms( 80,2) / 42098, 4.3880/
       data rms( 81,1), rms( 81,2) / 42100, 4.4300/
       data rms( 82,1), rms( 82,2) / 46104, 4.4370/
       data rms( 83,1), rms( 83,2) / 46106, 4.4670/
       data rms( 84,1), rms( 84,2) / 46108, 4.5240/
       data rms( 85,1), rms( 85,2) / 46110, 4.5900/
       data rms( 86,1), rms( 86,2) / 48110, 4.5780/
       data rms( 87,1), rms( 87,2) / 48112, 4.6080/
       data rms( 88,1), rms( 88,2) / 50112, 4.6205/
       data rms( 89,1), rms( 89,2) / 48114, 4.6305/
       data rms( 90,1), rms( 90,2) / 50114, 4.6020/
       data rms( 91,1), rms( 91,2) / 49115, 4.6460/
       data rms( 92,1), rms( 92,2) / 48116, 4.6390/
       data rms( 93,1), rms( 93,2) / 50116, 4.6240/
       data rms( 94,1), rms( 94,2) / 50117, 4.6250/
       data rms( 95,1), rms( 95,2) / 50118, 4.6630/
       data rms( 96,1), rms( 96,2) / 50119, 4.6390/
       data rms( 97,1), rms( 97,2) / 50120, 4.6430/
       data rms( 98,1), rms( 98,2) / 50122, 4.6580/
       data rms( 99,1), rms( 99,2) / 51122, 4.6300/
       data rms(100,1), rms(100,2) / 50124, 4.6807/
       data rms(101,1), rms(101,2) / 56138, 4.8360/
       data rms(102,1), rms(102,2) / 57139, 4.8500/
       data rms(103,1), rms(103,2) / 60142, 4.9253/
       data rms(104,1), rms(104,2) / 60144, 4.9260/
       data rms(105,1), rms(105,2) / 62144, 4.9470/
       data rms(106,1), rms(106,2) / 60146, 4.9815/
       data rms(107,1), rms(107,2) / 60148, 5.0020/
       data rms(108,1), rms(108,2) / 62148, 4.9890/
       data rms(109,1), rms(109,2) / 60150, 5.0037/
       data rms(110,1), rms(110,2) / 62150, 5.0450/
       data rms(111,1), rms(111,2) / 62152, 5.0947/
       data rms(112,1), rms(112,2) / 62154, 5.1259/
       data rms(113,1), rms(113,2) / 64154, 5.1240/
       data rms(114,1), rms(114,2) / 64156, 5.0680/
       data rms(115,1), rms(115,2) / 64158, 5.1720/
       data rms(116,1), rms(116,2) / 67165, 5.2100/
       data rms(117,1), rms(117,2) / 68166, 5.2593/
       data rms(118,1), rms(118,2) / 70174, 5.4100/
       data rms(119,1), rms(119,2) / 71175, 5.3700/
       data rms(120,1), rms(120,2) / 70176, 5.3790/
       data rms(121,1), rms(121,2) / 73181, 5.4800/
       data rms(122,1), rms(122,2) / 74184, 5.4200/
       data rms(123,1), rms(123,2) / 74186, 5.4000/
       data rms(124,1), rms(124,2) / 76192, 5.4130/
       data rms(125,1), rms(125,2) / 78196, 5.3800/
       data rms(126,1), rms(126,2) / 79197, 5.3000/
       data rms(127,1), rms(127,2) / 81203, 5.4630/
       data rms(128,1), rms(128,2) / 82204, 5.4790/
       data rms(129,1), rms(129,2) / 81205, 5.4745/
       data rms(130,1), rms(130,2) / 82206, 5.4963/
       data rms(131,1), rms(131,2) / 82207, 5.5050/
       data rms(132,1), rms(132,2) / 82208, 5.5016/
       data rms(133,1), rms(133,2) / 83209, 5.5167/
       data rms(134,1), rms(134,2) / 90232, 5.7087/
       data rms(135,1), rms(135,2) / 92238, 5.8470/

c----------------------------------------------------------------------

      do i = 1, irnm

         if( nint( z*1000 + a ) .eq. nint( rms(i,1) ) ) then

            radius = rms(i,2)
            return

         end if

      end do

c----------------------------------------------------------------------
c     if there is no data, a fitted values is used.
c     the formula is cited from NIM-B 152(1999)425-431
c----------------------------------------------------------------------

      radius = ( 0.84 * a**(1./3.) + 0.55 )

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      block data decayd
*                                                                      *
*       life time of decay particles                                   *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      common /dcyd/ rtdcy(50), kfdcy(50), ndcy

      data ndcy / 18 /
      data ( kfdcy(i), rtdcy(i), i = 1, 18 ) /
     &       211,   2.6029d-8,      ! pion+
     &       111,   1.0d-30,        ! pion0
     &      -211,   2.6029d-8,      ! pion-
     &       -13,   2.19703d-6,     ! muon+
     &        13,   2.19703d-6,     ! muon-
     &       321,   1.2371d-8,      ! kaon+
     &       311,   0.8922d-10,     ! kaon0
     &      -321,   1.2371d-8,      ! kaon+
     &       221,   1.0d-30,        ! eta
     &       331,   1.0d-30,        ! eta'
     &      -311,   0.8922d-10,     ! k0bar
     &      3122,   2.631d-10,      ! Lambda0
     &      3222,   0.799d-10,      ! Sigma+
     &      3212,   1.0d-30,        ! Sigma0
     &      3112,   1.479d-10,      ! Sigma-
     &      3322,   2.90d-10,       ! Xi0
     &      3312,   1.639d-10,      ! Xi-
     &      3334,   0.822d-10/      ! Omega-
      end
