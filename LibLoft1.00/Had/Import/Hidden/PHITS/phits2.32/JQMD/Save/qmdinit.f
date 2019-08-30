************************************************************************
*                                                                      *
      subroutine qmdint
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to initialize the input parameters for QMD              *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param01.inc'
      include 'param02.inc'

*-----------------------------------------------------------------------

      common /input1/ mstq1(mxpa1), parq1(mxpa1)

      common /const0/ idnta, idnpr, massta, masspr, mstapr, msprpr
      common /const1/ elab, rdist, bmin, bmax, ibch, ibin
      common /const2/ dt, ntmax, iprun, iprun0
      common /const3/ nfreq, nfrec, nfred

      common /const4/ plab, srtcm, ylabb, pincm
      common /const5/ pzpr, pxpr, rzpr, rxpr
      common /const6/ pzta, pxta, rzta, rxta
      common /const7/ betpr, gampr, prmas, radpr
      common /const8/ betta, gamta, tamas, radta

      common /vriab0/ massal, massba, mmeson
      common /vriab1/ b, llnow, ntnow
      common /rannum/ iseed, iseed0, iseed1

      common /impact/ bval(1000), ibnum(1000), bweight(1000), bdef

      common /poten0/ aaa, bbb, rpot, esymm
      common /poten1/ gamm, c0, c3, cs, cl, wl
      common /poten2/ t0, t3, rkk

      common /swich0/ icoul, irelcr, iavoid
      common /swich1/ ipot, insys, irkg, icolt
      common /swich2/ icfg, imany, icpus, idatm
      common /swich3/ ielst, jelst, kelst
      common /swich4/ ifin, ifout

      common /grndc0/ dsam, ddif, dsam2, ddif2
      common /grndc1/ cdp, c0p, c3p, csp, clp
      common /grndc2/ r00, r01, saa, rada, radb
      common /grndc3/ ipchs, mntry, dtg, fric

      common /gradu0/ c0g, c3g, csg, pag
      common /caldis/ c0w, c3w, clw, c0sw
      common /pauli0/ cpw, cph, cpc

      common /startf/ iday0,imon0,iyer0,ihor0,imin0,isec0

      common /coultr/ eccm, pzcc, rmax0, zeroz
      common /framtr/ betafr(0:2), gammfr(0:2)

      common /cood2b/ rha(nnn,nnn), rhe(nnn,nnn), rhc(nnn,nnn)
      common /cood2r/ rr2(nnn,nnn), rbij(nnn,nnn), pp2(nnn,nnn)

      common /fnamid/ fname(mxfnm), idf(mxfnm), ifnl(mxfnm)
      character       fname*80

      common /summ04/ ianal
      common /sdmcut/ sdmemin

      logical         exex

*-----------------------------------------------------------------------
*           output unit for error message
*-----------------------------------------------------------------------

                  ieo  = 6

*-----------------------------------------------------------------------

                  sdmemin = parq1(120)

*-----------------------------------------------------------------------
*           Basic Input [1-9]
*-----------------------------------------------------------------------

                  idnpr  = mstq1(1)
                  masspr = mstq1(2)
                  msprpr = mstq1(3)
                  idnta  = mstq1(4)
                  massta = mstq1(5)
                  mstapr = mstq1(6)
                  iprun  = mstq1(7)
                  ntmax  = mstq1(8)
                  insys  = mstq1(9)

                  elab   = parq1(1)
                  plab   = parq1(2)
                  bmin   = parq1(3)
                  bmax   = parq1(4)
                  dt     = parq1(5)

*-----------------------------------------------------------------------
*           Numerics and Control [10-29]
*-----------------------------------------------------------------------

                  iseed  = mstq1(10)
                  ibch   = mstq1(11)
                  ibin   = mstq1(12)
                  irkg   = mstq1(13)
                  imany  = mstq1(14)
                  ifout  = mstq1(15)
                  ifin   = mstq1(16)
                  ielst  = mstq1(17)

*-----------------------------------------------------------------------
*           Interaction [30-59]
*-----------------------------------------------------------------------

                  ipot   = mstq1(30)
                  icoul  = mstq1(31)
                  irelcr = mstq1(32)

*-----------------------------------------------------------------------

                  wl     = parq1(30)
                  rpot0  = parq1(31)
                  esymm  = parq1(32) * 0.001

*-----------------------------------------------------------------------
*           Collision [60-89]
*-----------------------------------------------------------------------

                  icolt  = mstq1(60)
                  iavoid = mstq1(61)

*-----------------------------------------------------------------------
*           Ground state [90-119]
*-----------------------------------------------------------------------

                  ipchs  = mstq1(90)
                  mntry  = mstq1(91)

*-----------------------------------------------------------------------

                  saa    = parq1(90)
                  r00    = parq1(91)
                  r01    = parq1(92)

                  rada   = parq1(93)
                  radb   = parq1(94)

                  dtg    = parq1(95)
                  fric   = parq1(96)

                  rdist  = parq1(97)

*-----------------------------------------------------------------------
*           Analysis [150-200]      detect cpu time, date and time.
*                                   the others are in sm_init
*-----------------------------------------------------------------------

                  icpus  = mstq1(169)
                  idatm  = mstq1(170)

                  ianal  = mstq1(175)

*-----------------------------------------------------------------------
*        Initial vlues of llnow, ntnow and iprun0
*-----------------------------------------------------------------------

                  llnow  = 0
                  ntnow  = 0
                  iprun0 = iprun

*-----------------------------------------------------------------------
*        Set the mass of projectile and target particle
*-----------------------------------------------------------------------

                  prmas  = ulmass( idnpr )
                  tamas  = ulmass( idnta )

*-----------------------------------------------------------------------
*        Set beam energy and momentum.
*-----------------------------------------------------------------------

               if( elab .gt. 0.0 .and. plab .le. 0.0 ) then

                     plab = sqrt( elab * ( 2.0 * prmas + elab ) )

               else if( elab .le. 0.0 .and. plab .gt. 0.0 ) then

                     elab = sqrt( prmas**2 + plab**2 ) - prmas

               else if( elab .gt. 0.0 .and. plab .gt. 0.0 ) then

                     plabp = sqrt( elab * ( 2.0 * prmas + elab ) )

                  if( abs( plab - plabp ) .gt. 0.000001 ) then

                     write(ieo,'('' Error: elab and plab'',
     &                           '' is mismatched'')')
                     write(ieo,'('' ====='')')
                     write(ieo,'(/
     &                '' Program is stopped at subroutine qmdint'')')

                     call parastop( 222 )

                  else

                     plab = plabp

                  end if

               else

                     elab = 0.0
                     plab = 0.0

               end if

*-----------------------------------------------------------------------
*        Set parameters according to the choice of frame.
*-----------------------------------------------------------------------

                  n1 = masspr
                  n2 = massta

*-----------------------------------------------------------------------
*         no target is error
*-----------------------------------------------------------------------

         if( n2 .eq. 0 ) then

                  write(ieo,'('' Error: Target is not specified'')')
                  write(ieo,'('' ====='')')
                  write(ieo,'(/
     &                '' Program is stopped at subroutine qmdint'')')

                  call parastop( 222 )

*-----------------------------------------------------------------------

         else if( n1 .gt. 0 ) then

*-----------------------------------------------------------------------

                  ee   = elab + prmas + tamas
                  srt  = sqrt( ee**2 - plab**2 )
                  pstn = pcmsr( srt, prmas, tamas )

                  ptot = plab * n1
                  etot = elab * n1 + prmas * n1 + tamas * n2
                  stot = sqrt( etot**2 - ptot**2 )
                  pstt =  pcmsr( stot, prmas*n1, tamas*n2 )

                  pzcc = pstt
                  eccm = stot - ( prmas * n1 + tamas * n2 )

*-----------------------------------------------------------------------
*           Events defined in the CM frame.
*-----------------------------------------------------------------------

            if( insys .eq. 1 ) then

                  p1   =  pstt / n1
                  p2   = -pstt / n2
                  e1   =  sqrt( prmas**2 + p1**2 )
                  e2   =  sqrt( tamas**2 + p2**2 )

                  ylab = 0.5 * log( ( etot + ptot ) / ( etot - ptot ) )

                  betacm = 0.0
                  gammcm = 1.0

                  betalb = - ptot / etot
                  gammlb =   etot / stot

                  betann = ( p1 + p2 ) / ( e1 + e2 )
                  gammnn = ( e1 + e2 )
     &                   / sqrt( ( e1 + e2 )**2 - ( p1 + p2 )**2 )

*-----------------------------------------------------------------------
*           Events defined in Lab ( fixed target ) frame.
*-----------------------------------------------------------------------

            else if( insys .eq. 0 ) then

                  p1   = plab
                  p2   = 0.0
                  e1   = sqrt( prmas**2 + p1**2 )
                  e2   = sqrt( tamas**2 + p2**2 )

                  ylab = 0.0

                  betacm = ( p1 * n1 + p2 * n2 )
     &                   / ( e1 * n1 + e2 * n2 )
                  gammcm = ( e1 * n1 + e2 * n2 )
     &                   / sqrt( ( e1 * n1 + e2 * n2 )**2
     &                         - ( p1 * n1 + p2 * n2 )**2 )

                  betalb = 0.0
                  gammlb = 1.0

                  betann = ( p1 + p2 ) / ( e1 + e2 )
                  gammnn = ( e1 + e2 )
     &                   / sqrt( ( e1 + e2 )**2 - ( p1 + p2 )**2 )

*-----------------------------------------------------------------------
*           Frame defined in NN CM.
*-----------------------------------------------------------------------

            else if( insys .eq. 2 ) then

                  p1   =  pstn
                  p2   = -pstn

                  e1   = sqrt( prmas**2 + p1**2 )
                  e2   = sqrt( tamas**2 + p2**2 )

                  ylab = 0.5 * log( ( ee + plab ) / ( ee - plab ) )

                  betacm = ( p1 * n1 + p2 * n2 )
     &                   / ( e1 * n1 + e2 * n2 )
                  gammcm = ( e1 * n1 + e2 * n2 )
     &                   / sqrt( ( e1 * n1 + e2 * n2 )**2
     &                         - ( p1 * n1 + p2 * n2 )**2 )

                  betalb = - plab / ee
                  gammlb =   ee / srt

                  betann = 0.0
                  gammnn = 1.0

*-----------------------------------------------------------------------
*           Unrecognize frame : Error
*-----------------------------------------------------------------------

            else

                  write(ieo,'('' Error: Unrecognized input'',
     &                        '' coordinate frame'')')
                  write(ieo,'('' ====='')')

                  write(ieo,'(/
     &                '' Program is stopped at subroutine qmdint'')')

                  call parastop( 222 )

            end if

*-----------------------------------------------------------------------
*           Set values
*-----------------------------------------------------------------------

                  srtcm  = srt
                  ylabb  = ylab
                  pincm  = pstn

                  pzpr   = p1
                  pxpr   = 0.0
                  betpr  = p1 / e1
                  gampr  = e1 / prmas

                  pzta   = p2
                  pxta   = 0.0
                  betta  = p2 / e2
                  gamta  = e2 / tamas

                  betafr(0) = betalb
                  gammfr(0) = gammlb
                  betafr(1) = betacm
                  gammfr(1) = gammcm
                  betafr(2) = betann
                  gammfr(2) = gammnn

*-----------------------------------------------------------------------
*        Target only
*-----------------------------------------------------------------------

         else if( n1 .eq. 0 ) then

                  elab   = 0.0
                  plab   = 0.0

                  srtcm  = 0.0
                  ylabb  = 0.0
                  pincm  = 0.0

                  pzpr   = 0.0
                  pxpr   = 0.0
                  betpr  = 0.0
                  gampr  = 1.0

                  pzta   = 0.0
                  pxta   = 0.0
                  betta  = 0.0
                  gamta  = 1.0

                  eccm   = 0.0
                  pzcc   = 0.0

                  betafr(0) = 0.0
                  gammfr(0) = 1.0
                  betafr(1) = 0.0
                  gammfr(1) = 1.0
                  betafr(2) = 0.0
                  gammfr(2) = 1.0

         end if

*-----------------------------------------------------------------------
*        massal, massba
*-----------------------------------------------------------------------

                  massal = masspr + massta

                  massba = 0.0
                  if( idnpr .eq. 0 ) massba = massba + masspr
                  if( idnta .eq. 0 ) massba = massba + massta


*-----------------------------------------------------------------------
*        Radius of Nucleus
*-----------------------------------------------------------------------

                  radpr = 0.0
                  radta = 0.0

                  if( n1 .gt. 2 ) radpr = r00 * float( n1 )**(1./3.)
                  if( n2 .gt. 2 ) radta = r00 * float( n2 )**(1./3.)

*-----------------------------------------------------------------------
*        Initial Distance : rmax0
*-----------------------------------------------------------------------
*
*        rdist ( D = -1.0 fm ) : initial distance
*
*            < 0.0: 
*                 elab < 200 MeV
*                       rmax0 = radta/gamta + radpr/gampr + 6.0
*                 200 MeV < elab < 1 GeV
*                       rmax0 = radta/gamta + radpr/gampr + 4.0
*                 elab > 1 GeV
*                       rmax0 = radta/gamta + radpr/gampr + 4.0
*                 elab > 10 GeV
*                       rmax0 = radta/gamta + radpr/gampr + 2.0
*            >=0.0: rmax0 = radta/gamma + radpr/gamma + rdist
*
*-----------------------------------------------------------------------

               if( rdist .lt. 0.1 ) then

                  if( elab .ge. 10.0 ) then

                     rminm = 2.0

                  else if( elab .ge. 1.0 ) then

                      rminm = 4.0

                  else if( elab .ge. 0.2 ) then

                     rminm = 4.0

                  else

                     rminm = 6.0

                  end if

               else

                     rminm = rdist

               end if

                  rmax0 = radpr / gampr + radta / gamta + rminm

*-----------------------------------------------------------------------
*        Shift distance of z and x.
*-----------------------------------------------------------------------

               rzpr  = -rmax0 * float( n2 ) / float( n2 + n1 )
               rxpr  =  0.0

               rzta  =  rmax0 * float( n1 ) / float( n2 + n1 )
               rxta  =  0.0

*-----------------------------------------------------------------------
*        shift distance for Lab. system
*-----------------------------------------------------------------------

            if( insys .eq. 0 ) then

               zeroz = rmax0 * float( n1 ) / float( n2 + n1 )
     &               + 16.0 - radpr - rmax0

            else

               zeroz = 0.0

            end if

*-----------------------------------------------------------------------
*        check of impact parameter
*-----------------------------------------------------------------------

            if( bmin .gt. bmax .or. bmin .lt. 0.0 ) then

               write(ieo,'('' Error: invalid input bmin or bmax'')')
               write(ieo,'('' ====='')')

               write(ieo,'(/
     &             '' Program is stopped at subroutine qmdint'')')

               call parastop( 222 )

            end if

*-----------------------------------------------------------------------
*        check of irkg
*-----------------------------------------------------------------------

            if( irkg .ne. 2 .and. irkg .ne. 4 ) irkg = 2


*-----------------------------------------------------------------------
*        Choice of the impact parameter bin
*-----------------------------------------------------------------------

                  if( ibch .lt. 0  .or.  ibch .gt. 2 ) ibch = 1

                  if( ifout .gt. 0 .and. ibch .eq. 2 ) ibch = 1
                  if( ifin  .gt. 0 .and. ibch .eq. 1 ) ibch = 2

*-----------------------------------------------------------------------
*           initialization of impact parameter and weight
*-----------------------------------------------------------------------

               if( ibch .eq. 0 ) then

                     ibin = 1
                     bdef = ( bmax - bmin ) / float(ibin)

                     ibnum(1)   = 0
                     bval(1)    = bmin + bdef / 2.0
                     bweight(1) = 2.0 * pi * bval(1) * bdef * 10.0

               else if( ibch .eq. 1 ) then

                     if( iprun .lt. ibin ) ibin = iprun

                     jprun  = iprun / ibin
                     iprun  = jprun * ibin
                     iprun0 = iprun

                     bdef = ( bmax - bmin ) / float(ibin)

                  do i = 1, ibin

                     ibnum(i)   = jprun
                     bval(i)    = bmin + bdef / 2.0 + float(i-1) * bdef
                     bweight(i) = 2.0 * pi * bval(i) * bdef * 10.0
     &                          * float(ibin)

                  end do

               else if( ibch .eq. 2 ) then

                     bdef = ( bmax - bmin ) / float(ibin)

                  do i = 1, ibin

                     ibnum(i)   = 0
                     bval(i)    = bmin + bdef / 2.0 + float(i-1) * bdef
                     bweight(i) = 2.0 * pi * bval(i) * bdef * 10.0

                  end do

               end if

*-----------------------------------------------------------------------
*        Choice of the judgement of elastic collision
*-----------------------------------------------------------------------

                  if( ielst .gt. 3 .or. ielst .lt. 0 ) ielst = 3

*-----------------------------------------------------------------------
*        QMD with only skyrme : Set potential parameters
*-----------------------------------------------------------------------

               if( ipot .eq. 1 ) then

                  rpot = 0.333333

               else if( ipot .eq. 2 ) then

                  rpot = 0.666667

               else if( ipot .eq. 3 ) then

                  rpot = 1.0

               else

                  write(ieo,'('' Error: invalid input mstq1(30)'')')
                  write(ieo,'('' ====='')')

                  write(ieo,'(/
     &                '' Program is stopped at subroutine qmdint'')')

                  call parastop( 222 )

               end if

               if( rpot0 .gt. 0.0 ) rpot = rpot0

*-----------------------------------------------------------------------

               gamm =  rpot + 1.0

               call kcomp(rpot,t3,t0,aaa,bbb,rkk)


*-----------------------------------------------------------------------
*        Local potentials.
*-----------------------------------------------------------------------

               c0 = aaa/(rho0*(4*pi*wl)**1.5*2.0)
               c3 = bbb/(rho0**gamm*(4.0*pi*wl)**(1.5*gamm)*(gamm+1.0))
               cs = esymm/(rho0*(4.0*pi*wl)**1.5*2.0)
               cl = ccoul/2.0 * icoul

*-----------------------------------------------------------------------
*        Pauli
*-----------------------------------------------------------------------

               cpw = 1.0 / 2.0 / wl
               cph = 2.0 * wl / hbc**2
               cpc = 4.0

*-----------------------------------------------------------------------
*        Ground State
*-----------------------------------------------------------------------

               dsam  = 1.5
               ddif  = 1.0

               dsam2 = dsam * dsam
               ddif2 = ddif * ddif

               cdp = 1.0 / ( 4.0 * pi * wl )**1.5

               c0p = c0 * 2.0
               c3p = c3 * ( gamm + 1.0 )
               csp = cs * 2.0
               clp = cl * 2.0

            if( masspr .eq. 0 ) then

               rzpr = 0.0

               bmin = 0.0
               bmax = 0.0

            end if

*-----------------------------------------------------------------------
*        gradu
*-----------------------------------------------------------------------

               c0g = - c0 / ( 2.0 * wl )
               c3g = - c3 / ( 4.0 * wl ) * gamm
               csg = - cs / ( 2.0 * wl )
               pag =   gamm - 1.0

*-----------------------------------------------------------------------
*        caldis
*-----------------------------------------------------------------------

               c0w  = 1.0 /   4.0 / wl
               c3w  = 1.0 /   4.0 / wl
               clw  = 2.0 / ( 4.0 * pi * wl )**0.5
               c0sw = sqrt(c0w)

*-----------------------------------------------------------------------
*        set diagonal elements of two-body quantities to be zero
*-----------------------------------------------------------------------

            do i = 1, nnn
*
               rha(i,i)  =   0.0
               rhe(i,i)  =   0.0
               rhc(i,i)  =   0.0
               rr2(i,i)  =   0.0
               pp2(i,i)  =   0.0
              rbij(i,i)  =   0.0
*
            end do

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine kcomp(rpot,t3,t0,aaa,bbb,rkk)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 26                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to set the parameter of Skyrme force and                *
*              calculate incompressibility in the case of simple       *
*              Skyrme interaction.                                     *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              rpot        pawer of density dependent                  *
*              t0, t3,     coefficients of Skyrme                      *
*              aaa, bbb    coefficients of Skyrme                      *
*              rkk         incompressibility                           *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param01.inc'
      include 'param02.inc'

*-----------------------------------------------------------------------
*        Calculate nuclear matter properties.
*-----------------------------------------------------------------------

               ebin  = ebinm * 0.001
               pfer  = hbc * ( 3./2. * pi**2 * rho0 )**(1./3.)
               efer  = pfer**2 / 2. / rmass

*-----------------------------------------------------------------------
*        t0, t3 and aaa, bbb
*-----------------------------------------------------------------------

               t3   =  8./3./rpot/rho0**(1.+rpot)*(efer/5.-ebin)
               t0   = -16./15.*efer/rho0 - (1.+rpot)*t3*rho0**rpot

               aaa  =  3./4.*t0*rho0
               bbb  =  3./8.*t3*(2.+rpot)*rho0**(1.+rpot)

*-----------------------------------------------------------------------
*        Incompressibility
*-----------------------------------------------------------------------

               rkk  =  6./5.*efer+9./4.*t0*rho0
     &              + 9./8.*(1.+rpot)*(2.+3.*rpot)*t3*rho0**(1.+rpot)


*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      function ulmass(kf)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 26                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to give the mass of the partilce with id = kf           *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              kf  : particle kf id                                    *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param01.inc'
      include 'param02.inc'

*-----------------------------------------------------------------------

            if( kf .eq. 0 ) then

                  ulmass = rmass

            else if( kf .eq. 2112 .or. kf .eq. 2212 ) then

                  ulmass = rmass

            else if( kf .eq.  211 .or.
     &               kf .eq. -211 .or.
     &               kf .eq.  111 ) then

                  ulmass = pmass

            else

                  ulmass = 0.0

            end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      function erf(x)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to calculate the eror function                          *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              x           : x value for error function                *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      parameter(  dmax = 5.0 )

*-----------------------------------------------------------------------

      dimension      a(6)

      data a / 0.0705230784, 0.0422820123,
     &         0.0092705272, 0.0001520143,
     &         0.0002765672, 0.0000430638 /

*-----------------------------------------------------------------------

               dnm = 1.0
               xx  = 1.0

            do i = 1, 6

               xx  = xx * x
               dnm = dnm + xx * a(i)

            end do

            if( dnm .gt. dmax ) then

               erf = 1.0

            else

               erf = 1.0 - 1.0 / dnm**16

            end if

*-----------------------------------------------------------------------

      return
      end


