************************************************************************
*                                                                      *
      subroutine bertin(jcasc,ityp,eein,mmas,mchg)
*                                                                      *
*                                                                      *
*       control routine of bertini and isobert                         *
*       modified by K.Niita on 2005/08/15                              *
*                                                                      *
*        call subroutine : bert and isbert                             *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       jcasc   : =1 : bertini, =2 : isobert, =3 : JAM                 *
*       ityp    : particle type of projectile                          *
*       eein    : energy of projectile (MeV)                           *
*       mmas    : mass of target                                       *
*       mchg    : charge of target                                     *
*                                                                      *
*---- in common -------------------------------------------------------*
*                                                                      *
*        nclst   : total number of out going particles and nuclei      *
*                                                                      *
*        iclust(nclst)                                                 *
*                                                                      *
*                i = 0, nucleus                                        *
*                  = 1, proton                                         *
*                  = 2, neutron                                        *
*                  = 3, pion                                           *
*                  = 4, photon                                         *
*                  = 5, kaon                                           *
*                  = 6, muon                                           *
*                  = 7, others                                         *
*                                                                      *
*        jclust(i,nclst)                                               *
*                                                                      *
*                i = 0, angular momentum                               *
*                  = 1, proton number                                  *
*                  = 2, neutron number                                 *
*                  = 3, ip, see below                                  *
*                  = 4,                                                *
*                  = 5, charge                                         *
*                  = 6, baryon number                                  *
*                  = 7, kf code                                        *
*                                                                      *
*        qclust(i,nclst)                                               *
*                                                                      *
*                i = 0, impact parameter                               *
*                  = 1, px (GeV/c)                                     *
*                  = 2, py (GeV/c)                                     *
*                  = 3, pz (GeV/c)                                     *
*                  = 4, etot = sqrt( p**2 + rm**2 ) (GeV)              *
*                  = 5, rest mass (GeV)                                *
*                  = 6, excitation energy (MeV)                        *
*                  = 7, kinetic energy (MeV)                           *
*                  = 8, weight change                                  *
*                  = 9, delay time                                     *
*                  = 10, x-displace                                    *
*                  = 11, y-displace                                    *
*                  = 12, z-displace                                    *
*                                                                      *
*        numpat(i) : total number of out going particles or nuclei     *
*                                                                      *
*                i =  0, nuclei                                        *
*                  =  1, proton                                        *
*                  =  2, neutron                                       *
*                  =  3, pi+                                           *
*                  =  4, pi0                                           *
*                  =  5, pi-                                           *
*                  =  6, mu+                                           *
*                  =  7, mu-                                           *
*                  =  8, K+                                            *
*                  =  9, K0                                            *
*                  = 10, K-                                            *
*                                                                      *
*                  = 11, other particles                               *
*                                                                      *
*                  = 12, electron                                      *
*                  = 13, positron                                      *
*                  = 14, photon                                        *
*                                                                      *
*                  = 15, deuteron                                      *
*                  = 16, triton                                        *
*                  = 17, 3He                                           *
*                  = 18, Alpha                                         *
*                  = 19, residual nucleus                              *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include '../JQMD/param00.inc'
c      #include "param00.inc" ; dont use this

      parameter ( pi  = 3.1415926535898d0 )

*-----------------------------------------------------------------------

      common /bparm/  andt,jevap,npidk
      common /qparm/  ielas,icasc,iqstep,lvlopt,igamma
      common /preeq/  npcle, nhole, efermi, atar, ztar

*-----------------------------------------------------------------------

      common/rtcom/tapcrs(6600)
      common/rtco1/tapcr1(352),tapcr2(352)

      common /cugnon/ icugn

*-----------------------------------------------------------------------

      common /clustf/ nclst, iclust(nnn)
      common /clustg/ jclust(0:7,nnn), qclust(0:12,nnn)
      common /clustp/ rumpat(0:20), numpat(0:20)

*-----------------------------------------------------------------------

      dimension rmspat(20)

      data rmspat/ 0.93827, 0.93958, 0.1396, 0.1350, 0.1396, 0.0, 0.0,
     &             0.4936, 0.4977, 0.4936, 10*0.0/

      dimension finput(7)

      dimension kind(nomp),ep(nomp),alpha(nomp),beta(nomp),gam(nomp)

      data initbe /0/

*-----------------------------------------------------------------------
*     initialization of cross section
*-----------------------------------------------------------------------
c////////////
c      write(0,*) 'bertin(jcasc,ityp,eein,mmas,mchg'
c      write(0,*) jcasc,ityp,eein,mmas,mchg
c      write(0,*) ' initbe=', initbe
c  mainly n, and p  1MeV< E < GeV
c         pi+ of energy 1.45GeV, 1 MeV.... came
c         pi-  //        1.75 //, 365.... MeV
c         pi0 can come ? but never seen.
c      
c////////////
         if( initbe .eq. 0 ) then

            initbe = initbe + 1

            call bertcin

         end if

*-----------------------------------------------------------------------
*        initial values
*-----------------------------------------------------------------------

               finput(1) = dble(mmas)
               finput(2) = dble(mchg)
               finput(3) = eein
               finput(4) = 0.0d0
               finput(5) = 1.d+0
               finput(6) = andt
               finput(7) = dble(ityp-1)

*-----------------------------------------------------------------------

               nneut = 0
               nprot = 0

               npipo = 0
               nping = 0
               npine = 0

               sume = 0.0

               poutx = 0.0
               pouty = 0.0
               poutz = 0.0

*-----------------------------------------------------------------------

               eppin = eein / 1000.0
               rmsin = rmspat(ityp)

            if( ityp .eq. 1 ) then

               masim = mmas + 1
               mchim = mchg + 1
               einad = 0.0

            else if( ityp .eq. 2 ) then

               masim = mmas + 1
               mchim = mchg
               einad = 0.0

            else if( ityp .eq. 3 ) then

               masim = mmas
               mchim = mchg + 1
               einad = 139.9d+0

            else if( ityp .eq. 4 ) then

               masim = mmas
               mchim = mchg
               einad = 135.1d+0

            else if( ityp .eq. 5 ) then

               masim = mmas
               mchim = mchg - 1
               einad = 139.9d+0

            else

               goto 340

            end if

*-----------------------------------------------------------------------
*        Bertini or Isobert
*-----------------------------------------------------------------------

         if( jcasc .eq. 1 ) then

            call bert(finput,nopart,nomp,kind,ep,alpha,beta,gam)

         else if( jcasc .eq. 2 ) then

            call isbert(finput,nopart,nomp,kind,ep,alpha,beta,gam)

         end if

*-----------------------------------------------------------------------
*        pseudo collision
*-----------------------------------------------------------------------

               nclst = 0

         if( nopart .lt. 0 ) then

               nclst = -1

               return

         end if

*-----------------------------------------------------------------------
*        real collision
*-----------------------------------------------------------------------

         if( nopart .gt. 0 ) then

            do n = 1, nopart

                  nclst = nclst + 1

                  lk = kind(n) + 1

               if( lk .eq. 1 ) then

                  nprot = nprot + 1

                  kf    = 2212
                  ibary = 1
                  ipid  = 1
                  ippad = 1
                  ipprt = 1
                  ipneu = 0
                  ipchg = 1

               else if( lk .eq. 2 ) then

                  nneut = nneut + 1

                  kf    = 2112
                  ibary = 1
                  ipid  = 2
                  ippad = 2
                  ipprt = 0
                  ipneu = 1
                  ipchg = 0

               else if( lk .eq. 3 ) then

                  npipo = npipo + 1

                  kf    = 211
                  ibary = 0
                  ipid  = 3
                  ippad = 3
                  ipprt = 0
                  ipneu = 0
                  ipchg = 1

               else if( lk .eq. 4 ) then

                  npine = npine + 1

                  kf    = 111
                  ibary = 0
                  ipid  = 3
                  ippad = 4
                  ipprt = 0
                  ipneu = 0
                  ipchg = 0

               else if( lk .eq. 5 ) then

                  nping = nping + 1

                  kf    = -211
                  ibary = 0
                  ipid  = 3
                  ippad = 5
                  ipprt = 0
                  ipneu = 0
                  ipchg = -1

               else

                  goto 340

               end if

*-----------------------------------------------------------------------
*           enhancement of pion kinetic energy
*-----------------------------------------------------------------------

c              if( lk .ge. 3 .and. lk .le. 5 .and.
c    &             ep(n) .lt. 200.0 ) then
c
c                 aa = 60.0
c                 bb = 50.0
c                 cc = 50.0
c
c                 ep(n) = ep(n) + aa * ep(n) / bb / 2.0
c    &                  * exp( - ( ep(n) - bb )**2 / 2.0 / cc**2 )
c
c              end if

*-----------------------------------------------------------------------

                  sume = sume + ep(n)

                  epp = ep(n) / 1000.0
                  rms = rmspat(lk)

                  pouta = sqrt( epp**2 + 2.0 * epp * rms )

                  pxrv  = pouta * alpha(n)
                  pyrv  = pouta * beta(n)
                  pzrv  = pouta * gam(n)

                  poutx = poutx + pxrv
                  pouty = pouty + pyrv
                  poutz = poutz + pzrv

*-----------------------------------------------------------------------
*        booking of outgoing particles
*-----------------------------------------------------------------------

                  iclust(nclst)    = ipid

                  jclust(0,nclst)  = 0
                  jclust(1,nclst)  = ipprt
                  jclust(2,nclst)  = ipneu
                  jclust(3,nclst)  = ippad
                  jclust(4,nclst)  = 0
                  jclust(5,nclst)  = ipchg
                  jclust(6,nclst)  = ibary
                  jclust(7,nclst)  = kf

                  qclust(0,nclst)  = -1.0
                  qclust(1,nclst)  = pxrv
                  qclust(2,nclst)  = pyrv
                  qclust(3,nclst)  = pzrv
                  qclust(4,nclst)  = epp + rms
                  qclust(5,nclst)  = rms
                  qclust(6,nclst)  = 0.0
                  qclust(7,nclst)  = epp * 1000.
                  qclust(8,nclst)  = 1.0
                  qclust(9,nclst)  = 0.0
                  qclust(10,nclst) = 0.0d0
                  qclust(11,nclst) = 0.0d0
                  qclust(12,nclst) = 0.0d0

            end do

         end if

*-----------------------------------------------------------------------
*              booking of produced particles
*-----------------------------------------------------------------------

                  npcle  = npcle - nprot - nneut

                  numpat(1)  = nprot
                  numpat(2)  = nneut
                  numpat(3)  = npipo
                  numpat(4)  = npine
                  numpat(5)  = nping

                  rumpat(1)  = nprot
                  rumpat(2)  = nneut
                  rumpat(3)  = npipo
                  rumpat(4)  = npine
                  rumpat(5)  = nping

*-----------------------------------------------------------------------
*           residual nucleus
*-----------------------------------------------------------------------

                  noutn = nprot + nneut
                  noutc = nprot + npipo - nping

                  masrs = masim - noutn
                  mchrs = mchim - noutc

            if( masim .ge. noutn .and. mchim .ge. noutc .and.
     &          masrs .gt. 0 ) then

                  nclst = nclst + 1

                  piabs = sqrt( eppin**2 + 2.0 * eppin * rmsin )

                  presx = - poutx
                  presy = - pouty
                  presz = - poutz + piabs

                  pabst = sqrt( presx**2 + presy**2 + presz**2 )
                  rsmas = 0.9385d0 * dble( masrs )

                  etota = sqrt( pabst**2 + rsmas**2 )

                  erres = ( etota - rsmas ) * 1000.0

                  exres = eein + einad
     &                  - sume
     &                  - 139.9d+0 * ( npipo + nping )
     &                  - 135.1d+0 * npine
     &                  - erres
     &                  - bindeg(mchg,mmas-mchg)
     &                  + bindeg(mchrs,masrs-mchrs)

                  exres = max( 0.0d0, exres )

*-----------------------------------------------------------------------
*           renormalization by neglecting px and py
*-----------------------------------------------------------------------

c                 preno = abs( presz )
c
c                 presx = presx * preno / pabst
c                 presy = presy * preno / pabst
c                 presz = presz * preno / pabst
c
c                 pabst = sqrt( presx**2 + presy**2 + presz**2 )
c                 etota = sqrt( presz**2 + rsmas**2 )

*-----------------------------------------------------------------------
*           booking of the residual nucleus
*-----------------------------------------------------------------------

                  numpat(0) = 1
                  rumpat(0) = 1

                  iclust(nclst)    = 0

                  jclust(0,nclst)  = 0
                  jclust(1,nclst)  = mchrs
                  jclust(2,nclst)  = masrs - mchrs
                  jclust(3,nclst)  = 19
                  jclust(4,nclst)  = 0
                  jclust(5,nclst)  = mchrs
                  jclust(6,nclst)  = masrs
                  jclust(7,nclst)  = mchrs * 1000000 + masrs

                  qclust(0,nclst)  = 0.0
                  qclust(1,nclst)  = presx
                  qclust(2,nclst)  = presy
                  qclust(3,nclst)  = presz
                  qclust(4,nclst)  = etota
                  qclust(5,nclst)  = rsmas
                  qclust(6,nclst)  = exres
                  qclust(7,nclst)  = ( etota - rsmas ) * 1000.
                  qclust(8,nclst)  = 1.0
                  qclust(9,nclst)  = 0.0
                  qclust(10,nclst) = 0.0d0
                  qclust(11,nclst) = 0.0d0
                  qclust(12,nclst) = 0.0d0

            else

                  numpat(0) = 0
                  rumpat(0) = 0

            end if

*-----------------------------------------------------------------------
*           x-y rotation by randum
*-----------------------------------------------------------------------

                  theta = 2.0d0 * pi * rn(0)

               do i = 1, nclst

                  px = qclust(1,i)
                  py = qclust(2,i)

                  qclust(1,i)  = px * cos( theta ) - py * sin( theta )
                  qclust(2,i)  = px * sin( theta ) + py * cos( theta )

               end do

*-----------------------------------------------------------------------

      return

*-----------------------------------------------------------------------
*        error in betini
*-----------------------------------------------------------------------

  340 continue

               write(6,1002) ityp,lk
 1002          format(/' *** error message from s.bertin ***'
     &         /' invalid condition of itype or lk was found.'
     &         /' ityp =',i5,'  lk =',i5)
               call parastop( 844 )

*-----------------------------------------------------------------------

      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine bertcin
*                                                                      *
*       initialization of the cross sections                           *
*       for Cugnon parametrization                                     *
*       last modified by K.Niita on 08/02/2000                         *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       icugn   : =1 old Cugnon, =2 new Cugnon                         *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      common/rtcom/tapcrs(6600)
      common/rtco1/tapcr1(352),tapcr2(352)

      common /cugnon/ icugn

*-----------------------------------------------------------------------

            if( icugn .eq. 1 ) then

               do i = 1, 352

                  tapcrs(i+5991) = tapcr1(i)

               end do

            else if( icugn .eq. 2 ) then

               do i = 1, 352

                  tapcrs(i+5991) = tapcr2(i)

               end do

            end if

*-----------------------------------------------------------------------

      return
      end subroutine

