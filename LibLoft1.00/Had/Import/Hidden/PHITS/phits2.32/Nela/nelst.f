#include "ZcheckPHITS.h"
************************************************************************
*                                                                      *
      subroutine nelst(ityp,ktyp,ata,atz,epin)

*                                                                      *
*       control of elastic collision                                   *
*       last modified by K.Niita on 09/02/2000                         *
*                                                                      *
*     input:                                                           *
*                                                                      *
*        ityp     : type of incident particle, proton or neutron only  *
*        ktyp     : kf code of incident particle                       *
*        ata      : target mass number                                 *
*        atz      : target charge                                      *
*        epin     : kinetic energy of incident particle in lab (MeV)   *
*                                                                      *
*     output:                                                          *
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
*                  = 8, weight change change                           *
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
*                  = 11, other particles                               *
*                  = 14, gammma                                        *
*                  = 15, deuteron                                      *
*                  = 16, triton                                        *
*                  = 17, 3He                                           *
*                  = 18, Alpha                                         *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

#include "param00.inc"

*-----------------------------------------------------------------------

      common /clustf/ nclst, iclust(nnn)
      common /clustg/ jclust(0:7,nnn), qclust(0:12,nnn)
      common /clustp/ rumpat(0:20), numpat(0:20)

*-----------------------------------------------------------------------

      dimension rmspat(20)
      data rmspat/ 0.93827, 0.93958, 0.1396, 0.1350, 0.1396, 0.0, 0.0,
     &             0.4936, 0.4977, 0.4936, 10*0.0/

*-----------------------------------------------------------------------
*        proton or neutron
*-----------------------------------------------------------------------
c///////////
#if defined  (CHECKPHITS)
      write(0,*) ' in nelst (ityp,ktyp,ata,atz,epin)'
      write(0,*) ityp,ktyp,ata,atz,epin
#endif
c//////////////////
               if( ityp .gt. 2 ) return

*-----------------------------------------------------------------------
*        determine csthcm : call seldsd
*        csthcm   : cosine of scattering angle in cm
*-----------------------------------------------------------------------

                  ian = 1

                  call seldsd(ityp,epin,ata,atz,csthcm,ian)

*-----------------------------------------------------------------------
*        initial values
*-----------------------------------------------------------------------

                  sinecm = sqrt( 1.0 - csthcm**2 )

                  mast = nint( ata )
                  mcht = nint( atz )

                  rmst = ata * 0.93895
                  rmsp = rmspat(ityp)

                  ekin = epin / 1000.0

                  etop = ekin + rmsp
                  pinp = sqrt( ekin**2 + 2.0 * rmsp * ekin )

                  etot = rmst

                  eall = etop + etot

*-----------------------------------------------------------------------
*        Lorentz factor
*-----------------------------------------------------------------------

                  betc = pinp / eall
                  gamc = eall / sqrt( eall**2 - pinp**2 )

*-----------------------------------------------------------------------
*        Lorentz transform
*-----------------------------------------------------------------------

                  pinc = pinp * gamc - betc * gamc * etop
                  etpc = sqrt( pinc**2 + rmsp**2 )
                  ettc = sqrt( pinc**2 + rmst**2 )

                  picz = pinc * csthcm
                  pixy = pinc * sinecm

                  thet = unirn(dummy) * 3.1415926536 * 2.0

                  csph =  cos(thet)
                  snph =  sin(thet)

                  pipx =  pixy * csph
                  pipy =  pixy * snph
                  pipz =  picz * gamc + betc * gamc * etpc

                  pipl =  sqrt( pipx**2 + pipy**2 + pipz**2 )
                  etpl =  sqrt( pipl**2 + rmsp**2 )

                  pitx = -pixy * csph
                  pity = -pixy * snph
                  pitz = -picz * gamc + betc * gamc * ettc

                  pitl =  sqrt( pitx**2 + pity**2 + pitz**2 )
                  ettl =  sqrt( pitl**2 + rmst**2 )

*-----------------------------------------------------------------------
*        booking of outgoing particles
*-----------------------------------------------------------------------

                  nclst = 1

                  qclust(0,nclst)  = -1.0
                  qclust(1,nclst)  = pipx
                  qclust(2,nclst)  = pipy
                  qclust(3,nclst)  = pipz
                  qclust(4,nclst)  = etpl
                  qclust(5,nclst)  = rmsp
                  qclust(6,nclst)  = 0.0
                  qclust(7,nclst)  = ( etpl - rmsp ) * 1000.
                  qclust(8,nclst)  = 1.0
                  qclust(9,nclst)  = 0.0
                  qclust(10,nclst) = 0.0d0
                  qclust(11,nclst) = 0.0d0
                  qclust(12,nclst) = 0.0d0

                  iclust(nclst)    = ityp

               if( ityp .eq. 1 ) then

                  jclust(1,nclst) = 1
                  jclust(2,nclst) = 0
                  jclust(3,nclst) = 1
                  jclust(5,nclst) = 1

               else if( ityp .eq. 2 ) then

                  jclust(1,nclst) = 0
                  jclust(2,nclst) = 1
                  jclust(3,nclst) = 2
                  jclust(5,nclst) = 0

               end if

                  jclust(0,nclst) = 0
                  jclust(4,nclst) = 0
                  jclust(6,nclst) = 1
                  jclust(7,nclst) = ktyp

                  numpat(1)  = jclust(1,nclst)
                  numpat(2)  = jclust(2,nclst)

                  rumpat(1)  = jclust(1,nclst)
                  rumpat(2)  = jclust(2,nclst)

*-----------------------------------------------------------------------
*           booking of the target
*-----------------------------------------------------------------------

                  numpat(0) = 1
                  rumpat(0) = 1

                  nclst = nclst + 1

                  iclust(nclst) = 0

                  jclust(0,nclst)  = 0
                  jclust(1,nclst)  = mcht
                  jclust(2,nclst)  = mast - mcht
                  jclust(3,nclst)  = 19
                  jclust(4,nclst)  = 0
                  jclust(5,nclst)  = mcht
                  jclust(6,nclst)  = mast
                  jclust(7,nclst)  = mcht * 1000000 + mast

                  qclust(0,nclst)  = 0.0
                  qclust(1,nclst)  = pitx
                  qclust(2,nclst)  = pity
                  qclust(3,nclst)  = pitz
                  qclust(4,nclst)  = ettl
                  qclust(5,nclst)  = rmst
                  qclust(6,nclst)  = 0.0
                  qclust(7,nclst)  = ( ettl - rmst ) * 1000.
                  qclust(8,nclst)  = 1.0
                  qclust(9,nclst)  = 0.0
                  qclust(10,nclst) = 0.0d0
                  qclust(11,nclst) = 0.0d0
                  qclust(12,nclst) = 0.0d0

*-----------------------------------------------------------------------
c//////////////
#if defined (CHECKPHITS)
      write(0,*) 'on ret form nelst:  nclst=',nclst
      write(0,*) '#   kf   Px      Py     Pz (MeV)     KE (MeV)'//
     * '        p#    n#    Z'
      do i = 1, nclst
         write(0,'(i3, i9,4g12.4,3i4)')
     *     i, jclust(7,i), qclust(1:3,i)*1000. ,  qclust(7,i),
     *     jclust(1,i), jclust(2,i), jclust(5,i)
*        jclust    = 1, proton number                                  *
*                  = 2, neutron number                                 *
*                  = 3, ip, see below                                  *
*                  = 4,                                                *
*                  = 5, charge                                         *
*                  = 6, baryon number                                  *
*                  = 7, kf code                                        *
      enddo
#endif
c////////////
      return
      end

