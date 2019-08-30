************************************************************************
*                                                                      *
      subroutine jqmdin(ityp,ktyp,eein,mmas,mchg,bmax0)
*                                                                      *
*       control routine of JQMD                                        *
*       modified by K.Niita on 2005/08/15                              *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*        icfg = 1, 2, 4                                                *
*                                                                      *
*              input       : content                     ; variables   *
*                                                                      *
*              'proj'      : projectile                  ; mstq1(1)    *
*                                                        ; mstq1(2)    *
*                                                        ; mstq1(3)    *
*              'targ'      : target                      ; mstq1(4)    *
*                                                        ; mstq1(5)    *
*                                                        ; mstq1(6)    *
*              'event'     : number of events            ; mstq1(7)    *
*              'tstep'     : total number of time step   ; mstq1(8)    *
*              'frame'     : reference frame             ; mstq1(9)    *
*                                                                      *
*                                                                      *
*              'win'       : incident energy or momentum ; parq1(1)    *
*                                                        ; parq1(2)    *
*              'bmin'      : minimum impact parameter    ; parq1(3)    *
*              'bmax'      : maximum impact parameter    ; parq1(4)    *
*              'dt'        : time step                   ; parq1(5)    *
*                                                                      *
*              'fname(i)'  : file name                   ; fname(i)    *
*              'mstq1(i)'  : integer parameters          ; mstq1(i)    *
*              'parq1(i)'  : real parameters             ; parq1(i)    *
*                                                                      *
*     output: in cldist                                                *
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
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param01.inc'

      parameter ( pi  = 3.1415926535898d0 )

*-----------------------------------------------------------------------

      common /swich2/ icfg, imany, icpus, idatm
      common /input1/ mstq1(mxpa1), parq1(mxpa1)

      common /swich3/ ielst, jelst, kelst

      common /clustf/ nclst, iclust(nnn)
      common /clustg/ jclust(0:7,nnn),  qclust(0:12,nnn)

      common /framtr/ betafr(0:2), gammfr(0:2)
      common /crshi/  bplus, icrhi, ijudg, imadj, iqmax

*-----------------------------------------------------------------------
*     input data and initialization
*-----------------------------------------------------------------------

               ieo  = 6
               ierr = 0

               mstq1(1) = 0
               mstq1(2) = ibryf(ityp,ktyp)
               mstq1(3) = ichgf(ityp,ktyp)

               msp = mstq1(3)
               msn = mstq1(2) - mstq1(3)

               if( ( msp .eq. 0 .and. msn .gt. 1 ) .or.
     &             ( msp .gt. 1 .and. msn .le. 0 ) ) then

                  write(ieo,'(/
     &                '' **** Collision is skipped by JQMD ****''/
     &                '' projectile charge = '',i4/
     &                ''            mass   = '',i4)') msp, msp+msn

                  nclst = -2
                  return

               end if

               mstq1(4) = 0
               mstq1(5) = mmas
               mstq1(6) = mchg

               if( ( mchg .eq. 0 .and. mmas-mchg .gt. 1 ) .or.
     &             ( mchg .gt. 1 .and. mmas-mchg .le. 0 ) ) then

                  write(ieo,'(/
     &                '' **** Collision is skipped by JQMD ****''/
     &                ''     target charge = '',i4/
     &                ''            mass   = '',i4)') mchg, mmas

                  nclst = -2
                  return

               end if

*-----------------------------------------------------------------------
*           bmax for QMD calculation
*           ielst, moving frame adjust
*-----------------------------------------------------------------------

               ap = mstq1(2)
               at = mstq1(5)

c              bmax = bmax0

               bmax = 1.15 * ( ap**(1./3.) + at**(1./3.) )
     &              - 0.4 + bplus

c              bmax = 1.2 * ( ap**(1./3.) + at**(1./3.) )
c    &              - 0.5 + bplus

c              bmax = 1.189 * ( ap**(1./3.) + at**(1./3.) )
c    &              - 0.9612 + bplus

*-----------------------------------------------------------------------

               mstq1(190) = imadj    ! 0:no, 1:moving frame adjust
               mstq1(191) = ijudg    ! 0:no, 1:new ielst judge

               mstq1(17)  = 1        ! ielst for judge

*-----------------------------------------------------------------------

               parq1(1)   = eein / dble( mstq1(2) ) / 1000.0d0
               parq1(2)   = -1.0
               parq1(3)   = 0.0
               parq1(4)   = bmax
               mstq1(11)  = 0

               mstq1(7)   = 1
               mstq1(8)   = iqmax
               parq1(5)   = 1.0

               mstq1(90)  = 1          ! 0:no, 1:ground energy adjust

               mstq1(9)   = imadj + 1  ! 0: lab,  1: cm,  2: nn 
               mstq1(32)  = 1          ! 0:no, 1:rel.correction

*-----------------------------------------------------------------------
*           nucleon-nucleus collisions
*-----------------------------------------------------------------------

            if( mstq1(2) .eq. 1 .or. mstq1(5) .eq. 1 ) then

               mstq1(9)   = 0
               mstq1(17)  = 3        ! ielst for judge
               mstq1(190) = 0        ! 0:no, 1:moving frame adjust

               bmax = 1.2 * ( ap**(1./3.) + at**(1./3.) )
     &              + bplus

               parq1(4)   = bmax

            end if

*-----------------------------------------------------------------------

               parq1(120) = 1.0      ! sdmemin

               icfg = 4

               call qmdint

*-----------------------------------------------------------------------
*           QMD Events
*-----------------------------------------------------------------------

               call qmdevent(ierr)

                  if( ierr .ne. 0 ) goto 999

*-----------------------------------------------------------------------
*           elastic judge and frame trnasform to lab system
*-----------------------------------------------------------------------

               call qmdjudge

                  if( kelst .eq. 1 ) goto 999

                     ifrm = 0

               do i = 1, nclst

                     px = qclust(1,i)
                     py = qclust(2,i)
                     pz = qclust(3,i)
                     et = qclust(4,i)
                     rm = qclust(5,i)

                     gamm = gammfr(ifrm)
                     beta = betafr(ifrm)

                     pz = pz * gamm - beta * gamm * et
                     et = sqrt( px**2 + py**2 + pz**2 + rm**2 )

                     qclust(3,i)  = pz
                     qclust(4,i)  = et
                     qclust(7,i)  = ( et - rm ) * 1000.

               end do

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
*     error or nothing happened
*-----------------------------------------------------------------------

  999 continue

                  nclst = -1

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine qmdevent(ierr)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 26                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to simulate one event by QMD                            *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      common /const1/ elab, rdist, bmin, bmax, ibch, ibin
      common /const2/ dt, ntmax, iprun, iprun0
      common /const3/ nfreq, nfrec, nfred
      common /vriab1/ b, llnow, ntnow
      common /swich1/ ipot, insys, irkg, icolt

      common /impact/ bval(1000), ibnum(1000), bweight(1000), bdef

      common /rannum/ iseed, iseed0, iseed1
      common /coln01/ iccoll

      data ibsf /0/
      save ibsf

*-----------------------------------------------------------------------
*        Initialization of one QMD event
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*           event number, initial time, and initial randum seed
*           total collision flag
*-----------------------------------------------------------------------

                  ierr =  0

                  ntnow  = 0
                  iseed1 = iseed
                  iccoll = 0

*-----------------------------------------------------------------------
*           impact parameter and multi run control
*-----------------------------------------------------------------------

               if( ibch .eq. 0 ) then

                     b = sqrt( max( 0.0d0, 
     &                   bmin**2 + ( bmax**2 - bmin**2 ) * rn(0) ) )

               else if( ibch .eq. 1 ) then

                     ibsf = ibsf + 1

                  if( ibsf .gt. ibin ) then

                     ibsf = 1

                  end if

                     b = bval(ibsf) - bdef / 2.0 + bdef * rn(0)

               end if

*-----------------------------------------------------------------------
*           make ground state and boost
*-----------------------------------------------------------------------

                  call ground(ierr)

                     if( ierr .ne. 0 ) return

                  call rboost(ierr)

                     if( ierr .ne. 0 ) return

*-----------------------------------------------------------------------
*        Time Evolution
*-----------------------------------------------------------------------

         do 100 nt = 1, ntmax

               ntnow = nt

*-----------------------------------------------------------------------
*           time integration
*-----------------------------------------------------------------------

               if( irkg .eq. 2 ) then

                  call rk12(dt)

               else if( irkg .eq. 4 ) then

                  call rkg4(dt)

               end if

*-----------------------------------------------------------------------
*           collision term
*-----------------------------------------------------------------------

               if( icolt .eq. 1 ) then

                  call pionem(dt)
                  call relcol
                  call pionab

               end if

*-----------------------------------------------------------------------

  100    continue

*-----------------------------------------------------------------------
*        Final pion decay and Final analysis of clusters
*-----------------------------------------------------------------------

                  call fpidecay
                  call cldist

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine qmdjudge
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 26                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to determin weight of the event, [qmdfac] and           *
*              to judge the elastic or inelastic reaction type and     *
*              to sum up the reaction cross section                    *
*                                                                      *
*                                                                      *
*        Variables: in common block /swich3/                           *
*                                                                      *
*              ielst       : input flag, jelst < ielst: elastic        *
*              jelst       : elastic or inelastic flag                 *
*                     = 0  : elastic without collision                 *
*                     = 1  : elastic with collision                    *
*                     = 2  : inelastic without collision               *
*                     = 3  : inelastic with collision                  *
*              kelst       : 1 -> elastic, 0-> inelastic               *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param01.inc'

*-----------------------------------------------------------------------

      common /const0/ idnta, idnpr, massta, masspr, mstapr, msprpr
      common /const1/ elab, rdist, bmin, bmax, ibch, ibin

      common /vriab3/ qmdfac, sdmfac

      common /clustf/ nclst, iclust(nnn)
      common /clustg/ jclust(0:7,nnn),  qclust(0:12,nnn)

      common /coln01/ iccoll
      common /swich3/ ielst, jelst, kelst
      common /sdmcut/ sdmemin

      common /impact/ bval(1000), ibnum(1000), bweight(1000), bdef

      common /summ03/ rcross(4), ireac(0:10)
      common /input1/ mstq1(mxpa1), parq1(mxpa1)

*-----------------------------------------------------------------------

                     iselc = mstq1(191)

*-----------------------------------------------------------------------
*        judgement of elastic or inelastic reaction
*-----------------------------------------------------------------------

                     iels = 1

            if( nclst .eq. 2 ) then

               if( ( jclust(1,1) .eq. mstapr .and.
     &               jclust(2,1) .eq. massta - mstapr .and.
     &               jclust(1,2) .eq. msprpr .and.
     &               jclust(2,2) .eq. masspr - msprpr )     .or.
     &             ( jclust(1,2) .eq. mstapr .and.
     &               jclust(2,2) .eq. massta - mstapr .and.
     &               jclust(1,1) .eq. msprpr .and.
     &               jclust(2,1) .eq. masspr - msprpr ) )   then

               if( iselc .eq. 0 ) then

                     if( qclust(6,1) .lt. sdmemin .and.
     &                   qclust(6,2) .lt. sdmemin  ) then

                        iels = 0

                     end if

               else

                  if( massta .gt. 1 .and. masspr .gt. 1 .and.
     &                iccoll .eq. 0 ) then

                        remin1 = 0.3d0 * ( jclust(1,1) + jclust(2,1) )
                        remin2 = 0.3d0 * ( jclust(1,2) + jclust(2,2) )

                     if( qclust(6,1) .lt. remin1 .and.
     &                   qclust(6,2) .lt. remin2  ) then

                        iels = 0

                     end if

                  else

                     if( qclust(6,1) .lt. sdmemin .and.
     &                   qclust(6,2) .lt. sdmemin  ) then

                        iels = 0

                     end if

                  end if

               end if

               end if

            end if

*-----------------------------------------------------------------------

            if( iccoll .eq. 0 ) then

               if( iels .eq. 0 ) then

                  jelst = 0

               else

                  jelst = 2

               end if

            else

               if( iels .eq. 0 ) then

                  jelst = 1

               else

                  jelst = 3

               end if

            end if

*-----------------------------------------------------------------------
*        elastic flag
*-----------------------------------------------------------------------

            if( jelst .ge. ielst ) then

               kelst = 0

            else

               kelst = 1

            end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine trfram(px,py,pz,et,rm,ifrm)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to determine energy and momentum by Lorentz transform   *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              px, py, pz     : momentum                               *
*              et             : energy                                 *
*              rm             : rest mass                              *
*              ifrm           : 0-> to lab                             *
*                               1-> to cm                              *
*                               2-> to n-n                             *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      common /framtr/ betafr(0:2), gammfr(0:2)

*-----------------------------------------------------------------------

         if( ifrm .lt. 0 .or. ifrm .gt. 2 ) then

            write(*,*) ' **** Error at [trfram], unrecognized frame '
            stop 999

         end if

*-----------------------------------------------------------------------

            gamm = gammfr(ifrm)
            beta = betafr(ifrm)

*-----------------------------------------------------------------------

            pz = pz * gamm - beta * gamm * et
            et = sqrt( px**2 + py**2 + pz**2 + rm**2 )

*-----------------------------------------------------------------------

      return
      end


