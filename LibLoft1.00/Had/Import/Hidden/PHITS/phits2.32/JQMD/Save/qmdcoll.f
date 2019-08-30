************************************************************************
*                                                                      *
*        PART 5: Collision Term                                        *
*                                                                      *
*                                                                      *
*   List of subprograms in rough order of relevance with main purpose  *
*      ( s = subroutine, f = function, b = block data, e = entry )     *
*                                                                      *
*                                                                      *
*  s  relcol    to calculate the kinematics between two particles      *
*  s  crosw     to determine collisoin channel and final state         *
*  s  pionem    to calculate the decay of delta or N*                  *
*  s  pionab    to calculate pion absorption                           *
*  s  fpidecay  to calculate final decay of the resonances             *
*  s  resmas    to calcurate the mass and width of Delta and N*        *
*  f  pideno    to calculate the correction factor                     *
*               of inverse cross section                               *
*  f  soo       to calculate delta, N* cross section                   *
*                                                                      *
************************************************************************
*                                                                      *
************************************************************************
*                                                                      *
      subroutine relcol
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to calculate the kinematics in a collision              *
*              between two particles                                   *
*                                                                      *
*                                                                      *
*        Comments : counting variables                                 *
*                                                                      *
*              lcoll( 1) : all collisions                              *
*                                                                      *
*              lcoll( 2) : N  - DR  elastic                            *
*              lcoll( 3) : P  - N   elastic                            *
*              lcoll( 4) : DR - DR  elastic                            *
*                                                                      *
*              lcoll( 5) : Puli-blocked collisions                     *
*              lcoll( 6) : energetically forbidden collision           *
*                                                                      *
*              lcoll( 7) : elastic collision                           *
*                                                                      *
*              lcoll( 8) : N + N -> N + D                              *
*              lcoll( 9) : N + D -> N + N                              *
*              lcoll(10) : N + N -> N + R                              *
*              lcoll(11) : N + R -> N + N                              *
*              lcoll(12) : N + N -> D + D                              *
*              lcoll(13) : D + D -> N + N                              *
*              lcoll(14) : N + D -> D + D                              *
*              lcoll(15) : N + R -> D + R                              *
*              lcoll(16) : N + D -> R + D                              *
*              lcoll(17) : N + R -> R + R                              *
*              lcoll(18) : D + D -> N + D                              *
*              lcoll(19) : D + R -> N + R                              *
*              lcoll(20) : R + D -> N + D                              *
*              lcoll(21) : R + R -> N + R                              *
*                                                                      *
*              lcoll(22) : D + N -> pi                                 *
*              lcoll(23) : R + N -> pi                                 *
*              lcoll(24) : R + D -> pi                                 *
*                                                                      *
*              lcoll(25) : N + pi -> D                                 *
*              lcoll(26) : N + pi -> R                                 *
*              lcoll(27) : D + pi -> R                                 *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param01.inc'
      include 'param02.inc'

*-----------------------------------------------------------------------

      parameter     ( deltar = 4.0 )
      parameter     ( sig0 =  55.0, bcmax0 = 1.323142 )
      parameter     ( sig1 = 200.0, bcmax1 = 2.523 )

*-----------------------------------------------------------------------
*
*     deltar        : maximum spatial distance for which a collision
*                     still can occur only for below 200 MeV
*
*     bcmax         : maximum impact parameter
*                     bcmax0 for NN, bcmax1 for the other
*
*     sig           : corresponding cross section
*                     sig0 for NN, sig1 for the other
*
*-----------------------------------------------------------------------

      common /const1/ elab, rdist, bmin, bmax, ibch, ibin
      common /const2/ dt, ntmax, iprun, iprun0

      common /vriab0/ massal, massba, mmeson

      common /coodrp/ r(5,nnn),  p(5,nnn)
      common /cood2r/ rr2(nnn,nnn), rbij(nnn,nnn), pp2(nnn,nnn)
      common /coodid/ ichg(nnn), inuc(nnn), ibry(nnn), inds(nnn),
     &                inun(nnn), iavd(nnn), ihis(nnn)

      common /swich0/ icoul, irelcr, iavoid
      common /swich1/ ipot, insys, irkg, icolt

      common /coln00/ lcoll(30)
      common /coln01/ iccoll

*-----------------------------------------------------------------------

      dimension       pcm(3), beta(3)

*-----------------------------------------------------------------------
*     initialization of counting variables
*-----------------------------------------------------------------------

            do i = 1, 21

               lcoll(i) = 0

            end do

*-----------------------------------------------------------------------

      do 800 i1 = 2, massba

*-----------------------------------------------------------------------

                  px1  = p(1,i1)
                  py1  = p(2,i1)
                  pz1  = p(3,i1)

                  e1   = p(4,i1)
                  em1  = p(5,i1)

                  iz1  = ichg(i1)
                  id1  = inds(i1)
                  inc1 = inuc(i1)

*-----------------------------------------------------------------------

         do 600 i2 = 1, i1 - 1

*-----------------------------------------------------------------------
*           avoid first collisions within the same nucleus
*-----------------------------------------------------------------------

               if( iavd(i1) * iavd(i2) .eq. iavoid ) goto 600

*-----------------------------------------------------------------------
*           avoid second collisions for the same pairs
*-----------------------------------------------------------------------

               if( inun(i1) .eq. i2 .and.
     &             inun(i2) .eq. i1 ) goto 600

*-----------------------------------------------------------------------
*           the following prescription (deltar) is not covariant,
*           but useful to save cpu time especially below 200 MeV/u
*-----------------------------------------------------------------------

               if( elab .le. 5.0 .and.
     &             rr2(i1,i2) .gt. deltar**2 ) goto 600

*-----------------------------------------------------------------------
*           now particles are close enough to each other
*-----------------------------------------------------------------------

                  px2  = p(1,i2)
                  py2  = p(2,i2)
                  pz2  = p(3,i2)

                  e2   = p(4,i2)
                  em2  = p(5,i2)

                  iz2  = ichg(i2)
                  id2  = inds(i2)
                  inc2 = inuc(i2)

                  s    = (  e1 +  e2 )**2
     &                 - ( px1 + px2 )**2
     &                 - ( py1 + py2 )**2
     &                 - ( pz1 + pz2 )**2

                  srt  = sqrt(s)

*-----------------------------------------------------------------------
*           low energy cutoff and max cross section
*           ( 1.8966  =  em1 + em2 + 0.02 GeV ; nucleon case )
*-----------------------------------------------------------------------

               if( em1 .lt. 0.94 .and. em2 .lt. 0.94 ) then

                  cutoff =  em1 + em2 + 0.02

                  bcmax  = bcmax0
                  sig    = sig0

               else

                  cutoff =  em1 + em2

                  bcmax  = bcmax1
                  sig    = sig1

               end if

               if ( srt .lt. cutoff ) goto 600

*-----------------------------------------------------------------------
*           are their impact parameter small enough ?
*-----------------------------------------------------------------------

                  dx   = r(1,i1) - r(1,i2)
                  dy   = r(2,i1) - r(2,i2)
                  dz   = r(3,i1) - r(3,i2)
                  rsq  = dx**2 + dy**2 + dz**2

                  p12  = e1 * e2 - px1 * px2 - py1 * py2 - pz1 * pz2
                  p1dr = px1 * dx + py1 * dy + pz1 * dz
                  p2dr = px2 * dx + py2 * dy + pz2 * dz
                  a12  = 1.0 - ( em1 * em2 / p12 ) ** 2
                  b12  = p1dr / em1 - p2dr * em1 / p12
                  c12  = rsq + ( p1dr / em1 )**2
                  brel = sqrt( abs(c12 - b12**2/a12) )

               if( brel .gt. bcmax ) goto 600

*-----------------------------------------------------------------------
*           average time-shift of the collision in the fixed frame
*           will particles get closest point in this time interval ?
*-----------------------------------------------------------------------

                  b21    =   - p2dr / em2 + p1dr * em2 / p12
                  t1     = (   p1dr / em1 - b12 / a12 ) * e1 / em1
                  t2     = ( - p2dr / em2 - b21 / a12 ) * e2 / em2

               if( abs( t1 + t2 ) .gt. dt ) goto 600

*-----------------------------------------------------------------------
*           lorentz-transformation in i1-i2-c.m. system
*-----------------------------------------------------------------------

                  etot12  =    e1 +  e2
                  beta(1) = ( px1 + px2 ) / etot12
                  beta(2) = ( py1 + py2 ) / etot12
                  beta(3) = ( pz1 + pz2 ) / etot12
                  betasq  = beta(1)**2 + beta(2)**2 + beta(3)**2
                  gamma   = 1.0 / sqrt( 1.0 - betasq )

*-----------------------------------------------------------------------
*           transformation of momenta ( p1c(i) = - p2c(i) )
*-----------------------------------------------------------------------

                  p1beta = px1 * beta(1) + py1 * beta(2) + pz1 * beta(3)
                  transf = gamma
     &                   * ( gamma * p1beta / ( gamma + 1.0 ) - e1 )
                  pcm(1) = beta(1) * transf + px1
                  pcm(2) = beta(2) * transf + py1
                  pcm(3) = beta(3) * transf + pz1
                  pcm2   = pcm(1)**2 + pcm(2)**2 + pcm(3)**2
                  prcm   = sqrt( pcm2 )

               if(prcm .le. 0.00001) goto 600

*-----------------------------------------------------------------------
*           calculate cross section and collide two particles
*-----------------------------------------------------------------------


               call crosw(ichannel,sig,cutoff,
     &                    pcm,prcm,srt,beta,gamma,i1,i2)


*-----------------------------------------------------------------------
*              ichannel : channel of this collision
*-----------------------------------------------------------------------
*
*                 =  0 ; nothing has happened
*                 =  1 ; elastic N-N collision
*                 =  2 ; N + N -> N + D
*                 =  3 ; N + D -> N + N
*                 =  4 ; N + N -> N + R
*                 =  5 ; N + R -> N + N
*                 =  6 ; N + N -> D + D
*                 =  7 ; D + D -> N + N
*                 =  8 ; N + D -> D + D
*                 =  9 ; D + D -> N + D
*                 = 10 ; N + R -> D + R
*                 = 11 ; D + R -> N + R
*                 = 12 ; N + D -> R + D
*                 = 13 ; R + D -> N + D
*                 = 14 ; N + R -> R + R
*                 = 15 ; R + R -> N + R
*                 = 99 ; energetically forbidden collision
*
*-----------------------------------------------------------------------

               if( ichannel .eq. 0 ) goto 600

                  lcoll(1) = lcoll(1) +1

*-----------------------------------------------------------------------
*           a collision has taken place
*-----------------------------------------------------------------------

                  ntag  = 0

               if( ichannel .eq. 99 ) then

                  ntag  = 1
                  lcoll(6) = lcoll(6) +1

               end if

*-----------------------------------------------------------------------
*           check of pauli-blocking
*-----------------------------------------------------------------------

               if( ntag .eq. 0 ) call pauli( i1, ntag, phase )
               if( ntag .eq. 0 ) call pauli( i2, ntag, phase )

*-----------------------------------------------------------------------
*              pauli blocked case, reset variables
*-----------------------------------------------------------------------

            if( ntag .eq. 1 ) then

                  lcoll(5) = lcoll(5) + 1

                  p(1,i1)  = px1
                  p(2,i1)  = py1
                  p(3,i1)  = pz1

                  p(5,i1)  = em1
                  p(4,i1)  = e1
                  inds(i1) = id1
                  ichg(i1) = iz1
                  inuc(i1) = inc1

                  p(1,i2)  = px2
                  p(2,i2)  = py2
                  p(3,i2)  = pz2

                  p(5,i2)  = em2
                  p(4,i2)  = e2
                  inds(i2) = id2
                  ichg(i2) = iz2
                  inuc(i2) = inc2

                  call caldis2(i1,i2)

*-----------------------------------------------------------------------
*              collision is realy happened
*-----------------------------------------------------------------------

            else

                  ihis1 = 1
                  ihis2 = 1

               if( ihis(i1) .lt. 0 ) ihis1 = -1
               if( ihis(i2) .lt. 0 ) ihis2 = -1

               if( ichannel .gt. 1 ) then

                  ihis1 = -1
                  ihis2 = -1

               end if

                  ihis0    = abs(ihis(i1)) + abs(ihis(i2)) + 1
                  ihis(i1) = ihis1 * ihis0
                  ihis(i2) = ihis2 * ihis0

                  inun(i1) = i2
                  inun(i2) = i1

                  px1  = p(1,i1)
                  py1  = p(2,i1)
                  pz1  = p(3,i1)

                  em1  = p(5,i1)
                  e1   = p(4,i1)

                  iz1  = ichg(i1)
                  id1  = inds(i1)
                  inc1 = inuc(i1)

                  iavd(i1) = iavd(i1) + iavd(i1) / abs(iavd(i1))
                  iavd(i2) = iavd(i2) + iavd(i2) / abs(iavd(i2))

*-----------------------------------------------------------------------
*                 counting of the collision
*-----------------------------------------------------------------------

                     iccoll = iccoll + 1

               if( ichannel .eq. 1 ) then

                     lcoll(7) = lcoll(7) + 1

                  if( ( em1 .lt. 0.94 .and. em2 .lt. 0.94 ) .and.
     &                  ichg(i1) + ichg(i2) .eq. 1 )
     &               lcoll(3) = lcoll(3) + 1

                  if( ( em1 .gt. 0.94 .and. em2 .lt. 0.94 ) .or.
     &                ( em1 .lt. 0.94 .and. em2 .gt. 0.94 ) )
     &               lcoll(2) = lcoll(2) + 1

                  if( em1 .gt. 0.94 .and. em2 .gt. 0.94 )
     &               lcoll(4) = lcoll(4) + 1

               end if

                  if( ichannel .eq.  2 ) lcoll( 8) = lcoll( 8) + 1
                  if( ichannel .eq.  3 ) lcoll( 9) = lcoll( 9) + 1
                  if( ichannel .eq.  4 ) lcoll(10) = lcoll(10) + 1
                  if( ichannel .eq.  5 ) lcoll(11) = lcoll(11) + 1
                  if( ichannel .eq.  6 ) lcoll(12) = lcoll(12) + 1
                  if( ichannel .eq.  7 ) lcoll(13) = lcoll(13) + 1
                  if( ichannel .eq.  8 ) lcoll(14) = lcoll(14) + 1
                  if( ichannel .eq.  9 ) lcoll(18) = lcoll(18) + 1
                  if( ichannel .eq. 10 ) lcoll(15) = lcoll(15) + 1
                  if( ichannel .eq. 11 ) lcoll(19) = lcoll(19) + 1
                  if( ichannel .eq. 12 ) lcoll(16) = lcoll(16) + 1
                  if( ichannel .eq. 13 ) lcoll(20) = lcoll(20) + 1
                  if( ichannel .eq. 14 ) lcoll(17) = lcoll(17) + 1
                  if( ichannel .eq. 15 ) lcoll(21) = lcoll(21) + 1

*-----------------------------------------------------------------------

            end if

*-----------------------------------------------------------------------

  600     continue

  800   continue

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine crosw(ichannel,sig,cutoff,
     &                 pcm,prcm,srt,beta,gamma,i1,i2)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to determine collisoin channel and final state          *
*                                                                      *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              ichannel    : channel information back                  *
*              sig         : max cross section at cutoff energy        *
*              cutoff      : cutoff energy                             *
*                            em1 + em2 + 0.02 GeV for N + N            *
*                            em1 + em2            for the others       *
*              pcm(3)      : momentum coordinates of one particle      *
*                            in cm frame                               *
*              prcm        : sqrt(pcm)                                 *
*              srt         : sqrt of s                                 *
*              beta(3)     : beta of cm frame to reference frame       *
*              gamma       : gamma of above beta(3)                    *
*              i1,i2       : identificator of particle 1 and 2         *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param01.inc'
      include 'param02.inc'

*-----------------------------------------------------------------------

      parameter     ( epse = 0.0001 )

*-----------------------------------------------------------------------

      common /poten0/ aaa, bbb, rpot, esymm

      common /coodrp/ r(5,nnn),  p(5,nnn)
      common /coodid/ ichg(nnn), inuc(nnn), ibry(nnn), inds(nnn),
     &                inun(nnn), iavd(nnn), ihis(nnn)

      common /swich1/ ipot, insys, irkg, icolt

*-----------------------------------------------------------------------

      dimension       pcm(3), beta(3)

*-----------------------------------------------------------------------
*           choice of channels
*-----------------------------------------------------------------------

*              (1) ratio of NN -> DD  to  NN -> NR

                        fcdd = 0.0

                        fnnr = 1.5 * ( 1.0 - fcdd )
                        fndd = 1.5 *   fcdd


*              (2) ND -> DD  (2-2), (4-2)

                        inddd = 1

*              (3) ND -> RD  (2-3), (5-2)

                        indrd = 1

*              (4) NR -> DR  (3-2), (5-1)

                        inrdr = 1

*              (5) NR -> RR  (3-3), (6-1)

                        inrrr = 1

*-----------------------------------------------------------------------

                  em1 = p(5,i1)
                  id1 = inds(i1)
                  iz1 = ichg(i1)

                  em2 = p(5,i2)
                  id2 = inds(i2)
                  iz2 = ichg(i2)

*-----------------------------------------------------------------------

                  rm2 = rmass**2
                  pr  = prcm
                  c2  = pcm(3) / pr

                  csrt = srt - cutoff

                  pri  = prcm
                  prf  = sqrt( 0.25 * srt**2 - rm2 )

                  asrt = srt - em1 - em2
                  pra  = prcm

*-----------------------------------------------------------------------
*              random number
*-----------------------------------------------------------------------

                  x1  = rn(0)
                  x2  = rn(0)

*-----------------------------------------------------------------------
*           elastic cross section;       modified Cugnon
*-----------------------------------------------------------------------

            if( iz1 .eq. iz2 ) then

               if( csrt .lt. 0.4286 ) then

                  sigel = 35.0 / ( 1. + csrt * 100.0 )  +  20.0

               else

                  sigel = ( - atan( ( csrt - 0.4286 ) * 1.5 - 0.8 )
     &                  *   2. / pi + 1.0 ) * 9.65 + 7.0

               end if

            else

               if( csrt .lt. 0.4286 ) then

                  sigel = 28.0 / ( 1. + csrt * 100.0 )  +  27.0

               else

                  sigel = ( - atan( ( csrt - 0.4286 ) * 1.5 - 0.8 )
     &                  *   2. / pi + 1.0 ) * 12.34 + 10.0

               end if

            end if

*-----------------------------------------------------------------------

                  sigtt = sigel

                  ichannel = 1


*======================================================================*

*     Start of Inelastic Collisions,   at least one pion production

      if( x1 .gt. sigel / sig ) then

*======================================================================*


            ichannel = 0

            if( srt .le. 2.0 * rmass + pmass ) return


*======================================================================*

*        [1] Nucleon - Nucleon  Inelastic Collsions

*                        (1-1) N + N -> N + D
*                        (1-2) N + N -> N + R
*                        (1-3) N + N -> D + D

         if( id1 .eq. 1 .and. id2 .eq. 1 ) then

*======================================================================*

                  signd = 0.0
                  signr = 0.0
                  sigdd = 0.0

*-----------------------------------------------------------------------

*           (1-1) N + N -> N + D

*-----------------------------------------------------------------------

                  z11   = soo(1,srt)
                  z10   = soo(2,srt)

               if( iz1 .eq. iz2 ) then

                  signd = 2.0  * z11 + z10

               else

                  signd = z11 + 0.5 * z10

               end if

*-----------------------------------------------------------------------

            if( x1 .lt. ( sigel + signd ) / sig ) then

                     ichannel = 2

                     idd1 = 2
                     idd2 = 1

                     izz1 = iz1
                     izz2 = iz2

               if( x2 .gt. 0.5 ) then

                     izz1 = iz2
                     izz2 = iz1

               end if

                     dem1 = rmass + pmass
                     dem2 = rmass

                  call resmas(1,idd1,idd2,srt,dem1,dem2,gam0,fqrq)


               if( ( iz1 + iz2 .ne. 1 ) .and.
     &             ( x1 .lt. ( sigel + 0.5 * z11 + z10 ) / sig ) ) then

                     izz1 = 3 * iz1 - 1
                     izz2 = iz2 - 2 * iz1 + 1

               end if

               goto 100

            end if

*-----------------------------------------------------------------------

*           (1-2) N + N -> N + R

*-----------------------------------------------------------------------

                  z01 = soo(3,srt)

                  signr = fnnr * z01

*-----------------------------------------------------------------------

            if( x1 .lt. ( sigel + signd + signr ) / sig ) then

                     ichannel = 4

                     idd1 = 3
                     idd2 = 1

                     izz1 = iz1
                     izz2 = iz2

               if( x2 .gt. 0.5 ) then

                     izz1 = iz2
                     izz2 = iz1

               end if

                     dem1 = rmass + pmass
                     dem2 = rmass

                  call resmas(1,idd1,idd2,srt,dem1,dem2,gam0,fqrq)


               goto 100

            end if


            if( srt .le. 2.0 * rmass + 2.0 * pmass ) return


*-----------------------------------------------------------------------

*           (1-3) N + N -> D + D

*-----------------------------------------------------------------------

                  sigdd = fndd * z01

*-----------------------------------------------------------------------

            if( x1 .lt. ( sigel + signd + signr + sigdd ) / sig ) then

                     ichannel = 6

                     idd1 = 2
                     idd2 = 2

                     izz1 = iz1
                     izz2 = iz2

                     dem1 = rmass + pmass
                     dem2 = rmass + pmass

                  call resmas(1,idd1,idd1,srt,dem1,dem2,gam0,fqrq)


               if( iz1 + iz2 .eq. 1) then

                  if( x1 .lt. ( sigel + signd + signr
     &                          + 7./10. * sigdd ) / sig ) then

                     izz1 =   2
                     izz2 = - 1

                  end if

               else

                  if( x1 .lt. ( sigel + signd + signr
     &                          + 3./7. * sigdd ) / sig ) then

                     izz1 = ( iz1 + iz2 ) / 2 + 1
                     izz2 = ( iz1 + iz2 ) / 2 - 1

                  end if

               end if

               goto 100

            end if

*-----------------------------------------------------------------------

         return

         end if

*-----------------------------------------------------------------------

*======================================================================*

*        [2] Nucleon - Delta ;  Inelastic Collisions

*                        (2-1) N + D -> N + N
*                        (2-2) N + D -> D + D
*                        (2-3) N + D -> R + D

         if( id1 + id2 .eq. 3 ) then

*======================================================================*

                  sigdn = 0.0
                  sigdd = 0.0
                  sigdr = 0.0

*-----------------------------------------------------------------------

*           (2-1) N + D -> N + N

            if( ( iz1 + iz2 + 3 ) / 3 .eq. 1 ) then

*-----------------------------------------------------------------------

                  signd = 0.0
                  z11   = soo(1,srt)
                  z10   = soo(2,srt)

               if( iz1 + iz2 .eq. 1 ) then

                  signd = 0.5 * z11 + 0.25 * z10
                  sigiv = 2.0

               end if

               if( iz1 .eq. iz2 ) then

                  signd = 1.5 * z11
                  sigiv = 4.0

               end if

               if( ( iz1 .ne. iz2 ) .and. ( iz1 + iz2 .ne. 1 ) ) then

                  signd = 0.5 * z11 + z10
                  sigiv = 4.0

               end if

                  pidn  = pideno(2,rmass+pmass,rmass,srt)

                  sigdn = signd * prf**2 / pri**2 / pidn / sigiv


*-----------------------------------------------------------------------

               if( x1 .lt. ( sigel + sigdn ) / sig ) then

                     ichannel = 3

                     idd1 = 1
                     idd2 = 1

                     izz1 = iz1
                     izz2 = iz2

                     dem1 = rmass
                     dem2 = rmass

                  if( ( iz1 + iz2 .eq. 0 ) .or.
     &                ( iz1 + iz2 .eq. 2 ) ) then

                     izz1 = ( iz1 + iz2 ) / 2
                     izz2 = ( iz1 + iz2 ) / 2

                  end if

                  goto 100

               end if

            end if

            if( srt .lt. 2.0 * rmass + 2.0 * pmass ) return

*-----------------------------------------------------------------------

*           (2-2) N + D -> D + D

            if( inddd .ne. 0 ) then

*-----------------------------------------------------------------------

                  srtd  = srt - ( em1 + em2 - 2.0 * rmass )

                  z11   = soo(1,srtd)
                  z10   = soo(2,srtd)

               if( ( iz1 + iz2 .eq. 3 ).or.( iz1 + iz2 .eq. -1 ) ) then

                  sigdd = 1.5 * z11

               else

                  sigdd = 2.0 * z11 + z10

               end if

*-----------------------------------------------------------------------

               if( x1 .lt. ( sigel + sigdn + sigdd ) / sig ) then

                     ichannel = 8

                     idd1 = 2
                     idd2 = 2

                  if( id1 .eq. 1 ) then

                     izz1 = iz1
                     izz2 = iz2

                  else

                     izz1 = iz2
                     izz2 = iz1

                  end if

                     dem1 = rmass + pmass
                     dem2 = rmass + pmass

                  call resmas(1,idd1,idd1,srt,dem1,dem2,gam0,fqrq)


               if( x1 .gt. ( sigel + sigdn + 1.5 * z11 ) / sig ) then

                  if( iz1 + iz2 .eq. 1 ) then

                     izz1 =  3 * iz1 - 1
                     izz2 =  iz2 - 2 * iz1 +1

                  else if( iz1 .eq. iz2 ) then

                     izz1 = ( iz1 + iz2 ) / 2 + 1
                     izz2 = ( iz1 + iz2 ) / 2 - 1

                  else

                     izz1 = ( iz1 + iz2 ) / 2
                     izz2 = ( iz1 + iz2 ) / 2

                  end if

               end if

                  goto 100

               end if

            end if

*-----------------------------------------------------------------------

*           (2-3) N + D -> R + D

            if( indrd .ne. 0 ) then

*-----------------------------------------------------------------------

                  z01   = soo(3,srtd)

                  sigdr = 1.5 * z01

*-----------------------------------------------------------------------

               if( x1 .lt.
     &           ( sigel + sigdn + sigdd + sigdr ) / sig ) then

                     ichannel = 12

                     idd1 = 3
                     idd2 = 2

                  if( id1 .eq. 1 ) then

                     izz1 = iz1
                     izz2 = iz2

                  else

                     izz1 = iz2
                     izz2 = iz1

                  end if

                     dem1 = rmass + pmass
                     dem2 = rmass + pmass

                  call resmas(1,idd1,idd2,srt,dem1,dem2,gam0,fqrq)

                  goto 100

               end if

            end if

*-----------------------------------------------------------------------

            return

         end if

*-----------------------------------------------------------------------

*======================================================================*

*        [3] Nucleon - Nstar  ;  Inelastic Collisions

*                        (3-1) N + R -> N + N
*                        (3-2) N + R -> D + R
*                        (3-3) N + R -> R + R

         if( ( id1 + id2 .eq. 4 ) .and. ( id1 .ne. 2 ) ) then

*======================================================================*

                  sigrn = 0.0
                  sigdr = 0.0
                  sigrr = 0.0

*-----------------------------------------------------------------------

*           (3-1) N + R -> N + N

*-----------------------------------------------------------------------

                  z01   = soo(3,srt)

                  signr = fnnr * z01

                  pidn  = pideno(3,rmass+pmass,rmass,srt)

               if( iz1 .eq. iz2 ) then

                  sigiv = 2.0

               else

                  sigiv = 1.0

               end if


                  sigrn = signr * prf**2 / pri**2 / pidn / sigiv

*-----------------------------------------------------------------------

               if( x1 .lt. ( sigel + sigrn ) / sig ) then

                     ichannel = 5

                     idd1 = 1
                     idd2 = 1

                     izz1 = iz1
                     izz2 = iz2

                     dem1 = rmass
                     dem2 = rmass

                  goto 100

               end if

               if( srt .le. 2.0 * rmass + 2.0 * pmass ) return

*-----------------------------------------------------------------------

*           (3-2) N + R -> D + R

            if( inrdr .ne. 0 ) then

*-----------------------------------------------------------------------

                  srtd  = srt - ( em1 + em2 - 2.0 * rmass )

                  z11   = soo(1,srtd)
                  z10   = soo(2,srtd)

                  sigdr = 2.0 * z11 + z10

*-----------------------------------------------------------------------

               if( x1 .lt. ( sigel + sigrn + sigdr ) / sig ) then

                     ichannel = 10

                     idd1 = 2
                     idd2 = 3

                  if( id1 .eq. 1 ) then

                     izz1 = iz1
                     izz2 = iz2

                  else

                     izz1 = iz2
                     izz2 = iz1

                  end if

                     dem1 = rmass + pmass
                     dem2 = rmass + pmass

                  call resmas(1,idd2,idd1,srt,dem2,dem1,gam0,fqrq)


               if( x1 .gt. ( sigel + sigrn + 1.5 * z11 ) / sig ) then

                  if( iz1 + iz2 .eq. 1 ) then

                     izzd = izz1
                     izz1 = izz2
                     izz2 = izzd

                  else 

                     izz1 = 3 * ( iz1 + iz2 ) / 2 -1
                     izz2 = iz1 + iz2 - izz1

                  end if

               end if

               goto 100

               end if

            end if

*-----------------------------------------------------------------------

*           (3-3) N + R -> R + R

            if( inrrr .ne. 0 ) then
*-----------------------------------------------------------------------

                  z01   = soo(3,srtd)

                  sigrr = 1.5 * z01

*-----------------------------------------------------------------------

               if( x1 .lt.
     &           ( sigel + sigrn + sigdr + sigrr ) / sig ) then

                     ichannel = 14

                     idd1 = 3
                     idd2 = 3

                  if( id1 .eq. 1 ) then

                     izz1 = iz1
                     izz2 = iz2

                  else

                     izz1 = iz2
                     izz2 = iz1

                  end if

                     dem1 = rmass + pmass
                     dem2 = rmass + pmass

                  call resmas(1,idd1,idd2,srt,dem1,dem2,gam0,fqrq)

                  goto 100

               end if

            end if

*-----------------------------------------------------------------------

            return

         end if

*-----------------------------------------------------------------------

*======================================================================*

*        [4] Delta - Delta  ;  Inelastic Collisions

*                        (4-1) D + D -> N + N
*                        (4-2) D + D -> D + N

         if( ( id1 .eq. 2 ) .and. ( id2 .eq. 2 ) ) then

*======================================================================*

                  sigdnn = 0.0
                  sigdnd = 0.0

*-----------------------------------------------------------------------

*           (4-1) D + D -> N + N

            if( ( iz1 + iz2 + 3 ) / 3 .eq. 1 .and.
     &            fndd .gt. 0.0 ) then

*-----------------------------------------------------------------------

                  z01    = soo(3,srt)
                  sigdd  = fndd * z01

               if( ( iz1 + iz2 .eq. 1 ) .or. ( iz1 .eq. iz2 ) ) then

                  sigiv = 4.0

               else 

                  sigiv = 8.0

               end if

                  sigdnn = sigdd * prf**2 / pri**2 / sigiv

*-----------------------------------------------------------------------

               if( x1 .lt. ( sigel + sigdnn ) / sig ) then

                     ichannel = 7

                     idd1 = 1
                     idd2 = 1

                     dem1 = rmass
                     dem2 = rmass

                  if( iz1 + iz2 .eq. 1 ) then

                     izz1 = 1
                     izz2 = 0

                  else

                     izz1 = ( iz1 + iz2 ) / 2
                     izz2 = ( iz1 + iz2 ) / 2

                  end if

                  goto 100

               end if

            end if

*-----------------------------------------------------------------------

*           (4-2) D + D -> D + N

            if( inddd .ne. 0 ) then

*-----------------------------------------------------------------------


               if( ( iz1 + iz2 .ge. -1 ) .and.
     &             ( iz1 + iz2 .le. 3 ) ) then

                     dem1 = rmass

                  if( x2 .gt. 0.5 ) then

                     dem2 = em2

                  else

                     dem2 = em1

                  end if

                     srtd  = srt - ( dem2 - rmass )

                     z11   = soo(1,srtd)
                     z10   = soo(2,srtd)

                     prfd  = sqrt( ( srt**2 - dem1**2 - dem2**2 )**2
     &                    - 4.0 * ( dem1 * dem2 )**2 ) / ( 2.0 * srt )


                  if( ( ( iz1 + iz2 .eq. 3 ) .or.
     &                  ( iz1 + iz2 .eq. -1 ) ) .or.
     &                ( ( iz1 + iz2 .eq. 1 ) .and.
     &                  ( iz1 .eq. 1 .or. iz2 .eq. 1 ) ) ) then

                     signdd = 3.0 * z11
                     icc    = 1

                  else if( iz1 + iz2 .eq. 1 ) then

                     signdd = z11 + 2.0 * z10
                     icc    = 2

                  else

                     signdd = 2.0 * z11 + z10
                     icc    = 3

                  end if

                     pidn   = pideno(2,rmass+pmass,dem2,srt)

                  if( iz1 .eq. iz2 ) then

                     sigiv = 1.0

                  else

                     sigiv = 2.0

                  end if
*
                     sigdnd = signdd * prfd**2 / pri**2 / pidn /sigiv

*-----------------------------------------------------------------------

                  if( x1 .lt. ( sigel + sigdnn + sigdnd ) / sig ) then

                        ichannel = 9

                        idd1 = 1
                        idd2 = 2

                        izz1 =  iz1
                        izz2 =  iz2


                     if( ( icc .eq. 2 ) .or.
     &                 ( ( icc .eq. 3 ) .and. 
     &                   ( x1 .gt.
     &                   ( sigel + sigdnn + 1.5 * z11 ) / sig ) ) ) then

                        if( iz1 + iz2 .eq. 1 ) then

                           if( x2 .gt. 0.5 ) then

                              izz1 =  1
                              izz2 =  0

                           else

                              izz1 =  0
                              izz2 =  1

                           end if

                        else if( iz1 .eq. iz2 ) then

                              izz1 = ( iz1 + iz2 ) / 2 + 1
                              izz2 = ( iz1 + iz2 ) / 2 - 1

                        else

                              izz1 = ( iz1 + iz2 ) / 2
                              izz2 = ( iz1 + iz2 ) / 2

                        end if

                     end if

                     if( ( izz1 + 2 ) / 2 .ne. 1 ) then

                           izzd = izz1
                           izz1 = izz2
                           izz2 = izzd

                     end if

                        goto 100

                  end if

               end if

            end if

*-----------------------------------------------------------------------

            return

         end if

*-----------------------------------------------------------------------

*======================================================================*

*        [5] Delta - Nstar  ;  Inelastic Collisions

*                        (5-1) D + R -> N + R
*                        (5-2) R + D -> N + D

         if( ( id1 .eq. 2 ) .and. ( id2 .eq. 3 ) .or.
     &       ( id1 .eq. 3 ) .and. ( id2 .eq. 2 ) ) then

*======================================================================*

                  sigdnr = 0.0
                  sigrnd = 0.0

*-----------------------------------------------------------------------

*           (5-1) D + R -> N + R

            if( inrdr .ne. 0 ) then

*-----------------------------------------------------------------------


               if( ( iz1 + iz2 + 3 ) / 3 .eq. 1 ) then

                     dem1 = rmass

                  if( id1 .eq. 2 ) then

                     dem2 = em2
                     izz1 = iz1
                     izz2 = iz2

                  else

                     dem2 = em1
                     izz1 = iz2
                     izz2 = iz1

                  end if


                     srtd  = srt - ( dem2 - rmass )

                     z11   = soo(1,srtd)
                     z10   = soo(2,srtd)

                     prfd  = sqrt( ( srt**2 - dem1**2 - dem2**2 )**2
     &                     - 4.0 * ( dem1 * dem2 )**2 ) / ( 2.0 * srt )


                  if( iz1 + iz2 .eq. 1 ) then

                     sigrdr = 2.0 * z11 + z10

                  else if( iz1 .eq. iz2 ) then

                     sigrdr = 1.5 * z11

                  else

                     sigrdr = 0.5 * z11 + z10

                  end if

                     pidn   = pideno(2,rmass+pmass,dem2,srt)

                     sigiv  = 2.0

                     sigdnr = sigrdr * prfd**2 / pri**2 / pidn / sigiv

*-----------------------------------------------------------------------

                  if( x1 .lt. ( sigel + sigdnr ) / sig ) then

                        ichannel = 11

                        idd1 = 1
                        idd2 = 3

                     if( ( iz1 + iz2 .eq. 1 ) .and.
     &                   ( x1 .gt. ( sigel + 1.5 * z11 ) / sig ) ) then

                        izzd = izz1
                        izz1 = izz2
                        izz2 = izzd

                     end if

                     if( ( iz1 + iz2 .ne. 1 ) .and.
     &                   ( iz1 .ne. iz2 ) ) then

                        izz1 = ( iz1 + iz2 ) / 2
                        izz2 = ( iz1 + iz2 ) / 2

                     end if

                     goto 100

                  end if

               end if

            end if

*-----------------------------------------------------------------------

*           (5-2) R + D -> N + D

            if( indrd .ne. 0 ) then

*-----------------------------------------------------------------------

                  if( id1 .eq. 3 ) then

                     dem1 = rmass
                     dem2 = em2

                     izz1 = iz1
                     izz2 = iz2

                  else

                     dem1 = rmass
                     dem2 = em1

                     izz1 = iz2
                     izz2 = iz1

                  end if


                     srtd  = srt - ( dem2 - rmass )

                     z01   = soo(3,srtd)

                     prfd  = sqrt( ( srt**2 - dem1**2 - dem2**2 )**2
     &                     - 4.0 * ( dem1 * dem2 )**2 ) / ( 2.0 * srt )

                     sigdrd = 1.5 * z01

                     pidn   = pideno(3,rmass+pmass,dem2,srt)

                     sigiv  = 1.0

                     sigrnd = sigdrd * prfd**2 / pri**2 / pidn / sigiv

*-----------------------------------------------------------------------

                  if( x1 .lt. ( sigel + sigdnr + sigrnd ) / sig ) then

                     ichannel = 13

                     idd1 = 1
                     idd2 = 2

                     goto 100

                  end if

            end if

*-----------------------------------------------------------------------

            return

         end if

*-----------------------------------------------------------------------

*======================================================================*

*        [6] Nstar - Nstar  ;  Inelastic Collisions

*                        (6-1) R + R -> N + R

         if( ( id1 .eq. 3 ) .and. ( id2 .eq. 3 ) ) then

*======================================================================*

                  sigrnr = 0.0

*-----------------------------------------------------------------------

*           (6-1) R + R -> N + R

            if( inrrr .ne. 0 ) then

*-----------------------------------------------------------------------

                  dem1 = rmass
                  dem2 = em2

                  srtd  = srt - ( dem2 - rmass )

                  z01   = soo(3,srtd)

                  prfd  = sqrt( ( srt**2 - dem1**2 - dem2**2 )**2
     &                    - 4.0 * ( dem1 * dem2 )**2 ) / ( 2.0 * srt )

                  sigrrr = 1.5 * z01

                  pidn   = pideno(3,rmass+pmass,dem2,srt)

                  if( iz1 .eq. iz2 ) then

                    sigiv = 0.5

                  else

                    sigiv = 1.0

                  end if

                  sigrnr = sigrrr * prfd**2 / pri**2 / pidn / sigiv

*-----------------------------------------------------------------------

               if( x1 .lt. ( sigel + sigrnr ) / sig ) then

                  ichannel = 15

                  izz1 = iz1
                  izz2 = iz2

                  idd1 = 1
                  idd2 = 3

                  goto 100

               end if

            end if


*-----------------------------------------------------------------------

            return

         end if

*-----------------------------------------------------------------------

*======================================================================*

         return

*======================================================================*


  100    continue

               if( x2 .gt. 0.5 ) then

                  em1 = dem1
                  em2 = dem2
                  iz1 = izz1
                  iz2 = izz2
                  id1 = idd1
                  id2 = idd2

               else

                  em1 = dem2
                  em2 = dem1
                  iz1 = izz2
                  iz2 = izz1
                  id1 = idd2
                  id2 = idd1

               end if


                  pr  = sqrt( ( srt**2 - em1**2 - em2**2 )**2
     &                - 4.0 * ( em1 * em2 )**2 ) / ( 2.0 * srt )


*======================================================================*

      end if

*======================================================================*


*-----------------------------------------------------------------------
*     Angular Distribution
*-----------------------------------------------------------------------

      if( ichannel .eq. 1 ) then

              ta  = -2.0 * pra**2
              x   = rn(0)

              plsq = srt**4 / ( 4.0 * rmass**2 ) - srt**2
              pl   = 0.0
              if( plsq .gt. 0.0 ) pl = sqrt( plsq )

*----------------------------------------------------------------------*
*       p + n collision
*----------------------------------------------------------------------*

        if( iz1 .ne. iz2 .or. id1 .ne. id2 ) then

          if( pl .le. 0.225 ) then

             c1 = 1.0 - 2.0 * x

          else

             if( pl .le. 0.6 ) then
               a = 6.2 * ( pl - 0.225 )/ 0.375
             else if( pl .le. 1.6 ) then
               a =  - 1.63 * pl + 7.16
             else if( pl .le. 2.0 ) then
               a = 5.5 * pl**8 / ( 7.7 + pl**8 )
             else
               a = 5.34 + 0.67 * ( pl - 2.0 )
             end if

             tt1 = log( ( 1.0 - x ) * exp( 2.0 * a * ta ) + x ) / a

             bprob = 1.0
cKN          if( pl .gt. 0.8 ) bprob = 0.64 / pl / pl
             if( pl .gt. 0.8 ) bprob = 0.64 / pl / pl

             bprob = bprob / ( 1. + bprob )

             if( rn(0) .gt. bprob ) then
               c1 = 1.0 - tt1 / ta
             else
               c1  = - 1.0 + tt1 / ta
             end if

             if( abs(c1) .gt. 1.0 ) c1 = 1.0 - 2.0 * x

          end if

*----------------------------------------------------------------------*
*       p + p collisions
*----------------------------------------------------------------------*

        else

          if( pl .le. 0.0 ) then

             c1 = 1.0 - 2.0 * X

          else 

             if( pl .le. 2.0 ) then
               a = 5.5 * pl**8 / ( 7.7 + pl**8 )
             else
               a = 5.34 + 0.67 * ( pl - 2.0 )
             end if

             tt1 = log( (1.-x) * exp(2.*a*ta) + x ) / a
             c1  = 1.0 - tt1 / ta

             if( abs(c1) .gt. 1.0 ) c1 = 2.0 * x - 1.0

          end if

        end if

*-----------------------------------------------------------------------
*        old Cugnon
*-----------------------------------------------------------------------

c              as  = ( 3.65 * asrt )**6
c              a   = 6.0 * as / (1.0 + as)
c              ta  = -2.0 * pra**2
c              x   = rn(0)
c              t1  = log( (1-x) * exp(2.*a*ta) + x )  /  a
c              c1  = 1.0 - t1/ta
c              if(abs(c1).gt.1.0) c1 = 2.0 * x - 1.0

*-----------------------------------------------------------------------

      else

*-----------------------------------------------------------------------

               srtd  = srt
               asrtd = asrt
               asrtn = asrt
               prad  = pra

*-----------------------------------------------------------------------

         if( ichannel .eq.  4 .or. ichannel .eq.  5 .or.
     &       ichannel .eq. 12 .or. ichannel .eq. 13 .or.
     &       ichannel .eq. 14 .or. ichannel .eq. 15 ) then

               srtd  = srt  - ( smass - dmass )
               asrtn = asrt - ( smass - dmass )

            if( asrtn .gt. 0.0 ) then

               asrtd = asrtn
               prad  = sqrt( ( srtd**2 - p(5,i1)**2 - p(5,i2)**2 )**2
     &               - 4.0 * ( p(5,i1) * p(5,i2) )**2 ) / ( 2.0 * srtd )

            end if

         end if

*-----------------------------------------------------------------------

         if( rn(0) .lt. 0.5 ) then

            if( asrtn .gt. 0.0 ) then

               as  = ( 3.65 * asrtd )**6
               a   = as / (1.0 + as) * srtd**4 * 0.14
               ta  = -2.0 * prad**2
               x   = rn(0)
               t1  = log( (1-x) * exp( max(-50.d0,2.*a*ta) ) + x ) / a
               c1  = 1.0 - t1/ta
               if(abs(c1).gt.1.0) c1 = 2.0 * x - 1.0

            else

               c1 = 2.0 * rn(0) - 1.0

            end if

         else

            if( srtd .lt. 2.14 ) then

                c1 = 2.0 * rn(0) - 1.0

            else

               if( srtd .gt. 2.4 ) then

                  b1 = 0.06
                  b3 = 0.4

               else

                  b1 = 29.0286 - 23.749  * srtd + 4.86549 * srtd**2
                  b3 =-30.3283 + 25.5257 * srtd - 5.30129 * srtd**2

               end if

                  pp3 = b1 / (3. * b3)
                  qq3 = 0.5 * (0.5 - rn(0)) / b3
                  pq3 = sqrt(qq3**2 + pp3**3)
                  uu  = (-qq3 + pq3 )**(1./3.)
                  vv  = ( qq3 + pq3 )**(1./3.)
                  c1  = uu - vv
                  if( abs(c1) .gt. 1. ) c1 = c1 / abs(c1)

            end if

         end if

      end if

*-----------------------------------------------------------------------
*     set the new momentum coordinates
*-----------------------------------------------------------------------

               t1  = 2.0 * pi * rn(0)

               if( pcm(1) .eq. 0.0 .and. pcm(2) .eq. 0.0 ) then

                 t2 = 0.0

               else

                 t2 = atan2(pcm(2),pcm(1))

               end if

               s1   = sqrt( 1.0 - c1**2 )
               s2  =  sqrt( 1.0 - c2**2 )
               ct1  = cos(t1)
               st1  = sin(t1)
               ct2  = cos(t2)
               st2  = sin(t2)
               ss   = c2 * s1 * ct1  +  s2 * c1

               pcm(1) = pr * ( ss*ct2 - s1*st1*st2 )
               pcm(2) = pr * ( ss*st2 + s1*st1*ct2 )
               pcm(3) = pr * ( c1*c2  - s1*s2 *ct1 )


*-----------------------------------------------------------------------
*     initial energy of two paticles
*-----------------------------------------------------------------------

               call epotall(epot)

               eini = epot + p(4,i1) + p(4,i2)
               etwo =        p(4,i1) + p(4,i2)

*-----------------------------------------------------------------------
*     charge and state identification
*-----------------------------------------------------------------------

         if( ichannel .gt. 1 ) then

               inds(i1) = id1
               inds(i2) = id2

            if( inds(i1).eq.1 ) then

               inuc(i1) = 1

            else

               inuc(i1) = 0

            end if

            if( inds(i2).eq.1 ) then

               inuc(i2) = 1

            else

               inuc(i2) = 0

            end if

               ichg(i1) = iz1
               ichg(i2) = iz2

               p(5,i1)  = em1
               p(5,i2)  = em2

         end if

*-----------------------------------------------------------------------

         do 7000 kk = 1, 4

*-----------------------------------------------------------------------
*           lorentz-transformation into reference frame
*-----------------------------------------------------------------------

            e1cm  = sqrt (em1**2 + pcm(1)**2 + pcm(2)**2 + pcm(3)**2)
            p1beta  = pcm(1)*beta(1) + pcm(2)*beta(2) + pcm(3)*beta(3)
            transf  = gamma * ( gamma * p1beta / (gamma + 1) + e1cm )

                  p(1,i1) = beta(1) * transf + pcm(1)
                  p(2,i1) = beta(2) * transf + pcm(2)
                  p(3,i1) = beta(3) * transf + pcm(3)

                  p(4,i1) = sqrt( p(5,i1)**2 + p(1,i1)**2
     &                          + p(2,i1)**2 + p(3,i1)**2 )


            e2cm  = sqrt (em2**2 + pcm(1)**2 + pcm(2)**2 + pcm(3)**2)
            transf  = gamma * (-gamma * p1beta / (gamma + 1.) + e2cm)

                  p(1,i2) = beta(1) * transf - pcm(1)
                  p(2,i2) = beta(2) * transf - pcm(2)
                  p(3,i2) = beta(3) * transf - pcm(3)

                  p(4,i2) = sqrt( p(5,i2)**2 + p(1,i2)**2
     &                          + p(2,i2)**2 + p(3,i2)**2 )

*-----------------------------------------------------------------------
*           check energy of two paticles
*-----------------------------------------------------------------------

            call caldis2(i1,i2)
            call epotall(epot)

            efin = epot + p(4,i1) + p(4,i2)

*-----------------------------------------------------------------------

            if( abs( eini - efin ) .lt. epse ) return

*-----------------------------------------------------------------------

               cona = ( eini - efin + etwo ) / gamma

               fac2 = 1. / ( 4.0 * cona**2 * pr**2 ) *
     &                ( ( cona**2 - ( em1**2 + em2**2 ) )**2 
     &                     - 4.0 * em1**2 * em2**2 )

               if( fac2 .gt. 0.0 ) then

                  fact   = sqrt( fac2 )
                  pcm(1) = pcm(1) * fact
                  pcm(2) = pcm(2) * fact
                  pcm(3) = pcm(3) * fact

               end if

 7000    continue

         ichannel = 99

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine pionem(dt)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to calculate the decay of delta or N*                   *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              dt          : time interval for decay                   *
*                                                                      *
*        Comments : counting variables                                 *
*                                                                      *
*              lcoll(22) : D + N -> pi                                 *
*              lcoll(23) : R + N -> pi                                 *
*              lcoll(24) : R + D -> pi                                 *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param01.inc'
      include 'param02.inc'

*-----------------------------------------------------------------------

      parameter     ( epse = 0.0001 )
      parameter     ( dtmax = 40.0 )
      parameter     ( prmin = 0.0001 )

*-----------------------------------------------------------------------

      common /vriab0/ massal, massba, mmeson

      common /coodrp/ r(5,nnn),  p(5,nnn)
      common /coodid/ ichg(nnn), inuc(nnn), ibry(nnn), inds(nnn),
     &                inun(nnn), iavd(nnn), ihis(nnn)

      common /coln00/ lcoll(30)

*-----------------------------------------------------------------------

      dimension       beta(3), pcm(3)

*-----------------------------------------------------------------------
*        output unit
*-----------------------------------------------------------------------

            ieo = 6

*-----------------------------------------------------------------------
*        zero set of counter
*-----------------------------------------------------------------------

                  lcoll(22)= 0
                  lcoll(23)= 0
                  lcoll(24)= 0

*-----------------------------------------------------------------------
*        factor of N* -> D + pi to N* -> N + pi
*-----------------------------------------------------------------------

                  frdp = 0.4

*-----------------------------------------------------------------------

      do 800 i  = 1, massba

            if( inds(i) .ne. 2 .and. inds(i) .ne. 3 ) goto 800

*-----------------------------------------------------------------------

                  px1 = p(1,i)
                  py1 = p(2,i)
                  pz1 = p(3,i)
                  e1  = p(4,i)
                  em1 = p(5,i)

                  iz1 = ichg(i)
                  in1 = inuc(i)
                  id1 = inds(i)
                  iu1 = inun(i)

                  beta(1)  = - px1 / e1
                  beta(2)  = - py1 / e1
                  beta(3)  = - pz1 / e1
                  betasq   =   beta(1)**2 + beta(2)**2 + beta(3)**2
                  gamma    =   1.0 / sqrt( 1.0 - betasq )

*-----------------------------------------------------------------------
*           can the Delta or N* decay in this time step ?
*-----------------------------------------------------------------------

                  dt0  = dt / gamma

                  dem1 = em1
                  dem2 = 0.0

                  call resmas(0,id1,1,em1,dem1,dem2,gam,fqrq)

                  w    = exp( - dt0 * gam / hbc )


            if( dt .le. dtmax .and. rn(0) .lt. w )  goto 800

*-----------------------------------------------------------------------

               if( id1 .eq. 2 ) then

*                    ---------  Delta -> Nucleon + pi ---------

                     dnmass = rmass
                     inid   = 1
                     idid   = 1

               else if( id1 .eq. 3 ) then


                  if( ( em1 .gt. rmass + 2.0 * pmass ) .and.
     &                ( rn(0) .lt. frdp )) then

*                    ---------  N* -> Delta + pi ---------

                     dem1 = rmass + pmass
                     dem2 = pmass

                     call resmas(1,2,1,em1,dem1,dem2,gam0,fqrq)

                     dnmass = dem1
                     inid   = 0
                     idid   = 2

                  else

*                    ---------  N* -> Nucleon + pi ---------

                     dnmass = rmass
                     inid   = 1
                     idid   = 1

                  end if

               end if

*-----------------------------------------------------------------------
*           now the Delta or N* may decay
*-----------------------------------------------------------------------

                  j = massal + 1


                  if( j .gt. nnn ) then

                     write(ieo,'('' Error: too many pions,'',
     &                           '' use larger nnn value.'')')
                     write(ieo,'('' ====='')')

                     call parastop( 222 )

                  end if


*-----------------------------------------------------------------------
*           check potential energy
*-----------------------------------------------------------------------

                  call epotall(epot)

                  eini = epot + e1
                  etwo =        e1

*-----------------------------------------------------------------------
*           charge and state identification
*-----------------------------------------------------------------------

                        ichg(j) = 0

                        xx = rn(0)

               if( inds(i) .eq. 2 ) then

                  if( ichg(i) .eq. 2 ) then

                        ichg(i) =  1
                        ichg(j) =  1

                  else if( ichg(i) .eq. -1 ) then

                        ichg(i) =  0
                        ichg(j) = -1

                  else

                     if( xx .gt. 0.66666667 ) then

                        ichg(j) =  2 * ichg(i) - 1
                        ichg(i) =  1 - ichg(i)

                     end if

                  end if

               else if( inid .eq. 1 ) then

                     if( xx .gt. 0.33333 ) then

                        ichg(j) =  2 * ichg(i) - 1
                        ichg(i) =  1 - ichg(i)

                     end if

               else

                  if( ichg(i) .eq. 1 ) then

                     if( xx .lt. 0.5 ) then

                        ichg(i) =  2
                        ichg(j) = -1

                     else if( xx .lt. 0.66666667 ) then

                        ichg(i) =  0
                        ichg(j) =  1

                     else

                        ichg(i) =  1
                        ichg(j) =  0

                     end if

                  else

                     if( xx .lt. 0.5 ) then

                        ichg(i) = -1
                        ichg(j) =  1

                     else if( xx .lt. 0.66666667 ) then

                        ichg(i) =  1
                        ichg(j) = -1

                     else

                        ichg(i) =  0
                        ichg(j) =  0

                     end if

                  end if

               end if

*-----------------------------------------------------------------------

                        inuc(i) = inid
                        inds(i) = idid
                        inun(i) = j
                        p(5,i)  = dnmass

                        inuc(j) = 0
                        inun(j) = i
                        inds(j) = 4
                        p(5,j)  = pmass

                        massal  = massal + 1

*-----------------------------------------------------------------------
*           the outgoing nucleon momentum : pcm(i) in the Delta cms
*           and determine the pion position
*-----------------------------------------------------------------------

                  ntag = 1

                  pr = max( prmin,
     &                     .25 * ( em1**2 - dnmass**2 - pmass**2 )**2
     &                                    - dnmass**2 * pmass**2 )

                  pr = sqrt( pr ) / em1

                  rr = 10.0

               do while( ( rr .lt. 0.001 ) .or. ( rr .gt. 1.0 ) )

                  xx = 1. - 2. * rn(0)
                  yy = 1. - 2. * rn(0)
                  zz = 1. - 2. * rn(0)
                  rr = sqrt( xx**2 + yy**2 + zz**2 )

               end do

                  pcm(1)  = pr * xx / rr
                  pcm(2)  = pr * yy / rr
                  pcm(3)  = pr * zz / rr

                  r(1,j) = r(1,i)
                  r(2,j) = r(2,i)
                  r(3,j) = r(3,i)

*-----------------------------------------------------------------------

         do 7000 kk = 1, 4

*-----------------------------------------------------------------------
*           Lorentz-transformation into reference frame
*-----------------------------------------------------------------------

                  e1cm   = sqrt( dnmass**2 + pcm(1)**2
     &                                     + pcm(2)**2
     &                                     + pcm(3)**2 )

                  p1beta = pcm(1) * beta(1)
     &                   + pcm(2) * beta(2)
     &                   + pcm(3) * beta(3)

                  transf = gamma * ( gamma * p1beta
     &                   / ( gamma + 1 ) - e1cm )

                  p(1,i) = beta(1) * transf + pcm(1)
                  p(2,i) = beta(2) * transf + pcm(2)
                  p(3,i) = beta(3) * transf + pcm(3)

                  p(4,i) = sqrt( p(5,i)**2 + p(1,i)**2
     &                         + p(2,i)**2 + p(3,i)**2 )


                  e2cm   = sqrt( pmass**2 + pcm(1)**2
     &                                    + pcm(2)**2
     &                                    + pcm(3)**2 )

                  transf = gamma * ( -gamma * p1beta 
     &                   / ( gamma + 1 ) - e2cm )

                  p(1,j) = beta(1) * transf - pcm(1)
                  p(2,j) = beta(2) * transf - pcm(2)
                  p(3,j) = beta(3) * transf - pcm(3)

                  p(4,j) = sqrt( p(5,j)**2 + p(1,j)**2
     &                         + p(2,j)**2 + p(3,j)**2 )

*-----------------------------------------------------------------------
*   CHECK ENERGY OF TWO PATICLES
*-----------------------------------------------------------------------

                  call caldis2(i,j)
                  call epotall(epot)

                  efin = epot + p(4,i) + p(4,j)

*-----------------------------------------------------------------------

            if( abs( eini - efin ) .gt. epse ) then

                  cona = ( eini - efin + etwo ) / gamma

                  fac2 = 1. / ( 4.0 * cona**2 * pr**2 )
     &                 * ( ( cona**2 - ( p(5,i)**2 + p(5,j)**2 ) )**2
     &                     - 4.0 * p(5,i)**2 * p(5,j)**2 )

               if( fac2 .gt. 0.0 ) then

                  fact = sqrt( fac2 )

                  pcm(1) = pcm(1) * fact
                  pcm(2) = pcm(2) * fact
                  pcm(3) = pcm(3) * fact

               else

                  if( dt .gt. dtmax ) ntag = 0
                  goto 9000

               end if

            else

                  ntag = 0
                  goto 9000

            end if

 7000    continue

*-----------------------------------------------------------------------
*           check final pauli-blocking of out going nucleon
*-----------------------------------------------------------------------

 9000          if( ntag .eq. 0 .and. dt .lt. dtmax )
     &               call pauli(i,ntag,phase)

*-----------------------------------------------------------------------
*              it does not decay, reset variables
*-----------------------------------------------------------------------

               if( ntag .eq. 1 ) then

                  p(1,i)  = px1 
                  p(2,i)  = py1
                  p(3,i)  = pz1
                  p(4,i)  = e1 
                  p(5,i)  = em1

                  ichg(i) = iz1
                  inuc(i) = in1
                  inun(i) = iu1
                  inds(i) = id1

                  ichg(j) = 0
                  inuc(j) = 0
                  inun(j) = 0
                  inds(j) = 0

                  massal  = massal - 1

                  call caldis2(i,j)

*-----------------------------------------------------------------------
*              Decay, count the summary variables
*-----------------------------------------------------------------------

               else

                     iavd(i) = iavd(i) + iavd(i)/abs(iavd(i))

                     ihis(i) = - abs(ihis(i)) - 1

                  if( id1 .eq. 2 ) ihis(j) =   abs(ihis(i))
                  if( id1 .eq. 3 ) ihis(j) = - abs(ihis(i))

                  if( id1 .eq. 2 ) then

                     lcoll(22) = lcoll(22) + 1

                  else if( inid .eq. 1 ) then

                     lcoll(23) = lcoll(23) + 1

                  else if( inid .eq. 0 ) then

                     lcoll(24) = lcoll(24) + 1

                  end if

               end if

*-----------------------------------------------------------------------

  800  continue

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine pionab
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to calculate pion absorption                            *
*                                                                      *
*        Comments : counting variables                                 *
*                                                                      *
*              lcoll(25) : N + pi -> D                                 *
*              lcoll(26) : N + pi -> R                                 *
*              lcoll(27) : D + pi -> R                                 *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param01.inc'
      include 'param02.inc'

*-----------------------------------------------------------------------

      parameter     ( epse = 0.0001 )

*-----------------------------------------------------------------------

      parameter     ( delpi  = 4.0 )

      parameter     ( pirr   = 0.977, pir0  = 2.07)
      parameter     ( pir1   = 2.52,  pir2  = 1.49)

*-----------------------------------------------------------------------
*
*     delpi         : maximum spatial distance for which a absorption
*                     still can occur only for 200 MeV
*
*     bcmax         : maximum impact parameter for
*                     pirr = 0.977 fm corresponds  to  30 mb
*                     pir0 = 2.07  fm corresponds  to 135 mb
*                     pir1 = 2.52  fm corresponds  to 200 mb
*                     pir2 = 1.49  fm corresponds  to  70 mb
*
*-----------------------------------------------------------------------

      common /const1/ elab, rdist, bmin, bmax, ibch, ibin
      common /const2/ dt, ntmax, iprun, iprun0

      common /vriab0/ massal, massba, mmeson

      common /coodrp/ r(5,nnn),  p(5,nnn)
      common /coodid/ ichg(nnn), inuc(nnn), ibry(nnn), inds(nnn),
     &                inun(nnn), iavd(nnn), ihis(nnn)
      common /cood2b/ rha(nnn,nnn), rhe(nnn,nnn), rhc(nnn,nnn)
      common /cood2r/ rr2(nnn,nnn), rbij(nnn,nnn), pp2(nnn,nnn)

      common /coln00/ lcoll(30)

*-----------------------------------------------------------------------
*        zero set of counter
*-----------------------------------------------------------------------

               lcoll(25) = 0
               lcoll(26) = 0
               lcoll(27) = 0

      if( massal .eq. massba ) return

*-----------------------------------------------------------------------

      do 800 i1  = massba + 1, massal

                  px1 = p(1,i1)
                  py1 = p(2,i1)
                  pz1 = p(3,i1)

                  e1  = p(4,i1)
                  em1 = p(5,i1)

                  iz1 = ichg(i1)
                  in1 = inuc(i1)
                  id1 = inds(i1)
                  iu1 = inun(i1)
                  ia1 = iavd(i1)
                  ih1 = ihis(i1)

*-----------------------------------------------------------------------
*        look for an absorbent nucleon
*-----------------------------------------------------------------------

         do 600 i2 = 1, massba

*-----------------------------------------------------------------------
*           avoid second collisions for the same pairs
*-----------------------------------------------------------------------

               if( i2 .eq. inun(i1) .and.
     &             i1 .eq. inun(i2) ) goto 600

*-----------------------------------------------------------------------
*           look only for nucleons or Delta + pion with proper charge
*-----------------------------------------------------------------------

               if( inds(i2) .eq. 3 ) goto 600

               if( inds(i2) .eq. 2 .and.
     &           ( iz1 + ichg(i2) + 2 ) / 2 .ne. 1 ) goto 600

*-----------------------------------------------------------------------
*           delpi is appled only for 5000 MeV
*-----------------------------------------------------------------------

               if( elab .le. 5.0 .and.
     &             rr2(i1,i2) .gt. delpi**2 ) goto 600

*-----------------------------------------------------------------------
*           check impact parameter and collision time
*-----------------------------------------------------------------------

                  px2 = p(1,i2)
                  py2 = p(2,i2)
                  pz2 = p(3,i2)

                  e2  = p(4,i2)
                  em2 = p(5,i2)

                  iz2 = ichg(i2)
                  in2 = inuc(i2)
                  id2 = inds(i2)
                  iu2 = inun(i2)
                  ia2 = iavd(i2)
                  ih2 = ihis(i2)

                  s   = (e1+e2)**2 - (px1+px2)**2 - (py1+py2)**2
     &                                             - (pz1+pz2)**2
                  srt = sqrt(s)

*-----------------------------------------------------------------------
*           impact parameter
*-----------------------------------------------------------------------

                  dx   = r(1,i1) - r(1,i2)
                  dy   = r(2,i1) - r(2,i2)
                  dz   = r(3,i1) - r(3,i2)
                  rsq  = dx**2 + dy**2 + dz**2

                  p12  = e1 * e2 - px1 * px2 - py1 * py2 - pz1 * pz2
                  p1dr = px1 * dx + py1 * dy + pz1 * dz
                  p2dr = px2 * dx + py2 * dy + pz2 * dz
                  a12  = 1.0 - ( em1 * em2 / p12 ) ** 2
                  b12  = p1dr / em1 - p2dr * em1 / p12
                  c12  = rsq + ( p1dr / em1 )**2

                  brel = sqrt( abs(c12 - b12**2/a12) )

*-----------------------------------------------------------------------
*           bcmax of maximum cross section
*-----------------------------------------------------------------------

                  izzz = ( iz2 + iz1 + 2 ) / 2

                  pind = 0.0
                  pird = 0.0

               if( inds(i2) .eq. 1 ) then

                  pind = pir0**2

                  if( izzz .ne. 1 )   pind = pir1**2

                  if( ( izzz .eq. 1 ) .and.
     &                ( iz1  .ne. 0 ) ) pind = pir2**2

                  if( izzz .eq. 1 ) pird = pirr**2

               else

                  pird = pirr**2

               end if

                  pirn = sqrt( pind + pird )

*-----------------------------------------------------------------------
*           is their impact parameter small enough ?
*-----------------------------------------------------------------------

               if( brel .gt. pirn ) goto 600

*-----------------------------------------------------------------------
*           average time-shift of the collision in the fixed frame
*           will particles get closest point in this time interval ?
*-----------------------------------------------------------------------

                  b21 =   - p2dr / em2 + p1dr * em2 / p12
                  t1  = (   p1dr / em1 - b12 / a12 ) * e1 / em1
                  t2  = ( - p2dr / em2 - b21 / a12 ) * e2 / em2

               if ( abs( t1 + t2 ) .gt. dt )  goto 600

*-----------------------------------------------------------------------
*           now the pion may be absorbed in this time step
*           Check the cross section with phase space factor
*-----------------------------------------------------------------------

                  xx = rn(0)
                  yy = 0.0
                  zz = 0.0

              if( inds(i2) .eq. 1 ) then

                  dem1 = srt
                  dem2 = 0.0

                  call resmas(0,2,1,srt,dem1,dem2,gamdl,fqrq)

                  yy = fqrq 
     &               / ( 4.0 * ( srt - dmass )**2 / gamdl**2 + 1.0 )
                  yy = ( pind * yy ) / ( pind + pird )

              end if

              if( izzz .eq. 1 ) then

                  dem1 = srt
                  dem2 = 0.0

                  call resmas(0,3,1,srt,dem1,dem2,gamrs,fqrq)

                  zz  = fqrq
     &                / ( 4. * ( srt - smass )**2 / gamrs**2 + 1.0 )
                  zz = ( pird * zz ) / ( pind + pird )

              end if


               if( inds(i2) .eq. 1 .and. xx .gt. yy + zz ) goto 600
               if( inds(i2) .eq. 2 .and. xx .gt.      zz ) goto 600

*-----------------------------------------------------------------------
*           check potential energy
*-----------------------------------------------------------------------

                  call epotall(epot)
                  eini = epot + e1 + e2

*-----------------------------------------------------------------------
*           set variables
*-----------------------------------------------------------------------

               if( izzz .eq. 1 .and. xx .le. zz ) then

                  inds(i2) = 3

               else

                  inds(i2) = 2

               end if

                  p(1,i2)  = px1 + px2
                  p(2,i2)  = py1 + py2
                  p(3,i2)  = pz1 + pz2
                  p(5,i2)  = srt

                  pps      = p(1,i2)**2 + p(2,i2)**2 + p(3,i2)**2
                  p(4,i2)  = sqrt( p(5,i2)**2 + pps )

                  ichg(i2) = iz2 + iz1
                  inuc(i2) = 0

                  inun(i2) = 0
                  iavd(i2) = iavd(i2) + iavd(i2)/abs(iavd(i2))
                  ihis(i2) = - abs(ihis(i2)) - abs(ihis(i1)) - 1

                  ichg(i1) = 0
                  inds(i1) = 0
                  inun(i1) = 0
                  iavd(i1) = 0
                  ihis(i1) = 0

*-----------------------------------------------------------------------

         do 7000 kk = 1, 4

*-----------------------------------------------------------------------
*           check final energy
*-----------------------------------------------------------------------

                  call caldis2(i2,i1)
                  call epotall(epot)

                  efin = epot + p(4,i2)

*-----------------------------------------------------------------------
*           pion is absorbed, counting varables
*-----------------------------------------------------------------------

            if( abs( eini - efin ) .lt. epse ) then

                  if( inds(i2) .eq. 2 ) then

                        lcoll(25) = lcoll(25) + 1

                  else if( inds(i2) .eq. 3 ) then

                     if( id2 .eq. 1 ) then

                        lcoll(26) = lcoll(26) + 1

                     else if( id2 .eq. 2 ) then

                        lcoll(27) = lcoll(27) + 1

                     end if

                  end if

                  goto 800

*-----------------------------------------------------------------------

            else

                  enew = eini - epot
                  srt  = sqrt( enew**2 - pps )

               if( srt .gt. rmass + pmass ) then

                  p(5,i2)  = srt
                  p(4,i2) = enew

               else

                  goto 8000

               end if

            end if

 7000    continue

*-----------------------------------------------------------------------
*              reset variables
*-----------------------------------------------------------------------

 8000             p(1,i2)  = px2
                  p(2,i2)  = py2
                  p(3,i2)  = pz2

                  p(5,i2)  = em2
                  p(4,i2)  = e2 

                  ichg(i2) = iz2
                  inuc(i2) = in2
                  inds(i2) = id2
                  inun(i2) = iu2
                  iavd(i2) = ia2
                  ihis(i2) = ih2

                  ichg(i1) = iz1
                  inuc(i1) = in1
                  inds(i1) = id1
                  inun(i1) = iu1
                  iavd(i1) = ia1
                  ihis(i1) = ih1

                  call caldis2(i2,i1)

*-----------------------------------------------------------------------

  600   continue

  800 continue

*-----------------------------------------------------------------------
*        renumbering of pion part
*-----------------------------------------------------------------------

                     inpion    = 0

         do 900 i = massba + 1, massal

            if( inds(i) .eq. 0 ) then

               do 950 j = massal, i + 1, -1

                  if( inds(j) .eq. 4 ) then

                     r(1,i)    =   r(1,j)
                     r(2,i)    =   r(2,j)
                     r(3,i)    =   r(3,j)
                     p(1,i)    =   p(1,j)
                     p(2,i)    =   p(2,j)
                     p(3,i)    =   p(3,j)
                     p(4,i)    =   p(4,j)
                     p(5,i)    =   p(5,j)
                     ichg(i)   =   ichg(j)
                     inds(i)   =   inds(j)
                     ihis(i)   =   ihis(j)
                     inds(j)   =   0
                     ihis(j)   =   0

                  do 100 k = 1, massal

                    rbij(i,k)  =   rbij(j,k)
                     rr2(i,k)  =   rr2 (j,k)
                     pp2(i,k)  =   pp2 (j,k)
                     rha(i,k)  =   rha (j,k)
                     rhe(i,k)  =   rhe (j,k)
                     rhc(i,k)  =   rhc (j,k)

                    rbij(k,i)  = - rbij(i,k)
                    rr2 (k,i)  =   rr2 (i,k)
                    pp2 (k,i)  =   pp2 (i,k)
                     rha(k,i)  =   rha (i,k)
                     rhe(k,i)  =   rhe (i,k)
                     rhc(k,i)  =   rhc (i,k)

  100             continue

                     inpion    = inpion + 1

                     goto 900

                  end if

  950          continue

            else

                     inpion = inpion + 1

            end if

  900    continue


                     massal = massba + inpion

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine fpidecay
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to calculate final decay of the resonances              *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      common /coln00/ lcoll(30)
      common /swich1/ ipot, insys, irkg, icolt

*-----------------------------------------------------------------------
*        final pion decay
*-----------------------------------------------------------------------

            if( icolt .eq. 1 ) then

                     timepi = 1000.

                  call pionem(timepi)

                     ncol22 = lcoll(22)
                     ncol23 = lcoll(23)
                     ncol24 = lcoll(24)

                  call pionem(timepi)

                     lcoll(22) = ncol22 + lcoll(22)
                     lcoll(23) = ncol23 + lcoll(23)
                     lcoll(24) = ncol24 + lcoll(24)

            end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine resmas(ipw,isd1,isd2,srt,dem1,dem2,gam0,fqrq)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to calcurate the mass and width of Delta and N*         *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              ipw = 0     : only width                                *
*              ipw = 1     : determine the mass                        *
*                                                                      *
*              isd1 = 2    : for Delta                                 *
*              isd1 = 3    : for N*(1440)                              *
*                                                                      *
*              isd2 = 1    : mass2 = dem2                              *
*              isd2 = 2    : for Delta                                 *
*              isd2 = 3    : for N*(1440)                              *
*                                                                      *
*              srt         : sqrt(s)                                   *
*                                                                      *
*              dem1        : minimum mass of resonance 1               *
*                            mass of resonance                         *
*                                                                      *
*              dem2        : minimum mass of resonance 2               *
*                            mass of resonance                         *
*                                                                      *
*              gam0        : width                                     *
*                                                                      *
*              fqrq        : phase space factor for absorption         *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param01.inc'
      include 'param02.inc'

*-----------------------------------------------------------------------

      dimension      rmss(3), gamr(3),
     &               bet2(3), qqr2(3)

      dimension      isdm(3)

*-----------------------------------------------------------------------

      data           rmss / 0.938, 1.232, 1.440 /
      data           gamr / 0.000, 0.110, 0.200 /
      data           bet2 / 1.000, 0.090, 0.274 /
      data           qqr2 / 1.000, 0.051936, 0.158125 /
      data           isdm / 0, 1, 1 /

*-----------------------------------------------------------------------

            gam0   = 0.0
            fqrq   = 0.0

            rm2    = rmass**2
            pm2    = pmass**2

*-----------------------------------------------------------------------
*     calculate only the width
*-----------------------------------------------------------------------

         if( ipw .eq. 0 ) then

            qq2  = 0.25 * ( ( dem1**2 - rm2 + pm2 ) / dem1 )**2 - pm2

            if( qq2 .le. 0.0 ) return

            form  = ( 1. + qqr2(isd1) / bet2(isd1) )
     &            / ( 1. + qq2        / bet2(isd1) )

            gam0  =   sqrt( qq2 / qqr2(isd1) )**3
     &              * rmss(isd1) / dem1 * gamr(isd1) * form**2

            fqrq  = qqr2(isd1) / qq2

            return

         end if


*-----------------------------------------------------------------------
*     maximum phase space factor
*-----------------------------------------------------------------------

            dem10 = dem1
            dem20 = dem2

            p02 = ( srt - dem10 - dem20 )
     &          * ( srt + dem10 + dem20 )
     &          * ( srt - dem10 + dem20 )
     &          * ( srt + dem10 - dem20 )

            if( p02 .le. 0.0 ) p02 = 0.0001

*-----------------------------------------------------------------------
*     determine the mass
*-----------------------------------------------------------------------

            idem = 0

  200       idem = idem + 1

            if( idem .gt. 500 )  return

            dem1 = dem10 + ( srt - dem10 - dem20 ) * rn(0)
            dem2 = dem20 + ( srt - dem1  - dem20 ) * rn(0) * isdm(isd2)


*-----------------------------------------------------------------------
*     check of phase space factor
*-----------------------------------------------------------------------

            pr2 = ( srt - dem1 - dem2 )
     &          * ( srt + dem1 + dem2 )
     &          * ( srt - dem1 + dem2 )
     &          * ( srt + dem1 - dem2 )

            gmcn =  pr2 / p02

            if( rn(0) .gt. gmcn ) goto 200

*-----------------------------------------------------------------------
*     calculate the width of Delta or N* by Moniz
*     check whether this mass satisfies the bright wigner distribution
*-----------------------------------------------------------------------

               demax1 = srt - dem20

               qq2  = 0.25 * ( ( dem1**2 - rm2 + pm2 ) / dem1 )**2 - pm2

               if( qq2 .le. 0.0 ) return

               form  = ( 1. + qqr2(isd1) / bet2(isd1) )
     &               / ( 1. + qq2        / bet2(isd1) )

               gam2  = ( sqrt( qq2 / qqr2(isd1) )**3
     &                 * rmss(isd1) / dem1 * gamr(isd1) * form**2 )**2


               bwtop = 1.0

            if( demax1 .lt. rmss(isd1) ) then

               bwtop = 0.25 * gam2 
     &               / ( ( rmss(isd1) - demax1 )**2 + 0.25 * gam2 )

            end if

               fmcn  = 0.25 * gam2 
     &               / ( ( rmss(isd1) - dem1   )**2 + 0.25 * gam2 )

            if( rn(0) * bwtop .gt. fmcn )  goto 200

*-----------------------------------------------------------------------
*     for second resonances
*-----------------------------------------------------------------------

         if( isd2 .ge. 1 ) then

               demax2 = srt - dem1

               qq2  = 0.25 * ( ( dem2**2 - rm2 + pm2 ) / dem2 )**2 - pm2

               if( qq2 .le. 0.0 ) return

               form  = ( 1. + qqr2(isd2) / bet2(isd2) )
     &               / ( 1. + qq2        / bet2(isd2) )

               gam2  = ( sqrt( qq2 / qqr2(isd2) )**3
     &                 * rmss(isd2) / dem2 * gamr(isd2) * form**2 )**2

               bwtop = 1.0

            if( demax2 .lt. rmss(isd2) ) then

               bwtop = 0.25 * gam2 
     &               / ( ( rmss(isd2) - demax2 )**2 + 0.25 * gam2 )

            end if

               fmcn  = 0.25 * gam2 
     &               / ( ( rmss(isd2) - dem2   )**2 + 0.25 * gam2 )

            if( rn(0) * bwtop .gt. fmcn ) goto 200

         end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      function pideno(isd1,dem1,dem2,srt)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to calculate the correction factor                      *
*              of inverse cross section                                *
*                                                                      *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              isd1 = 2    ; for Delta                                 *
*              isd1 = 3    ; for N*(1440)                              *
*                                                                      *
*              dem1        ; minimum mass of resonance 1               *
*              dem2        ; fixed mass of 2 ( mass of partner )       *
*                                                                      *
*              srt         ; sqrt(s)                                   *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param01.inc'
      include 'param02.inc'

*-----------------------------------------------------------------------

      parameter     ( prmin = 0.0001 )

*-----------------------------------------------------------------------

      dimension      rmss(3), gamr(3),
     &               bet2(3), qqr2(3)

*-----------------------------------------------------------------------

      data           rmss / 0.938, 1.232, 1.440 /
      data           gamr / 0.000, 0.110, 0.200 /
      data           bet2 / 1.000, 0.090, 0.274 /
      data           qqr2 / 1.000, 0.051936, 0.158125 /

*-----------------------------------------------------------------------

               rm2    = rmass**2
               pm2    = pmass**2

               dlmsin = dem1
               dlmsfn = srt    - dem2
               dlmsdf = dlmsfn - dlmsin

*-----------------------------------------------------------------------

            if( dlmsfn .lt. 2.0 * rmss(isd1) - dem1 ) then

               defi = dlmsdf / 20.0
               nsimp = 20

            else

               defi = ( rmss(isd1) - dem1 ) / 5.0
               nsimp = nint( dlmsdf / defi )

            end if

*-----------------------------------------------------------------------

            seki = 0.0

         do 10 i = 0, nsimp

            defj = defi

            if( i .eq. 0 .or. i .eq. nsimp ) defj = defi / 2.0

            dlmas = defi * float(i)  + dlmsin

            qq2   = max( prmin, 0.25 *
     &                 ( ( dlmas**2 - rm2 + pm2 ) / dlmas )**2 - pm2 )

            form  = ( 1. + qqr2(isd1) / bet2(isd1) )
     &            / ( 1. + qq2        / bet2(isd1) )

            gam   = 0.5 * ( sqrt( qq2 / qqr2(isd1) )**3
     &              * rmss(isd1) / dlmas * gamr(isd1) * form**2 )

            fm    =  gam / ( ( dlmas - rmss(isd1) )**2 + gam**2 ) / pi

            seki  = seki + fm * defj

   10    continue

            pideno = seki

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
       function soo(is,srt)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to calculate delta, N* cross section                    *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              srt         : sqrt of s                                 *
*                                                                      *
*              IS = 1      : S11                                       *
*              IS = 2      : S10                                       *
*              IS = 3      : S01                                       *
*              IS = 4      : S10D                                      *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param01.inc'
      include 'param02.inc'

*-----------------------------------------------------------------------

      real          mn,mpi,md,m0,m02,mavs,lm0,lden

*-----------------------------------------------------------------------

              soo   = 0.0
              mn    = 938.9
              mpi   = 138.0
              md    = 1877.8
              hbarc = 197.33

*-----------------------------------------------------------------------

              hbc2  = 10. * hbarc**2
              s     = ( srt * 1000. )**2
              p2    = s/4.-mn * mn
              fterm = pi * hbc2 / 2.0 / p2

*-----------------------------------------------------------------------

           if( is .eq. 1 ) then

              alfa = 3.0
              beta = 0.9
              m0   = 1188.
              gam  = 99.02
              lm0  = 1220.
              gam0 = 120.

           end if

           if( is .eq. 2 ) then

              alfa = 14.00
              beta = -0.3
              m0   = 1245.
              gam  = 120.0
              lm0  = 1220.
              gam0 = 120.

           end if

           if( is .eq. 3 ) then

              alfa = 23
              beta = 1.5
              m0   = 1472.
              gam  = 300.0
              lm0  = 1430.
              gam0 = 200.

           end if

              s0    = ( mn + m0 )**2
              p02   = s0 / 4. - mn * mn
              m02   = m0 * m0
              gam2  = gam * gam
              zpl   = 2. / gam0 * ( sqrt(s) - mn - lm0 )
              zpl2  = zpl * zpl
              zmin  = 2. / gam0 * ( mn + mpi - lm0 )
              zmin2 = zmin * zmin
              avss  = atan(zpl) - atan(zmin)

              if( abs(avss) .lt. 1.e-5 ) return

              mavs  = 1. / avss
              mavs  = mavs * gam0 / 4.
     &              * log( ( 1.0 + zpl2 ) / ( 1.0 + zmin2 ) )
              mavs  = mavs + lm0
              sstar = mavs * mavs

              pr2   = ( s - ( mn - mavs )**2 )
     &              * ( s - ( mn + mavs )**2 ) / 4. / s

              q2    = ( sstar - ( mn - mpi )**2 )
     &              * ( sstar - ( mn + mpi )**2 ) / 4. / sstar

              q02   = ( m02 - ( mn - mpi )**2 )
     &              * ( m02 - ( mn + mpi )**2 ) / 4. / m02

              if( pr2 .le. 0.0. or.
     &            q2  .le. 0.0 .or.
     &            q02 .le. 0.0 .or.
     &            p02 .le. 0.0 ) return

              pr = sqrt(pr2)
              q  = sqrt(q2)
              q0 = sqrt(q02)
              p0 = sqrt(p02)

              soo = fterm * alfa * ( pr / p0 )**beta
              den = ( sstar - m02 )**2 + m02 * gam2
              soo = soo * m02 * gam2 * ( q / q0 )**3. / den


           if( is .eq. 2 ) then

              alfa = 6.03
              beta = 1.7
              m0   = 1203.
              gam  = 134.3

              s0   = ( mn + m0 )**2
              p02  = s0 / 4. - mn * mn
              gam2 = gam * gam
              m02  = m0 * m0
              spin = ( sqrt(s) - mn )**2

              pr2  = ( s - ( md - mpi )**2)
     &             * ( s - ( md + mpi )**2) / 4. / s

              if( pr2 .le. 0.0 .or. p02 .le. 0.0 ) return

              sood = fterm * alfa * ( sqrt ( pr2 / p02 ) )**beta
              lden = ( spin - m02 )**2 + m02 * gam2
              sood = sood * m02 * gam2 / lden

              soo  = soo + sood

           end if

*-----------------------------------------------------------------------

      return
      end


