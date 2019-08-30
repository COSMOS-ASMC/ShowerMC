************************************************************************
*                                                                      *
*        PART 4: Mean Field Part                                       *
*                                                                      *
*                                                                      *
*   List of subprograms in rough order of relevance with main purpose  *
*      ( s = subroutine, f = function, b = block data, e = entry )     *
*                                                                      *
*                                                                      *
*  s  caldisa   to calculate the two-body quantities for all particle  *
*  s  caldis2   to calculate the two-body quantities for i1 and i2     *
*  s  gradu     to calculate the graduate of e.o.m  for rkg            *
*  s  rk12      to propagate particles by 2th order Runge-Kutta-Gill   *
*  s  rkg4      to propagate particles by 2th order Runge-Kutta-Gill   *
*  s  pauli     to calculate the Pauli blocking factor                 *
*  s  epotall   to calculate total potential energy                    *
*  s  etotal    to calculate total energy of cluster                   *
*  s  pcmcl     to calculate energy and angular momentum of cluster    *
*  s  cldist    to determine nuclear cluster                           *
*                                                                      *
*                                                                      *
************************************************************************
*                                                                      *
************************************************************************
*                                                                      *
      subroutine caldisa
*                                                                      *
*                                                                      *
*        Last Revised:     2004 09 23                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to calculate the two-body quantities for all particle   *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param01.inc'
      include 'param02.inc'

*-----------------------------------------------------------------------

      common /vriab0/ massal, massba, mmeson

      common /coodrp/ r(5,nnn),  p(5,nnn)
      common /coodid/ ichg(nnn), inuc(nnn), ibry(nnn), inds(nnn),
     &                inun(nnn), iavd(nnn), ihis(nnn)

      common /poten1/ gamm, c0, c3, cs, cl, wl

      common /cood2b/ rha(nnn,nnn), rhe(nnn,nnn), rhc(nnn,nnn)
      common /cood2r/ rr2(nnn,nnn), rbij(nnn,nnn), pp2(nnn,nnn)

      common /swich0/ icoul, irelcr, iavoid
      common /swich1/ ipot, insys, irkg, icolt

      common /caldis/ c0w, c3w, clw, c0sw

*-----------------------------------------------------------------------

      do 100 j = 2, massal
            if( inds(j) .eq. 0 ) goto 100

      do 110 i = 1, j - 1
            if( inds(i) .eq. 0 ) goto 110

*-----------------------------------------------------------------------

            rbrb = 0.0
            bij2 = 0.0
            rij2 = 0.0
            pij2 = 0.0

            eij  = p(4,i) + p(4,j)

         do l = 1, 3

            rij  =   r(l,i) - r(l,j)
            pij  =   p(l,i) - p(l,j)
            bij  = ( p(l,i) + p(l,j) ) /  eij

            rbrb = rbrb + rij * bij
            bij2 = bij2 + bij * bij
            rij2 = rij2 + rij * rij
            pij2 = pij2 + pij * pij

         end do

            rbrb      = irelcr * rbrb
            gij2      = 1. / ( 1. - bij2 )

            rr2(i,j)  = rij2 + gij2 * rbrb * rbrb
            rr2(j,i)  = rr2(i,j)

            rbij(i,j) = gij2 * rbrb
            rbij(j,i) = - rbij(i,j)

            pp2(i,j)  = pij2
     &                + irelcr * ( -( p(4,i) - p(4,j) )**2
     &                + gij2 * ( ( p(5,i)**2 - p(5,j)**2 ) / eij )**2 )
*
            pp2(j,i)  = pp2(i,j)

*-----------------------------------------------------------------------
*        Gauss term
*-----------------------------------------------------------------------

               expa1 = - rr2(i,j) * c0w

            if( expa1 .gt. epsx ) then

               rh1 = exp( expa1 )

            else

               rh1 = 0.0

            end if

               rha(i,j) = ibry(i) * ibry(j) * rh1
               rha(j,i) = rha(i,j)

*-----------------------------------------------------------------------
*        Coulomb
*-----------------------------------------------------------------------

               rrs2  = rr2(i,j) + epscl
               rrs   = sqrt( rrs2 )

              erfij  = erf( rrs * c0sw ) / rrs

            rhe(i,j) = ichg(i) * ichg(j) * erfij
            rhe(j,i) = rhe (i,j)

            rhc(i,j) = ichg(i) * ichg(j)
     &               * ( - erfij + clw * rh1 ) / rrs2
            rhc(j,i) = rhc (i,j)

*-----------------------------------------------------------------------

  110 continue
  100 continue

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine caldis2(i1,i2)
*                                                                      *
*                                                                      *
*        Last Revised:     2004 09 23                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to calculate the two-body quantities for i1 and i2      *
*                                                                      *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param01.inc'
      include 'param02.inc'

*-----------------------------------------------------------------------

      common /vriab0/ massal, massba, mmeson

      common /coodrp/ r(5,nnn),  p(5,nnn)
      common /coodid/ ichg(nnn), inuc(nnn), ibry(nnn), inds(nnn),
     &                inun(nnn), iavd(nnn), ihis(nnn)

      common /poten1/ gamm, c0, c3, cs, cl, wl

      common /cood2b/ rha(nnn,nnn), rhe(nnn,nnn), rhc(nnn,nnn)
      common /cood2r/ rr2(nnn,nnn), rbij(nnn,nnn), pp2(nnn,nnn)

      common /swich0/ icoul, irelcr, iavoid
      common /swich1/ ipot, insys, irkg, icolt

      common /caldis/ c0w, c3w, clw, c0sw

*-----------------------------------------------------------------------

         idd = i2 - i1

         if( i1 .eq. i2 ) idd = 1

      do 100 j = i1, i2, idd
            if( inds(j) .eq. 0 ) goto 100
      do 110 i = 1, massal
            if( inds(i) .eq. 0 ) goto 110

            if( i .eq. j ) goto 110

*-----------------------------------------------------------------------

            rbrb = 0.0
            bij2 = 0.0
            rij2 = 0.0
            pij2 = 0.0

            eij  = p(4,i) + p(4,j)

         do l = 1, 3

            rij  =   r(l,i) - r(l,j)
            pij  =   p(l,i) - p(l,j)
            bij  = ( p(l,i) + p(l,j) ) /  eij

            rbrb = rbrb + rij * bij
            bij2 = bij2 + bij * bij
            rij2 = rij2 + rij * rij
            pij2 = pij2 + pij * pij

         end do

            rbrb      = irelcr * rbrb
            gij2      = 1. / ( 1. - bij2 )

            rr2(i,j)  = rij2 + gij2 * rbrb * rbrb
            rr2(j,i)  = rr2(i,j)

            rbij(i,j) = gij2 * rbrb
            rbij(j,i) = - rbij(i,j)

            pp2(i,j)  = pij2
     &                + irelcr * ( -( p(4,i) - p(4,j) )**2
     &                + gij2 * ( ( p(5,i)**2 - p(5,j)**2 ) / eij )**2 )
*
            pp2(j,i)  = pp2(i,j)

*-----------------------------------------------------------------------
*        Gauss term
*-----------------------------------------------------------------------

               expa1 = - rr2(i,j) * c0w

            if( expa1 .gt. epsx ) then

               rh1 = exp( expa1 )

            else

               rh1 = 0.0

            end if

               rha(i,j) = ibry(i) * ibry(j) * rh1
               rha(j,i) = rha(i,j)

*-----------------------------------------------------------------------
*        Coulomb
*-----------------------------------------------------------------------

               rrs2  = rr2(i,j) + epscl
               rrs   = sqrt( rrs2 )

              erfij  = erf( rrs * c0sw ) / rrs

            rhe(i,j) = ichg(i) * ichg(j) * erfij
            rhe(j,i) = rhe (i,j)

            rhc(i,j) = ichg(i) * ichg(j)
     &               * ( - erfij + clw * rh1 ) / rrs2
            rhc(j,i) = rhc (i,j)

*-----------------------------------------------------------------------

  110 continue
  100 continue

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine gradu(ffr,ffp)
*                                                                      *
*                                                                      *
*        Last Revised:     2004 09 23                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to calculate the graduate of e.o.m  for rkg             *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              ffr         :  gradu of r                               *
*              ffp         :  gradu of p                               *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param01.inc'
      include 'param02.inc'

*-----------------------------------------------------------------------

      common /vriab0/ massal, massba, mmeson

      common /coodrp/ r(5,nnn),  p(5,nnn)
      common /coodid/ ichg(nnn), inuc(nnn), ibry(nnn), inds(nnn),
     &                inun(nnn), iavd(nnn), ihis(nnn)
      common /cood2b/ rha(nnn,nnn), rhe(nnn,nnn), rhc(nnn,nnn)
      common /cood2r/ rr2(nnn,nnn), rbij(nnn,nnn), pp2(nnn,nnn)

      common /poten1/ gamm, c0, c3, cs, cl, wl

      common /swich0/ icoul, irelcr, iavoid
      common /swich1/ ipot, insys, irkg, icolt

      common /gradu0/ c0g, c3g, csg, pag

*-----------------------------------------------------------------------

      dimension       ffr(3,nnn), ffp(3,nnn), rh3d(nnn)

*-----------------------------------------------------------------------

               do 100 i = 1, massal
                     if( inds(i) .eq. 0 ) goto 100

                     rho3 = 0.0

                  do 110 j = 1, massal
                     if( inds(j) .eq. 0 ) goto 110

                     rho3 = rho3 + rha(j,i)

  110             continue

                     rh3d(i) = rho3 ** pag

  100          continue

*-----------------------------------------------------------------------

         do 200 i = 1, massal
                  if( inds(i) .eq. 0 ) goto 200

                  bix = p(1,i) / p(4,i)
                  biy = p(2,i) / p(4,i)
                  biz = p(3,i) / p(4,i)

                  ffr(1,i) = bix
                  ffr(2,i) = biy
                  ffr(3,i) = biz

                  ffp(1,i) = 0.0
                  ffp(2,i) = 0.0
                  ffp(3,i) = 0.0

*-----------------------------------------------------------------------

            do 210 j = 1, massal
                  if( inds(j) .eq. 0 ) goto 210

                  eij  = p(4,i) + p(4,j)

                  ccpp =  c0g * rha(j,i)
     &                 +  c3g * rha(j,i) * ( rh3d(j) + rh3d(i) )
     &                 +  csg * rha(j,i) * inuc(j) * inuc(i)
     &                        * ( 1. - 2. * abs(ichg(j)-ichg(i)) )
     &                 +  cl  * rhc(j,i)

                  grbb = - rbij(j,i)
                  ccrr =   grbb * ccpp / eij

                  rijx =   r(1,i) - r(1,j)
                  rijy =   r(2,i) - r(2,j)
                  rijz =   r(3,i) - r(3,j)

                  bijx = ( p(1,i) + p(1,j) ) / eij
                  bijy = ( p(2,i) + p(2,j) ) / eij
                  bijz = ( p(3,i) + p(3,j) ) / eij

                  cijx =   bijx - bix
                  cijy =   bijy - biy
                  cijz =   bijz - biz

               ffr(1,i) = ffr(1,i) + 2.0 * ccrr * ( rijx + grbb * cijx )
               ffr(2,i) = ffr(2,i) + 2.0 * ccrr * ( rijy + grbb * cijy )
               ffr(3,i) = ffr(3,i) + 2.0 * ccrr * ( rijz + grbb * cijz )

               ffp(1,i) = ffp(1,i) - 2.0 * ccpp * ( rijx + grbb * bijx )
               ffp(2,i) = ffp(2,i) - 2.0 * ccpp * ( rijy + grbb * bijy )
               ffp(3,i) = ffp(3,i) - 2.0 * ccpp * ( rijz + grbb * bijz )

  210       continue

  200    continue

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine rk12(dt)
*                                                                      *
*                                                                      *
*        Last Revised:     2004 09 23                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to propagate particles by 2th order Runge-Kutta-Gill    *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              dt       : time step                                    *
*                                                                      *
*----------------------------------------------------------------------*
*        comment:                                                      *
*                                                                      *
*              second order Runge-Kutta-Gill method                    *
*                                                                      *
*              z(t+h) = z(t)                                           *
*                     + h*{ c1*f(z(t)) + c2*f( z(t)+c3*h*f(z(t) ) }    *
*                                                                      *
*              c2 = 0.5  => modified euler method i                    *
*              c2 = 1    => modified euler method ii                   *
*              c2 = 0.75 => heun method                                *
*                                                                      *
*              c2 must be 0< c2 <=1                                    *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param01.inc'
      include 'param02.inc'

*-----------------------------------------------------------------------

      parameter ( c2 = 1.0 )
      parameter ( c1 = 1.0 - c2 )
      parameter ( c3 = 1.0 / 2.0 / c2 )

*-----------------------------------------------------------------------

      common /vriab0/ massal, massba, mmeson
      common /coodrp/ r(5,nnn),  p(5,nnn)
      common /coodid/ ichg(nnn), inuc(nnn), ibry(nnn), inds(nnn),
     &                inun(nnn), iavd(nnn), ihis(nnn)

      dimension       f0r(3,nnn), d1r(3,nnn),
     &                f0p(3,nnn), d1p(3,nnn)

*-----------------------------------------------------------------------

               dt3 = dt * c3
               dt1 = dt * ( c1 - c3 )
               dt2 = dt * c2

*-----------------------------------------------------------------------

               call gradu(d1r,d1p)

         do 20 i = 1, massal
               if( inds(i) .eq. 0 ) goto 20

            do j = 1, 3

               r(j,i) = r(j,i) + dt3 * d1r(j,i)
               p(j,i) = p(j,i) + dt3 * d1p(j,i)

               f0r(j,i) = d1r(j,i)
               f0p(j,i) = d1p(j,i)

            end do

               p(4,i) = sqrt( p(5,i)**2 + p(1,i)**2
     &                      + p(2,i)**2 + p(3,i)**2 )

   20    continue

               call caldisa

*-----------------------------------------------------------------------

               call gradu(d1r,d1p)

         do 30 i = 1, massal
               if( inds(i) .eq. 0 ) goto 30

            do j = 1, 3

               r(j,i) = r(j,i) + f0r(j,i) * dt1 + d1r(j,i) * dt2
               p(j,i) = p(j,i) + f0p(j,i) * dt1 + d1p(j,i) * dt2

            end do

               p(4,i) = sqrt( p(5,i)**2 + p(1,i)**2
     &                      + p(2,i)**2 + p(3,i)**2 )

   30    continue

                call caldisa


      return
      end


************************************************************************
*                                                                      *
      subroutine rkg4(dt)
*                                                                      *
*                                                                      *
*        Last Revised:     2004 09 23                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to propagate particles by 4th order Runge-Kutta-Gill    *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              dt     -  time step                                     *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param01.inc'
      include 'param02.inc'

*-----------------------------------------------------------------------

      common /vriab0/ massal, massba, mmeson
      common /coodrp/ r(5,nnn),  p(5,nnn)

      common /coodid/ ichg(nnn), inuc(nnn), ibry(nnn), inds(nnn),
     &                inun(nnn), iavd(nnn), ihis(nnn)

      dimension       qqr(3,nnn), qqp(3,nnn)

      dimension       ggr(3,nnn), ccr(3,nnn), ffr(3,nnn),
     &                ggp(3,nnn), ccp(3,nnn), ffp(3,nnn)

      dimension       fac1(4), fac2(4), fac3(4)

      save            qqr, qqp

*-----------------------------------------------------------------------

      data fac1 / 0.5, 0.29289322, 1.70710678, 0.1666666667 /
      data fac2 / 2.0, 1.0, 1.0, 2.0 /
      data fac3 / 0.5, 0.29289322, 1.70710678, 0.5 /

      data qqr / nnn*0.0, nnn*0.0, nnn*0.0 /
      data qqp / nnn*0.0, nnn*0.0, nnn*0.0 /

*-----------------------------------------------------------------------

            do 30 ik = 1, 4

               call gradu(ffr,ffp)

               fa1 = fac1(ik)
               fa2 = fac2(ik)
               fa3 = fac3(ik)

            do 20 i = 1, massal
               if( inds(i) .eq. 0 ) goto 20

            do 10 j = 1, 3

               ccr(j,i) = dt * ffr(j,i)
               ccp(j,i) = dt * ffp(j,i)

               ggr(j,i) = fa1 * ( ccr(j,i) - fa2 * qqr(j,i) )
               ggp(j,i) = fa1 * ( ccp(j,i) - fa2 * qqp(j,i) )

               qqr(j,i) = qqr(j,i) + 3.0 * ggr(j,i) - fa3 * ccr(j,i)
               qqp(j,i) = qqp(j,i) + 3.0 * ggp(j,i) - fa3 * ccp(j,i)

               r(j,i) = r(j,i) + ggr(j,i)
               p(j,i) = p(j,i) + ggp(j,i)

   10       continue

               p(4,i) = sqrt( p(5,i)**2 + p(1,i)**2
     &                      + p(2,i)**2 + p(3,i)**2 )

   20       continue

               call caldisa

   30       continue

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine pauli(i,ntag,phase)
*                                                                      *
*                                                                      *
*        Last Revised:     2004 09 23                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to calculate the Pauli blocking factor                  *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              i           : number of particle                        *
*              ntag        : flag which tells if phase-space is        *
*                            Pauli-blocked                             *
*                                                                      *
*                            ntag =  0 => phase space open             *
*                            ntag =  1 => phase space blocked          *
*                                                                      *
*              phase       : phase space factor                        *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param01.inc'
      include 'param02.inc'

*-----------------------------------------------------------------------

      common /vriab0/ massal, massba, mmeson

      common /coodrp/ r(5,nnn),  p(5,nnn)
      common /coodid/ ichg(nnn), inuc(nnn), ibry(nnn), inds(nnn),
     &                inun(nnn), iavd(nnn), ihis(nnn)
      common /cood2r/ rr2(nnn,nnn), rbij(nnn,nnn), pp2(nnn,nnn)

      common /swich0/ icoul, irelcr, iavoid
      common /swich1/ ipot, insys, irkg, icolt

      common /poten1/ gamm, c0, c3, cs, cl, wl

      common /pauli0/ cpw, cph, cpc

*-----------------------------------------------------------------------

               ntag  = 0
               phase = 0.0

               if( inuc(i) .ne. 1 ) return

*-----------------------------------------------------------------------

         do j = 1, massal

            if( ichg(j) .eq. ichg(i) .and.
     &          inuc(j) .eq. 1 ) then

                     expa   = - rr2(i,j) * cpw

               if( expa .gt. epsx ) then

                     expa   = expa - pp2(i,j) * cph

                  if( expa .gt. epsx ) then

                     phase  = phase + exp( expa )

                  end if

               end if

            end if

         end do

*-----------------------------------------------------------------------

               phase  = ( phase - 1.0 ) * cpc

               if( phase .gt. rn(0) ) ntag = 1

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine epotall(epot)
*                                                                      *
*                                                                      *
*        Last Revised:     2004 09 23                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to calculate total potential energy                     *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              epot        : total potential energy                    *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param01.inc'
      include 'param02.inc'

*-----------------------------------------------------------------------

      common /vriab0/ massal, massba, mmeson

      common /coodrp/ r(5,nnn),  p(5,nnn)
      common /coodid/ ichg(nnn), inuc(nnn), ibry(nnn), inds(nnn),
     &                inun(nnn), iavd(nnn), ihis(nnn)
      common /cood2b/ rha(nnn,nnn), rhe(nnn,nnn), rhc(nnn,nnn)

      common /poten1/ gamm, c0, c3, cs, cl, wl

*-----------------------------------------------------------------------

      dimension       rhoa(nnn), rho3(nnn), rhos(nnn), rhoc(nnn)

*-----------------------------------------------------------------------

               do i = 1, massal

                  rhoa(i) = 0.0
                  rho3(i) = 0.0
                  rhos(i) = 0.0
                  rhoc(i) = 0.0

               end do

               do 100 i = 1, massal
                  if( inds(i) .eq. 0 ) goto 100

               do 110 j = 1, massal
                  if( inds(j) .eq. 0 ) goto 110

                  rhoa(i) = rhoa(i) + rha(j,i)
                  rhoc(i) = rhoc(i) + rhe(j,i)
                  rhos(i) = rhos(i) + rha(j,i) * inuc(j) * inuc(i)
     &                    * ( 1. - 2. * abs(ichg(j)-ichg(i)) )

  110          continue
  100          continue

               do i = 1, massal

                  rho3(i) = rhoa(i) ** gamm

               end do

*-----------------------------------------------------------------------

                  epot1 = 0.0
                  epot3 = 0.0
                  epots = 0.0
                  epotc = 0.0

               do i = 1, massal

                  epot1 = epot1 + rhoa(i)
                  epot3 = epot3 + rho3(i)
                  epots = epots + rhos(i)
                  epotc = epotc + rhoc(i)

               end do

*-----------------------------------------------------------------------

                  epot = c0 * epot1
     &                 + c3 * epot3
     &                 + cs * epots
     &                 + cl * epotc

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine etotal(nsel,it,rs,ps,ekin,epot,ebin,emas,epin,jj)
*                                                                      *
*                                                                      *
*        Last Revised:     2004 09 23                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to calculate total energy of cluster                    *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              nsel        : 1; cm system of cluster,                  *
*                            0; reference frame                        *
*              it          : assignment of considered particles        *
*              rs          : coordinate in cm system of cluster        *
*              ps          : momentum   in cm system of cluster        *
*              ekin        : kinetic energy   / nucleon                *
*              epot        : potential energy / nucleon                *
*              ebin        : binding energy   / nucleon                *
*              emas        : rest mass        / nucleon                *
*              epin        : pion energy      / nucleon                *
*              jj          : angular momentum of the fragment          *
*                            only for ic = 0                           *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param01.inc'
      include 'param02.inc'

*-----------------------------------------------------------------------

      common /poten1/ gamm, c0, c3, cs, cl, wl

      common /coodrp/ r(5,nnn),  p(5,nnn)
      common /coodid/ ichg(nnn), inuc(nnn), ibry(nnn), inds(nnn),
     &                inun(nnn), iavd(nnn), ihis(nnn)
      common /cood2b/ rha(nnn,nnn), rhe(nnn,nnn), rhc(nnn,nnn)

*-----------------------------------------------------------------------

      dimension       it(0:nnn)
      dimension       rhoa(nnn), rho3(nnn), rhos(nnn), rhoc(nnn)
      dimension       rs(3,nnn), ps(3,nnn),  es(nnn)

*-----------------------------------------------------------------------

                  ipion = 0

*-----------------------------------------------------------------------
*        not cluster
*-----------------------------------------------------------------------

            if( it(0) .eq. 0 ) then

                  ekin = 0.0
                  epot = 0.0
                  ebin = 0.0
                  emas = 0.0
                  epin = 0.0
                  jj   = 0

                  return

            else if( it(0) .eq. 1 ) then

                  ekin = p(4,it(1)) - p(5,it(1))
                  epot = 0.0
                  ebin = 0.0
                  emas = p(5,it(1))
                  jj   = 0

               if( inds(it(1)) .eq. 4 ) then

                  epin = p(4,it(1))

               else

                  epin = 0.0

               end if

               return

            end if

*-----------------------------------------------------------------------
*        potential energy
*-----------------------------------------------------------------------

            do ii = 1, it(0)

                  i = it(ii)

                  if( inds(i) .eq. 4 ) ipion = ipion + 1

                  rhoa(i) = 0.0
                  rho3(i) = 0.0
                  rhos(i) = 0.0
                  rhoc(i) = 0.0

               do jj = 1, it(0)

                  j = it(jj)

                  rhoa(i) = rhoa(i) + rha(j,i)
                  rhoc(i) = rhoc(i) + rhe(j,i)
                  rhos(i) = rhos(i) + rha(j,i) * inuc(j) * inuc(i)
     &                    * ( 1. - 2. * abs(ichg(j)-ichg(i)) )

               end do

                  rho3(i) = rhoa(i) ** gamm

                  es(i)   = p(4,i)

            end do

*-----------------------------------------------------------------------
*        Lorentz transform to the c.m. frame of the cluster
*-----------------------------------------------------------------------

                  jj = 0

                  if( nsel .eq. 1 )  call pcmcl(it,rs,ps,es,jj)

*-----------------------------------------------------------------------
*        potential, kinetic and pion energy
*-----------------------------------------------------------------------

                  epot1 = 0.0
                  epot3 = 0.0
                  epots = 0.0
                  epotc = 0.0

                  ekin  = 0.0
                  emas  = 0.0
                  epin  = 0.0

            do ii = 1, it(0)

                  i = it(ii)

               if( inds(i) .eq. 4 ) then

                  epin = epin + es(i)

               else

                  emas = emas + p(5,i)
                  ekin = ekin + es(i) - p(5,i)

               end if

                  epot1 = epot1 + rhoa(i)
                  epot3 = epot3 + rho3(i)
                  epots = epots + rhos(i)
                  epotc = epotc + rhoc(i)

            end do

*-----------------------------------------------------------------------

                  epot = c0 * epot1
     &                 + c3 * epot3
     &                 + cs * epots
     &                 + cl * epotc

                  emas =   emas          / float( it(0) - ipion )
                  ebin = ( ekin + epot ) / float( it(0) - ipion )
                  ekin =   ekin          / float( it(0) - ipion )
                  epot =   epot          / float( it(0) - ipion )
                  epin =   epin          / float( it(0) - ipion )

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine pcmcl(it,rs,ps,es,jj)
*                                                                      *
*                                                                      *
*        Last Revised:     2004 09 23                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to calculate the energy and angular momentum            *
*              in the rest frame of the cluster                        *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              it     - assignment of considered particles             *
*              rs     - cm coordinate of particle                      *
*              ps     - cm momentum of particle                        *
*              es     - cm energy of particle                          *
*              jj     - angular momentum of the fragment               *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param01.inc'
      include 'param02.inc'

*-----------------------------------------------------------------------

      common /coodrp/ r(5,nnn),  p(5,nnn)
      common /coodid/ ichg(nnn), inuc(nnn), ibry(nnn), inds(nnn),
     &                inun(nnn), iavd(nnn), ihis(nnn)

*-----------------------------------------------------------------------

      dimension       it(0:nnn)
      dimension       rs(3,nnn), ps(3,nnn),  es(nnn)

*-----------------------------------------------------------------------

                  betax = 0.0
                  betay = 0.0
                  betaz = 0.0

                  pcmx  = 0.0
                  pcmy  = 0.0
                  pcmz  = 0.0

                  esum  = 0.0

               do ii = 1, it(0)

                  i = it(ii)

                  pcmx  = pcmx + p(1,i)
                  pcmy  = pcmy + p(2,i)
                  pcmz  = pcmz + p(3,i)

                  esum  = esum + p(4,i)

               end do

                  gamma = esum
     &                  / sqrt( esum**2 - pcmx**2 - pcmy**2 - pcmz**2 )

                  betax = pcmx / esum
                  betay = pcmy / esum
                  betaz = pcmz / esum


                  pcm1 = 0
                  pcm2 = 0
                  pcm3 = 0

               do ii = 1, it(0)

                  i = it(ii)

                  transp =  gamma / ( gamma + 1.0 )
     &                   * ( p(1,i) * betax
     &                     + p(2,i) * betay
     &                     + p(3,i) * betaz )


                  ps(1,i) = p(1,i) - betax * transp
                  ps(2,i) = p(2,i) - betay * transp
                  ps(3,i) = p(3,i) - betaz * transp

                  pcm1 = pcm1 + ps(1,i)
                  pcm2 = pcm2 + ps(2,i)
                  pcm3 = pcm3 + ps(3,i)

               end do

                  pcm1 = pcm1 / float(it(0))
                  pcm2 = pcm2 / float(it(0))
                  pcm3 = pcm3 / float(it(0))


               do ii = 1, it(0)

                  i = it(ii)

                  ps(1,i) = ps(1,i) - pcm1
                  ps(2,i) = ps(2,i) - pcm2
                  ps(3,i) = ps(3,i) - pcm3

               end do


                  tmass = 0.0
                  cmrx  = 0.0
                  cmry  = 0.0
                  cmrz  = 0.0


               do ii = 1, it(0)

                  i = it(ii)

                  transr =  gamma *  gamma / ( gamma + 1.0 )
     &                   * ( r(1,i) * betax
     &                     + r(2,i) * betay
     &                     + r(3,i) * betaz )

                  es(i) = sqrt( p(5,i)**2
     &                        + ps(1,i)**2
     &                        + ps(2,i)**2
     &                        + ps(3,i)**2 )

                  rs(1,i) = r(1,i) + betax * transr
                  rs(2,i) = r(2,i) + betay * transr
                  rs(3,i) = r(3,i) + betaz * transr

                  cmrx = cmrx + rs(1,i) * es(i)
                  cmry = cmry + rs(2,i) * es(i)
                  cmrz = cmrz + rs(3,i) * es(i)

                  tmass = tmass + es(i)

               end do

                  cmrx = cmrx / tmass
                  cmry = cmry / tmass
                  cmrz = cmrz / tmass


               do ii = 1, it(0)

                  i = it(ii)

                  rs(1,i) = rs(1,i) - cmrx
                  rs(2,i) = rs(2,i) - cmry
                  rs(3,i) = rs(3,i) - cmrz

               end do


*-----------------------------------------------------------------------
*        Angluar momentum
*-----------------------------------------------------------------------

                  rlx = 0.0
                  rly = 0.0
                  rlz = 0.0

               do ii = 1, it(0)

                  i = it(ii)

                  rlx = rlx + rs(2,i) * ps(3,i)
     &                      - rs(3,i) * ps(2,i)
                  rly = rly + rs(3,i) * ps(1,i)
     &                      - rs(1,i) * ps(3,i)
                  rlz = rlz + rs(1,i) * ps(2,i)
     &                      - rs(2,i) * ps(1,i)

               end do

                  jj = nint( 
     &                 sqrt( rlx**2 + rly**2 + rlz**2 ) / hbc )


*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine cldist
*                                                                      *
*                                                                      *
*        Last Revised:     2004 09 23                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to determine nuclear cluster.                           *
*                                                                      *
*        Variables:  in common                                         *
*                                                                      *
*              common /clusti/ itc(0:nnn,0:nnn)                        *
*              common /clustf/ nclst, iclust(nnn),                     *
*              common /clustg/ jclust(0:7,nnn),  qclust(0:12,nnn)      *
*                                                                      *
*              itc(0,0)    : total number of clusters plus particles   *
*              itc(i,0)    : number  of baryon in i-th cluster         *
*              itc(i,j)    : j-th id of baryon in i-th cluster         *
*                                                                      *
*              nclst       : total number of clusters plus particles   *
*                                                                      *
*              iclust(i)                                               *
*                                                                      *
*                   i = 0, nucleus                                     *
*                     = 1, proton                                      *
*                     = 2, neutron                                     *
*                     = 3, pion                                        *
*                     = 4, Gamma                                       *
*                     = 5, kaon                                        *
*                     = 6, muon                                        *
*                     = 7, others                                      *
*                                                                      *
*              jclust(k,i)                                             *
*                                                                      *
*                   i = 0, angular momentum                            *
*                     = 1, proton number                               *
*                     = 2, neutron number                              *
*                     = 3, ip, see below                               *
*                     = 4, status of the particle 0: real, <0 : dead   *
*                     = 5, charge                                      *
*                     = 6, baryon number                               *
*                     = 7, kf code                                     *
*                                                                      *
*              qclust(k,i)                                             *
*                                                                      *
*                   i = 0, impact parameter                            *
*                     = 1, px (GeV/c)                                  *
*                     = 2, py (GeV/c)                                  *
*                     = 3, pz (GeV/c)                                  *
*                     = 4, etot = sqrt( p**2 + rm**2 ) (GeV)           *
*                     = 5, rest mass (GeV)                             *
*                     = 6, excitation energy (MeV)                     *
*                     = 7, kinetic energy (MeV)                        *
*                     = 8, weight change                               *
*                     = 9, delay time                                  *
*                     = 10, x-displace                                 *
*                     = 11, y-displace                                 *
*                     = 12, z-displace                                 *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param01.inc'
      include 'param02.inc'

*-----------------------------------------------------------------------

      common /vriab0/ massal, massba, mmeson
      common /vriab1/ b, llnow, ntnow

      common /cood2b/ rha(nnn,nnn), rhe(nnn,nnn), rhc(nnn,nnn)
      common /cood2r/ rr2(nnn,nnn), rbij(nnn,nnn), pp2(nnn,nnn)
      common /coodrp/ r(5,nnn),  p(5,nnn)
      common /coodid/ ichg(nnn), inuc(nnn), ibry(nnn), inds(nnn),
     &                inun(nnn), iavd(nnn), ihis(nnn)

      common /poten1/ gamm, c0, c3, cs, cl, wl

      common /clusti/ itc(0:nnn,0:nnn)
      common /clustf/ nclst, iclust(nnn)
      common /clustg/ jclust(0:7,nnn),  qclust(0:12,nnn)

      common /clustp/ rumpat(0:20), numpat(0:20)

*-----------------------------------------------------------------------

      common/jamevnt1/v(5,nnn),k(11,nnn)

*-----------------------------------------------------------------------

      dimension       mascl(nnn), num(nnn)
      dimension       isort(nnn), isorti(nnn)
      dimension       rhoa(nnn)
      dimension       it(0:nnn), rs(3,nnn), ps(3,nnn)

*-----------------------------------------------------------------------
*        Cluster distance
*-----------------------------------------------------------------------

               cpf2 = ( 1.5 * pi**2
     &              * ( 4.0 * pi * wl )**(-1.5) )**(2./3.)
     &              * hbc ** 2

               rcc2 = rclds ** 2

*-----------------------------------------------------------------------
*        Calculate overlap of the wave packets
*-----------------------------------------------------------------------

            do i = 1, massal

                     rhoa(i) = 0.0d00

               if( inds(i) .ne. 0 ) then

                  do j = 1, massal

                     rhoa(i) = rhoa(i) + rha(i,j)

                  end do

                     rhoa(i) = ( rhoa(i) + 1.0 ) ** (1./3.)

               end if

            end do

*-----------------------------------------------------------------------
*     identification of the cluster
*-----------------------------------------------------------------------

               do i = 1, nnn

                  mascl(i) = 1
                  num(i)   = i

               end do

                  nclst = 1
                  ichek = 1

         do 100 i = 1, massal - 1

                  j1 = ichek + 1
                  id1 = inds( num(i) )

            do 110 j = j1, massal

                  id2 = inds( num(j) )

                  rdist2 = rr2( num(i), num(j) )
                  pdist2 = pp2( num(i), num(j) )
                  pcc2   = cpf2
     &                   * ( rhoa( num(i) ) + rhoa( num(j) ) )**2

               if( rdist2 .lt. rcc2 .and. pdist2 .lt. pcc2 .and.
     &             id1 .eq. 1 .and. id2 .eq. 1 ) then

                  ibuf             = num( ichek + 1 )
                  num( ichek + 1 ) = num( j )
                  num( j )         = ibuf
                  ichek            = ichek + 1
                  mascl( nclst )   = mascl( nclst ) + 1

               end if

  110       continue

               if( ichek .eq. i ) then

                  nclst = nclst + 1
                  ichek = ichek + 1

               end if

  100    continue

*-----------------------------------------------------------------------
*        sort for summary
*-----------------------------------------------------------------------

               do i = 1, nclst

                  isort(i) = i

               end do

  510    continue

                  nexch = 0

            do i = 1, nclst - 1

               if( mascl(isort(i)) .lt. mascl(isort(i+1)) ) then

                  nexch      = nexch + 1
                  ibuf       = isort(i)
                  isort(i)   = isort(i+1)
                  isort(i+1) = ibuf

               end if

            end do

            if( nexch .ne. 0 ) goto 510

               do i = 1, nclst

                  isorti(isort(i)) = i

               end do

                  itc(0,0) = nclst

                  inum = 0

            do i = 1, nclst

                  itc(isorti(i),0) = mascl(i)

               do j = 1, mascl(i)

                  inum = inum +1
                  itc(isorti(i),j) = num(inum)

               end do

            end do


*-----------------------------------------------------------------------
*     Loop over clusters
*     Summarry of the cluster
*-----------------------------------------------------------------------

                  mnucl = 0
                  mprot = 0
                  mneut = 0
                  mpipo = 0
                  mpine = 0
                  mping = 0
                  mmupo = 0
                  mmune = 0
                  mkapo = 0
                  mkane = 0
                  mkang = 0
                  mothe = 0
                  mgamm = 0

                  ii = 0

*-----------------------------------------------------------------------

      do 1000 i = 1, nclst

*-----------------------------------------------------------------------

                     mclst = itc(i,0)

                     nchpa = 0

                     nprot = 0
                     nnuet = 0
                     ndelt = 0
                     nstar = 0

                     jj    = 0
                     jc    = 0
                     js    = 0

*-----------------------------------------------------------------------
*     This is cluster
*-----------------------------------------------------------------------

      if( mclst .gt. 1 ) then

                     it(0) = mclst

                     ipcst = 0

                     mnucl = mnucl + 1

*-----------------------------------------------------------------------
*           Momentum of cluster.
*-----------------------------------------------------------------------

                     pclx = 0.0
                     pcly = 0.0
                     pclz = 0.0

*-----------------------------------------------------------------------
*           Loop over mass of one cluster,
*           determine property of cluster.
*-----------------------------------------------------------------------

            do j = 1, mclst

                     it(j) = itc(i,j)
                     jp    = it(j)

                     pclx  = pclx + p(1,jp)
                     pcly  = pcly + p(2,jp)
                     pclz  = pclz + p(3,jp)

                  if( ichg(jp) .eq. 1 ) then

                     nprot = nprot + 1

                  else if( ichg(jp) .eq. 0 ) then

                     nnuet = nnuet + 1

                  end if

                     nchpa = nchpa + ichg(jp)

            end do

*-----------------------------------------------------------------------
*           Calculate binding energy of cluster.
*           Exciation energy and mass of the cluster.
*-----------------------------------------------------------------------

         if( nprot .gt. 0 .and. nnuet .gt. 0 ) then

               call etotal(1,it,rs,ps,ekin,epot,ebin,emas,epin,jj)

                     texc = max( 0.d0, ebin * float( nprot + nnuet )
     &                         + bindeg( nprot, nnuet ) / 1000.0 )

                     tmas = ( emas + ebin ) * it(0)
                     pabs = sqrt( pclx**2 + pcly**2 + pclz**2 )
                     etot = sqrt( tmas**2
     &                           + pclx**2 + pcly**2 + pclz**2 )

                     kf = nprot * 1000000 + nprot + nnuet
                     ibary = nprot + nnuet

*-----------------------------------------------------------------------

                     ii = ii + 1

                     iclust(ii)    = ipcst

                     jclust(0,ii)  = jj
                     jclust(1,ii)  = nprot
                     jclust(2,ii)  = nnuet
                     jclust(3,ii)  = 19
                     jclust(4,ii)  = 0
                     jclust(5,ii)  = nchpa
                     jclust(6,ii)  = ibary
                     jclust(7,ii)  = kf

                     qclust(0,ii)  = b
                     qclust(1,ii)  = pclx
                     qclust(2,ii)  = pcly
                     qclust(3,ii)  = pclz
                     qclust(4,ii)  = etot
                     qclust(5,ii)  = tmas
                     qclust(6,ii)  = texc * 1000.0
                     qclust(7,ii)  = ( etot - tmas ) * 1000.0
                     qclust(8,ii)  = 1.0
                     qclust(9,ii)  = 0.0
                     qclust(10,ii) = 0.0d0
                     qclust(11,ii) = 0.0d0
                     qclust(12,ii) = 0.0d0

*-----------------------------------------------------------------------
*        strange nucleus
*-----------------------------------------------------------------------

         else

            do j = 1, mclst

                     ip = itc(i,j)

                  if( ichg(ip) .eq. 1 ) then

                     nprot = 1
                     nnuet = 0
                     ipcst = 1
                     ippad = 1
                     kf    = 2212

                     mprot = mprot + 1

                  else if( ichg(ip) .eq. 0 ) then

                     nprot = 0
                     nnuet = 1
                     ipcst = 2
                     ippad = 2
                     kf    = 2112

                     mneut = mneut + 1

                  end if

                     pclx = p(1,ip)
                     pcly = p(2,ip)
                     pclz = p(3,ip)
                     pabs = sqrt( pclx**2 + pcly**2 + pclz**2 )

                     etot = p(4,ip)
                     tmas = p(5,ip)

                     texc = 0.0

                     ii = ii + 1

                     iclust(ii)    = ipcst

                     jclust(0,ii)  = jj
                     jclust(1,ii)  = nprot
                     jclust(2,ii)  = nnuet
                     jclust(3,ii)  = ippad
                     jclust(4,ii)  = 0
                     jclust(5,ii)  = nprot
                     jclust(6,ii)  = 1
                     jclust(7,ii)  = kf

                     qclust(0,ii)  = b
                     qclust(1,ii)  = pclx
                     qclust(2,ii)  = pcly
                     qclust(3,ii)  = pclz
                     qclust(4,ii)  = etot
                     qclust(5,ii)  = tmas
                     qclust(6,ii)  = texc * 1000.0
                     qclust(7,ii)  = ( etot - tmas ) * 1000.0
                     qclust(8,ii)  = 1.0
                     qclust(9,ii)  = 0.0
                     qclust(10,ii) = 0.0d0
                     qclust(11,ii) = 0.0d0
                     qclust(12,ii) = 0.0d0

            end do

         end if

*-----------------------------------------------------------------------
*     This is a baryon or meson
*-----------------------------------------------------------------------

      else

*-----------------------------------------------------------------------

                     ip  = itc(i,1)

               if( inds(ip) .eq. 1 ) then

                  if( ichg(ip) .eq. 1 ) then

                     nprot = 1
                     nnuet = 0
                     ipcst = 1
                     ippad = 1
                     ibary = 1
                     kf    = 2212

                     mprot = mprot + 1

                  else if( ichg(ip) .eq. 0 ) then

                     nprot = 0
                     nnuet = 1
                     ipcst = 2
                     ippad = 2
                     ibary = 1
                     kf    = 2112

                     mneut = mneut + 1

                  end if

               else if( inds(ip) .eq. 2 ) then

                     ndelt = 1
                     ippad = 11
                     ibary = 1
                     nprot = 0
                     nnuet = 0

                  if( ichg(ip) .eq. 2 ) then

                     kf    = 2224

                  else if( ichg(ip) .eq.  1 ) then

                     kf    = 2214

                  else if( ichg(ip) .eq.  0 ) then

                     kf    = 2114

                  else if( ichg(ip) .eq. -1 ) then

                     kf    = 1114

                  end if

                     ipcst = 7
                     mothe = mothe + 1

               else if( inds(ip) .eq. 3 ) then

                     nstar = 1
                     ippad = 11
                     ibary = 1
                     nprot = 0
                     nnuet = 0

                  if( ichg(ip) .eq. 1 ) then

                     kf    = 12212

                  else if( ichg(ip) .eq. 0 ) then

                     kf    = 12112

                  end if

                     ipcst = 7
                     mothe = mothe + 1

               else if( inds(ip) .eq. 4 ) then

                  if( ichg(ip) .eq. 1 ) then

                     kf    = 211
                     ippad = 3
                     mpipo = mpipo + 1

                  else if( ichg(ip) .eq.  0 ) then

                     kf    = 111
                     ippad = 4
                     mpine = mpine + 1

                  else if( ichg(ip) .eq. -1 ) then

                     kf    = -211
                     ippad = 5
                     mping = mping + 1

                  end if

                     ipcst = 3
                     ibary = 0
                     nprot = 0
                     nnuet = 0

*-----------------------------------------------------------------------

               else if( inds(ip) .eq. 5 ) then

                     kf    = k(2,ip)
                     ibary = k(9,ip) / 3
                     nprot = 0
                     nnuet = 0

                  if( kf .eq. 22 ) then

                        ipcst = 4
                        ippad = 14

                        mgamm = mgamm + 1

                  else if( kf .eq. -321 .or.
     &                     kf .eq.  311 .or.
     &                     kf .eq.  321 ) then

                        ipcst = 5

                     if( kf .eq. 321 ) then

                        mkapo = mkapo + 1
                        ippad = 8

                     else if( kf .eq. 311 ) then

                        mkane = mkane + 1
                        ippad = 9

                     else if( kf .eq. -321 ) then

                        mkang = mkang + 1
                        ippad = 10

                     end if

                  else if( kf .eq. -13 .or.
     &                     kf .eq.  13 ) then

                        ipcst = 6

                     if( kf .eq. 13 ) then

                        mmupo = mmupo + 1
                        ippad = 6

                     else if( kf .eq. -13 ) then

                        mmune = mmune + 1
                        ippad = 7

                     end if

                  else

                     mothe = mothe + 1

                        ipcst = 7
                        ippad = 11

                  end if

*-----------------------------------------------------------------------

               end if

*-----------------------------------------------------------------------

                     pclx = p(1,ip)
                     pcly = p(2,ip)
                     pclz = p(3,ip)
                     pabs = sqrt( pclx**2 + pcly**2 + pclz**2 )

                     etot = p(4,ip)
                     tmas = p(5,ip)

                     texc = 0.0

                     nchpa = ichg(ip)

*-----------------------------------------------------------------------

                     ii = ii + 1

                     iclust(ii)    = ipcst

                     jclust(0,ii)  = jj
                     jclust(1,ii)  = nprot
                     jclust(2,ii)  = nnuet
                     jclust(3,ii)  = ippad
                     jclust(4,ii)  = 0
                     jclust(5,ii)  = nchpa
                     jclust(6,ii)  = ibary
                     jclust(7,ii)  = kf

                     qclust(0,ii)  = b
                     qclust(1,ii)  = pclx
                     qclust(2,ii)  = pcly
                     qclust(3,ii)  = pclz
                     qclust(4,ii)  = etot
                     qclust(5,ii)  = tmas
                     qclust(6,ii)  = texc * 1000.0
                     qclust(7,ii)  = ( etot - tmas ) * 1000.0
                     qclust(8,ii)  = 1.0
                     qclust(9,ii)  = 0.0
                     qclust(10,ii) = 0.0d0
                     qclust(11,ii) = 0.0d0
                     qclust(12,ii) = 0.0d0

*-----------------------------------------------------------------------

      end if

 1000 continue

                     nclst = ii

*-----------------------------------------------------------------------

                     numpat(0)  = mnucl
                     numpat(1)  = mprot
                     numpat(2)  = mneut
                     numpat(3)  = mpipo
                     numpat(4)  = mpine
                     numpat(5)  = mping
                     numpat(6)  = mmupo
                     numpat(7)  = mmune
                     numpat(8)  = mkapo
                     numpat(9)  = mkane
                     numpat(10) = mkang
                     numpat(11) = mothe
                     numpat(14) = mgamm

                     rumpat(0)  = mnucl
                     rumpat(1)  = mprot
                     rumpat(2)  = mneut
                     rumpat(3)  = mpipo
                     rumpat(4)  = mpine
                     rumpat(5)  = mping
                     rumpat(6)  = mmupo
                     rumpat(7)  = mmune
                     rumpat(8)  = mkapo
                     rumpat(9)  = mkane
                     rumpat(10) = mkang
                     rumpat(11) = mothe
                     rumpat(14) = mgamm

*-----------------------------------------------------------------------

      return
      end


