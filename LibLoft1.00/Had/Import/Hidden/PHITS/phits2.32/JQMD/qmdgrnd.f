************************************************************************
*                                                                      *
*        PART 3: Ground State and Boost                                *
*                                                                      *
*                                                                      *
*   List of subprograms in rough order of relevance with main purpose  *
*      ( s = subroutine, f = function, b = block data, e = entry )     *
*                                                                      *
*                                                                      *
*  s  ground    to make the ground state                               *
*  s  packqmd   to make the ground state by random packing method      *
*  s  gcmang    to kill cm motion and angular momentum of the nucleus  *
*  s  rboost    to boost the ground state                              *
*  s  bcoul     to determine initial position and momentum             *
*                  according to coulomb trajectory                     *
*                                                                      *
*                                                                      *
************************************************************************
*                                                                      *
************************************************************************
*                                                                      *
      subroutine ground(ierr)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              make the goround state of target and projectile         *
*              for QMD mode using the random packing method            *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param01.inc'
      include 'param02.inc'

*-----------------------------------------------------------------------

      common /const0/ idnta, idnpr, massta, masspr, mstapr, msprpr

      common /const7/ betpr, gampr, prmas, radpr
      common /const8/ betta, gamta, tamas, radta

      common /vriab0/ massal, massba, mmeson

      common /coodrp/ r(5,nnn),  p(5,nnn)
      common /coodid/ ichg(nnn), inuc(nnn), ibry(nnn), inds(nnn),
     &                inun(nnn), iavd(nnn), ihis(nnn)

      dimension       r0(5,nnn),  p0(5,nnn)
      dimension       ichg0(nnn), inuc0(nnn), ibry0(nnn), inds0(nnn),
     &                inun0(nnn), iavd0(nnn), ihis0(nnn)


*-----------------------------------------------------------------------


            ierr = 0
*-----------------------------------------------------------------------
*     Loop over target and projectile
*-----------------------------------------------------------------------


      do 10000 itp = 1, 2


*-----------------------------------------------------------------------
*           Set target and projectile quantities
*-----------------------------------------------------------------------

            if( itp .eq. 1 ) then

                  massal = massta
                  msprgr = mstapr
                  massba = massal

                  ichtp  = msprgr
                  tpmas  = tamas
                  gambb  = gamta
                  idntp  = idnta
                  iavtp  = 1

                  iofset = 0

            else if( itp .eq. 2 ) then

                  massal = masspr
                  msprgr = msprpr
                  massba = massal

                  ichtp  = msprgr
                  tpmas  = prmas
                  gambb  = gampr
                  idntp  = idnpr
                  iavtp  = -1

                  iofset = massta

            end if

*-----------------------------------------------------------------------
*           Nucleon or Pion
*-----------------------------------------------------------------------

               if( idntp .eq. 0 ) then

                  inutp = 1
                  ibrtp = 1
                  indtp = 1

               else if( idntp .ne. 0 ) then

                  inutp = 0
                  ibrtp = 0
                  indtp = 4

               end if

*-----------------------------------------------------------------------

         if( massal .eq. 0 ) goto 10000

*-----------------------------------------------------------------------
*        Set all coordinate of each particle
*-----------------------------------------------------------------------

         do 500 i = 1, massal

*-----------------------------------------------------------------------
*              Charge of particle
*-----------------------------------------------------------------------

               if( massal .gt. 1 ) then

                  if( i .le. msprgr ) then

                     ichtp = 1

                  else

                     ichtp = 0

                  end if

               end if

*-----------------------------------------------------------------------
*              Set variables
*-----------------------------------------------------------------------

                  do m = 1, 3

                     r(m,i)  = 0.0
                     p(m,i)  = 0.0

                  end do

                     p(5,i)  = tpmas
                     p(4,i)  = p(5,i)

                     ichg(i) = ichtp
                     inuc(i) = inutp
                     ibry(i) = ibrtp
                     inds(i) = indtp
                     inun(i) = 0
                     iavd(i) = iavtp
                     ihis(i) = 0

*-----------------------------------------------------------------------

  500    continue


*-----------------------------------------------------------------------
*        Set positions and momenta of nucleons
*-----------------------------------------------------------------------


               call packqmd(massal,msprgr,gambb,ierr)

                  if( ierr .ne. 0 ) return

*-----------------------------------------------------------------------
*        Now ground state was successfully created !
*        Store these value
*-----------------------------------------------------------------------

            do i = 1, massal

                  j = i + iofset

               do m = 1, 3

                  r0(m,j) = r(m,i)
                  p0(m,j) = p(m,i)

               end do

                  p0(4,j) = p(4,i)
                  p0(5,j) = p(5,i)

                  ichg0(j) = ichg(i)
                  inuc0(j) = inuc(i)
                  ibry0(j) = ibry(i)
                  inds0(j) = inds(i)
                  inun0(j) = inun(i)
                  iavd0(j) = iavd(i)
                  ihis0(j) = ihis(i)

            end do

*-----------------------------------------------------------------------
*     End loop over proj. and targ.
*-----------------------------------------------------------------------

10000 continue

*-----------------------------------------------------------------------
*     reset total mass number and meson number
*-----------------------------------------------------------------------

                  massal = massta + masspr
                  massba = massal
                  mmeson = 0

            if( massta .eq. 1 .and. idnta .ne. 0 ) then

                  mmeson = mmeson + 1
                  massba = massba - 1

            end if

            if( masspr .eq. 1 .and. idnpr .ne. 0 ) then

                  mmeson = mmeson + 1
                  massba = massba - 1

            end if

*-----------------------------------------------------------------------
*        reset r(m,i), p(m,i) and ichg(i), inuc(i), ibry(i), inds(i)
*        iavd(i), and set inun(i), ihis(i)
*-----------------------------------------------------------------------

            do i = 1, massal

               do m = 1, 3

                  r(m,i)  = r0(m,i)
                  p(m,i)  = p0(m,i)

               end do

                  p(4,i)  = p0(4,i)
                  p(5,i)  = p0(5,i)

                  ichg(i) = ichg0(i)
                  inuc(i) = inuc0(i)
                  ibry(i) = ibry0(i)
                  inds(i) = inds0(i)
                  inun(i) = inun0(i)
                  iavd(i) = iavd0(i)
                  ihis(i) = ihis0(i)

            end do

*-----------------------------------------------------------------------
*        set remaining parameter to be zero
*-----------------------------------------------------------------------

            do i = massal + 1, nnn

               do m = 1, 3

                  r(m,i)  = 0.0
                  p(m,i)  = 0.0

               end do

                  p(4,i)  = 0.0
                  p(5,i)  = 0.0

                  ichg(i) = 0
                  inuc(i) = 0
                  ibry(i) = 0
                  inds(i) = 0
                  iavd(i) = 0
                  inun(i) = 0
                  ihis(i) = 0

            end do


*-----------------------------------------------------------------------

      return
       end subroutine



************************************************************************
*                                                                      *
      subroutine packqmd(massal,msprgr,gambb0,ierr)
*                                                                      *
*                                                                      *
*        Last Revised:     2004 03 04                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to sample positions and momenta of nucleus              *
*              according to the Woods-Saxon type distribution.         *
*              The binding energy is cheked with the mass formula.     *
*              If the binding energy per nucleaon lies                 *
*              within E_bin=+-0.5MeV, we adopt this configuration.     *
*              Pauli principle is also checked.                        *
*                                                                      *
*              If ipchs ( mstq1(90) ) is 1, the binding                *
*              energy is adjusted.                                     *
*              If ipchs ( mstq1(190) ) is 1, the binding               *
*              energy is adjusted in moving frame                      *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              massal      : total mass                                *
*              msprgr      : number of  proton                         *
*              gambb0      : gamma factor of nucleus                   *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param01.inc'
      include 'param02.inc'

*-----------------------------------------------------------------------

      parameter     ( edepth = 0.0 )

      parameter     ( epse = 0.000001 )

*-----------------------------------------------------------------------

      common /const2/ dt, ntmax, iprun, iprun0

      common /vriab1/ b, llnow, ntnow

      common /coodrp/ r(5,nnn),  p(5,nnn)
      common /coodid/ ichg(nnn), inuc(nnn), ibry(nnn), inds(nnn),
     &                inun(nnn), iavd(nnn), ihis(nnn)

      common /cood2b/ rha(nnn,nnn), rhe(nnn,nnn), rhc(nnn,nnn)
      common /cood2r/ rr2(nnn,nnn), rbij(nnn,nnn), pp2(nnn,nnn)

      common /poten1/ gamm, c0, c3, cs, cl, wl

      common /grndc0/ dsam, ddif, dsam2, ddif2
      common /grndc1/ cdp, c0p, c3p, csp, clp
      common /grndc2/ r00, r01, saa, rada, radb
      common /grndc3/ ipchs, mntry, dtg, fric

      common /pauli0/ cpw, cph, cpc

      common /rannum/ iseed, iseed0, iseed1
      common /input1/ mstq1(mxpa1), parq1(mxpa1)

*-----------------------------------------------------------------------

      dimension       rhol(nnn), dpot(nnn)
      dimension       rhoa(nnn), rhos(nnn), rhoc(nnn)
      dimension       phase0(nnn), phase(nnn)

      dimension       d1r(3,nnn), d1p(3,nnn)

*-----------------------------------------------------------------------

            ierr = 0
            ieo  = 0
c !!!!!!!!!!!!!!!!!
c            write(*,*) ' now ground'
c            write(*,*) ' mntry=', mntry
c            call epcurrent
c!!!!!!!!!!!!!!!!

*-----------------------------------------------------------------------
*           Maxmum number of baryons for QMD
*-----------------------------------------------------------------------

               if( massal .le. 1 ) then

                  return

               else if( massal .gt. nnn ) then

                  write(ieo,'('' Error: nnn is smaller then massal'')')
                  write(ieo,'('' ====='')')

                  call parastop( 222 )

               end if

*-----------------------------------------------------------------------
*     Woods-Saxon radius (rt00), cut-off radius (radm)
*     and Maxmum value (rmax).
*-----------------------------------------------------------------------

               rad0 = r00  * float( massal )**(1./3.)
               rt00 = rad0 - r01
               radm = rad0 - rada * ( gamm - 1.0 ) + radb

               rmax = 1.0 / ( 1.0 + exp( - rt00 / saa ) )

*-----------------------------------------------------------------------
*     Binding energy to be fitted from mass formula.
*-----------------------------------------------------------------------

               ebini  = - bindeg( msprgr, massal - msprgr )
     &                / float( massal )
               ebin00 =   ebini * 0.001

            if( massal .gt. 4 ) then

               ebin0  = ( ebini - 1.5 ) * 0.001
               ebin1  = ( ebini + 1.5 ) * 0.001

            else

               ebin0  = ( ebini - 2.5 ) * 0.001
               ebin1  = ( ebini + 2.5 ) * 0.001

            end if

*-----------------------------------------------------------------------
*     (1) Sampling position of nucleons
*-----------------------------------------------------------------------

               ntry = 0

 1000    continue

               ntry = ntry + 1

               if( ntry .gt. mntry / 2 ) then

                  ebin0  = ( ebini - 5.0 ) * 0.001
                  ebin1  = ( ebini + 5.0 ) * 0.001

               end if

               if( ntry .gt. mntry - 5 ) then

                  ebin0  = ( ebini - 10.0 ) * 0.001
                  ebin1  = ( ebini + 10.0 ) * 0.001

               end if

               if( ntry .gt. mntry ) then

                  write(ieo,'('' Error: Ground state cannot be'',
     &                        '' created. Try again with another'',
     &                        '' parameters.'')')
                  write(ieo,'('' ====='')')

                  iprun = llnow

                  write(ieo,'(/
     &                '' **** Run is skipped by GROUND ****''/
     &                ''      charge = '',i4/
     &                ''      mass   = '',i4)') msprgr, massal

                  ierr = 1
                  return

               end if

*-----------------------------------------------------------------------
*              First particle
*-----------------------------------------------------------------------

                  rwod = -1.0

            do while( rn(0) * rmax .gt. rwod )

                  rsqr = 10.0

               do while( rsqr .gt. 1.0 )

                  rx = 1.0 - 2.0 * rn(0)
                  ry = 1.0 - 2.0 * rn(0)
                  rz = 1.0 - 2.0 * rn(0)

                  rsqr = rx*rx + ry*ry + rz*rz

               end do

                  rrr  = radm * sqrt(rsqr)
                  rwod = 1.0 / ( 1.0 + exp( ( rrr - rt00 ) / saa ) )

            end do

                  r(1,1) = radm * rx
                  r(2,1) = radm * ry
                  r(3,1) = radm * rz

*-----------------------------------------------------------------------
*              the other particles i > 2
*-----------------------------------------------------------------------

         do 100 i = 2, massal

                  ntry1 = 0

  200       continue

                  ntry1 = ntry1 + 1

                  if( ntry1 .gt. mntry ) goto 1000

*-----------------------------------------------------------------------

                  rwod = -1.0

            do while( rn(0) * rmax .gt. rwod )

                  rsqr = 10.0

               do while( rsqr .gt. 1.0 )

                  rx = 1.0 - 2.0 * rn(0)
                  ry = 1.0 - 2.0 * rn(0)
                  rz = 1.0 - 2.0 * rn(0)

                  rsqr = rx*rx + ry*ry + rz*rz

               end do

                  rrr  = radm * sqrt(rsqr)
                  rwod = 1.0 / ( 1.0 + exp( ( rrr - rt00 ) / saa ) )

            end do

                  r(1,i) = radm * rx
                  r(2,i) = radm * ry
                  r(3,i) = radm * rz

*-----------------------------------------------------------------------
*           Check if distance is large enough.
*-----------------------------------------------------------------------

            do j = 1, i - 1

                  rdis2 = ( r(1,i) - r(1,j) )**2
     &                  + ( r(2,i) - r(2,j) )**2
     &                  + ( r(3,i) - r(3,j) )**2

               if( ichg(i) .eq. ichg(j) ) then

                  dmin2 = dsam2

               else

                  dmin2 = ddif2

               end if

                  if( rdis2. lt. dmin2 ) goto 200

            end do

*-----------------------------------------------------------------------

 100     continue

*-----------------------------------------------------------------------
*        Two-body quantities, local density rhol(i),
*        depth of potential dpot(i), and total potential energy
*-----------------------------------------------------------------------

                  call caldisa

               do i = 1, massal

                  rhoa(i) = 0.0
                  rhos(i) = 0.0
                  rhoc(i) = 0.0

               end do

               do i = 1, massal
               do j = 1, massal

                  rhoa(i) = rhoa(i) + rha(j,i)

                  rhos(i) = rhos(i) + rha(j,i) * inuc(j) * inuc(i)
     &                    * ( 1.0 - 2.0 * abs( ichg(j) - ichg(i) ) )

                  rhoc(i) = rhoc(i) + rhe(j,i)

               end do
               end do


               do i = 1, massal

                  rhol(i) = cdp * ( rhoa(i) + 1.0 )

                  dpot(i) = c0p * rhoa(i)
     &                    + c3p * rhoa(i)**gamm
     &                    + csp * rhos(i)
     &                    + clp * rhoc(i)

               end do

*-----------------------------------------------------------------------
*        (2) Determine momentum of nucelons
*-----------------------------------------------------------------------

               do i = 1, massal

                  phase0(i) = 0.0

               end do

*-----------------------------------------------------------------------
*           First particle
*-----------------------------------------------------------------------

                  pfm = hbc * ( 3.0 / 2.0 * pi**2 * rhol(1) )**(1./3.)

               if( massal .gt. 10 .and. ebini .gt. -5.5 ) then

                  pfm = pfm
     &                * ( 1.0 + 0.2 * sqrt( abs( 8.0 + ebini ) / 8.0 ) )

               end if

*-----------------------------------------------------------------------

                  ntry2 = 0

  705       continue

                  ntry2 = ntry2 + 1

               if( ntry2 .gt. mntry ) goto 1000

                  psqr = 10.0

               do while( psqr .gt. 1.0 )

                  px = 1.0 - 2.0 * rn(0)
                  py = 1.0 - 2.0 * rn(0)
                  pz = 1.0 - 2.0 * rn(0)

                  psqr = px*px + py*py + pz*pz

               end do

                  p(1,1) = pfm * px
                  p(2,1) = pfm * py
                  p(3,1) = pfm * pz
                  p(4,1) = sqrt( psqr * pfm**2 + p(5,1)**2 )

                  ekini = p(4,1) - p(5,1)

               if( ekini + dpot(1) .gt. edepth ) goto 705

*-----------------------------------------------------------------------
*           the other particles i > 2
*-----------------------------------------------------------------------

         do 300 i = 2, massal

                  pfm = hbc * ( 3.0 / 2.0 * pi**2 * rhol(i) )**(1./3.)

               if( massal .gt. 10 .and. ebini .gt. -5.5 ) then

                  pfm = pfm
     &                * ( 1.0 + 0.2 * sqrt( abs( 8.0 + ebini ) / 8.0 ) )

               end if

*-----------------------------------------------------------------------

                  ntry3 = 0

  710       continue

                  ntry3 = ntry3 + 1

               if( ntry3 .gt. mntry ) goto 1000


  700       continue

                  psqr = 10.0

               do while( psqr .gt. 1.0 )

                  px = 1.0 - 2.0 * rn(0)
                  py = 1.0 - 2.0 * rn(0)
                  pz = 1.0 - 2.0 * rn(0)

                  psqr = px*px + py*py + pz*pz

               end do

                  p(1,i) = pfm * px
                  p(2,i) = pfm * py
                  p(3,i) = pfm * pz
                  p(4,i) = sqrt( psqr * pfm**2 + p(5,i)**2 )

                  ekini = p(4,i) - p(5,i)

               if( ekini + dpot(i) .gt. edepth ) goto 700

*-----------------------------------------------------------------------
*              Check Pauli principle
*-----------------------------------------------------------------------

                        phase(i) = 0.0

            do j = 1, i - 1

                        phase(j) = 0.0

               if( ichg(i) .eq. ichg(j) ) then

                        expa = - rr2(i,j) * cpw

                  if( expa .gt. epsx ) then

                        pdist2 = ( p(1,i) - p(1,j) )**2
     &                         + ( p(2,i) - p(2,j) )**2
     &                         + ( p(3,i) - p(3,j) )**2

                        pdist2 = pdist2 * cph
                        expa   = expa - pdist2

                     if( expa .gt. epsx ) then

                        phase(j)  = exp( expa )

                        if( phase(j) * cpc .gt. 0.2 ) goto 710

                        if(( phase0(j) + phase(j) ) * cpc .gt. 0.5 )
     &                                                goto 710

                        phase(i)  = phase(i) + phase(j)

                        if( phase(i) * cpc .gt. 0.3 ) goto 710

                     end if

                  end if

               end if

            end do

                        phase0(i) = phase(i)

            do j = 1, i-1

                        phase0(j) = phase0(j) + phase(j)

            end do


*-----------------------------------------------------------------------

  300    continue

*-----------------------------------------------------------------------
*        Shift nucleus to thier CM frame and kill angular momentum
*-----------------------------------------------------------------------

               call gcmang

*-----------------------------------------------------------------------
*        Is the binding energy O.K. ?
*-----------------------------------------------------------------------

               ekinal = 0.0

            do i = 1, massal

               ekinal = ekinal + p(4,i) - p(5,i)

            end do

               call caldisa
               call epotall(epotal)

               ebinal = ( epotal + ekinal ) / float( massal )

*-----------------------------------------------------------------------
*        Random Packing :
*        ipchs = 0
*-----------------------------------------------------------------------

            if( ebinal .lt. ebin0 .or.
     &          ebinal .gt. ebin1 ) goto 1000

c!!!!         if( ipchs .eq. 0 ) then
          if( ipchs .eq. 0 ) return
c!!!!!!!!!!!!!!!!!!!!!!!!!!
c            write(*,'(a, 3i5,2i3,1p,g12.3)') 
c     *       '#try= ', ntry, ntry1, mntry, 
c     *        massal, msprgr, gambb0
c             return
c          endif
c !!!!!!!!!!!!!!!!!!

*-----------------------------------------------------------------------
*        Adjust by frictional cooling or heating
*        ipchs = 1
*        mstq1(190) = 0: normal, 1: adjust in moving frame
*-----------------------------------------------------------------------

         if( mstq1(190) .eq. 0 ) goto 123

*-----------------------------------------------------------------------
*        adjust in moving frame
*-----------------------------------------------------------------------

               gambb = min( gambb0, 1.3d0 )

               gamaa = ( gambb - 1.0 ) * 0.055
     &               * ( dble(massal)**(1./3.) * 0.152 + 0.1 )

               nnmax = 120

*-----------------------------------------------------------------------
*           over cooling in moving frame
*-----------------------------------------------------------------------

               ekinal = 0.0

            do i = 1, massal

               r(3,i) = r(3,i) / gambb
               p(3,i) = p(3,i) * gambb
               p(4,i) = sqrt( p(5,i)**2
     &                      + p(1,i)**2
     &                      + p(2,i)**2
     &                      + p(3,i)**2 )
               ekinal = ekinal + p(4,i) - p(5,i)

            end do

               call caldisa

               call epotall(epotal)

               ebinal = ( epotal + ekinal ) / float( massal )

                  dtc  =  1.0
                  frg  = -0.1
                  rdf0 =  0.5

                  edif0 = ebinal - ebin00 + gamaa

               if( edif0 .gt. 0.0 ) then

                  cfrc =   frg

               else

                  cfrc = - frg

               end if

                  ifrc = 1

*-----------------------------------------------------------------------

            do k = 1, nnmax

               edif = ebinal - ebin00 + gamaa

               if( abs( edif ) .lt. epse ) goto 351

               if( edif .lt. 0.0 ) then

                  jfrc =  1

               else

                  jfrc = -1

               end if

               if( jfrc .ne. ifrc ) then

                  cfrc = - rdf0 * cfrc
                  dtc  =   rdf0 * dtc

               end if

               if( jfrc .eq. ifrc .and.
     &             abs(edif) .gt. abs(edif0) ) then

                  cfrc = - rdf0 * cfrc
                  dtc  =   rdf0 * dtc

               end if

                  ifrc  = jfrc
                  edif0 = edif

                  call gradu(d1r,d1p)

               do i = 1, massal
               do j = 1, 3

                  r(j,i) = r(j,i) + dtc * ( d1r(j,i) - cfrc * d1p(j,i) )
                  p(j,i) = p(j,i) + dtc * ( d1p(j,i) + cfrc * d1r(j,i) )

               end do
               end do

                  ekinal = 0.0

               do i = 1, massal

                  p(4,i) = sqrt( p(5,i)**2
     &                         + p(1,i)**2 + p(2,i)**2 + p(3,i)**2 )
                  ekinal = ekinal + p(4,i) - p(5,i)

               end do

                  call caldisa
                  call epotall(epotal)

                  ebinal = ( epotal + ekinal ) / float( massal )

            end do

*-----------------------------------------------------------------------

  351       continue

*-----------------------------------------------------------------------
*           heating up and adjust in moving frame
*-----------------------------------------------------------------------

                  dtc  =  1.0
                  frg  = -0.1
                  rdf0 =  0.5

                  edif0 = ebinal - ebin00

               if( edif0 .gt. 0.0 ) then

                  cfrc =   frg

               else

                  cfrc = - frg

               end if

                  ifrc = 1

*-----------------------------------------------------------------------

            do k = 1, nnmax

               edif = ebinal - ebin00

               if( abs( edif ) .lt. epse ) goto 352

               if( edif .lt. 0.0 ) then

                  jfrc =  1

               else

                  jfrc = -1

               end if

               if( jfrc .ne. ifrc ) then

                  cfrc = - rdf0 * cfrc
                  dtc  =   rdf0 * dtc

               end if

               if( jfrc .eq. ifrc .and.
     &             abs(edif) .gt. abs(edif0) ) then

                  cfrc = - rdf0 * cfrc
                  dtc  =   rdf0 * dtc

               end if

                  ifrc  = jfrc
                  edif0 = edif

                  call gradu(d1r,d1p)

               do i = 1, massal
               do j = 1, 3

                  r(j,i) = r(j,i) + dtc * ( d1r(j,i) - cfrc * d1p(j,i) )
                  p(j,i) = p(j,i) + dtc * ( d1p(j,i) + cfrc * d1r(j,i) )

               end do
               end do

                  ekinal = 0.0

               do i = 1, massal

                  p(4,i) = sqrt( p(5,i)**2
     &                         + p(1,i)**2 + p(2,i)**2 + p(3,i)**2 )
                  ekinal = ekinal + p(4,i) - p(5,i)

               end do

                  call caldisa
                  call epotall(epotal)

                  ebinal = ( epotal + ekinal ) / float( massal )

            end do

*-----------------------------------------------------------------------

  352       continue

*-----------------------------------------------------------------------
*           go back to the rest frame
*-----------------------------------------------------------------------

               ekinal = 0.0

            do i = 1, massal

               r(3,i) = r(3,i) * gambb
               p(3,i) = p(3,i) / gambb
               p(4,i) = sqrt( p(5,i)**2
     &                      + p(1,i)**2
     &                      + p(2,i)**2
     &                      + p(3,i)**2 )
               ekinal = ekinal + p(4,i) - p(5,i)

            end do

               call caldisa

               call epotall(epotal)

               ebinal = ( epotal + ekinal ) / float( massal )

*-----------------------------------------------------------------------
*        normal adjust in rest frame
*-----------------------------------------------------------------------

  123    continue

                  dtc  =  1.0
                  frg  = -0.1
                  rdf0 =  0.5

                  edif0 = ebinal - ebin00

               if( edif0 .gt. 0.0 ) then

                  cfrc =   frg

               else

                  cfrc = - frg

               end if

                  ifrc = 1

            do k = 1, 200

               edif = ebinal - ebin00

               if( ntry .lt. mntry / 2 ) then
                  if( abs( edif ) .lt. epse ) goto 353
               else if( ntry .lt. mntry - 10 ) then
                  if( abs( edif ) .lt. epse * 1000.0 ) goto 353
               end if

               if( edif .lt. 0.0 ) then

                  jfrc =  1

               else

                  jfrc = -1

               end if

               if( jfrc .ne. ifrc ) then

                  cfrc = - rdf0 * cfrc
                  dtc  =   rdf0 * dtc

               end if

               if( jfrc .eq. ifrc .and.
     &             abs(edif) .gt. abs(edif0) ) then

                  cfrc = - rdf0 * cfrc
                  dtc  =   rdf0 * dtc

               end if

                  ifrc  = jfrc
                  edif0 = edif

                  call gradu(d1r,d1p)

               do i = 1, massal
               do j = 1, 3

                  r(j,i) = r(j,i) + dtc * ( d1r(j,i) - cfrc * d1p(j,i) )
                  p(j,i) = p(j,i) + dtc * ( d1p(j,i) + cfrc * d1r(j,i) )

               end do
               end do

                  ekinal = 0.0

               do i = 1, massal

                  p(4,i) = sqrt( p(5,i)**2
     &                         + p(1,i)**2 + p(2,i)**2 + p(3,i)**2 )
                  ekinal = ekinal + p(4,i) - p(5,i)

               end do

                  call caldisa
                  call epotall(epotal)

                  ebinal = ( epotal + ekinal ) / float( massal )

            end do

c                 write(ieo,'(/
c    &                '' **** Fault in Adjust Ground Energy ****'')')
c                 write(ieo,'(
c    &                ''      Edif = '',e13.4,'' MeV'')') edif*1000.0

            goto 1000

  353       continue

*-----------------------------------------------------------------------
c !!!!!!!!!!!!!!!!!!
c            write(*,'(a, 3i5,2i3,1p,g12.3)') 
c     *       '#try= ', ntry, ntry1, mntry, 
c     *        massal, msprgr, gambb0
c !!!!!!!!!!!!!!!!!!!
      return
       end subroutine


************************************************************************
*                                                                      *
      subroutine gcmang
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to kill cm motion and angular momentum of the nucleus   *
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

*-----------------------------------------------------------------------

      dimension       sekp(3), sekr(3)
      dimension       rll(3), rr(3,3), ss(3,3), opl(3)

*-----------------------------------------------------------------------
*        Move to cm system
*-----------------------------------------------------------------------


            do j = 1, 3

               sekp(j) = 0.0
               sekr(j) = 0.0

            end do

               sekm = 0.0

         do i = 1, massal

            do j = 1, 3

               sekp(j) = sekp(j) + p(j,i)
               sekr(j) = sekr(j) + r(j,i) * p(5,i)

            end do

               sekm = sekm + p(5,i)

         end do

            do j = 1, 3

               sekp(j) = sekp(j) / float(massal)
               sekr(j) = sekr(j) / sekm

            end do


         do i = 1, massal

               sek = 0.0

            do j = 1, 3

               p(j,i) = p(j,i) - sekp(j)
               r(j,i) = r(j,i) - sekr(j)

               sek = sek + p(j,i)**2

            end do

               p(4,i) = sqrt( sek + p(5,i)**2 )

         end do


*-----------------------------------------------------------------------
*        kill the angular momentum
*-----------------------------------------------------------------------

         if( massal .gt. 2 ) then

*-----------------------------------------------------------------------

                  rll(1) = 0.0
                  rll(2) = 0.0
                  rll(3) = 0.0

               do i = 1, massal

                  rll(1) = rll(1) + r(2,i) * p(3,i)
     &                            - r(3,i) * p(2,i)
                  rll(2) = rll(2) + r(3,i) * p(1,i)
     &                            - r(1,i) * p(3,i)
                  rll(3) = rll(3) + r(1,i) * p(2,i)
     &                            - r(2,i) * p(1,i)

               end do

                  rr(1,1) = 0.0
                  rr(2,1) = 0.0
                  rr(3,1) = 0.0
                  rr(1,2) = 0.0
                  rr(2,2) = 0.0
                  rr(3,2) = 0.0
                  rr(1,3) = 0.0
                  rr(2,3) = 0.0
                  rr(3,3) = 0.0

               do i = 1, massal

                  rr(1,1) = rr(1,1) + r(2,i)**2 + r(3,i)**2
                  rr(2,1) = rr(2,1) - r(2,i) * r(1,i)
                  rr(3,1) = rr(3,1) - r(3,i) * r(1,i)
                  rr(1,2) = rr(1,2) - r(1,i) * r(2,i)
                  rr(2,2) = rr(2,2) + r(1,i)**2 + r(3,i)**2
                  rr(3,2) = rr(3,2) - r(3,i) * r(2,i)
                  rr(1,3) = rr(1,3) - r(1,i) * r(3,i)
                  rr(2,3) = rr(2,3) - r(2,i) * r(3,i)
                  rr(3,3) = rr(3,3) + r(1,i)**2 + r(2,i)**2

               end do

                  ss(1,1) = 1.0
                  ss(2,1) = 0.0
                  ss(3,1) = 0.0
                  ss(1,2) = 0.0
                  ss(2,2) = 1.0
                  ss(3,2) = 0.0
                  ss(1,3) = 0.0
                  ss(2,3) = 0.0
                  ss(3,3) = 1.0

               do i = 1, 3

                     ap = rr(i,i)

                  do m = 1, 3

                     rr(i,m) = rr(i,m) / ap
                     ss(i,m) = ss(i,m) / ap

                  end do

                  do j = 1, 3

                     if( j .ne. i ) then

                           ap = rr(j,i)

                        do m = 1, 3

                           rr(j,m) = rr(j,m) - ap * rr(i,m)
                           ss(j,m) = ss(j,m) - ap * ss(i,m)

                        end do

                     end if

                  end do

               end do

               do i = 1, 3

                     opl(i) = 0

                  do j = 1, 3

                     opl(i) = opl(i) + ss(i,j) * rll(j)

                  end do

               end do

               do i = 1, massal

                 p(1,i) = p(1,i) + r(2,i) * opl(3) - r(3,i) * opl(2)
                 p(2,i) = p(2,i) + r(3,i) * opl(1) - r(1,i) * opl(3)
                 p(3,i) = p(3,i) + r(1,i) * opl(2) - r(2,i) * opl(1)

                 p(4,i) = sqrt( p(5,i)**2
     &                        + p(1,i)**2
     &                        + p(2,i)**2
     &                        + p(3,i)**2 )

              end do

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------

      return
       end subroutine


************************************************************************
*                                                                      *
      subroutine rboost(ierr)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to boost the ground state                               *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param01.inc'
      include 'param02.inc'

*-----------------------------------------------------------------------

      common /const0/ idnta, idnpr, massta, masspr, mstapr, msprpr

      common /const5/ pzpr, pxpr, rzpr, rxpr
      common /const6/ pzta, pxta, rzta, rxta
      common /const7/ betpr, gampr, prmas, radpr
      common /const8/ betta, gamta, tamas, radta

      common /coodrp/ r(5,nnn),  p(5,nnn)
      common /coodid/ ichg(nnn), inuc(nnn), ibry(nnn), inds(nnn),
     &                inun(nnn), iavd(nnn), ihis(nnn)

*-----------------------------------------------------------------------

            ierr = 0

*-----------------------------------------------------------------------
*        Coulomb Trajectory
*-----------------------------------------------------------------------

            call bcoul(ierr)

*-----------------------------------------------------------------------
*        boost the initial ground state
*-----------------------------------------------------------------------

         do i = 1, massta + masspr

            if( i .le. massta ) then

               rxb  = rxta
               rzb  = rzta
               pxb  = pxta
               pzb  = pzta
               gmb  = gamta

            else

               rxb  = rxpr
               rzb  = rzpr
               pxb  = pxpr
               pzb  = pzpr
               gmb  = gampr

            end if

               r(1,i) = r(1,i) + rxb
               r(2,i) = r(2,i)
               r(3,i) = r(3,i) / gmb + rzb

               p(1,i) = p(1,i) + pxb
               p(2,i) = p(2,i)
               p(3,i) = p(3,i) * gmb + pzb

               p(4,i) = sqrt( p(5,i)**2
     &                      + p(1,i)**2
     &                      + p(2,i)**2
     &                      + p(3,i)**2 )

         end do

*-----------------------------------------------------------------------

               call caldisa

*-----------------------------------------------------------------------

      return
       end subroutine


************************************************************************
*                                                                      *
      subroutine bcoul(ierr)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to determine initial position according to b            *
*              and initial position and momentum according to          *
*              coulomb trajectory (non-relativistic)                   *
*              only for elab < 5.0 GeV                                 *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param01.inc'
      include 'param02.inc'

*-----------------------------------------------------------------------

      common /const0/ idnta, idnpr, massta, masspr, mstapr, msprpr
      common /const1/ elab, rdist, bmin, bmax, ibch, ibin
      common /const2/ dt, ntmax, iprun, iprun0

      common /const5/ pzpr, pxpr, rzpr, rxpr
      common /const6/ pzta, pxta, rzta, rxta
      common /const7/ betpr, gampr, prmas, radpr
      common /const8/ betta, gamta, tamas, radta

      common /coultr/ eccm, pzcc, rmax0, zeroz
      common /framtr/ betafr(0:2), gammfr(0:2)

      common /vriab0/ massal, massba, mmeson
      common /vriab1/ b, llnow, ntnow

      common /swich0/ icoul, irelcr, iavoid
      common /swich1/ ipot, insys, irkg, icolt

      common /rannum/ iseed, iseed0, iseed1

*-----------------------------------------------------------------------

            ierr = 0
            ieo  = 0

*-----------------------------------------------------------------------

            betacm = betafr(1)
            gammcm = gammfr(1)

*-----------------------------------------------------------------------
*     without Coulomb trajectory, for no Coulomb or Elab > 5GeV
*-----------------------------------------------------------------------

            if( b .eq. 0.0 .or.
     &          elab .gt. 5.0 .or. 
     &          msprpr * mstapr * icoul .eq. 0 ) then

                  rxpr =   b / 2.0
                  rxta = - b / 2.0

                  rzpr = rzpr - zeroz
                  rzta = rzta - zeroz

                  return

            end if

*-----------------------------------------------------------------------
*     with Coulomb trajectory
*-----------------------------------------------------------------------

         irmax = 0

  100    continue

            rmax   = sqrt( rmax0**2 + b**2 )

            pcca   = 1.0 - float( msprpr * mstapr ) * ccoul
     &             / eccm / rmax - ( b / rmax )**2

*-----------------------------------------------------------------------

         if( pcca .le. 0.0 ) then

                  write(ieo,'(/'' Warning: rdis ( initial distance )'',
     &                         '' is too small for Coulomb'',
     &                         '' trajectory, we increased it.''/
     &                         '' Now it is '',g14.5,'' fm'')')
     &                            rmax0 * 1.2
                  write(ieo,'( '' ======='')')

            irmax = irmax + 1

               if( irmax .gt. 20 ) goto 999

            rmax0 = rmax0 * 1.2

            goto 100

         end if

            pccf   = sqrt(pcca)

*-----------------------------------------------------------------------

               aas    = 2.0 * eccm * b
     &                / float( msprpr * mstapr ) / ccoul
               bbs    = 1.0 / sqrt( 1.0 + aas**2 )
               aas1   = ( 1.0 + aas * b / rmax ) * bbs

            if( 1.0 - aas1**2 .le. 0.0 .or.
     &          1.0 - bbs**2  .le. 0.0 ) then

               cost   = 1.0
               sint   = 0.0

            else

               aat1   = aas1 / sqrt( 1.0 - aas1**2 )
               aat2   = bbs  / sqrt( 1.0 - bbs**2 )

               thet1  = atan( aat1 )
               thet2  = atan( aat2 )

            end if

               theta  = thet1 - thet2

               cost   = cos( theta )
               sint   = sin( theta )

*-----------------------------------------------------------------------
*        Coulomb trajectory in classical CM system
*-----------------------------------------------------------------------

            rzpr   = - rmax * cost * tamas * massta
     &               / ( prmas * masspr + tamas * massta )

            rzta   =   rmax * cost * prmas * masspr
     &               / ( prmas * masspr + tamas * massta )

            rxpr   =   rmax / 2.0 * sint
            rxta   = - rxpr

            pzpc   =   pzcc * ( cost * pccf + sint * b / rmax)
            pxpr   =   pzcc * (-sint * pccf + cost * b / rmax)

            pztc   = - pzpc
            pxta   = - pxpr

            epc    =   sqrt( pzpc**2 + pxpr**2
     &             +   ( prmas * masspr )**2 )

            etc    =   sqrt( pztc**2 + pxta**2
     &             +   ( tamas * massta )**2 )

*-----------------------------------------------------------------------
*        Transformation from CM to reference frame
*-----------------------------------------------------------------------

            pzpr   = pzpc + betacm * gammcm
     &             * ( gammcm / ( 1. + gammcm ) * pzpc * betacm + epc )

            pzta   = pztc + betacm * gammcm
     &             * ( gammcm / ( 1. + gammcm ) * pztc * betacm + etc )

            epr    = gammcm * ( epc + betacm * pzpc )
            eta    = gammcm * ( etc + betacm * pztc )


*-----------------------------------------------------------------------
*        Initial position, momentum and beta, gamma for Coulomb Tr.
*-----------------------------------------------------------------------

            betpr  = pzpr / epr
            betta  = pzta / eta

            gampr  = epr / ( prmas * masspr )
            gamta  = eta / ( tamas * massta )

            pzta   = pzta / float( massta )
            pxta   = pxta / float( massta )

            pzpr   = pzpr / float( masspr )
            pxpr   = pxpr / float( masspr )

            rzpr   = rzpr - zeroz
            rzta   = rzta - zeroz

*-----------------------------------------------------------------------

      return

*-----------------------------------------------------------------------

  999       continue

                  write(ieo,'(/'' Erorr: rdis ( initial distance )'',
     &                         '' is too small for Coulomb'',
     &                         '' trajectory.'',
     &                         '' Now it is '',g14.5,'' fm'')') rmax0
                  write(ieo,'( '' ====='')')

                  iprun = llnow

                  write(ieo,'(/
     &                '' **** Run is skipped by BCOUL ****''/
     &                ''      at Event = '',i6/
     &                ''         iseed = '',i12)') iprun, iseed1

                  ierr = 1
                  return

*-----------------------------------------------------------------------

       end subroutine


