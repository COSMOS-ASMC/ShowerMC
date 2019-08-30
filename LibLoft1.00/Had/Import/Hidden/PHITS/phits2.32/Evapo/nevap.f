************************************************************************
*                                                                      *
      subroutine nevap(ipos)
*                                                                      *
*       main control routine of Nuclear Evaporation                    *
*       modified by K.Niita on 2005/12/22                              *
*                                                                      *
*     input:                                                           *
*                                                                      *
*        ipos  : = 0, call from ovly12 and 13 ,=1 from tally           *
*                = 2, call from sctneut                                *
*            /// always 0                                              *
*     output:                                                          *
*                                                                      *
*---- in common -------------------------------------------------------*
*                                                                      *
*        nclsts   : total number of out going particles and nuclei     *
*                                                                      *
*        iclusts(nclsts)                                               *
*        
*                i = 0, nucleus                                        *
*                  = 1, proton                                         *
*                  = 2, neutron                                        *
*                  = 3, pion                                           *
*                  = 4, photon                                         *
*                  = 5, kaon                                           *
*                  = 6, muon                                           *
*                  = 7, others                                         *
*                                                                      *
*        jclusts(i,nclsts)                                             *
*                                                                      *
*                i = 0, angular momentum                               *
*                  = 1, proton number                                  *
*                  = 2, neutron number                                 *
*                  = 3, ip, see below                                  *
*                  = 4, status of the particle 0: real, <0 : dead      *
*                  = 5, charge                                         *
*                  = 6, baryon number                                  *
*                  = 7, kf code                                        *
*                                                                      *
*        qclusts(i,nclsts)                                             *
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
*        numpam(i) : after evaporation                                 *
*                                                                      *
*        numpal(i) : event of final = cascade + evaoparation           *
*        rumpal(i) : weight                                            *
*                  : total number of out going particles or nuclei     *
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
*        kdecay(4) = 0 : no fission                                    *
*                  = 1 : with fission                                  *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param02.inc'

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp
      common /qparm/  ielas,icasc,iqstep,lvlopt,igamma
      common /bparm/  andt,jevap,npidk

*-----------------------------------------------------------------------

      common /clustf/ nclst, iclust(nnn)
      common /clustg/ jclust(0:7,nnn),  qclust(0:12,nnn)
      common /clustp/ rumpat(0:20), numpat(0:20)

*-----------------------------------------------------------------------

      common /clusts/ nclust, kclust(3,nnn)
      common /clustu/ lclust(0:7,nnn), sclust(0:12,nnn)
      common /clustm/ rumpam(0:20), numpam(0:20)
      common /clustv/ kdecay(4)
      common /mathzn/ mathz, mathn, jcoll, kcoll

*-----------------------------------------------------------------------

      common /clustt/ nclsts, iclusts(nnn)
      common /clustw/ jclusts(0:7,nnn),  qclusts(0:12,nnn)
      common /clustl/ rumpal(0:20), numpal(0:20)

*-----------------------------------------------------------------------

      common /kcount/ rncnt(40), rnint(40)
      common /wtsave/ oldwt

*-----------------------------------------------------------------------

      dimension epht(200)

*-----------------------------------------------------------------------
*        do loop for the clusters from CASCAD
*-----------------------------------------------------------------------
c/////////////
c      write(0,*) ' now in nevap; nclst, jevap=', nclst, jevap
c///////////////
               if( ipos .ne. 2 ) nclsts = 0

                  numnut = 0

               do i = 1, 4

                  kdecay(i) = 0

               end do

               do i = 0, 20

                     numpam(i) = 0
                     rumpam(i) = 0.0d0

                  if( ipos .ne. 2 ) then

                     numpal(i) = 0
                     rumpal(i) = 0.0d0

                  end if

               end do

*-----------------------------------------------------------------------

!///      if( nclst .gt. 0 ) then

      do 100 i = 1, nclst

                  ex = qclust(6,i)
*-----------------------------------------------------------------------
*        statistical decay of cluster
*-----------------------------------------------------------------------

         if( iclust(i) .eq. 0 .and. jevap .gt. 0 ) then

*-----------------------------------------------------------------------

                  iz = jclust(1,i)
                  in = jclust(2,i)

               if( iz .lt. 0 .or. in .lt. 0 .or.
     &             iz + in .eq. 0 ) goto 100

!/////                  if( ipos .eq. 0 .or. ipos .eq. 2 ) call cputime(9)

                  jj = jclust(0,i)
                  bi = qclust(0,i)

                  px = qclust(1,i)
                  py = qclust(2,i)
                  pz = qclust(3,i)
                  et = qclust(4,i)
                  rm = qclust(5,i)
                  ex = qclust(6,i)
                  ek = qclust(7,i)
                  wt = qclust(8,i)

                  pt = sqrt( px**2 + py**2 + pz**2 )

                  ex0 = ex
                  er0 = ek

*-----------------------------------------------------------------------
*           statistical decay, call sdmexec
*-----------------------------------------------------------------------
!           iqstep  is always 6
!            if( iqstep .eq. 4 ) then
!
!               call sdmexec(iz,in,jj,ex,px,py,pz,wt,ierr)
!
!
!            else if( iqstep .eq. 2 .or. iqstep .eq. 3 .or.
!     &               iqstep .eq. 5 ) then
!
!               call erupin(iz,in,ex,px,py,pz,pt,et,rm,wt)
!
!
!            else if( iqstep .eq. 6 ) then

               call gemexec(iz,in,ex,px,py,pz,pt,et,rm,wt,ierr)

!            end if

*-----------------------------------------------------------------------
*           new booking
*-----------------------------------------------------------------------

            do j = 1, nclust

               if( kclust(1,j) .ge. 100 ) then

                     nclsts = nclsts + 1

                  if( lclust(1,j) .eq. 0 .and.
     &                lclust(2,j) .eq. 0 ) then

                     kf    = 22
                     ibary = 0
                     ipid  = 4
                     ippad = 14

                  else if( lclust(1,j) .eq. 1 .and.
     &                     lclust(2,j) .eq. 0 ) then

                     kf    = 2212
                     ibary = 1
                     ipid  = 1
                     ippad = 1

                  else if( lclust(1,j) .eq. 0 .and.
     &                     lclust(2,j) .eq. 1 ) then

                     kf    = 2112
                     ibary = 1
                     ipid  = 2
                     ippad = 2

                  else

                     ipid  = 0
                     ibary = lclust(1,j) + lclust(2,j)

                     kf    = lclust(1,j) * 1000000
     &                     + lclust(1,j) + lclust(2,j)

                     if( lclust(1,j) .eq. 1 .and.
     &                   lclust(2,j) .eq. 1 ) then

                        ippad  = 15

                     else if( lclust(1,j) .eq. 1 .and.
     &                        lclust(2,j) .eq. 2 ) then

                        ippad  = 16

                     else if( lclust(1,j) .eq. 2 .and.
     &                        lclust(2,j) .eq. 1 ) then

                        ippad  = 17

                     else if( lclust(1,j) .eq. 2 .and.
     &                        lclust(2,j) .eq. 2 ) then

                        ippad  = 18

                     else

                        ippad  = 19

                     end if

                  end if

                     iclusts(nclsts)   = ipid

                     jclusts(0,nclsts) = lclust(0,j)
                     jclusts(1,nclsts) = lclust(1,j)
                     jclusts(2,nclsts) = lclust(2,j)
                     jclusts(3,nclsts) = ippad
                     jclusts(4,nclsts) = 0
                     jclusts(5,nclsts) = lclust(1,j)
                     jclusts(6,nclsts) = ibary
                     jclusts(7,nclsts) = kf

                     qclusts(0,nclsts) = bi


                  do l = 1, 12

                     qclusts(l,nclsts) = sclust(l,j)

                  end do

                     if( qclusts(6,nclsts) .lt. 0.0 )
     &                   qclusts(6,nclsts) = 0.0

                     numpam(ippad) = numpam(ippad) + 1
                     rumpam(ippad) = rumpam(ippad) + sclust(8,j)

               end if

            end do
         
*-----------------------------------------------------------------------

!///               if( ipos .eq. 0 .or. ipos .eq. 2 ) then

!/////                  rncnt(9) = rncnt(9) + 1.0
!/////                  call cputime(9)

!////               end if

*-----------------------------------------------------------------------
*        particle or without evaporation : pass through
*-----------------------------------------------------------------------

         else
c        else if( iclust(i) .gt. 0 .or. jevap .le. 0 ) then

*-----------------------------------------------------------------------

               if( iclust(i) .eq. 0 ) then

                     if( jclust(1,i) .eq. 1 .and.
     &                   jclust(2,i) .eq. 1 ) then

                        ippad  = 15

                     else if( jclust(1,i) .eq. 1 .and.
     &                        jclust(2,i) .eq. 2 ) then

                        ippad  = 16

                     else if( jclust(1,i) .eq. 2 .and.
     &                        jclust(2,i) .eq. 1 ) then

                        ippad  = 17

                     else if( jclust(1,i) .eq. 2 .and.
     &                        jclust(2,i) .eq. 2 ) then

                        ippad  = 18

                     else

                        ippad  = 19

                     end if

                        jclust(3,i) = ippad

               end if

*-----------------------------------------------------------------------

                     nclsts = nclsts + 1

                     iclusts(nclsts)   = iclust(i)

                  do j = 0, 7

                     jclusts(j,nclsts) = jclust(j,i)

                  end do

                  do j = 0, 12

                     qclusts(j,nclsts) = qclust(j,i)

                  end do

                     ippad = jclust(3,i)

                     numpam(ippad) = numpam(ippad) + 1
                     rumpam(ippad) = rumpam(ippad) + qclust(8,i)

         end if

*-----------------------------------------------------------------------

  100 continue

!///      end if

*-----------------------------------------------------------------------
*        gamma production from deexcited residual nuclei
*-----------------------------------------------------------------------

      if( igamma .ne. 0 .and. nclsts .gt. 0 ) then
!////////// igamma is 0
         write(0,*) ' igamma=',igamma, ' in cnevap; strange' 
         stop
!/////////////
         do i = 1, nclsts

            if( iclusts(i) .eq. 0 .and.
     &          qclusts(6,i) .gt. 0.0 .and.
     &          jclusts(1,i) .gt. 0 .and.
     &          jclusts(2,i) .gt. 0 ) then

                     iz  = jclusts(1,i)
                     in  = jclusts(2,i)

                     pxn = qclusts(1,i)
                     pyn = qclusts(2,i)
                     pzn = qclusts(3,i)
                     ern = qclusts(4,i)
                     emn = qclusts(5,i)
                     exn = qclusts(6,i)
                     ekn = qclusts(7,i)
                     wt  = qclusts(8,i)

!/////                  if( ipos .eq. 0 .or. ipos .eq. 2 ) call cputime(2)

                     apr = dble( iz + in )
                     zpr = dble( iz )
                     exi = exn

!///                  call dexgam(apr,zpr,exi,npht,epht)

                  if( ipos .eq. 0 .or. ipos .eq. 2 ) then

!////                     rncnt(2) = rncnt(2) + 1.0
!///                     call cputime(2)

                  end if

               if( npht .gt. 0 ) then

                     pxg = 0.d0
                     pyg = 0.d0
                     pzg = 0.d0

                  do j = 1, npht

                     pr   = epht(j) / 1000.0

                     cos1 = 1.0 - 2.0 * rn(0)
                     sin1 = sqrt( 1.0 - cos1**2 )
                     phi1 = 2.0 * pi * rn(0)

                     pxr = pr * sin1 * cos(phi1)
                     pyr = pr * sin1 * sin(phi1)
                     pzr = pr * cos1

                     pcs = pxn * pxr + pyn * pyr + pzn * pzr

                     trans = ( pcs / ( ern + emn ) + pr ) / emn

                     pxl = pxr + trans * pxn
                     pyl = pyr + trans * pyn
                     pzl = pzr + trans * pzn

                     egl = sqrt( pxl**2 + pyl**2 + pzl**2 )

                     nclsts = nclsts + 1
                     iclusts(nclsts) = 4

                     jclusts(0,nclsts)  = 0
                     jclusts(1,nclsts)  = 0
                     jclusts(2,nclsts)  = 0
                     jclusts(3,nclsts)  = 14
                     jclusts(4,nclsts)  = 0
                     jclusts(5,nclsts)  = 0
                     jclusts(6,nclsts)  = 0
                     jclusts(7,nclsts)  = 22

                     qclusts(1,nclsts)  = pxl
                     qclusts(2,nclsts)  = pyl
                     qclusts(3,nclsts)  = pzl
                     qclusts(4,nclsts)  = egl
                     qclusts(5,nclsts)  = 0.0d0
                     qclusts(6,nclsts)  = 0.0d0
                     qclusts(7,nclsts)  = egl * 1000.
                     qclusts(8,nclsts)  = wt
                     qclusts(9,nclsts)  = 0.0
                     qclusts(10,nclsts) = 0.0d0
                     qclusts(11,nclsts) = 0.0d0
                     qclusts(12,nclsts) = 0.0d0

                     numpam(14) = numpam(14) + 1
                     rumpam(14) = rumpam(14) + wt

                     pxg = pxg + pxl
                     pyg = pyg + pyl
                     pzg = pzg + pzl

                  end do

                     plmx = pxn - pxg
                     plmy = pyn - pyg
                     plmz = pzn - pzg

                     pabs = plmx**2 + plmy**2 + plmz**2
                     etot = sqrt( pabs + emn**2 )
                     elab = etot - emn

                     qclusts(1,i) = plmx
                     qclusts(2,i) = plmy
                     qclusts(3,i) = plmz
                     qclusts(4,i) = etot
                     qclusts(5,i) = emn
                     qclusts(6,i) = max( 0.0d0, exi )
                     qclusts(7,i) = elab * 1000.d0

               end if

            end if

         end do

      end if

*-----------------------------------------------------------------------
*     total number of particles and nuclei
*-----------------------------------------------------------------------

            if( ipos .ne. 2 ) then

               do i = 0, 19

                  numpal(i) = numpam(i)
                  rumpal(i) = rumpam(i) * oldwt

               end do

               if( numpal(ityp) .gt. 0 .and. 
     &           ( jcoll .eq. 8 .or. jcoll .eq. 11 ) ) then

                  numpal(ityp) = numpal(ityp) - 1
                  rumpal(ityp) = rumpal(ityp) - oldwt

               end if

               if( kdecay(4) .gt. 0 ) then

                  numpal(20) = 1
                  rumpal(20) = oldwt

               end if

            end if

*-----------------------------------------------------------------------
!///////////
!            write(0,*)
!     *     ' return nevap; nclst, nclsts=', nclst, nclsts
!///////
c////////
!      write(*,*) 'eva ', nclsts,  nclst
!      do i =1, nclsts
!         write(*,'(a, 4i5, 1p,g12.3)')
!     *   'eva ',  iclusts(i), jclusts(1,i), jclusts(2,i),
!     *                   jclusts(5,i), qclusts(7,i)
!         if(i <= nclst ) then 
!            write(*,'(4x, 4i5, 1p,g12.3)')
!     *          iclust(i), jclust(1,i), jclust(2,i),
!     *                   jclust(5,i), qclust(7,i)
!         endif  
!      enddo
!c//////////
      return
      end subroutine


