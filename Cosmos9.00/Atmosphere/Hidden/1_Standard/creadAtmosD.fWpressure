      subroutine creadAtmosD
!
!   creadAtmosD; read atmosphere data
!   catmosConsts:  compute some basic constants
!  
!       creadAtmosD:
!  read segmented atmosphere data consisting of 
!      height temperature pressure density       

!       catmosConsts:
!      a,  b,  d0, cum, H
! where
!     a, b; const shown below
!       d0: see the formula below.
!     cumd: amount of atmospher above the node.
!        H: scale height at the node
!     all in mks unit.
!
!     The scale height is approximated by a
!   number of stright lines as a function of height. The data in
!   stdatmos1.d gives height, temperatur, etc at  each nodal point.
!
!         The scale height, H, is expressed by H = H0 + a(z-z0)
!                                                = kT/mg      
!   in each region.
!   We neglect height dependence of gravitational accelleration g,
!   and the average mass of  air molecules, m.   
!   Since the data table gives T(z)= T0 + b(z-z0) at the nodal points, 
!   we can first get b,
!   and then a by a = dH/dz = k/mg * b.  H at a nodal point, z,  
!   is obtained as H(z) =kT(z)/mg.
!
!   The density off a nodal point is given by
!
!              rho = rho0 * (1+ a(z-z0)/H(z0))**(-1-1/a)      (a != 0)
!                  = rho0 * exp(- (z-z0)/H)            (a =0; hence H is const)
!
!   (We employ H(z0) as the scale height in the segment)
!   The amount of air between  given heights, z1 and z2  is by
!
!              d = d0 *(fd(z1) - fd(z2))  where
!         
!            fd(z) = (1+ a(z-z0)/H(z0))**(-1/a)                   (a != 0)
!
!                  =  exp(-(z-z0)/H )                             (a = 0)
!
!   where  d0 = rho0*H(z0)
!
!   If z1=z0,  d becomes
!
!              d= d0 ( 1 - fd(z2))
!
!
!
!
      use modAtmosDef
      implicit none
#include  "Zmanagerp.h"
#include  "ZmediaLoft.h"
      
      integer  ios, nodes, icon, nf

      character*150 msg
      character(100):: line
      character*60 adata(maxnodes)
      character(20):: matter
      integer i, locs
!          adata is the same one as
!          Data/Atmos/stdatmos1.d
!      To avoid reading it as default, it is stored here.
       data (adata(i), i = 1, 15)/
     1 '-400	290.75	1062.2	1.2790 Air',
     2 '0	288.15	1013.25	1.2250',
     3 '3.e3	268.659	7.0121e2 9.0925e-1',
     4 '6.e3	249.187	4.7217e2 6.6011e-1',
     5 '11.1e3	216.65	223.46	0.35932',
     6 '20.0e3	216.65	55.293	8.891e-2',
     7 '32.2e3	228.756	8.6314	1.3145e-2',
     8 '47.4e3	270.65	1.1022	1.4187e-3',
     9 '51.0e3	270.65	7.0458e-1 9.0696e-4',
     a '72.0e3	214.263	3.8362e-2 6.2374e-5',
     b '86.0e3	186.87	3.7388e-3 6.958e-6',
     c '91.0e3	186.87	1.5381e-3 2.860e-6',
     d '110.0e3	240.0	7.1042e-5 9.708e-8',
     e '130.0e3	469.27	1.2505e-5 8.152e-9',
     f '160.0e3	696.29	3.0359e-6 1.233e-9'/

       data (adata(i), i=16, 24)/
     1 '250.0e3	941.33	2.4767e-7 6.073e-11',
     2 '300.0e3	976.01	8.7704e-8 1.916e-11',
     3 '400.0e3 995.83	1.4518e-8 2.803e-12',
     4 '500.0e3	999.24	3.0236e-9 5.215e-13',
     5 '600.0e3	999.85	8.2130e-10 1.137e-13',
     6 '700.0e3	999.97	3.1908e-10 3.070e-14',
     7 '800.0e3	999.99	1.7036e-10 1.136e-14',
     8 '900.0e3	1000.0	1.0873e-10 5.759e-15',
     9 '1.0e6	1000.0	7.5138e-11 3.561e-15'/
       data adata(25)/' '/

!     nodeal points are from low height to high height       
!        read basic data
      if(AtmosFile .ne. '     ') then
         call copenf(TempDev, AtmosFile, icon)
         if(icon .ne. 0) stop 9999
         call cskipComment(TempDev, icon)
         if(icon .ne. 0) stop 9999
      endif
      nodes = 0
!     atmos%matter(:) = "  "
      matter = " "
      NoOfMedia = 0      
      do while(.true.)
         if(AtmosFile .ne. '     ') then
            read(TempDev, '(a)',end=10, iostat=ios) line
         else
            line = adata(nodes+1)
         endif
         call kcountFields(line, nf)
         if( nf > 0 ) then
            nodes = nodes + 1
            if( nf == 4 ) then
               read(line, *)
     *           atmos%z(nodes), atmos%T(nodes),
     *              atmos%P(nodes), atmos%rho(nodes)
               locs = 0
            elseif( nf == 5 ) then
               read(line, *)
     *           atmos%z(nodes), atmos%T(nodes),
!     *           atmos%P(nodes), atmos%rho(nodes), atmos%matter(nodes)
     *              atmos%P(nodes), atmos%rho(nodes), matter
               locs = index(matter,"*")
            else
               write(0,*) ' error in adata or input data creadAtmosD.fh'
               write(0,*) ' data: ', line
               write(0,*) ' nf=', nf
               stop
            endif
         else
            goto 10
         endif
         if(locs == 0 ) then
            atmos%matter(nodes) = matter
            atmos%rhoc(nodes) = 1.0
         elseif(locs > 1 ) then
            atmos%matter(nodes) = matter(1:locs-1)
            read(matter(locs+1:), *) atmos%rhoc(nodes)
         else
            write(0,*) ' ivalide media spec=', matter
            stop
         endif
!         if(ios .ne. 0) then
!            call cerrorMsg(
!     *       'something wrong about atmosphere data file:',1)
!            call cerrorMsg(AtmosFile, 0)
!        endif
      enddo
 10   continue

      atmos%nodes = nodes
         !         assume default is Air
      if( atmos%matter(1) == "  " )  atmos%matter(1) = "Air"
      do i = 2, nodes
         if( atmos%matter(i) == " " )
     *            atmos%matter(i) = atmos%matter(i-1)
      enddo
!     Count # of uniq media (normally 1; only Air)
!     and read media data (by eprdMFfile) to be ready for EM
!      process sampling
      call cgetMedia

      
      if(AtmosFile .ne. '  ') then
         close(TempDev)
         write(msg, *) 
     *   "Atmosphere data has been read: # of nodal points =",
     *    nodes
      else
         write(msg, *) 
     *   "Atmosphere data stored internally is read nodes=",
     *    nodes
      endif
      call cerrorMsg(msg, 1)
      end   subroutine creadAtmosD

      subroutine cqAtmosModel(modelno)
      implicit none
      integer,intent(out):: modelno
      modelno = 1
      end subroutine cqAtmosModel

