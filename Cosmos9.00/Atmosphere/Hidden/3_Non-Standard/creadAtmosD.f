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
      character(20):: matter
      integer i, locs

!     nodeal points are from low height to high height       
!        read basic data
      if(AtmosFile  == '     ') then
         AtmosFile =
     *     "$COSMOSTOP/Atmosphere/AtmosModel/withIceH2O.d"
      endif
      call copenf(TempDev, AtmosFile, icon)
      if(icon .ne. 0) stop 9999
      call cskipComment(TempDev, icon)
      if(icon .ne. 0) stop 9999

      nodes = 0
!     atmos%matter(:) = "  "
      matter = " "
      NoOfMedia = 0      
      do while(.true.)
         read(TempDev, '(a)',end=10, iostat=ios) line
!!!!!!!!!!
!         write(0,*) ' line: ', trim(line)
!!!!!!!!!!!!!!!         
         call kcountFields(line, nf)
!!!!!!!!!!
!         write(0,*) ' nf =', nf
!!!!!!!!!!!!         
         if( nf > 0 ) then
            nodes = nodes + 1
            if( nf == 3 ) then
               read(line, *)
     *           atmos%z(nodes), atmos%T(nodes), atmos%rho(nodes)
!     *              atmos%P(nodes), 
               locs = 0
            elseif( nf == 4 ) then
               read(line, *)
     *              atmos%z(nodes), atmos%T(nodes),
     *                   atmos%rho(nodes), matter
!     *              atmos%P(nodes), atmos%rho(nodes), matter
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
      modelno = 3
      end subroutine cqAtmosModel
