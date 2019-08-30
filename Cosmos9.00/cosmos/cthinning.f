      subroutine cthinning(Tracks, n, iTrack, nout)
! 
!       The user may replace this routine by his/her own  thinning method.
!       For 'n' tracks in 'Tracks', change weight if necessary or remove it.
!       If a track is removed, the user must move the remaining tracks
!       to the upper position in Tracks.
!       This will be done simply as
!          nout = 0
!          do i = 1, n
!            Examine Tracks(i)
!            if it is accepted, (do neceeary weight change)
!                nout++
!                store it in Tracks(nout)
!            else do nothing
!          enddo
! 
!   The standard thinning routine is supploed as csetThinwgt
!   A paremeter in $Hparam, EthinRatio, is used as follows.
!   EthinRation(1) to (4) are used as (Ethin, MaxWgt) for e/g
!   and for mu/had.  If EthiRatio(3) and (4) are not given
!   (1)/10 and (2)/10 are used.
!   Others are hard wired and can be chaged in csetThinwgt
!   (see  top part of csetThinwgt.  They  are in the lines
!   above ^^^^^^^^^^^^^^^^^^^^^^^^^.
! 

      implicit none
#include  "Zcode.h"
#include  "Ztrack.h"
! #include  "Zmagfield.h"
#include  "Ztrackp.h"
#include  "Ztrackv.h"
#include  "Zincidentv.h"
      integer n     ! input. no. of produced particles and stored in Tracks
      type(track):: Tracks(n)  ! input and outut.  track info.
      type(track):: iTrack  ! input. incident particle track info. for the coll.
      integer  nout  !  output.  number of particles accepted in the
                     !     thinning
!  
!       
      integer i
      type(track):: aTrack
      integer icon

      nout = 0
      do i = 1, n
         aTrack=Tracks(i)
         call csetThinwgt(iTrack, aTrack, icon)
         if(icon .eq. 0 )  then
            nout = nout +1
            Tracks(nout)=aTrack
         endif
      enddo
      end

      subroutine csetThinwgt(iTrack, aTrack, icon)
      use modEMcontrol
      implicit none
#include  "Zcode.h"
#include  "Ztrack.h"
! #include  "Zmagfield.h"      
#include  "Ztrackp.h"
#include  "Ztrackv.h"
#include  "Zincidentv.h"
#include  "Zelemagp.h"
#include  "Zobs.h"
#include  "Zobsp.h"
#include  "Zobsv.h"
!
! **************** if you want to limit the max wwight exactly 
!                  give 1 to the next line
!                  else give 0.
!                  For exact max limit, you  have to give larger 
!                  weight in EthinRatio(2) and EthinRatio(4) for
!                  defending cpu time increase. (~10 times) 
!            If you give  0 then  if the weight is > the max weight,
!            no more  thinning is tried. However, if the weight is < max,
!            thinning is tried and resultant weight may be
!            larger than max limit.
!            If 1, the weight managed so that it never exceeds
!            the max.  0 has probably better performance.
#define EXACTWGT 0
!           Define `far` by DETAILFAR.  (particle location is far
!           from the axis or not.)
!            0:  far is judged at  depeth1<depth<depth2 for e/g
!                and mu/h 
!            1:  far is judged at depeth1<depth<depth2 for e/g
!                              at all depths < depth2 for mu/h
#define DETAILFAR 0

      real*8 big  ! if you don't want to control the thinnig
                  ! by the distance of the current particle
                  ! from the shower axis, give big to rfar 
      parameter (big=1.d20)
!   ******************************************************************
!   **************************** fix the following *******************
!  
      real*8 depth1/4000.d0/  ! between depth1 and depth2, check if
      real*8 depth2/8750.d0/  ! a ptcl is far from the core. If  so 
             ! we employ lesser thinning or no thinning.
             ! depth2  should be the last obs. depth where lateral
             ! information is taken. 
             ! the unit is kg/m2 (devide by 10 for g/cm2)
      real*8 rfar       ! r>rfar is regared as "far from axis" (m).
                 ! if want to skip distance dependent thinning
                 !  set this to be "big" (see below)
      real*8 rfar2      ! rfar**2
      parameter (rfar =20.0, rfar2= rfar*rfar)
      real*8 deepfactor/10./  ! at depth  > depht2, stronger thinning
               !      by factor 'deepfactor' than standard thinning
               !      specified by EthinRatio.  
               !  Ethin and Weight are multiplied by this factor.
      real*8 farfactor/0.01/  ! thinning factor is weekened by this factor
                  ! than standard one at far points.  
                  !  Ethin and Weight are multiplied by this factor.
!  
!   
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!  




      type(track):: iTrack   ! input  parent particle
      type(track):: aTrack  ! input   child particle
!                            ! output.    wgt
      integer icon           ! output.  0 if this is not tobe discarded
                             !          1 if this is to be discarded
      real*8 u, p 
      real*8 iergpn, aergpn
      logical dothin
      real*8 dd, rhoE
      real*8 cvh2den
      real*8 Z1, E0
      integer ii, jj  ! ii is for incident and jj is for child; for e/g 1,  for mu/had 3
      data Z1/-1./, E0/-1./
      save Z1, E0, rhoE, dd
      type(coord):: oxyz, axisxyz
      real*8 len, h, relax, dist
      logical far
!  
!          
      if(IncidentCopy%p%code .eq. kphoton  .and.
     *    HowPhotoP > 0 )  then
!     *    PhotoProd )  then
!            photon  primary and muon is interested
!            so thinsamling must be carfull
!            we apply thinning only if current depth is
!            > 120 g/cm2 from the first col. point.
         if(Z1 .ne.  Zfirst%pos%depth  .or.
     *      E0 .ne.  IncidentCopy%p%fm%p(4) ) then
            Z1 = Zfirst%pos%depth 
            E0 = IncidentCopy%p%fm%p(4) 
            if( LpmEffect ) then
               rhoE=cvh2den( Zfirst%pos%height )* 1.e-3 * 
     *             IncidentCopy%p%fm%p(4) 
               if(rhoE .lt. 1.e6)  then
                  dd = 200.
               else
                  dd =min( 200.* sqrt(rhoE/1.e6), 1000.d0)
               endif
            else
               dd = 200.
            endif
         endif

         if(  MagBrem .eq. 2 .and.   Z1 .lt. 1.e-6 ) then
            dothin = iTrack%pos%depth/Zfirst%vec%coszenith
     *           .gt. 300.
         else
            dothin=( iTrack%pos%depth-Zfirst%pos%depth)/
     *           Zfirst%vec%coszenith .gt. (1000.+ dd)
         endif
      else
         dothin = .true.
      endif
!  
!  
!          
!  
!     
#if DETAILFAR == 0
      if( dothin .and. rfar .lt. big .and.
     *    iTrack%pos%depth .lt. depth2 .and.
     *    iTrack%pos%depth .gt. depth1  ) then
#elif DETAILFAR == 1
      if( dothin .and. rfar .lt. big .and.
     *    iTrack%pos%depth .lt. depth2 .and.
     *    (iTrack%p%code .gt. kelec  .or.
     *     iTrack%pos%depth .gt. depth1 ) ) then
#endif
!            compute distance form the shower axis
        h = iTrack%pos%height - ObsSites(aTrack%where)%pos%height
        len = h / (-AngleAtObsCopy%r(3))
         axisxyz%r(1:3) = ObsSites(aTrack%where)%pos%xyz%r(1:3) -
     *                      len*DcAtObsXyz%r(1:3)
!         axisxyz%r(2) = ObsSites(aTrack%where)%pos%xyz%y -
!     *                      len*DcAtObsXyz%r(2)
!         axisxyz%z = ObsSites(aTrack%where)%pos%xyz%z -
!     *                      len*DcAtObsXyz%r(3)
!              
!ok             write(0,*) ' h=',h, ' len =',len, ' cos=',
!ok     *               -AngleAtObsCopy.r(3)
!ok             write(0,*) ' axis at  depth =',
!ok     *               ObsSites(aTrack.where).pos.depth
!ok             write(0,*) ' is  x,y,z=', ObsSites(aTrack.where).pos.xyz.x,
!ok     *             ObsSites(aTrack.where).pos.xyz.y,
!ok     *             ObsSites(aTrack.where).pos.xyz.z
!ok             write(0,*) ' ptcl pos=',iTrack.pos.xyz.x,
!ok     *           iTrack.pos.xyz.y, iTrack.pos.xyz.z
!ok             write(0,*) ' dz=',DcAtObsXyz.r(3)
!ok             write(0,*) ' axisxyz.x, y, z=', axisxyz.x, 
!ok     *                   axisxyz.y, axisxyz.z

         if(ObsPlane .eq. horizontal) then
            call cxyz2det(axisxyz,
     *                 aTrack%pos%xyz, oxyz)
         elseif(ObsPlane .eq. perpendicular) then
            call cxyz2prim(axisxyz,
     *                 aTrack%pos%xyz, oxyz)
         endif
         dist = sqrt( oxyz%r(1)**2+ oxyz%r(2)**2)
         far= dist  .gt. rfar
      else
         far=.false.
      endif
!           ////////////////

      if(dothin) then 
         if(iTrack%pos%depth .gt. depth2) then
            relax = deepfactor
         elseif(far) then
             relax = farfactor
!            if( dist .lt. 10.) then
!               relax = farfactor
!            elseif ( dist .lt. 20.) then
!               relax = farfactor*0.3
!            elseif ( dist .lt. 100.) then
!               relax = farfactor*0.1
!            else
!               relax = farfactor*0.03
!            endif
         else
            relax = 1.
         endif
         iergpn = iTrack%p%fm%p(4)
!cc               kinetic or total  energy ?
!cc           by 100 TeV e- pirmary  case, total seems better.
!cc         if( iTrack.p.code .eq. knuc .and. 
!cc     *       iTrack.p.subcode .ne. antip) then
!cc            iergpn = iTrack.p.fm.p(4) - iTrack.p.mass
        if(iTrack%p%code .eq. kgnuc) then
!            iergpn = (iergpn-iTrack.p.mass)/iTrack.p.subcode
            iergpn = iergpn/iTrack%p%subcode
        endif

!         if(aTrack.p.code .eq. knuc .and.
!     *      aTrack.p.subcode .ne. antip) then
!            aergpn = aTrack.p.fm.p(4) - aTrack.p.mass
!         elseif(aTrack.p.code .eq. kgnuc) then
         aergpn = aTrack%p%fm%p(4) 
         if(aTrack%p%code .eq. kgnuc) then
!            aergpn = (aergpn-aTrack.p.mass) / aTrack.p.subcode
            aergpn = aergpn / aTrack%p%subcode
         endif
!        ----------------
         if( aTrack%p%code .le. kelec ) then
            jj = 1
         else
            jj = 3
         endif
         if( iTrack%p%code .le. kelec ) then
            ii = 1
         else
            ii = 3
         endif
!      
!      
         if(iergpn .gt. Ethin(ii)*relax ) then
            if(aergpn .gt. Ethin(jj)*relax) then
!                    Both   Ei, Ec> Ethin1; no thinning
               icon = 0
               aTrack%wgt = iTrack%wgt
!         elseif(aergpn .gt. Ethin(2)) then 
            else
#if EXACTWGT == 1
               p = aergpn/(Ethin(jj)*relax)
               if( iTrack%wgt/p .lt. Ethin(jj+1)*relax) then
#else
               if( aTrack%wgt .lt. Ethin(jj+1)*relax) then
                  p = aergpn/(Ethin(jj)*relax)
#endif

                  call rndc(u)
                  if(u .lt. p)  then
                     icon = 0
                     aTrack%wgt = iTrack%wgt / p
                  else
                     icon = 1
                  endif
               else
                  icon = 0
                  aTrack%wgt = iTrack%wgt 
               endif
            endif
         else
            if(aergpn .gt. Ethin(jj)*relax) then
               aTrack%wgt = iTrack%wgt
               icon = 0
            else
#if EXACTWGT == 1
               p = aergpn/iergpn
               if( iTrack%wgt/p .lt. Ethin(jj+1)*relax ) then
#else  
               if( aTrack%wgt .lt. Ethin(jj+1)*relax ) then
                  p = aergpn/iergpn
#endif

                  call rndc(u)
                  if(u .lt. p ) then
                     icon = 0
                     aTrack%wgt = iTrack%wgt / p
                  else
                     icon = 1
                  endif
               else
                  icon =0
                  aTrack%wgt = iTrack%wgt
               endif
            endif
         endif
      else
         icon = 0
         aTrack%wgt = iTrack%wgt
      endif

      end
