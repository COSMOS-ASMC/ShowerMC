      subroutine ciniTracking( incident )
!     this is called after primary is fixed.
      use modAtmosDef
      implicit none

#include "Zmanagerp.h"
#include "Ztrack.h"
! #include "Zmagfield.h"
#include "Ztrackp.h"
#include "Ztrackv.h"
#include "Zobs.h"
#include "Zobsp.h"
#include "Zobsv.h"
!  #include "Zatmos.h"
#include "Zincidentv.h"


      external cxyz2prim, cxyz2det
      integer klena, leng, icon
      real*8  zAngleSave, r1, r2, cosx, lengx
      real*8 clenbetween2h, cnewcos
      character*80 tracefile
      type(track)::incident
      type(coord)::xyz
      integer i
      data zAngleSave/-1.d30/
      save zAngleSave
      

      Zfirst%pos%depth = 0.
      Zsave = -1.d30
#if LABELING > 0
      Labelcounter = 0
#endif
! 
      if(Job .eq. 'flesh') then
         call csetEmin(Generate2, KEminObs2(1), KEmin2, KEminCas2)
         call csetEmin(Generate, KEminObs(1), KEmin, KEminCas)
      else
         call csetEmin(Generate, KEminObs(1), KEmin, KEminCas)
         KEmin2 = KEmin
         KEminCas2 = KEminCas
      endif


      call cprimxyz   ! compute primary x,y,z axis in 'xyz' system
      if((Trace .gt. 0 .and. Trace .lt. 60) .or. 
     *    (Trace .gt. 100 .and. Trace .lt. 160) ) then
!            default trace output is requested. Open disk file
         write(tracefile, *) TraceDir(1:klena(TraceDir))//'/trace',
     *         EventNo
         call kseblk(tracefile, ' ', leng)
         call copenfw(TraceDev, tracefile(1:leng), icon)
         if(icon .ne. 0) then
            call cerrorMsg('tracefile couldnot be opened',0)
         endif
      elseif(Trace .gt. 60 .and. Trace .lt. 100 .or.
     *       Trace .gt. 160 .and. Trace .lt. 200) then
         call cputCerenkovS     ! cerenkov output header for each event
      endif
!
      call ciniASObs
      if(abs(ObsPlane) .eq. horizontal  .or. 
     *        ObsPlane .eq. spherical) then
         call csetObsZ(cxyz2det)
      elseif(abs(ObsPlane) .eq. perpendicular ) then
         call csetObsZ(cxyz2prim)
!          inicident.where can be fixed here. not  in cmkIncidnt
         call cxyz2prim(ObsSites(NoOfSites)%pos%xyz, 
     *                  incident%pos%xyz, xyz)         
         do i = 1, NoOfSites
            if(xyz%r(3) .gt. ObsSites(i)%zpl ) then
               incident%where = i
               goto 222
            endif
         enddo
         incident%where = NoOfSites + 1
         IncidentCopy%where = incident%where
 222     continue
      elseif(ObsPlane .eq.  notUsed ) then
!           Nothing to do. observation is at height =
!           BorderHeighL. (for neutrino)
      else
         call cerrorMsg('ObsPlane value is wrong', 0)
      endif
!         check some parameters
      call  cexamParam
!         if one dimensional and large angle, use table method
!         for lenth <--> thickness conversion.  make table for that
      UseTbl = OneDim .eq. 3 .or.
     *  ( OneDim .eq. 2  .and. 
     *    abs( AngleAtObsCopy%r(3) ) .lt. 0.5 )
      UseTbl = UseTbl .and.  .not. Upgoing
      if(UseTbl .and. zAngleSave .ne. AngleAtObsCopy%r(3) ) then

         Htop =  atmos%z(atmos%nodes)    !  2016/Jun/20
         Hbase = atmos%z(1)              !
!         Hbase = HeightList(NoOfSites) - 0.3d0  ! make little bit smaller
!         Htop =100.d3           ! old value was 30d3

         r1 = Htop + Eradius
         r2 = Hbase + Eradius
         lengx = clenbetween2h(r2, r1, -AngleAtObsCopy%r(3))
         cosx = cnewcos(r2, -AngleAtObsCopy%r(3), lengx)
         call cl2tTbl(Htop, Hbase, cosx, -AngleAtObsCopy%r(3),
     *     LenStep, 
     *     LenTbl, HeightTbl, CosTbl,  ThickTbl, maxl2t, NumStep)
         zAngleSave = AngleAtObsCopy%r(3)
      endif
!          fix B at the starting point
      if( HowGeomag .le. 2 .or. HowGeomag .eq. 31 ) then
         call cgetMagF( YearOfGeomag, incident%pos%xyz, Mag, icon )
         call ctransMagTo('xyz', incident%pos%xyz,  Mag, Mag)
      else
         Mag = MagfieldXYZ
      endif

      end
!     --------- compute primary system, x,y,z.  
!               The deepest detector is the reference.
      subroutine cprimxyz
      implicit none
#include "Ztrack.h"
! #include "Zmagfield.h"      
#include "Zincidentv.h"
#include "Zobs.h"
#include "Zobsv.h"
!           all are unit vector in 'xyz' system.
!
!        X = Z x V  ( Z  is primary direction; V is vertical axis )
!          =  V x(-Z)
!          =  V x DcAtObsXyz = DetZaxis.r(1) DcAtObsXyz
!        Y = Z x X  = X x (-Z)
!          = X x DcAtObsXyz
!              
!      compute Xprimary, Yprimary, Zprimary in 'xyz' system.

      real*8  temp
!      
      Zprimary%r(:) = - DcAtObsXyz%r(:)
      call cvecProd(Zprimary, DetZaxis, Xprimary)
!         see if Zprimary // DetZaxis; if so reset Xprimary
      temp= sum( Xprimary%r(:)**2 )
      if(temp .lt. 1.e-12) then
         Xprimary = DetXaxis
      else
         temp =sqrt(temp)
!                ifort v14 dies here if (:) is used  
         Xprimary%r(1) = Xprimary%r(1)/temp         
         Xprimary%r(2) = Xprimary%r(2)/temp         
         Xprimary%r(3) = Xprimary%r(3)/temp         
      endif
      call cvecProd(Zprimary, Xprimary, Yprimary)
      Txyz2prim(1,:) = Xprimary%r(:)
      Txyz2prim(2,:) = Yprimary%r(:)
      Txyz2prim(3,:) = Zprimary%r(:)
      Tprim2xyz(:,:) = transpose(Txyz2prim(:,:))
       
      end
!     ******************************
      subroutine csetObsZ(converter)
!     *****************************
!      compute Z of observation plane 
!      
      implicit none
#include "Zcoord.h"
#include "Zpos.h"
#include "Zmagfield.h"
#include "Zobs.h"
#include "Zobsp.h"
#include "Zobsv.h"

      external converter  ! cxyz2prim or cxyz2det
!                           depending on |ObsPlane| =2, 1
!                           or  cxyz2det for ObsPlane=3
      integer i
      type(coord)::temp
!        compute z from the base detector.
      do i = 1, NoOfSites
         call converter(ObsSites(NoOfSites)%pos%xyz, 
     *     ObsSites(i)%pos%xyz,
     *     temp)
         ObsSites(i)%zpl = temp%r(3)
      enddo
!            this is for A.S
      do i = 1, NoOfASSites
         call converter(ASObsSites(NoOfASSites)%pos%xyz, 
     *     ASObsSites(i)%pos%xyz,
     *     temp)
         ASObsSites(i)%zpl = temp%r(3)
      enddo
      end
!      **************************
      subroutine cqFirstID(depth)
!         inquire the first interaction depth
      implicit none

#include "Ztrack.h"
! #include "Zmagfield.h"     
#include "Ztrackv.h"
      real*8 depth   ! output. to get the first I%D
      
      depth = Zfirst%pos%depth   ! vertical depth
      end
!      **************************
      subroutine cqFirstIPI(ptrack)
!         inquire the complete first interaction point info. 
!         as a track of the primary.
      implicit none
      
#include "Ztrack.h"
! #include "Zmagfield.h"
#include "Ztrackv.h"
      type(track)::ptrack   ! output.    if ptrack%pos%depth=0.
                              !  no interaction has been occurred
      
      ptrack = Zfirst
      end
      subroutine cqFirstColMedia(A, Z, xs)
!        retrns first col. element and Xsecion on it
      implicit none
#include "Ztrack.h"
! #include "Zmagfield.h"      
#include "Ztrackv.h"
      integer,intent(out):: A
      integer,intent(out):: Z
      real(8),intent(out):: xs
      A = FirstColA
      Z = FirstColZ
      xs = FirstColXs
      end
! --------------------------------
      subroutine ciniASObs
      implicit none

#include "Zcoord.h"
#include "Zpos.h"
#include "Zmagfield.h"
#include "Zobs.h"
#include "Zobsp.h"
#include "Zobsv.h"

      integer i

      do i = 1, NoOfASSites
         ASObsSites(i)%esize = 0.
         ASObsSites(i)%age = 0.
      enddo
      end
!     ***********************
      subroutine  csetEmin(gen, eminob,  emin, emCas)
!     ***********************
      implicit none
#include  "Zmanagerp.h"
#include  "Zcode.h"
#include  "Ztrack.h"
! #include "Zmagfield.h"      
#include  "Ztrackp.h"
#include  "Ztrackv.h"
#include  "Zincidentv.h"
#include  "Zheavyc.h"
#include  "Zheavyp.h"
!
      character*(*)  gen ! input. a copy of Generate. 
      real*8  eminob   !   input. minim observational energy (kientic)
      real*8  emin     !  output. Minimum energy to be followed for non e-g
                       !  ptcls.
      real*8  emCas    !  output. Minimum energy to be followed for e-g
!
!
      real*8  ergpn    ! energy / nucleon
      logical cas      ! 
      logical obas     ! 

!
      cas = index(gen, 'em') .gt. 0 
      obas = index(gen, 'as') .gt. 0 .or.
     *   index(gen, 'lat') .gt. 0


!          energy / nucleon for heavy
      if(IncidentCopy%p%code .eq. kgnuc) then
         ergpn =
     *        IncidentCopy%p%fm%p(4)/IncidentCopy%p%subcode
      elseif(IncidentCopy%p%code .ge. kalfa .and.
     *     IncidentCopy%p%code .le.  khvymax) then
!           will  not come here
         ergpn =
     *        IncidentCopy%p%fm%p(4)/Code2massN(IncidentCopy%p%code)
      else
         ergpn = IncidentCopy%p%fm%p(4)
      endif
      if(obas ) then
          ! wait until electon  energy becomes < EasWait for AS
          ! generation. 
         if( WaitRatio == 0. ) then !            default
            if( IncidentCopy%p%code == kphoton .or.    ! 
     *           IncidentCopy%p%code == kelec ) then   
               if( ergpn <= 100.d0) then ! E1ry < 100 GeV/n
                  EasWait = ergpn ! E1ry is too small for hyb AS: no flucutuation if electron 1ry.
               elseif( ergpn < 1.0d4 ) then ! E1ry < 10 TeV/n
                  EasWait = 100.d0 !  , Hyb. AS  generation  when Ee becomes < 100 GeV.
               else
                  EasWait = 0.01d0* ergpn ! > 10TeV. when e- energy becomes < 1/100 of E1ry, Hyb. AS calc
                                       ! is tried.
               endif
            else
               if( ergpn <= 500.d0) then ! E1ry < 100 GeV/n
                  EasWait = ergpn ! E1ry is too small for hyb AS
               elseif( ergpn < 5.0d3 ) then ! E1ry < 5 TeV/n
                  EasWait = 500.d0 !  Hyb. AS  generation  when Ee becomes < 500 GeV.
               else
                  EasWait = 0.1d0* ergpn ! > 5TeV. when elec energy becomes < 1/10 of E1ry, Hyb. AS calc
                                       ! is tried.  (Note:  Ee > 0.1 ergpn is rare)
               endif
            endif
         elseif( WaitRatio > 0.d0 ) then
                     ! input value is forced. no check for > 1. 
            EasWait =ergpn * WaitRatio
         else
            call cerrorMsg('WaitRaio value(<0) invalid ',0)
         endif
      endif


      emCas = 1.d30

      if(obas) then
         emCas = EasWait
      endif
      if(cas) then
         emCas =min(eminob, emCas)
      endif
      emin = eminob
      if(obas) then
         emin = min(emin, max(RatioToE0* ergpn, 1.d0) )
      endif 

      if(ThinSampling) then
!           thin sampling
         if(EthinRatio(1) .lt. 0.) then
            Ethin(1) = abs(EthinRatio(1))
         else
            Ethin(1) = EthinRatio(1)*ergpn
         endif
         if(EthinRatio(3) .lt. 0.) then
            Ethin(3) = abs(EthinRatio(3))
         elseif(EthinRatio(3) .eq. 0.) then
            Ethin(3) = Ethin(1)/10.
         else
            Ethin(3) = EthinRatio(3)*ergpn
         endif


!         if(EthinRatio(2) .lt. 0.) then
!            Ethin(2) = abs(EthinRatio(2))
!         else
!            Ethin(2) = EthinRatio(2)*ergpn
!         endif
!             maximum weight above which no thinning
          Ethin(2) =  EthinRatio(2)
          if( EthinRatio(4) .eq. 0.) then
             Ethin(4) =  Ethin(2)/10.
          else
             Ethin(4) = EthinRatio(4)
          endif
!///////////////////
!          write(0, *) ' ****** Ethin1,2=', Ethin
!          write(0, *) ' ergpn=', ergpn, ' EthinRatio=',EthinRatio
!////////////////

      else
           Ethin(1) = 0.
           Ethin(2) = 0.
           Ethin(3) = 0.
           Ethin(4) = 0.
      endif
      end
      subroutine cexamParam
      implicit none
#include "Ztrackp.h"
#include "Zelemagp.h"

      end
