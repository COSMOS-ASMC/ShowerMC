!     Suppose an incident (1ry; say p) starting from the border of the
!     world (P0) and makes 1st hadronic collision at a gvien point (Pc).
!     Form P0 to PC,  we suppress hadronic collisions but make  all other
!     interactions take place as in normal simulation. When the 1ry
!     reaches Pc,a hadronic collision is forced to take place.  After
!     that, the simulation proceeds in the normaly way.
!     This program helps to do such type of simulaitons.
!
!     How to do:
!     We sample a 1ry of which starting point is Pc in a normal way.
!     Besides this, we specify FreeC =F in epicsfile.
!     In the old model, after the 1ry is fixed, a hadronic collision
!     is forced to  take place at Pc and the simulation continues 
!     disregarding the processes form P0 to Pc.
!     A) In the current model, P0 is computed by using direction cosine of
!        the 1ry and Pc with information of the world.  Then, the 1ry is
!        reset to start from P0. The place for this is just before calling
!        epgen in subroutine s1set of sepics.f
!     There are two essensial variables:
!     V1ry (integer)
!     Dist2colp  (distance from P0 to Pc in cm).
!
!      The default V1ry value 0 is made to be 1 here.  and Dist2colp
!     is computed.   
!     B) Then, simulation starts.  From P0 to Pc we must suppress hadronic
!       collsion.  So the 1ry from P0 to Pc is called virtual 1ry;
!     this is judged by looking at V1ry value(=1).
!     During virtual 1ry tracking, simulation is performed as if
!     normal 1ry is treated.  Distances for various physical processes
!     to take place are sampled and  the shortest distance one is
!     selected.  If the path is too long, truncation is made and
!     the 1ry is moved that distance. The shortest one may be employed
!     for case 1,2,3 below (x: current 1ry pos.  B: border of the
!     current component. * is the place the sampled  process to occur.) (Note:
!     if 1'), path is truncated at B and the process is not employed) 
!       
!    1)  x----*--B---Pc       1')  x-------B--*--Pc
!    2)  x----*--Pc----B
!    3)  x-------Pc--*-B      
!     The current  Dist2colp is from x to Pc.  When
!     the particle is moved, Dist2colp is updated everytime.
!     For case 1), LengthToB (from x to B) is < Dist2colp. 
!
!     For case 1) and 2), if the process sampled is not collision,
!     the process is simulated. If it is collision, it is neglected
!     and the 1ry is moved to *  (only dE/dx may be considered).
!     If case 3) happens, independently of the process,
!     we should make collsion at Pc after moving
!     1ry from x to Pc. So the 1ry is moved to Pc.
!     V1ry=2 is set at this moment.
!     During the tracking with V1ry=1, a particle (e-) may be  produced
!     by, say, knock-on, and then this is tracked. If V1ry is still 1,
!     Dist2colp may be affected, so we change V1ry=-1 until the incident
!     track appears again (V1ry=1 is resotred)  
      
      module modV1ry
      implicit none
!               the interaction to be forced at the
!     given point: select one of the list for each
!     1ry.
!     (hadint is for  p, n, pi, K, He... Fe..)
!     (muint   ;  muon)
!     (gint    :  gamma ray)
!     (eint    :  electron)      
      character(6):: hadint = "coll"  ! "decay", "coll"
      character(6):: muint = "nuci"   ! "brem", "pair", "nuci"
      character(6):: gint = "pair" ! "comp", "photop", "pair"
      character(6):: eint = "brem" ! "brem"
      logical:: FollowV1ry=.true.  ! if false. no virtual 1ry
!      (i.e no 1ry from P0 to PC).  This is the same  as
!     old FreeC =f.    Can be fixed by "epicsfile"
      
!     Other interactions might occure  between P0 and Pc but
!     the intraction specified by hadint etc is supprssed
!     during the viursual 1ry state.      
!      
      
      real(8),save:: Dist2colp
      real(8),save:: lengthToB
      integer,save:: V1ry
      
      contains
      subroutine epSetupV1ry(icon)
      implicit none
#include "Zmedia.h"      
#include "ZepTrackp.h"
#include "ZepTrackv.h"
#include "Zcode.h"
      integer,intent(out):: icon ! 0--> ok. else ng
      
      integer:: nw, ncp, n1ry, na
       type(epTrack)::  InciTrack ! original incident 

      character(16):: worldstruc
      character(12):: mediaName

       type(epTrack):: aTrack
      real(8):: vol(3), el
      integer:: kcon

      call epqworld(nw, worldstruc)  ! inquire # of world, 
! world struc.

      if( nw == 0 ) then
         write(0,*) ' For FreeC=F, world must be given'
         stop
      else
         call epqstn(n1ry)    ! # of stacked ptcls == # of 1ry here
         if( n1ry /= 1 ) then
            write(0,*) ' FreeC=F  but # of incidents=', n1ry, ' !=1'
            write(0,*) ' You must set only 1 incident in this case'
            stop
         endif
!     To get current incident,     
!     next call is too early.( epqinc is not yet ready)
        !         call epqinc( InciTrack ) !
!     Since this routine is to make a
!     virtual 1ry, some delicate treatment is needed
!     and this is to get track info. from the stack area
!
         call epgetTrack(1, InciTrack, kcon)
         if(kcon /= 0 ) then
            write(0,*) ' 1ry has not been stacked '
            stop
         endif

!     at this moment, cTrack, Cn has not been set so
!     we can use  them freely
         Cn = InciTrack%cn
         call epqCn2MediaName(Cn, mediaName)
         if( mediaName == "sp" .or. mediaName == "hollow") then
            write(0,*)
     *      'FreeC=F but madium @ input collision point is'
            write(0,*) trim( mediaName ), ' invalid'
            stop
         endif
!!!!!
!     if( InciTrack.p.code == kelec ) then
         if(.not. FollowV1ry ) then
              ! old style use.
            V1ry = 1
            Dist2colp = 0.
            icon = 0
            return              !!!!!!!!!!!!!!!!
         endif   
!!!!!!!!!!
         
         cTrack = InciTrack
!          invert the direction
         cTrack%w%x = -cTrack%w%x
         cTrack%w%y = -cTrack%w%y
         cTrack%w%z = -cTrack%w%z
         call epqncp(ncp)       ! Total # of comp.  ->world         
         call epbndry2(ncp, el, kcon) ! get length to the world
                   ! border.  (el)
         if( kcon /= 0 ) then
            write(0,*)
     *      ' Strange: V1ry with given col. point has no '
            write(0,*) ' crossing point with world '
            stop
         endif

         Dist2colp = el-Epsleng/2
! execept for the  starting point and Cn, everything should
!     be the same so change the  starting point
!              in local coord. is
         InciTrack%pos%x = InciTrack%pos%x + cTrack%w%x*Dist2colp
         InciTrack%pos%y = InciTrack%pos%y + cTrack%w%y*Dist2colp
         InciTrack%pos%z = InciTrack%pos%z + cTrack%w%z*Dist2colp
!     convert to world
         call epl2wTrack(InciTrack, InciTrack)
!            clear stack
         call epempty
!
!     With   InciTrack.cn = ncp (-->world),  we may inclined to
!     use next call:
!       call eppush(InciTrack)  
!     Howver, this   is N.G:
!     If the box_w is tight and closely contact to
!     an inner component, the starting point could be already
!     inside that compo. but the point could be judged inside the
!     box_w (which is normally "sp"). So energy loss, knock-on
!     etc will not happen. So we ask to find the component in which
!     the point exists by the next call and stack it.
         call epputTrack(InciTrack)
         V1ry = 1   ! virtual 1ry tracking starts
      endif
      icon = 0
      end subroutine epSetupV1ry
      
      
      subroutine epqmyinc(type)
      implicit none
#include "Zcode.h"
#include "ZepTrack.h"
      character(1),intent(out)::type
       type(epTrack):: inci
      call epqinc(inci)
      if( inci%p%code == knuc .or. inci%p%code == kgnuc ) then
         type = "h"
      elseif( inci%p%code == kpion .or. inci%p%code == kkaon ) then
         type = "h"
      elseif( inci%p%code == kphoton)  then
         type = "g"
      elseif(inci%p%code == kmuon ) then
         type = "m"
      elseif (inci%p%code == kelec ) then
         type = "e"
      else
         type = "o"
      endif
      end subroutine epqmyinc
      end module  modV1ry

      subroutine epSetIntType(ptype, itype)
      use modV1ry
      implicit none
      character(1),intent(in):: ptype ! ptcl type.
            ! one of h, g, e, m
      character(*),intent(in):: itype ! intracton type
!           to be forced to occur at Pc. 
!
      if( ptype == "h" ) then
         hadint = itype
      elseif( ptype == "g" ) then
         gint = itype
      elseif( ptype == "e" ) then
         eint = itype
      elseif( ptype == "m" ) then
         muint = itype
      else
         write(0,*) ' ptype=',ptype, ' to epSetIntType '
         write(0,*) ' in epSetupV1ry is invalid'
         stop
      endif
      write(0,*) 'epSetIntType is called for V1ry'
      write(0,*) ' ptype =', ptype, ' itype=',itype,' sepcified'
      end subroutine epSetIntType

      subroutine epSetFollowV1ry(fl)
      use modV1ry
      implicit none
      logical,intent(in):: fl
      FollowV1ry = fl
      end subroutine epSetFollowV1ry

      subroutine epManageInciflag(aTrack, bTrack)
      implicit none
#include  "Zmedia.h"            
#include  "ZepTrackv.h"
#include  "ZepStack.h"
#include  "ZsepManager.h"
#include  "Zcode.h"
       type(epTrack)::  aTrack ! input ptcl to be stacked
       type(epTrack)::  bTrack ! copy of aTrack
! after particle production, it is rather difficult
!     to judge  which particle is the incident.
!     For practical purpose, we regard a particle
!     the incdient if it satisfies the next.
! 1)    Non e-/e+  particle makes knock-on
!     Then, the non e-/e+ after knock-on keeps the inciflag
! 2)    Muon makes an inteaction, then surviving muon keeps
!     inciflag.  (muon direct production, if any, is disregarded)
! 3)  e+ makes knock-on.  Then e+ keeps inciflag
! 4)  e- makes knock-on.  Higher energy one keeps inciflag           
! 5)  e+,e-. brems.  e+e- keep inciflag
!     6)  gamma makes compton scattering.  gamma keeps inciflag
!     7)
      if( Move%Track%p%code == kelec ) then
         if( Move%proc == "brem" ) then
            if(bTrack%p%code == kphoton ) then
               bTrack%inciflag = 0 !     clear inciflag
            else
               bTrack%inciflag = 1
            endif
         elseif( Move%proc == "knoc" ) then
            if(Move%Track%p%charge == 1) then
                     ! e+
               if( bTrack%p%charge == -1 ) then
                  bTrack%inciflag = 0
               else
                  bTrack%inciflag = 1
               endif
            else
                  ! e-
               if( bTrack%p%fm%p(4)*2  <
     *              Move%Track%p%fm%p(4) + Move%Track%p%mass ) then
                  bTrack%inciflag = 0
               else
                  bTrack%inciflag = 1
               endif
            endif
         else
            bTrack%inciflag = 0
         endif
      elseif( Move%Track%p%code == kphoton ) then
         if( Move%proc == "comp" ) then
            if( bTrack%p%code == kelec ) then
               bTrack%inciflag = 0
            else
               bTrack%inciflag  =1
            endif
         else
            bTrack%inciflag = 0
         endif
      elseif( Move%Track%p%code == kmuon ) then
         if( bTrack%p%code /= kmuon ) then
            bTrack%inciflag = 0
         else
            bTrack%inciflag = 1
         endif
      else                      ! hadorons 
         if(Move%proc == "knoc" ) then
            if( bTrack%p%code == kelec ) then
               bTrack%inciflag = 0
            else
               bTrack%inciflag = 1
            endif
         elseif( Move%proc == "elas") then
!  "elas" is not yet used
         else
            bTrack%inciflag = 0
         endif
      endif
#if defined (OLDDEFINITION)
      if(   Move%Track%p%code /= kelec 
     *     .and. Move%proc == "knoc" 
     *     .and. bTrack%p%code /= kelec ) then
!  keep the current inciflag 
      else
            !  clear inciflag
         bTrack%inciflag = 0
      endif
#endif      
      end subroutine epManageInciflag
            !    nexts are for printing debuggin info. 
      subroutine epPrTrackInfo(id, aTrack)
      implicit none
#include "ZepTrack.h"
      character(*),intent(in):: id
       type(epTrack)::   aTrack
      call epPrPtclInfo(' ', aTrack%p)
      write(0,*) ' Cn=', aTrack%cn, ' inciflag=', aTrack%inciflag
      write(0,*) ' pos=',aTrack%pos%x, aTrack%pos%y, aTrack%pos%z
      write(0,*) ' dir=',aTrack%w%x, aTrack%w%y, aTrack%w%z
      end subroutine epPrTrackInfo

      subroutine epPrPtclInfo(id, aPtcl)
      implicit none
#include "ZepTrack.h"
      character(*),intent(in):: id
       type(ptcl)::   aPtcl
      if( trim(id) /= "") then
         write(0,*) '----- ',trim(id), ' ---- ptcl info'
      endif
      write(0,*) 'code ,subc, chg=',
     *     aPtcl%code, aPtcl%subcode, aPtcl%charge
      write(0,*) ' mass =', aPtcl%mass
      write(0,*) 'p(:) =', aPtcl%fm%p(1:4)
      end subroutine epPrPtclInfo

      subroutine eptruncV1ry
      !  trucation or force interaction
      use modV1ry
      implicit none
#include  "ZmediaLoft.h"                  
#include  "ZepTrackv.h"
#include  "ZepTrackp.h"
#include  "Zcode.h"
      
      if( Dist2colp <= lengthToB ) then
         if( Move%dl >= Dist2colp )  then
! case 3)  Be ready for specified interaction after moving
!             Dist2colp            
            Move%dl = Dist2colp 
            Move%dx = Move%dl/Media(MediaNo)%gtocm *
     *           Media(MediaNo)%rhoc
            Move%dt = Move%dx/ Media(MediaNo)%X0g
            if(Move%Track%p%code == knuc
     *        .or. Move%Track%p%code== kgnuc) then
               Move%proc = hadint
            elseif( Move%Track%p%code == kphoton ) then
               Move%proc = gint
            elseif( Move%Track%p%code == kelec ) then 
               Move%proc = eint
            elseif( Move%Track%p%code == kmuon ) then
               Move%proc = muint
            elseif(Move%Track%p%code == kpion
     *        .or. Move%Track%p%code== kkaon) then
               Move%proc = hadint
            else
               write(0,*) ' At present particle code ',
     *         Move%Track%p%code, ' cannot be used with FreeC=f'
               stop
            endif
            V1ry = 0      !from now,  do everything as if normal one
            if( Move%proc == hadint )  then
                ! hadron int. uses Cosmos info. so we must call this
               call epResetProcNoForV1ry   
            endif
         else
! case 2
!            if( Move%proc == hadint .or.
!     *           Move%proc == gint .or.
!     *           Move%proc == eint .or.
!     *           Move%proc == muint  ) then
            if( Move%proc /= "knoc" ) then
               Move%proc = "xx"
               Move%Trunc = .true. ! though could be already truncated.
            else
            !      write(0,*)
     *      !         'A) V1ry=1: proc =',Move.proc, ' accepted'
            endif
         endif
      else
!     case 1)
!         if( Move%proc == hadint .or.
!     *        Move%proc == gint .or.
!     *        Move%proc == eint .or.
!     *        Move%proc == muint  ) then
         if( Move%proc /= "knoc" ) then
            Move%proc = "yy"
            Move%Trunc = .true.
         else
                  !!!!!!!!
!            write(0,*)
!    *            'B) V1ry=1: proc =',Move.proc, ' accepted'
         endif
      endif
      
      end subroutine eptruncV1ry

      subroutine epForceV1ryInt
      use modV1ry
      use modSetIntInf
      implicit none
#include  "Zglobalc.h"
!#include  "Ztrack.h"
!#include  "Ztrackv.h"
      integer i
      character(1):: type
!           force the v1ry interaction
      call epqmyinc(type)
      if( type == "h" ) then
         do i = 1, NumberOfInte
            if( IntInfArray(i)%process == hadint ) then
               IntInfArray(i)%thickness = 0.
               exit
            endif
         enddo
      elseif( type == "g" ) then
         do i = 1, NumberOfInte
            if( IntInfArray(i)%process == gint ) then
               IntInfArray(i)%thickness = 0.
               exit
            endif
         enddo
      elseif( type == "m" ) then
         do i = 1, NumberOfInte
            if( IntInfArray(i)%process == muint ) then
               IntInfArray(i)%thickness = 0.
               exit
            endif
         enddo
      elseif( type == "e" ) then
         do i = 1, NumberOfInte
            if( IntInfArray(i)%process == eint ) then
               IntInfArray(i)%thickness = 0.
               exit
            endif
         enddo
      else
!     no specifcation
         write(0,*) ' V1ry type = ', type
!     stop
      endif
      end subroutine epForceV1ryInt
#define NORMALFirstInt
#if defined (NORMALFirstInt)
!     regards every 1ry interaction except knock-on 
!     as  first interaction      
      subroutine epSeeIf1stInt
      use modV1ry
      implicit none
#include  "Zmedia.h"
#include  "ZepTrackv.h"
#include  "ZepTrackp.h"
#include  "Zcode.h"
#include  "Zevhnv.h"

      integer:: inela
      if( Move%proc /= "knoc" ) then
         if( cTrack%inciflag == 1 ) then
            FirstInt = cTrack%pos
            FirstIntTrack =cTrack
            Proc1 = Move%proc

               ! hadron
            if(Move%proc == 'coll') then
                ! in jam we should introduce 'elas' for proc
                ! in future, for elastic scattering and avoid
                ! next 
               if(ActiveMdl == 'jam') then
                 !     if hadron proj, inela could be 0 or 1
                  call cjamElaInfo(1, inela)
               else
                  inela = 1
               endif
               if( inela == 1) then
                  FirstC=.false.
                  call epsaveFirstCol ! only for nuclear coll.
!     saved info canbe obtained by
!     #include "Zptcl.h"
!          integer,parameter::n =...
!          type(ptcl):: ptcls(n)
!          call  epqFirstColPtcls(ptcls, n, m)
!     then,
!          m-particls  are put in ptcls 
               endif
            else
               FirstC=.false.
            endif
         endif
      endif
      end subroutine epSeeIf1stInt
#else
!     regards 1ry 1st interaction as first interaction if
!     it is specified by hadint, gint, eint,  muint
      subroutine epSeeIf1stInt
      use modV1ry
      implicit none
#include  "Zmedia.h"      
#include  "ZepTrackv.h"
#include  "ZepTrackp.h"
#include  "Zcode.h"
#include  "Zevhnv.h"

      integer:: inela
      
      if( cTrack%inciflag == 1 ) then
         if(Incident%p%code .eq. kelec ) then
            if( Move%proc == eint ) then
               FirstInt = cTrack%pos
               FirstIntTrack =cTrack
               Proc1 = Move%proc
               FirstC=.false.
            endif
         elseif( Incident%p%code .eq. kphoton ) then
            if( Move%proc == gint ) then
!                only pair is employed
               FirstInt = cTrack%pos
               FirstIntTrack =cTrack
               Proc1 = Move%proc
               FirstC=.false.
            endif                  
         elseif(cTrack%p%code .eq. kmuon ) then
            if( Move%proc == muint  ) then
               FirstInt = cTrack%pos
               FirstIntTrack =cTrack
               Proc1 = Move%proc
               FirstC=.false.
            endif
         else
               ! hadron
            if(Move%proc ==  hadint ) then
               if(Move%proc == 'coll') then
                ! in jam we should introduce 'elas' for proc
                ! in future, for elastic scattering and avoid
                ! next 
                  if(ActiveMdl == 'jam') then
!     if hadron proj, inela could be 0 or 1
                     call cjamElaInfo(1, inela)
                  else
                     inela = 1
                  endif
               else
                  inela = 1
               endif
               if( inela == 1 ) then
                  FirstInt = cTrack%pos
                  FirstIntTrack = cTrack
                  Proc1 = Move%proc
                  FirstC=.false.
                  call epsaveFirstCol ! only for nuclear coll. 
               endif
            endif
         endif
      endif
      end subroutine epSeeIf1stInt
#endif
