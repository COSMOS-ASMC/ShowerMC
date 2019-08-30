      subroutine  epLightSeeNextComp(cnx, surfN, icon)
      use modepLightPty
      implicit none
#include "ZepTrackv.h"
#include "ZepTrackp.h"
#include "Zcnfig.h"


      integer,intent(in)::cnx    ! comp # of next comp.
      integer,intent(in)::surfN  ! surface # ; in the current comp.  
      integer,intent(out)::icon  ! 0 reflection;  1 refraction; 
                                 ! -1 absorbed or go out 
      

      ! ----------------------------

      real(8)::u   ! uniform random #
      integer::lcompno, mno, ptyno
       type(epDirec)::  normal
      real(8)::wrapperN  ! next component's refraction index.
                  !  it is not wrapper but we can use the same
                  !  routine for wrapper case by this def.
      real(8)::fzf  ! fuzzy 1 sigma (deg).
      real(8)::refp, refs, ref   ! polarization
       type(epDirec)::  dirp  !in.  original direction 


      lcompno = Det%cmp(cnx)%LightCompNo   ! light comp #

      if( lcompno > 0 ) then
              ! 
         mno =  Det%Cn2media(cnx)
         wrapperN =Media(mno)%n    ! this is actully not a wrpper but we 
                ! may go parallel to the wapper case
         if( wrapperN == 0. ) then
            icon = 0
            return              ! **** should be sensor
         endif
             ! get refractive index of the comp. cnx as a funtion of w.l
!///////////////
!         call Lcompchk('next', lcompno)
!         call Lcompchk('nextp',  Lcomp(lcompno)%comInfoNo )
!         call Lcompchk('nextpp',  cPtyNo)
!///////////////
         call  epLightwl2N(
     *        cTrack%wl,
     *        Lcomp( lcompno),
     *        comInfo( Lcomp(lcompno)%comInfoNo ),
     *        Media(mno)%gasF, 
     *        wrapperN)

                   ! get normal vector at the surface  of current comp.            
         call epLightNormalVec(Cn, surfN, normal)

                   !  get reflection ratio
         call epLightReflecRatio( Lcomp( cLcompNo )%refracN,
     *       wrapperN,  normal,  Move%Track%w, refp, refs)

      !  refp: reflection prob. of parallel polarization wave 
      !  refs: reflection prob. of transverse polarizaton wave
      !           fix which polization. 
         call rndc(u)
         if( u .lt. Move%Track%pol  ) then
            ref = refp
            Move%Track%pol = 1.
         else
            ref = refs
            Move%Track%pol = 0.
         endif
         call rndc(u)
          !   reflection or refraction
         if(u .lt. ref) then
          !   mirror reflection:  current dir in Move is reflected
          !   about normal vec. and stored in cTrack.w
            call epLightMirrorRef(normal, Move%Track%w, cTrack%w)
            icon = 0
         else
              !    refraction
            call epLightRefrac( Lcomp( cLcompNo )%refracN,  wrapperN,
     *          normal,    Move%Track%w, cTrack%w, icon)
           ! return value of icon; 0 --> no error 1--> suspected to
           ! be total reflection, 2--> input error ?
           ! we neglect errors (message should have issued in the
           ! above routine. 
           ! next icon=1   means refraction happend
            icon = 1
         endif

         dirp= cTrack%w
          ! fuzzify the angle if requested
         fzf =  comInfo( cPtyNo )%fuzzy(surfN) 
         if( icon  ==  0 ) then
                                ! reflection.  surfN is valid 
            if(fzf .gt. 0. ) then
         !           add some spread
               call epLightFuzzify( fzf,  normal,  dirp, cTrack%w )
            endif
         else
          !  refraction.  we have to know the surface # of the cnx comp
          !  but this is not possible if the current comp. is
          !  a matreshka. so we use current comp's surface infor fuzzy.
          !  may be questionable and may be skipped.
            if( fzf  > 0. ) then
               call epLightFuzzify(fzf,  normal,  dirp, cTrack%w )
            endif
         endif
      else
         icon = -1
      endif

      end
