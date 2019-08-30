      subroutine epLightCrossMatreshka(cnx, icon)
      use modepLightPty
      implicit none
#include "ZepTrackv.h"
#include "ZepTrackp.h"
#include "Zcnfig.h"
      integer,intent(inout):: cnx ! outside is this comp. #. 
                         ! if reflection happens, will become Cn
      integer,intent(out)::icon
!       current light hit matreshka


      real(8)::u   ! uniform random #
      integer::lcompno, mno, ptyno
      integer::jcon
       type(epDirec)::  normal
      real(8)::wrapperN  ! next component's refraction index.
                  !  it is not wrapper but we can use the same
                  !  routine for wrapper case by this def.
      real(8)::fzf  ! fuzzy 1 sigma (deg).
      real(8)::refp, refs, ref   ! polarization
       type(epDirec)::  dirp  !in.  original direction 
       type(epDirec)::  dirtemp
      integer::surfnx 

      lcompno = Det%cmp(cnx)%LightCompNo   ! light comp #

      if( lcompno > 0 ) then
         mno =  Det%Cn2media(cnx)   ! media no. of the next comp.
         wrapperN =Media(mno)%n    ! this is actully not a wrpper but we 
                ! may go parallel to the wapper case
         if( wrapperN == 0. ) then
            icon = 0
            write(0,*) ' Warning from epLightCrossMatreshka:'
            write(0,*)
     *       ' sensor? is contained inside another comp. (matreshka)'
            write(0,*) " it's stragne. light goes into there"
         else
             ! get refractive index of the comp. cnx as a funtion of w.l
            call  epLightwl2N(
     *        cTrack%wl,
     *        Lcomp( lcompno),
     *        comInfo( Lcomp(lcompno)%comInfoNo ),
     *        Media(mno)%gasF, 
     *        wrapperN)

            !  get surface # of matreshka
!         call epLightGetSurfN(cnx, Move.boundary,  surfnx)           
            call epLightGetSurfN(cnx, cTrack%pos,  surfnx)           
              ! get normal vector at the surface  of matreshka
            call epLightNormalVec(cnx, surfnx, normal)

                   !  get reflection ratio
            normal%x = - normal%x ! we have to revert the sign since normal vector
            normal%y = - normal%y ! is regarded as outgoing from the current comp. (not from matreshka)
            normal%z = - normal%z  

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
               jcon = 0
            else
              !    refraction
               call epLightRefrac( Lcomp( cLcompNo )%refracN,  wrapperN,
     *          normal,    Move%Track%w, cTrack%w, jcon)

           ! return value of jcon; 0 --> no error 1--> suspected to
           ! be total reflection, 2--> input error ?
           ! we neglect errors (message should have issued in the
           ! above routine. 
           ! next jcon=1   means refraction happend
          ! cTrack.w's  angle here is still local in the current comp
          !  convert it to cnx
               call epl2wd(Cn, cTrack%w,  dirtemp)
               call epw2ldm(cnx, dirtemp, cTrack%w, cTrack%p)
               jcon = 1
            endif

            dirp= cTrack%w
             ! fuzzify the angle if requested
            if( jcon  ==  0 ) then
              ! reflection.  from matreshka; we use surface info. of matreshka
              ! somewhat questionalbe to use surface info of matreshka.
               fzf= comInfo( Lcomp(lcompno)%comInfoNo )%fuzzy(surfnx) 
               if(fzf > 0.) then
         !           add some spread
                  call epLightFuzzify(fzf,  normal,  dirp, cTrack%w )
               endif
            else
         !   refracction.   can do same as above without doubt.
               fzf= comInfo( Lcomp(lcompno)%comInfoNo )%fuzzy(surfnx) 
               if(fzf > 0. ) then 
                  call epLightFuzzify(fzf,  normal,  dirp, cTrack%w )
               endif
            endif
            if(jcon .eq. 0 ) then
               cTrack%pos = Move%Track%pos
               cTrack%cn = Cn
               cnx = Cn
               Move%Cross = .false.
               icon = 0
            elseif(jcon .eq. 1) then
               icon = 0 
            else
               icon = 1
            endif
         endif
      else
         icon = 1
      endif
      end
