      subroutine  epLightAtBndry(  cnx,  icon )
        ! cTrack, Move and Cn are implicit  in/out param.
        !     Move.Track at x      x|*
        !                           |
        !                           |boundary
        !     cTrack's pos is at *
        ! In conclusion, icon = 0 if reflection or refraction (
        !   including simple pass-trough) happens.
        !   icon = 1, if light dies at the boundary by some reason.
        ! Other effect:  (-- means unchagned)
        !               cTrack        Move.Cross     cnx     Cn
        ! refraction.   new angle       --           --      --  
        ! reflection.   new angle       F            Cn      --
        !               pos=Move's  
        ! absorbed        --            F
        ! pass-thru       --            --           --      -- 

        !   ** angle and position are in the local coordinate   
        !
        !   -----------------cases shown above happens when
        !   wrapper is reflector:  reflection happens or
        !              no reflection (i.e, absorbed )
        !   wrapper is grease/thin transparent medium:
        !              reflection happens  or  refraction happens
        !   no wrapper and outside is "light" material & refraction inex
        !              is given  reflection or refraction takes place
        !   no wrapper and outside is  "light" material 
        !              but refraction index is 0 (Say Si)
        !              This should be a light sensor; however, if outside
        !              is a matreshka, it should not  be a sensor, so
        !              waring is issued.  In any case, the light
        !              passes thru the boundary.
        !   no wrapper and outside is not "light" material
        !              we don't follow the light  (absorbed). 
  
      use modepLightPty
      implicit none
#include "ZepTrackv.h"
#include "ZepTrackp.h"
#include "Zcnfig.h"
      integer,intent(inout):: cnx ! outside is this comp. #. 
                         ! if reflection happens, will become Cn
      integer,intent(out)::icon

      integer::jcon

      if(CrossMode == 0 ) then
         call epLightCrossSelfBoudary(cnx, jcon)
      elseif(CrossMode == 1 ) then
         call epLightCrossMatreshka(cnx, jcon)
      elseif(CrossMode == 2) then
         call epLightCrossPContainer(cnx, jcon)
      else
         write(0,*) ' CrossMode=', CrossMode, 'strange'
         write(0,*) ' detected at epLightAtBndry; cnx=',cnx
         stop
      endif
      icon = jcon
      end
      subroutine epLightGetSurfN(compno, aPos, surfn)
          ! get surface # where aTrack is on in component with
          ! compno.
      use modepLightMaxDef
      implicit none
#include "ZepTrackv.h"
#include "ZepTrackp.h"
#include "Zcnfig.h"

      integer,intent(in)::compno  ! component #
       type(epPos)::  aPos     ! input. track of which pos is examined
      integer,intent(out):: surfn  ! surface # obtained

      real(8)::vola(maxVolAttr)
      integer::na
              ! get volume attribute of current component #
      if(compno > Det%nct ) then
            !  compno is not in the ranage
         surfn = 0
      else
         call epqvolatr( compno, na, vola)
         if(na > maxVolAttr ) then
            write(0,*) ' too many vol. attributes for structure='
            write(0,*)  Det%cmp(Cn)%struc
            stop
         endif
         if( Det%cmp(compno)%struc(1:3) == "box" ) then
            call epLightGetBoxSurfN(vola(1), vola(2), vola(3),
     *           aPos%x, aPos%y, aPos%z,
     *           surfn )

         elseif( Det%cmp(compno)%struc(1:7) == "octagon" ) then
            call epLightGetOctagonSurfN( Det%cmp(compno), aPos, surfn )

         elseif( Det%cmp(compno)%struc(1:3) == "cyl" ) then
            call epLightGetCylSurfN( Det%cmp(compno), vola(1), vola(2),
     *                  aPos, surfn )

         elseif( Det%cmp(compno)%struc(1:4) == "ecyl" ) then
            call epLightGetEcylSurfN(vola(1), vola(2), vola(3),
     *           aPos%x, aPos%y, aPos%z,
     *           surfn )

         elseif( Det%cmp(compno)%struc(1:4) == "pipe" ) then
            call epLightGetPipeSurfN(Det%cmp(compno), vola, aPos, surfn)
         elseif( Det%cmp(compno)%struc(1:5) == "prism" ) then
            call epLightGetPrismSurfN(Det%cmp(compno),  aPos, surfn)
         else
            write(0,*) ' struc=',Det%cmp(compno)%struc, ' not yet'
            write(0,*) ' supported  in epLightReflection'
            stop
         endif
      endif
      end

      subroutine epLightReflector(surfN, jcon) 
        !  judge if the light is reflected or absorbed by the wrapper
        !  if reflection, judge mirror or diffusive and 
        !  fix the angle and put it in cTrack
        !  commmon to all structures at the given surface with surfN
      use modepLightPty
      implicit none
#include "ZepTrackv.h"


      integer,intent(in)::surfN  ! surface #
      integer,intent(out)::jcon  ! 0 reflection;  1 absorbed
      


      real(8)::u   ! uniform random #
      real(8)::reflectance  ! wrapper reflection coef.
       type(epDirec)::  normal
!////////
!      call Lcompchk(' A ', cLcompNo)
!      call Lcompchk(' Ap ', cPtyNo)
!////////
!      if( Lcomp( cLcompNo )%wrapperReflecAF == 0 ) then
      if( comInfo(cPtyNo)%wrapperReflecAF == 0 ) then
            !  reflectance indepenent of w.l
         reflectance = comInfo( cPtyNo )%wrapper(surfN)
      else
           !  get w.l dependent reflectance
!         call csampAFintp( Lcomp( cLcompNo )%wrapperReflecAF, 
         call csampAFintp( comInfo( cPtyNo )%wrapperReflecAF, 
     *               Move%Track%wl, reflectance)
      endif
      
      call rndc(u)
      if(u < reflectance) then
         ! reflection
         jcon = 0
          ! get normal vector at the surface             
         call epLightNormalVec(Cn, surfN, normal)
         call rndc(u)
         if( u < comInfo( cPtyNo)%mirror(surfN) ) then
            !  mirror reflection
            call epLightMirrorRef(normal, Move%Track%w, cTrack%w)
         else
            ! diffusive reflection
            call epLightDiffRef(normal, Move%Track%w, cTrack%w)
         endif
      else
         ! absorbed
         jcon = 1
      endif
      end


      subroutine epLightReflecOrRefrac(cnx, surfN, jcon)
        !  judge if the light is reflected or refracted
        !  commmon to all structures at the given surface with surfN
      use modepLightPty
      implicit none
#include "ZepTrackv.h"
#include "Zcnfig.h"
      integer,intent(in)::cnx    ! comp# of the next comp (may be void)
      integer,intent(in)::surfN  ! surface #
      integer,intent(out)::jcon  ! 0 reflection;  1 refraction 
      
      real(8)::u   ! uniform random #
      real(8)::wrapperN  ! wrapper refraction  index.
      real(8):: refp, refs, ref
      real(8):: fzf      ! fazzy angle 1 sigma deg
      integer::ptyno, compno
       type(epDirec)::  dirp, normal, dirtemp

!////////////
!      write(0,*) ' cLcompNo =', cLcompNo, ' cPtyNO=', cPtyNo,
!     *  ' N=',     comInfo( cPtyNo )%wrapper(surfN )
!      call Lcompchk(' B ', cLcompNo)
!      call Lcompchk(' Bp ', cPtyNo)
!////////

      if( comInfo( cPtyNo)%wrapperRefracAF == 0 ) then
                 ! w.l independent refraction index of wrapper
         wrapperN = comInfo( cPtyNo )%wrapper(surfN )
      else
          ! get refraction index of wrapper at this w.l
         call csampAFintp( comInfo( cPtyNo)%wrapperRefracAF, 
     *               Move%Track%wl, wrapperN)
      endif
!////////////
!      write(0,*) ' wrappeN=',wrapperN, Lcomp( cLcompNo )%refracN
!/////////
        ! get normal vector at the surface             
      call epLightNormalVec(Cn, surfN, normal)
        !  get reflection ratio     refracN --> wrapperN
!//////////
!      call Lcompchk(' C ', cLcompNo)
!////////

      call epLightReflecRatio( Lcomp( cLcompNo )%refracN,  wrapperN,
     *    normal,  Move%Track%w, refp, refs)

      !  refp: reflection prob. of parallel polarization wave 
      !  refs: reflection prob. of transverse polarizaton wave
      !           fix which polization. 
!////////////
!      write(0,*) ' refp=',refp, ' refs=',refs
!//////////////
      call rndc(u)
      if( u .lt. Move%Track%pol  ) then
         ref = refp
         Move%Track%pol = 1.
      else
         ref = refs
         Move%Track%pol = 0.
      endif
      call rndc(u)
         !      reflection or refraction
      if(u .lt. ref) then
         !     mirror reflection:  current dir in Move is reflected
         !     about normal vec. and stored in cTrack.w
         call epLightMirrorRef(normal, Move%Track%w, cTrack%w)
         jcon = 0
      else
         !   refraction
!//////////
!      call Lcompchk('D  ', cLcompNo)
!////////

         call epLightRefrac( Lcomp( cLcompNo )%refracN,  wrapperN,
     *    normal,    Move%Track%w, cTrack%w, jcon)

           ! return value of jcon; 0 --> no error 1--> suspected to
           ! be total reflection, 2--> input error ?
           ! we neglect errors (message should have issued in the
           ! above routine. 
           ! next jcon=1   means refraction happend
         !   cTrack's direction is still in the current comp. coord.
         !   so transorm to the next comp.  only angle need to be conv.
         call epl2wd(Cn, cTrack%w,  dirtemp)
!cc         call epw2ld(cnx, dirtemp, cTrack.w)
         call epw2ldm(cnx, dirtemp, cTrack%w, cTrack%p)
         jcon = 1
      endif

      dirp= cTrack%w
         ! fuzzify the angle if requested
      if(jcon == 0 ) then
         ! reflection.  easy to introduce randomness
         if( comInfo( cPtyNo )%fuzzy(surfN) >  0. ) then
         !           add some spread
            call epLightFuzzify(comInfo( cPtyNo )%fuzzy(surfN),
     *           normal,  dirp, cTrack%w )
         endif
      else
         ! jcon =1; refracion.
         ! we don't consider fuzzy factor in the wrapper
         ! since it is adjacent to light sensor if outside is
         ! not void
!         if(cnx <= Det.nct) then
!            call epLightGetSurfN(cnx, cTrack.pos,  surfN)
!            compno = Det.cmp(cnx).LightCompNo
!            if(compno > 0 ) then
!               ptyno = Lcomp( compno )%comInfoNo
!               fzf = comInfo( ptyno )%fuzzy(surfN)
!               if( fzf  > 0. ) then
!                                !        add some spread
!                  normal.x = -normal.x
!                  normal.y = -normal.y
!                  normal.z = -normal.z
!                  call epLightFuzzify(fzf,  normal,  dirp, cTrack.w )
!               endif
!            endif
!         endif
      endif

      end

      subroutine epLightReflecRatio(n1, n2, normal,  dir, refp,refs)
      implicit none
#include "ZepDirec.h"      
      !  media 1 --> 2; reflection or refraction
      real*8 n1  !  refraction. index of media 1
      real*8 n2  !  //                         2
       type(epDirec)::  normal !input normal unit vetctor at the surface 
                   ! (directed to the outward of the volume)
       type(epDirec)::  dir
!         ! input  direction cos of the light given in the 
!               canonical coordinate. 
      real*8 refp ! reflection ratio of light with perpendicular polarization
      real*8 refs ! //                            parallel //
!          perpendicular: perpendicular to the  incident plane (formed by
!          incdient and reflected light).
!          parallel: parallle to the boundary plane
!
      real*8 nn, ns 
      real*8 cost, sint2, temp  
!        l= sint*cosf
!        m = sint*sinf 
!        n = cost
!      
!      if( dir.z  .ge. 0.) then
!         write(0,*) 'error input to eprefDirCos; must n<0 ' 
!         write(0,*) 'n1, n2, l,m,nr=',
!     *        n1, n2, dir
!         stop
!      endif
!      cost = -dir.z
!      cost = abs(normal*dir)  ! inner prod.
      call cscalerProd(normal, dir, cost)
!////////
!      write(0,*) ' n1=',n1, ' n2=',n2
!      write(0,*) ' cos of light angle=',cost
!///////
!      sint2 = dir.x**2 + dir.y**2
      sint2 = 1. -cost**2
      nn = n2/n1
      ns = nn*nn
      temp = ns - sint2
!////////////
!      write(0,*) ' sint2=',sint2
!      write(0,*) ' n2/n1=',nn, ' nn^2=',ns
!      write(0,*) ' ns-sint2=',temp
!/////////
      if(temp .le. 0.) then
         refp = 1.
         refs = 1.
      else
         temp = sqrt(temp)
         refp = ( (ns *cost - temp)/(ns*cost +temp) )**2
         refs = ( (cost -temp)/(cost + temp) )**2
      endif
      end

      subroutine epLightNormalVec(cnx, surfN, normal)
         ! get normal vector at the surface (normal vector
         ! is directed to the outward (at the outside surface)
      use modepLightMaxDef
#include "Zcondc.h"
      implicit none
#include "ZepTrackv.h"
#include "ZepTrackp.h"
#include "Zcnfig.h"

      integer,intent(in):: cnx ! component #
      integer,intent(in):: surfN  ! surface # of the comp.
       type(epDirec)::   normal  ! output . normal vector of the current
            ! comp. at surface # surfN. directed to outwards of the vol.
      
         !    ------------------
       type(epDirec)::  dir(6)  ! normal vector at six surfaces of 
!             the canonical box
!             out going direction
      data dir(1:6)/
     * epDirec(0., 0., -1.), epDirec(1., 0., 0.), epDirec(0., -1.,  0.), 
     * epDirec(0., 1., 0.),  epDirec(-1., 0., 0.), epDirec(0., 0.,  1.)
     * /
      save dir

      real(8):: temp, cost, sint, tant, x, y
      integer:: na
      real(8):: vol(maxVolAttr)
      real(8):: a, b, c, h, d 
#if defined MATHLOUSY
      real(8),parameter::cos45=7.071067811865475244d-1
#else
      real(8),parameter::cos45=sqrt(2.d0)/2.d0 
#endif


      if(Det%cmp(cnx)%struc(1:3) == "box" ) then
         normal = dir(surfN)
      elseif(Det%cmp(cnx)%struc(1:7) == "octagon" ) then
         if(surfN <= 6) then
            normal = dir(surfN)
         else
!            call epqvolatr(cnx, na, vol)
!            b = vol(2)
!            c = vol(3)
!            d = vol(4)

            if(surfN == 7 ) then
               normal = epDirec(0.d0, -cos45,  cos45)
            elseif( surfN == 8) then
               normal = epDirec(0.d0,  -cos45, -cos45)
            elseif (surfN == 9 ) then
               normal = epDirec(0.d0,  cos45,  -cos45)
            elseif (surfN == 10 ) then
               normal = epDirec(0.d0, cos45,  cos45)
            else
               write(0,*) ' surfN=', surfN,' is  invalid '
               write(0,*) 'for octagon in epLightNormalVec'
               write(0,*) 'cn =', cnx
            endif
         endif
      elseif(Det%cmp(cnx)%struc(1:3) == "cyl" ) then
         if(surfN == 1) then
            normal = dir(1)
         elseif(surfN == 2 ) then
            normal = dir(6)
         else
            temp = sqrt( Move%Track%pos%x**2 + Move%Track%pos%y**2)
            normal%x = Move%Track%pos%x/temp
            normal%y = Move%Track%pos%y/temp
            normal%z = 0.
         endif
      elseif(Det%cmp(cnx)%struc(1:3) == "ecyl" ) then
         if(surfN == 1) then
            normal = dir(1)
         elseif(surfN == 2 ) then
            normal = dir(6)
         else
            x = Move%Track%pos%x
            y = Move%Track%pos%y
            call epqvolatr(cnx, na, vol)
            if( abs(y) <= Epslength2  ) then
               cost =0.
               sint =sign(1.d0, x) 
            else
               tant =-(vol(2)/vol(1))**2 * x/y
               cost = sign( 1.d0/sqrt(1.+tant**2), y)
               sint = sign(tant*cost, x)
            endif
            ! rotate -90 deg (tangential vec=(-cos,sin)
            normal%x = sint
            normal%y = cost
            normal%z = 0
         endif
      elseif (Det%cmp(cnx)%struc(1:4) == "pipe" ) then
         if(surfN == 1) then
            normal = dir(1)
         elseif(surfN == 2 ) then
            normal = dir(6)
         elseif(surfN == 3) then
           ! inner wall
            temp = sqrt( Move%Track%pos%x**2 + Move%Track%pos%y**2)
            normal%x = -Move%Track%pos%x/temp
            normal%y = -Move%Track%pos%y/temp
            normal%z = 0.
         elseif( surfN == 4 ) then
            ! outer wall
            temp = sqrt( Move%Track%pos%x**2 + Move%Track%pos%y**2)
            normal%x = Move%Track%pos%x/temp
            normal%y = Move%Track%pos%y/temp
            normal%z = 0.
         else
            write(0,*) ' surface #=', surfN, ' strage for pipp'
            stop
         endif
      elseif (Det%cmp(cnx)%struc(1:5) == "prism" ) then
!       surf #:   bottom 1:  \ 2:  / 5: x-z @y=0: 3; x-z@y=b; 4
         call epqvolatr(cnx, na, vol)
         a = vol(1)
         c = vol(3)
         h = vol(4)

         if(surfN == 1) then
            normal = epDirec(0.d0, 0.d0, sign(1.d0, -h))
         elseif(surfN == 2 ) then
!          \ is  ((a-c)/sqrt((a-c)**2 + h**2), 0,-h/sqrt(//))
!         normal to \:  (h/sqrt(//), 0, (a-c)/sqrt(//))
            temp = sqrt( (a-c)**2 + h**2)
            normal =
     *      epDirec(abs(h)/temp, 0.d0, (a-c)/temp*sign(1.d0,h))
         elseif(surfN == 3) then
            normal = dir(3)
         elseif( surfN == 4 ) then
            normal = dir(4)
         elseif( surfN == 5 ) then
!          / is  (c/sqrt(c**2 + h**2), 0, h/sqrt(//))
!         normal to /:  (-h/sqrt(//), 0, c/sqrt(//))
            temp = sqrt( c**2 + h**2 )
            normal =
     *      epDirec(-abs(h)/temp, 0.d0, c/temp*sign(1.d0,h))
         else
            write(0,*) ' surface #=', surfN, ' strage for prism'
            stop
         endif
         
      else
         write(0,*) 'structure ', Det%cmp(cnx)%struc(1:4), 
     *    ' is not supported for NormalVec'
         stop
      endif
      end

      
      subroutine epLightMirrorRef(normal, d1, d2)
      implicit none
#include "ZepDirec.h"      
       type(epDirec)::   normal  
         !   input. normal vector of the current surface (outward going)
       type(epDirec)::   d1
         !   input.  light direction cos in the local coord. of the current
         !           compo.
       type(epDirec)::   d2
         !   output.  light direction cos in the local coord. of the current
         !             compo. after mirror reflection.

      real(8)::a
  !          -d2 +  + d1 =a n  (n = normal ;  unit vector)
  !   so    (d1-d2)^2 = a^2 and    d2^2 =d1^2
  !   then   2d1^2 -2d1.d2 =a^2 ==> 2d1^2 -2d1.(d1-an)=a^2
  !          hence 2ad1.n=a^2 
  !   finally  a = 2d1.n  and  d2 = d1 - 2d1.n n
      call cscalerProd( d1, normal, a)
      d2%x = d1%x - 2*a* normal%x
      d2%y = d1%y - 2*a* normal%y
      d2%z = d1%z - 2*a* normal%z
      end

      subroutine epLightDiffRef(normal, d1, d2)
      implicit none
#include "ZepDirec.h"      
       type(epDirec)::   normal  
         !   input. normal vector of the current surface (outward going)
       type(epDirec)::   d1
         !   input.  light direction cos in the local coord. of the current
         !           compo.
         !  at present  NOT USED since isotropic
       type(epDirec)::   d2
         !   output.  light direction cos in the local coord. of the current
         !             compo. after diffusive reflection.


      real(8)::cost, sint, cs, sn

      call rndc(cost)
      sint = sqrt(1.d0 - cost**2)
      call kcossn(cs, sn)
      d2%x = sint*cs
      d2%y = sint*sn
      d2%z = -cost   ! since reflection, normal vector is -normal
         !               this can be equivalent to making -cost
         !  rotate this so that z axis is normal
      call eptransVect( normal, d2, d2)
      end


      subroutine epLightRefrac(n1, n2, Normal,  Dir1, Dir2, jcon)
!            can be tested by using Test/testLightRefrac.f 
!   in that case,   activate common below 
      implicit none
#include "ZepDirec.h"
      real(8),intent(in):: n1, n2 !  refraction indices of  media. n1-->n2 transmission
       type(epDirec)::  Normal ! input. unit normal vector at point on
                        ! the boundary  of two media.  direction is outgoing
       type(epDirec)::  Dir1   ! input.  direction cos of light inside comp.
       type(epDirec)::  Dir2   ! output. difracted light's dir. 
              !       **** this is the value in the local cooord. of current comp. *****
      integer,intent(out)::jcon ! =0 --> no error, =1 refraction cannot 
                    ! take place, probably total reflection occurs.
                    ! =2, Normal vector and Dir1 is physically strange
                    ! (angle between them >= 90 deg  or not normalized dir).

   !  -----------------------

       type(epDirec)::  Xnormal ! x axis of the system with z azis = normal
                               ! x perpendicular to  the plane formed by 
                               ! normal and dir1; defined in current comp.
!      common /fortest/ Xnormal  ! ****** activate this for test

       type(epDirec)::  d1, d2

!           sin1 n1 = sin2 n2
     
      real(8):: sin1, sin2,  cos2, cos1
      real(8):: cosf, sinf, fai, temp, maxcomp

      call cscalerProd( Dir1, Normal, cos1 )   ! 
      if(cos1 <= 0.  .or. abs(cos1) > 1.) then
         write(0,*) ' error : light direction is wrong '
         write(0,*) ' epLightRefrac: Normal vec= ', Normal
         write(0,*) '                 Light vec= ', Dir1
         write(0,*) ' cos angle between them is ', cos1
         jcon = 2
      else
         temp = 1. - cos1**2                    
         sin1 = sqrt(temp)
         sin2 = sin1*n1/n2      
         if( sin2 .gt. 1.) then
            write(0,*) ' error input to epLightRefrac'
            write(0,*) ' n1, n2, Normal, Dir1 =',
     *           n1, n2, Normal,  Dir1
            write(0,*)
     *        ' no refraction but total reflection should happen'
            jcon = 1
         else
            cos2 = sqrt(1.d0 - sin2**2)           
      !      get Xnormal
      !      perpendicular to the plane formed by Dir1 and Normal
            call cvecProd(Dir1, Normal,  Xnormal) 
      !      we can say Dir2 is defined in the xyz system.
      !      if Dir1 and Normal is parallel, we may choose 
      !      Xnormal = current x axis 
            call cscalerProd(Xnormal, Xnormal, temp) ! 1-Dir1.z^2 = sin1^2
            if( abs(temp) < 1.d-8) then
               call epgetNormal2vec(Normal, Xnormal)
            else         
               temp = sqrt(temp)
               Xnormal%x = Xnormal%x/temp
               Xnormal%y = Xnormal%y/temp
               Xnormal%z = Xnormal%z/temp
            endif
              ! system R: x is Xnormal, z is Normal
              ! vectors so far  defined in B system (local coord).
              ! get fai of Dir1 in system R, for this, convert
              ! Dir1 into R.
            call epitransVectZx(1, Normal, Xnormal, Dir1, d1) ! d1
              ! fai of d1 is the same as that of d2;  
            fai = atan2(d1%y, d1%x) 
            d2%z = cos2
            d2%x = sin2*cos(fai)
            d2%y = sin2*sin(fai)
               ! convert d2 into local coord.
            call eptransVectZx(1, Normal, Xnormal, d2, Dir2)
         endif
      endif
      end

      subroutine epLightFuzzify(sigma, normal, dirp, dirpo)
  !       add random factor to the direction by Gaussian noise
      implicit none
#include "ZepDirec.h"
      real*8 sigma  ! input sigma of the gaussian spread  (deg)
       type(epDirec)::  normal  ! normal vector at the surface.  may be
             ! outword or inword directed. If Sampled angel is 
             ! almost perpedicular to this normal vec. we resample
             !
       type(epDirec)::  dirp  !in.  original direction 
       type(epDirec)::  dirpo  !out.  sampled direction 

      real*8  teta, phi, hwhm, oa
       type(epDirec)::  dirp2
      integer i
      real*8 temp
!       get theta and phi in deg. of the direction cos


      i = 1
      call kdtoa(dirp%x, dirp%y, dirp%z, teta, phi)
      hwhm = sigma*1.177410     ! half width at half max.
      oa =min( hwhm * 3.d0, 90.d0) ! limit the spread by this opening angle.
      do while (.true.)
         call epgonSphere(i, hwhm, 1.d0, teta, phi, oa,
     *     dirp2)
         call cscalerProd(dirp2, normal, temp)
         if(  abs(temp) > 1.d-2) exit
              !  if angle is almost perpendicular to the normal vec.
              !  retry
         i = i + 1
         if( i > 100 ) then
            dirp2 = dirp
            exit
         endif
      enddo
      dirpo = dirp2
      end

      subroutine epLightNormalVecFE(posin, cnx, surfN, normal,icon)
         ! get normal vector at pos which is on the surface
         ! of component specified by component # cnx
         ! This is the front end
      use modepLightMaxDef
#include "Zcondc.h"
      implicit none
#include "ZepTrackv.h"
#include "ZepTrackp.h"
#include "Zcnfig.h"
!      real(8),intent(in):: pos(3)  ! position very close to 
!                   a surface of a compoenent
       type(epPos)::  posin
      integer,intent(in):: cnx ! the component is specified by 
                ! comp. #  cnx
      integer,intent(out):: surfN ! the surface where pos
              ! is has the surface # of this
      real(8),intent(out):: normal(3)  ! obtained normal
           !      vector directed to outword of the comp.
      integer,intent(out):: icon  !  0 pos surfN obtained
                     !     !  surfN cannot be obtained
      
      real(8):: pos(3)
      pos(1) =posin%x
      pos(2) =posin%y
      pos(3) =posin%z
! at present we keep old interface ; later we change
      call epLightGetSurfN(cnx, posin, surfN)
      call epLightNormalVec(cnx, surfN, normal)
      end      subroutine epLightNormalVecFE

      subroutine epLightGetSurfNFE(cnx, posin, surfN)
!     this is front end to use epLightGetSurfN
      implicit none
#include "ZepPos.h"
      integer,intent(in):: cnx ! comp. #
       type(epPos)::  posin   ! input position   

      integer,intent(out):: surfN  ! surface #
      real(8)::pos(3)

      pos(1) = posin%x
      pos(2) = posin%y
      pos(3) = posin%z
!          in many system,  posin could be
!        directly put in the next call,but
!        we use explicit array notation

!      call epLightGetSurfN(cn, pos, surfN)
!           to be replaced by above.
      call epLightGetSurfN(cnx, posin, surfN)
      end       subroutine epLightGetSurfNFE

      
