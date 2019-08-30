!             gaussian density beam
       subroutine epGaussb(hwhm, rcut, x, y)
       implicit none
       real*8 hwhm  ! input. half width at half max (cm)
       real*8 rcut  ! input. discard if radius is > rcut (cm)
       real*8 x, y  ! ouput.  sampled  position (cm)

       real*8 alg2, a, tmp, ome, u, r, cs, sn, omu
       
!               ln(2)
           data alg2/.69314718/
!             density is: f(r)=* exp(-(r/rh)**2 * ln(2))
!             distribution is f(r)*r*dr where rh is half
!             width at half
!             max of f(r) = hwhm
           a=alg2/hwhm**2
           tmp=a* rcut**2
           if(tmp .lt. 1.e-4) then
               ome= tmp*(1. - tmp/2)
           else
               ome= 1. -exp(-tmp)
           endif
           call rndc(u)
           tmp= u* ome
           if(tmp .lt. 1.e-4) then
                omu=tmp*(tmp - 1.)
           else
                omu=log(1. - tmp)
           endif
           r=sqrt(-omu/a)
!             sample (x,y)
           call kcossn(cs,sn)
           x=r*cs
           y=r*sn
       end
!     *********************************************************
      subroutine epgonSphere(ini, hwhm, rin, teta, phi, oa, pos)
      implicit none
#include "ZepTrackp.h"
#include "ZepPos.h"     
#include "Zptcl.h"
#include "Zglobalc.h"
!            generate a random point distributed on the
!            surface of a sphere with a Gaussian density.(Note below)
!            Points are distributed around
!            given polar angles (teta, phi) within a given opening angle
!            (oa).  

! Note:   Gaussian means that the beam has a Gaussian density if it is
!    projected to the plane which is tangent to the sphere at the
!    Gaussian center.
!    The rcut used for epGaussb is related to oa here as
!    rcut = r sin( oa ).
! ****** If hwhm is negative, uniform distribution is assumed instead of Gauss.
!
      integer ini               ! input
                      !  1-->  teta and phi are different from
                      !        previous call or this is the first call.
                      !  != 1 -->  teta, and phi are the same as
                      !        the previous call.
      real*8  hwhm      ! input.  half width at half maximum of the
                        !      Gaussian density beam.
              ! if negative, uniform beam is taken
      real*8  rin       ! input.  radius of the sphere
      real*8  teta      ! input.  polar angle in degree
      real*8  phi       ! input.  azimutal angle in degree
      real*8  oa        ! input.  opnening angle in degree. has meaning
!         if os <= 90.
      real(8),intent(out):: origpos(3)
       type(epPos)::  pos  ! output. an  obtained random point

       type(fmom),save::  xyz, xyz2
      real*8  a(4, 4), b(4, 4), ba(4, 4)

      real*8  rcut, x, y, z, z2, r, temp, u, fsin, fcos
      save ba

      r = rin 
      if(ini .eq. 1) then
         call cgetRotMat4(2, -teta*Torad, a)
         call cgetRotMat4(3, -phi*Torad, b)
         call cmultRotMat4(b, a, ba)
      endif
      rcut = r *sin(oa*Torad)
      if(hwhm .ge. 0.0) then
         call epGaussb(hwhm, rcut, x, y)
         z2= r**2 - (x**2 + y**2)
         if(z2 .lt. 0.) then
            call cerrorMsg(
     *           'opening angle for epgonShpere should be checked', 0)
         endif
         z = sqrt(z2)
      else
         call kcossn(fcos, fsin)
         call rndc(u)
         temp = rcut *  sqrt(u)
         x = temp * fcos
         y = temp * fsin 
         z =  sqrt(r**2 - x**2-y**2)
      endif
      xyz%p(1) = x
      xyz%p(2) = y
      xyz%p(3) = z
      xyz%p(4) = 1.
      call capplyRot4(ba, xyz, xyz2)
      pos%x = xyz2%p(1)
      pos%y = xyz2%p(2)
      pos%z = xyz2%p(3)
      return
!!!
      entry epgonSphere2(origpos)
      origpos(1:3) = xyz%p(1:3)
      end
!     **************************************************
      subroutine epuonSphere(ini, rin, teta, phi, oa, pos)
      implicit none
#include "ZepTrackp.h"
#include "ZepPos.h"     
#include "Zptcl.h"
#include "Zglobalc.h"
!          generate a random point uniformly distributed on the
!          surface of a sphere.  Points are distributed around
!          given polar angles (teta, phi) within a given opening angle
!          (oa). (Actually, the point is put little bit inside of
!          the sphere by an amount of EpsLeng so that it can be
!          judged on or inside the sphere safely even with some
!          numerical error). 
!    By uniform  is meant that the points are uniformly distributed on
!    the surface of the sphere but not on a projected plane.
!
      integer ini               ! input
                      !  1-->  teta and phi are different from
                      !        previous call or this is the first call.
                      !  != 1 -->  teta, and phi are the same as
                      !        the previous call.
      real*8  rin               ! input.  radius of the sphere
      real*8  teta              ! input.  polar angle in degree
      real*8  phi               ! input.  azimutal angle in degree
      real*8  oa                ! input.  opnening angle in degree
       type(epPos)::  pos        ! output. an  obtained random point

       type(fmom)::  xyz, xyz2
      real*8  a(4, 4), b(4, 4), ba(4, 4)
      real*8  u, r
      real*8 fcos,  fsin
      save ba

      r = rin 
      if(ini .eq. 1) then
         call cgetRotMat4(2, -teta*Torad, a)
         call cgetRotMat4(3, -phi*Torad, b)
         call cmultRotMat4(b, a, ba)
      endif

      call rndc(u)
      fcos = cos(oa*Torad)
      fcos = (1.d0- fcos) * u +  fcos
      fsin = sqrt(1.d0- fcos**2)
      call rndc(u)
      u = u*pi*2
      xyz%p(1) = r * (fsin * cos(u))
      xyz%p(2) = r * (fsin * sin(u))
      xyz%p(3) = r * fcos 
      xyz%p(4) = 1.
      call capplyRot4(ba, xyz, xyz2)
      pos%x = xyz2%p(1)
      pos%y = xyz2%p(2)
      pos%z = xyz2%p(3)
      end
      subroutine epuSphere(rin, dir,  pos)
      implicit none
#include "ZepTrackp.h"
#include "ZepDirec.h"
#include "ZepPos.h"     
!#include "Zptcl.h"
!#include "Zglobalc.h"
!          Note  diff. from epuonSphere
!          This one will
!          generates a random point on the  surface of a sphere.
!          Points are distributed on the projected plane (area is
!          pixR^2 where R is the radius of the sphere)
!          which is defined as a pependicular plane to the given
!          direction given by dir.  The position returned to the user
!          is little bit inside the sphere by an amount of EpsLeng.
!

      real*8  rin               ! input.  radius of the sphere
!                                 sphere center is assumed to be at the origin.
       type(epDirec)::  dir        ! input.  direction cosine
       type(epPos)::  pos        ! output. an  obtained random point
       type(epPos)::  xyz
      real*8 fcos,  fsin, r, u


      r = rin 

      call kcossn(fcos, fsin)
      call rndc(u)
      xyz%x = r *  sqrt(u)
      xyz%y = xyz%x * fsin 
      xyz%x = xyz%x * fcos
      xyz%z = r * sqrt(1.d0 -  u) 
      if(dir%z .gt. 0.) then 
         xyz%z= -xyz%z
      endif
      call ctransVectZ(dir, xyz, pos)
      end
!        ************************* real*8 data
       subroutine arprmr(vvalue, x)
        implicit none
        integer io
        character*(*) vvalue
        real*8 x
        read(vvalue, *)   x, x
        end
       subroutine arprmra(vvalue, x, n)
        implicit none
        integer io
        character*(*) vvalue
        integer n
        real*8 x(n)
        read(vvalue, *)   x, x
        end

!     *******************************
      subroutine epgetdatline(io, dat, j)
      implicit none
      integer io  ! input. read device logical number
      character*120  dat  ! output. data line
!       null lines or the lines which have two or more blank characters
!       at head is neglected
      integer j  ! input/output input. 0--> If  EOF,stop
                 !                     1--> If  Eof, 2 is returned.
                 !                          if not eof, 1 is unchanged   
      character*1 tab

      tab = char(9)

      dat ='  '
      do while ( dat(1:2) .eq. '  ' .or. dat(1:1) .eq. tab) 
         read(io, '(a)', end=80) dat
      enddo
!         check if "/" exists for safety
      if( index(dat, "/") .eq. 0 .and. j .eq. 0 ) then
         call cerrorMsg(dat, 1)
         call cerrorMsg(' has no "/" in the data line', 0)
      endif
      return
 80   continue
      if(j .eq.  0) then
         call cerrorMsg('EOF while reading input data', 1)
         call cerrorMsg(
     *   'Some of new parameters may be missing: see epicsfile'//
     *   ' or sepicsfile in FirstKiss', 0)
      else
         j=2
      endif
         
      end   
!     ************************* complex data
      subroutine arprmm(vvalue, c)
      implicit none
      character*(*) vvalue
      complex(8):: c
      read( vvalue, *)   c, c
      end
!     ************************ integer data
      subroutine arprmi(vvalue, i)
      implicit none
      character*(*) vvalue
      integer i
      read(vvalue, *)   i,i
      end
!        ************************* character data
      subroutine arprmc(vvalue, cha)
      implicit none
      character*(*) vvalue
      character*(*) cha
      read(vvalue, *)  cha, cha
      end
!        ***************************** logical data
      subroutine arprml(vvalue, logi)
      implicit none
      character*(*) vvalue
      logical logi
      read(vvalue, *)  logi, logi
      end     
!        ---------------------------------------------
      subroutine awprmr(io, vname, x)
      implicit none
      integer io
      character*(*) vname
      real*8  x
      
      write(io, *) ' ', vname,' ', x,' /'
      end
      subroutine awprmr2(io, vname, x, n)
      implicit none
      integer io
      integer n  ! arra size of x
      character*(*) vname
      real*8  x(n)
      
      write(io,*) ' ', vname,' ', x,' /'
      end
      subroutine awprmra(io, vname, x, n)
      implicit none
      integer io
      character*(*) vname
      integer n
      real*8  x(n)
      
      write(io, *) ' ', vname,' ', x,' /'
      end
      subroutine awprmm(io, vname, c)
      implicit none
      integer io
      character*(*) vname
      complex(8)::  c
      write(io,  *) ' ', vname,' ', c,' /'
      end
      subroutine awprmi(io, vname, i)
      implicit none
      integer io
      character*(*) vname
      integer i
      
      write(io,  *) ' ', vname,' ', i,' /'
      end
      subroutine awprmc(io, vname, cha)
      implicit none
      integer io
      character*(*) vname
      character*(*) cha
      integer klena
      character*2 qmk/" '"/             ! ' 
      if(klena(cha) .gt. 0) then
         write(io,  *) ' ', vname, qmk, cha(1:klena(cha)),
     *        qmk,' /'
      else
         write(io, *) ' ', vname, qmk, ' ', qmk, ' /'
      endif
      end
      subroutine awprml(io, vname, logi)
      implicit none
      integer io
      character*(*) vname
      logical  logi

      write(io,  *) ' ', vname,' ', logi,' /'
      end
      subroutine anamec(vn, dt)
      implicit none
!      
!            check variable name
!
      character*(*) vn, dt
      character*90 msg
      if(vn .ne.  dt) then
         write(msg,*) ' data name should be ',vn,
     *        ' but it is ',dt
         call cerrorMsg(msg, 0)
      endif
      end
!      ********************
      subroutine epqHookc(i, cv)
      implicit none
#include "ZepManager.h"
      integer i  ! input. i-th user defined char vairable is requested
      character*(*) cv  ! output. requested variable value

      integer klena
      if(i .le.  0 .or. i .gt. epHooks(1))  then
!         call cerrorMsg('out of range request to epqHookc',2)
         cv = ' '
      else
         cv = epHookc(i)(1:klena(epHookc(i)))
      endif
      end

      subroutine epqHooki(i, iv)
      implicit none
#include "ZepManager.h"
      integer i  ! input. i-th user defined integer vairable is requested
      integer iv  ! output. requested variable value


      if(i .le.  0 .or. i .gt. epHooks(2))  then
!         call cerrorMsg('out of range request to epqHooki',0)
         iv = -9999999
      else
         iv = epHooki(i)
      endif
      end


      subroutine epqHookr(i, rv)
      implicit none
#include "ZepManager.h"
      integer i  ! input. i-th user defined real vairable is requested
      real*8 rv  ! output. requested variable value


      if(i .le.  0 .or. i .gt. epHooks(3))  then
!         call cerrorMsg('out of range request to epqHookr',0)
         rv = -1.d-60
      else
         rv = epHookr(i)
      endif
      end
      logical function epgetParmN( io,  vname, vv )
!          get parameter variable name and given value(s)  from io
      implicit none
       integer io
       character*(*)  vname, vv  ! output
       integer linel
       parameter( linel = 100)
       character*(linel)  line
       integer loc, loc2
       vname = " "
       do while(.true.)
          read(io, '(a)', end=100 ) line
          if( line(1:1) .eq. " " .and. line(2:2) .ne. " ") then
             loc = index( line(3:linel), " ")  + 2
             vname = line(2:loc-1)
             loc2 = index( line, "/")
             if(loc2 .eq. 0 ) then
                write(0,* ) ' "/"  is missing in the input data file '
                write(0,*)  ' The line is: ', line
                stop 1234
             endif
             vv = line(loc+1:linel)  !  some data containes '/' such as '../../Media' so put all
                                     ! data.
             goto 50
          endif
       enddo
 50    continue
       epgetParmN = .true.
       return
 100   continue
       epgetParmN =.false.
       end
      subroutine epSeePMfile(file,  icon)
      implicit none
         !       see if file name contain + or - at head of the
         ! name 

      character(*)  file
      integer,intent(out)::icon !  0  -> no +/-
                                !  1  ->    +
                                !  -1 ->    -

      integer::j, n, i

      if( file(1:1) == "+"  ) then
         icon = 1
      elseif(file(1:1) == "-"  ) then
         icon = -1
      else
         n = len( trim(file))
         do i = n-1, 1, -1
            if( file(i:i) == "/" ) then
               if( file(i:i+1) == "/+") then
                  icon = 1
               elseif(file(i:i+1) == "/-") then
                  icon = -1
               else
                  icon = 0
               endif
               exit
            endif
         enddo
      endif
      end
      subroutine epskipcom(icon)
      implicitnone
#include  "ZsepManager.h"
      integer icon
       do while (.true.)
          if(Inp1ry > 0 ) then
             read(Ioprim, '(a)', end=100)  buf
          else
             read(Ioprim, end=100) buf
          endif
          if( buf(1:20) == "#-------------------") exit
       enddo
       icon = 0 
       return
 100   continue
       icon = 1
       end
      subroutine epqversion(cosv, epiv)
      implicit none
#include "ZepManager.h"

      character*8::cosv ! output  cosmos version such 7.58; left justified
      character*8::epiv ! output  epics  version such 9.99://
      character*64::EPICSTOP
      character*128::filen
      integer kgetenv2, icon

      if( kgetenv2("EPICSTOP", EPICSTOP) == 0 ) then
         write(0,*) " Environmental variable "
         write(0,*) "EPICSTOP is not defined in epqversion"
         stop
      endif
      filen = trim(EPICSTOP)//"/Version/version" 
      call copenf(iowk, filen, icon)
      read(iowk, '(a)') epiv
      close(iowk)
      call cqversion(cosv)
      end


