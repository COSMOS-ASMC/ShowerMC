!     prepare envirionmental variables in setupenv
!     at execution, give (x,y,z) to be converted from stdin 
!     converted data is put in stdout.
!
      program convcoord
      implicit none

#include  "Ztrack.h"
#include  "Ztrackp.h"
#include  "Ztrackv.h"
#include  "Zobs.h"
#include  "Zobsv.h"

      integer kgetenv2
      integer leng, icon
      character*120 input
      character*20  fromto
      real*8  x, y, z
      type(coord)::xyzin
      type(coord)::xyzout
!
      leng = kgetenv2("PARAM", input)
      write(0,*) ' leng =', leng
      call copenfw2(11, input(1:leng), 1, icon)
      if(icon .ne. 1) then
         write(0,*) input, icon
         write(0,*) ' cannnot be opened'
         stop 1234
      endif
      write(0,*) ' reading parameter file'
      call creadParam(11)
      close(11)
      write(0,*) ' parameter has been read'
      leng = kgetenv2("FROMTO", input)
      write(0,*) ' form to  is ', input
      fromto=input(1:leng)

      call cinitObs
      call cprintObs(0)
      do while(.true.)
         read(*,*, end=1000 ) x, y, z
         if(fromto .eq. 'xyz2det') then
            call csetCoord('xyz', x,y,z, xyzin)
!              to detector system
            call cxyz2det(ObsSites(NoOfSites).pos.xyz, 
     *      xyzin, xyzout)
         elseif( fromto .eq. 'det2xyz') then
            call csetCoord('det', x, y, z, xyzin)
            call cdet2xyz(ObsSites(NoOfSites).pos.xyz, 
     *      xyzin, xyzout)
         elseif( fromto .eq. 'xyz2prim') then
            call csetCoord('xyz', x, y, z, xyzin)
            call cxyz2prim(ObsSites(NoOfSites).pos.xyz,
     *      xyzin, xyzout)
!         elseif( fromto .eq. 'prim2zyz') then
!            call csetCoord('prim', x, y, z, xyzin)
!            call cprim2xyz(ObsSites(NoOfSites).pos.xyz,
!     *      xyzin, xyzout)
         else
            write(0,*) ' fromto=', fromto, ' not supported'
            stop 2222
         endif
         write(*,'(3g12.4)')  xyzout.x, xyzout.y, xyzout.z 
      enddo
 1000 continue
      end
!c#include "BlockData/cblkElemag.h"
