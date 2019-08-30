#include "Zcondc.h"
#include "BlockData/cblkGene.h"

c     Functions to describe the standard atmospher in the vertical direction.
c
c       cthick2den:   atmospheric depth to air density
c       cvh2den:      height to density
c       cvh2denp:     height to d rho/dz
c       cvh2den2p:     height to  d(d rho/dz)/dz
c       cthick2h:     atmospheric depth to height
c       cvh2temp:   h to temperature 
c       cvh2tthick:   height to thickness(depth)
c
      implicit none
c----      include 'Zstdatmos.h'
#include  "Zstdatmos.h"
#include "ZcosmosExt.h"
      real*8 z, den, cvh2den, t, temp

      real*8 cthick2h
      real*8 cvh2thick, denp, den2p
      real*8 cvh2denp, cvh2den2p, cvh2temp
c         den = 10.d0* ((ha-z)/hl)**hmi -
c     *                10.d0* ((ha-hc)/hl)**hmi 
c         write(0, *) " amount at z < hc=", den, " >hc =",
c     *             10.d0*exp( (hb-hc)/hn)
c        temp=  10.d0*exp( (hb-hc)/hn) - 10.d0* ((ha-hc)/hl)**hmi 
c         write(0, *) " const to be used in cvh2thick is", temp
      integer:: icon

      call copenNLf(11, "param",icon)
      if(icon .ne. 0) then
         call cerrorMsg('File param missing',0)
      endif
      call creadParam(11)
      close(11)
      call creadAtmosD
      write(0,*) "# h(m) rho(kg/m3) rho' rho'' thickness(kg/m2) T(K)"
      write(0,*) "   1       2       3    4       5              6" 
      do z= -500.d0, 100.d3, 100.d0
c         t = 10.d0**z
         den = cvh2den(z)
         denp = cvh2denp(z)
         den2p = cvh2den2p(z)
         t = cvh2thick(z)
         temp = cvh2temp(z)
         write(*, '(1p,6g14.4)')  z, den, denp, den2p, t,temp
      enddo
      end

