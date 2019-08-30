!>@file Cosmos7.655/Tracking/Geomag/cgeomag.f: These code is for calculating geomagnetic field. 
!>@author K.Kasahara
!>@date 1980/12/30 
!c          test crdGeomag
!      real*8 year
!      read(*,*) year
!      call crdGeomag('../../Data/Geomag/wmm', year)
!      end
!       **************************************************************
!       *
!       *  cgeomag: geomagnetic filed strength is obtained
!       *
!       **************************************************************
!  /usage/ call cgeomag(year, llh,  h, icon)
!   year: real*8. input.  such as 1990.5
!   llh:  /coord/ input.  position around the earth. 
!                         in 'llh' form is better. if not 'llh'
!                         conversion is done here.
!     h:  /magfield/. output.  magnetic field is set in
!                         the form of 'ned' (north, 
!                         east-down). The unit is T.
!  icon:  output. integer*4  0---> o.k
!                            1---> input parameter wrong.
!


      subroutine cgeomag(yearin, llh, h,  icon)
      use modAtmosDef
      implicit none
#include  "Zglobalc.h"
#include  "Zcoord.h"
#include  "Zmagfield.h"
!       #include  "Zearth.h"
       real*8 yearin
       type(coord)::llh
       type(magfield)::h
       integer icon

!     ************************************
      integer nmax
      parameter (nmax = 13)
      real*4 gnm(nmax, 0:nmax), hnm(nmax, nmax), year
      real*4 gnmd(nmax, 0:nmax), hnmd(nmax, nmax)
      integer nmx
      common /ZmagCoef/  gnm, gnmd, hnm, hnmd, year, nmx
!     ************************************      

!
!
       real*4 r, sumn, sume, sumd, t, cost, sint, x, tlonr, gmnc, 
     *       cosml, sinml, hmnc, temp
       real*8 gn, ge, gd
       real*4 ssumd, ssumn, ssume 
       integer m, n
       type(coord)::cdata
       real*4 kdpmnxsinn, kpmnxn, kpmnxisinn

!
!         check data type
       if(llh%sys .eq. 'llh') then
          cdata = llh
       else     ! convet to llh
          call ctransCoord2('llh', llh, cdata)
       endif   
!
       icon =  0
       if( abs(yearin-year) .gt. 7.) then
           icon=1
           call cerrorMsg(
     *     ' Year spec. for geomag data is invalid.', 1)
!       elseif(abs(cdata.lat) .gt. 90.) then
       elseif(abs(cdata%r(1)) .gt. 90.5d0) then
           icon=2
!       elseif(abs(cdata.long) .gt. 360.) then
       elseif(abs(cdata%r(2)) .gt. 360.5d0) then
           icon=2
!       elseif(cdata.h .gt. 5000.d3) then
       elseif(cdata%r(3) .gt. 100000.d3) then
           icon=1
!       elseif(cdata.h .lt. -3000.d3) then
       elseif(cdata%r(3) .lt. -3000.d3) then
           icon=1
       endif
       if(icon .ne. 2) then
!           r=1./( 1.+cdata.h/Eradius )
           r=1./( 1.+cdata%r(3)/Eradius )
           sumn=0.
           sume=0.
           sumd=0.
!           t=(90.-cdata.lat)*Torad
           t=(90.-cdata%r(1))*Torad
           cost=cos(t)
           sint=sin(t)
           x=cost
!           tlonr=cdata.long*Torad
           tlonr=cdata%r(2)*Torad
           do   n=1, nmx
              ssumn = 0.
              ssume = 0.
              ssumd = 0.
              do   m=0, n
                 gmnc= gnm(n, m)
                 cosml=cos(m*tlonr)
                 sinml=sin(m*tlonr)
                 if(m .gt. 0) then
                    hmnc=hnm(n,m) 
                    ssumn = ssumn+  (gmnc*cosml+hmnc*sinml)*
     *                kdpmnxsinn(m, n, sint, x)
                    ssume = ssume+  m *(-gmnc*sinml + hmnc*cosml)
     *                  *  kpmnxisinn(m, n, sint, x)
                    ssumd = ssumd + (gmnc*cosml+ hmnc*sinml)
     *                 * kpmnxn(m, n, sint, x) 
                 else
                    ssumn = ssumn + gmnc*kdpmnxsinn(m, n, sint, x)
                    ssumd = ssumd + gmnc*kpmnxn(m, n, sint, x)
                 endif
              enddo
              temp = r**(n+2)
              sumn = sumn + ssumn*temp
              sume = sume + ssume*temp
              sumd = sumd + ssumd*temp *(n + 1)
          enddo
!              original formula gives  data in nT.
!              north component
          gn = -sumn /1.e9  ! to T.
!              east  component
          ge = -sume /1.e9
!              down
          gd = -sumd /1.e9
          call csetMagField('ned', gn, ge, gd, h)
       endif
       if(icon .eq. 2) then
          call cerrorMsg('Geometrical input data wrong', 0)
       endif
       end
!>@brief Subroutine to read Field Model parameter from a file in $COSMOSTOP/Data/Geomag/ 
!>@param[in] filepath  File contain geomagnetic fileld model parameter
!>@param[in] yearin  calculation year
      subroutine crdGeomag(filepath, yearin)
      implicit none
!    #include "Zmanagerp.h"
      character*(*) filepath
      real*8 yearin
      integer klena
      character*65 msg


!     ************************************
      integer nmax
      parameter (nmax = 13)
      real*4 gnm(nmax, 0:nmax), hnm(nmax, nmax), year
      real*4 gnmd(nmax, 0:nmax), hnmd(nmax, nmax)
      integer nmx
      common /ZmagCoef/  gnm, gnmd, hnm, hnmd, year, nmx
!     ************************************      

      character*128 path
      character*1 NULL
      integer devn, icon, leng, kgetenv
      save devn
      data devn /12/

      NULL = char(0)

      if(filepath .eq. ' ') then
!        assume it is  Cosmos/Data/Geomag/igrf
!          get $COSMOSTOP
!     leng = kgetenv("COSMOSTOP"//NULL, path)
         leng = kgetenv("LIBLOFT"//NULL, path)         
         if(leng .eq. 0) then
            call cerrorMsg(
     *      'Env. variable "LIBLOFT" cannot be found', 0)
         endif
         path = path(1:leng)//"/Data/Geomag/igrf"
      else
         path = filepath
      endif
      call copenf(devn, path(1:klena(path)), icon)
      if(icon .eq. 1) then
         call cerrorMsg(path, 1)
         call cerrorMsg(
     *        'above file cannot be opened for geomagnetic data', 0)
      endif

      call cgmgIgrf(devn, path, yearin)
      end
!     *************************
!>@brief Read igrf data from COSMOSTOP/Data/Geomag/ 
!> Read coefficient of IGRF datafile called at subroutine cgeomag
!>@param[in] devn (Opened file's device number?)
!>@param[in] filepath filename which is scanned.
!>@param[in] yearin input year. such as 1990.5
      subroutine cgmgIgrf(devn, filepath, yearin)
      implicit none
!            assume at least two year data is included.      
      integer devn              ! input file dev. no.
      character*(*) filepath
      real*8  yearin ! input. year of geomag.

      integer i
      character*1024  msg


!     ************************************
      integer nmax
      parameter (nmax = 13)
      real*4 gnm(nmax, 0:nmax), hnm(nmax, nmax), year
      real*4 yearmin, yearmax
      real*4 gnmd(nmax, 0:nmax), hnmd(nmax, nmax)
      integer nmx
      common /ZmagCoef/  gnm, gnmd, hnm, hnmd, year, nmx
!     ************************************      
!           format from vesion Cosmos 7.55 is diff. from
!       old one;  the data format for old one cannot be
!       found in the web site of igrf. The format found
!       in the web site is 
!   g/h n m 1900.0 1905.0  ...   2005     2010      SV
!   g  1  0 -31543 -3146   ... -29554.63 -29496.5   11.4
!   g  1  1  -2298  -2298  ... -1669.05  -1585.9    16.7        
!   h  1  1   5922   5909  ...  5077.99   4945.1  -28.8
!   g  2  0   -677   -728  ... -2337.24  -2396.6  -11.3
!   ....
!   g 13 13      0      0  ...  -0.18     -0.3    0.0
!   h 13 13      0      0  ..   -0.82     -0.8    0.0
!
      integer klena, msglen

      integer nf ! # of  fields in the data
      integer n1, n2
      integer:: n11, n12, n21, n22, n, m
      integer:: targetf1, targetf2

      real(4):: year1, year2, sv, dy,  yearx

      rewind  devn

      
!       first find the appropriate year
      year = yearin
      msg = ' '
!     format changed from igrf12:  skip if first col. is # or c/s
      do
         read(devn, '(a)' ) msg
         if( msg(1:1) /= "#" .and. msg(1:3) /= "c/s") exit
      enddo
!     chagne for igrf12 ended.  igrf10 can be read also.
!         now msg should start with "g/h"

      call kcountFields(msg, nf)
      if(nf < 5) then
         write(0,*) ' input igrf data is invalid'
         write(0,*) ' the first line is ', msg
         write(0,*) ' detected in cgmgIgrf '
         stop 9999
      endif
      call kcgetaField(msg, 4, n1, n2)      

      if(n1 <= 0) then
         write(0,*) ' strange in cgmgIgrf'
         stop 9999
      endif

      read( msg(n1:n2), *)  year1
      yearmin = year1
      
      call kcgetaField(msg, nf-1, n1, n2)
      read( msg(n1:n2), *)  yearmax

      targetf2 = 5
      if(year < year1 ) then
         targetf1 = 4
         targetf2 = 5
         dy = year - year1
      else
         do i = 5, nf-1
            call kcgetaField(msg, i, n1, n2) ! 2nd year  to last year
            read(msg(n1:n2), *) year2
            if( year < year2 .and. year >= year1 ) then
               targetf1 = i - 1
               targetf2 = i
               dy = year - year1
               exit
            else
               year1 = year2
               if( i == nf-1 ) then
                  targetf1 = nf -1
                  targetf2 = nf -1
                  dy = year - year1
                  year2 = year1+5
               endif
            endif
         enddo
      endif
      if( abs(dy) > 11. ) then
         write(0,*) 'yearin to cgmgIgrf is out of range'
         write(0,*) 'yearin=', yearin
         stop 8888
      endif
!              read coeff. 
      nmx = 0
      do while (.true.)
         read(devn, '(a)', end=1000) msg
         call kcgetaField(msg, 2, n1, n2)
         read( msg(n1:n2), *)  n 
         call kcgetaField(msg, 3, n1, n2)
         read( msg(n1:n2), *)  m
         if( n > nmax .or.  m > nmax ) then
            write(0,*) ' too large m or n; m=', m, ' n=',n
            write(0,*) ' max should be ', nmax
            write(0,*) ' cgmgIgrf '
            stop 1111
         endif
         nmx = max(n,m,nmx)
         call kcgetaField(msg, 1, n1, n2)    !  g/h
         call kcgetaField(msg, targetf1, n11, n12)  ! coef. for year1
         call kcgetaField(msg, targetf2, n21, n22)  ! //            2  

            
         if( msg(n1:n2) == 'g' ) then
            read( msg(n11:n12), *) gnm(n,m) 
            read( msg(n21:n22), *) gnmd(n,m)
            if( year2 <= yearmax) then
               gnmd(n,m) = ( gnmd(n,m) - gnm(n,m) )/ (year2-year1)
            else
               call kcgetaField(msg, nf, n1, n2)
               read(msg(n1:n2), *) sv
               gnmd(n,m) = sv
            endif
         elseif(msg(n1:n2) == 'h' ) then
            read( msg(n11:n12), *)  hnm(n,m) 
            read( msg(n21:n22), *) hnmd(n,m)
            if( year2 <= yearmax ) then
               hnmd(n,m) = ( hnmd(n,m) - hnm(n,m) )/ (year2-year1)
            else
               call kcgetaField(msg, nf, n1, n2)
               read(msg(n1:n2), *) sv
               hnmd(n,m) = sv
            endif
         else
            write(0,*) ' igff data strange the line is' 
            write(0,*) msg
            stop 12345
         endif
      enddo
 1000 continue
      do n = 1, nmx
         do m=0,  n
            gnm(n, m) = gnm(n, m) + dy*gnmd(n, m)
            if(m .gt. 0) then
               hnm(n, m) = hnm(n,m) + dy*hnmd(n, m)
            endif
         enddo
      enddo


      write(0, *) 'Geomagnetic data has been read from', filepath
      write(0, *) ' year=',yearin, ' # of expansion terms=', nmx
      end


!     *************************
!>@brief Print out Geomag model coeff data
      subroutine cprGeomag
      implicit none
!        print current geomag coeff.


      integer m, n


!     ************************************
      integer nmax
      parameter (nmax = 13)
      real*4 gnm(nmax, 0:nmax), hnm(nmax, nmax), year
      real*4 gnmd(nmax, 0:nmax), hnmd(nmax, nmax)
      integer nmx
      common /ZmagCoef/  gnm, gnmd, hnm, hnmd, year, nmx
!     ************************************      

      write(*,*) 'year=',year, ' nmax=', nmx
      write(*,'("  m  n    gnm      hnm    dgnm    dhnm")')
      do n = 1, nmx
         do m =0, n
            if(m .eq. 0) then
               write(*,'(2i3,2f7.0,2f7.1)') m, n,
     *          gnm(n, m), 0., gnmd(n, m), 0.
            else
               write(*,'(2i3,2f7.0,2f7.1)') m, n,
     *                gnm(n, m), hnm(n,m), gnmd(n,m), hnmd(n,m)
            endif
         enddo
      enddo
      end


