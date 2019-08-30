#include "BlockData/cblkGene.h"
      include "ZlatfitBD.h"
      program latFit
      implicit none
      include "Zlatfit.h"

      real*8 xin(maxbin), yin(maxbin)
      character*128 buf 
      integer nbin0, icon, code, tcode
      integer count, status, cond
      nbin0 = 0
      count = NARGS()
      if(count .eq. 3) cond=1
      if(count .eq. 2) cond=0
      if( count .ne. 2 .and. count .ne. 3 ) then
         write(0,*) 'Usage: latFit..  code [c] < inputfile > outfile'
         write(0,*) '    code = particle code,',
     *   ' inputfile is output  from splitLat%csh '
         write(0,*) ' c:  if exists, write only coeff. on stdout '
         write(0,*) '     if omitted, coeff. is output  on stderr,',
     *              ' fitted (x,y) on stdout'
         stop 123
      endif
      call getarg(1, buf, status)
      read(buf, *) code
      tcode=code
!      if(code .le. 3) then
!         code = 1
!      else
!         code = 2
!          with negative c hadron can be treated same as e,g,mu
!           
!         code = 1
!      endif

      do while(.true.)
         read(*,*,end=100) xin(nbin0+1), yin(nbin0+1)
         nbin0 = nbin0 +1
      enddo
100   continue
      if(nbin0 .gt. maxbin) then
         write(0,*) ' too many lat data> ', maxbin
         stop
      endif

!      call copenfw2(minout, "./minout_@_#_%",   1,  icon)
      call copenfw2(minout, "/dev/null",   1,  icon)
!ccc      call copenfw2(minsave, "./minsave_@_#_%",  1,  icon)  !!!


      call fitlat0(cond, tcode,  code, xin, yin, nbin0)
      end

      subroutine fitlat0(cond, tcode,  code, xin, yin, nbin)
      implicit none
      include "Zlatfit.h"
      integer cond
      integer nbin, code, tcode
      real*8 xin(nbin), yin(nbin)
      real*8 xuse(nbin), yuse(nbin)
      real*8 prmout(nparam, 4, nregion)

      integer i 
      integer n1, n2, region
      real*8 xx, f, xb

!          fitting at region 
      do region=1, nregion
!         if(tcode .eq. 3 .and. region .eq. 4) then
!            pw = 0.5     ! %pw% <--this %pw%  is needed blanks/ = /
                         ! Used by mkLDD/Util/Lat/ 
!         else
!            pw = 0.5
!         endif

         n1  =0  
         do i = 1, nbin
            if(xin(i) .gt. x2(region) ) exit
            if(xin(i) .lt. x1(region) ) cycle
            n1 = n1 +1
            xuse(n1) = xin(i)
            yuse(n1) = yin(i)
         enddo
!////////////
!         write(0,*) ' region=',region, ' points=', n1
!         write(0,*) ' param(1,region)=', param(1, region)
!         write(0,*) ' param(2,region)=', param(2, region)
!         write(0,*) ' param(3,region)=', param(3, region)
!//////////
!                       region,  x,   y,   #, output param
         call fitlat1(region, code,  xuse, yuse, n1,
     *     param(1, code, region), prmout(1, code,region))

         if(cond .eq. 1) then
!                only coeff. is put on stdout
            write(*,'(5g12.4)')  prmout(1, code,region),
     *      prmout(2,code,  region),
     *      prmout(3, code, region),  prmout(4,code, region),
     *      maxdiff
!            if(code .eq. 2) then
!c                for hadron, region 1 is mssing so we repeat region 2 data
!               write(*,'(3g12.4)')
!     *          prmout(1, region), prmout(2, region), prmout(3, region)
!            endif
         else
!             coeff is put on stderr
            write(0,'(5g12.4)')  prmout(1,code, region),
     *       prmout(2, code, region),
     *      prmout(3, code,  region), prmout(4, code,region),
     *      maxdiff
         endif
         if(cond .eq. 0) then
!               to see fitted result (r, t) is put on stdout  
            xx = drx1( region)
            pw = prmout(4,code, region)
            do while ( xx .le. drx2(region) )
               f=prmout(1,code,region)/
     *         xx**(prmout(2,code, region) +
!     *              prmout(3,region)* log(xx) )
     *              prmout(3,code, region)* xx**pw )
               write(*,*) xx, f
               xx = xx*10.0**0.02
            enddo
         endif
      enddo
      end
!     ***************************************************
      subroutine fitlat1(region,  code, xin, yin, n, prmin, prmout )
      implicit none
      include "Zlatfit.h"
      integer region, code
      real*8  prmin(nparam), prmout(nparam)
      integer  n
      real*8 xin(n), yin(n)

      integer nlabel(nparam)
      character*10  pname(nparam)
      real*8 initval(nparam)
      real*8 step(nparam)

      data nlabel/ 1,  2, 3, 4/
      data pname/ 'p',  'q', 'r', 'pw'/
      data step/  1.,  0.001d0, 0.0001d0, 0.1d0/
      real*8 zero, one, three, four, five
      data zero,one,three,four, five / 0., 1., 3.,4., 5. /
      real*8 fval, xx
      integer i, ierflg

      external latfnc
      
!
!           in fortran mode, this must be called for a new fnc
!
      npoint = n
      do i = 1, npoint
         x(i) = xin(i)
         y(i) = yin(i)
!*************
         badindex(i)=i
!***************
      enddo

      do i = 1, nparam
         initval(i) = prmin(i)
      enddo

      call mninit( 5, minout, minsave)

      do  i= 1, nparam
!        nprm: a number given to a parameter: (label)
!        pnam: name of the parameer
!        vstrt: initial value of the parameter
!        stp:   initial step size of the //
!        next two: zero-->the parameter is not bounded (lower or upper)
!        ierflg: retrun value; cond code. 0--> ok

         call mnparm(nlabel(i), pname(i), initval(i), step(i),
     *     low(i, code,region), up(i, code,region), ierflg)

         if (ierflg .ne. 0)  then
            write (0,'(a,i3)')  ' unable to define parameter no.',i
            stop
         endif
      enddo
!
      call mnseti('lat as a function of core distance')
!       request fcn to read in (or generate random) data (iflag=1)
!            fcnk0: function to be minimuzed is calculated. also 
!              there are other funcitons
!            one is the  argument to fcnk0.  seems to be converted to
!            integer inside.
!            1 number of argument in one  (one could be array)
!           ierflf: ouptut. 0-->ok
!            0: no external function is used in fcnk0
      limit = 0.
      call mnexcm(latfnc, 'call fcn', one ,1,ierflg, 0)
!        fix the  3,4,5-th parameters,  
!      call mnexcm(timefnc,'fix', fixlist ,3, ierflg,0)
!       print minumum things   
      call mnexcm(latfnc,'set print', zero ,1,ierflg,0)
!                use migrad method for minimization
!                with default condtions
      call mnexcm(latfnc,'migrad', zero ,0,ierflg,0)
!                analysis of errors for all parameters
      call mnexcm(latfnc,'minos', zero ,0,ierflg,0)


      if(region .eq. 4 ) then
!             if max diff is < 10% no more trial
!                log(1.1)**2 = 0.009
         if(maxdiff .gt. 0.01) then
            badindex(maxindex)= -maxindex
            call mnexcm(latfnc,'migrad', zero ,0,ierflg,0)
!                analysis of errors for all parameters
            call mnexcm(latfnc,'minos', zero ,0,ierflg,0)
!             if there is still 20 % diff.. remove it
            if(maxdiff .gt. 0.033) then
!               write(0,*) ' maxdiff=',maxdiff,
!     *          ' idx=',maxindex
               badindex(maxindex)= -maxindex
               call mnexcm(latfnc,'migrad', zero ,0,ierflg,0)
!                analysis of errors for all parameters
               call mnexcm(latfnc,'minos', zero ,0,ierflg,0)
            endif
         endif
      endif

!
!                 call fcn with 3. i.e, ouput etc.  
      call mnexcm(latfnc,'call fcn', three , 1,ierflg, 0)

      do i = 1, nparam
         prmout(i) = oparam(i)
      enddo

      end
