#include "BlockData/cblkGene.h"
      include "ZfitBD.h"
      program timeFit
      implicit none
      include "Zfit.h"

      real*8 xin(maxbin), yin(maxbin)
      character*128 buf 
      integer nbin0, icon, code
      integer count, status, cond
      nbin0 = 0
      count = NARGS()
      if(count .eq. 3) cond=1
      if(count .eq. 2) cond=0
      if( count .ne. 2 .and. count .ne. 3 ) then
         write(0,*) 'Usage: timeFit..  code [c] < inputfile > outfile'
         write(0,*) '    code = particle code,',
     *   ' inputfile is output  from mkRvsT.csh '
         write(0,*) ' c:  if exists, write only coeff. on stdout '
         write(0,*) '     if omitted, coeff. is output  on stderr,',
     *              ' fitted (x,y) on stdout'
         stop
      endif
      call getarg(1, buf, status)
      read(buf, *) code
      if(code .le. 3) then
         code = 1
      else
c         code = 2
c          with negative c hadron can be treated same as e,g,mu
c           
         code = 1
      endif

      do while(.true.)
         read(*,*,end=100) xin(nbin0+1), yin(nbin0+1)
         nbin0 = nbin0 +1
      enddo
100   continue
      if(nbin0 .gt. maxbin) then
         write(0,*) ' too many time data> ', maxbin
         stop
      endif

c      call copenfw2(minout, "./minout_@_#_%",   1,  icon)
      call copenfw2(minout, "/dev/null",   1,  icon)
cccc      call copenfw2(minsave, "./minsave_@_#_%",  1,  icon)  !!!


      call fittime0(cond, code, xin, yin, nbin0)
      end

      subroutine fittime0(cond, code, xin, yin, nbin)
      implicit none
      include "Zfit.h"
      integer cond
      integer nbin, code
      real*8 xin(nbin), yin(nbin)
      real*8 xuse(nbin), yuse(nbin)
      real*8 prmout(nparam, nregion)

      integer i 
      integer n1, n2, region
      real*8 xx, f, xb

c          fitting at region 
      do region=1, nregion
c              for hadron fitting is done only at region2
         if(code .eq. 2 .and. region .eq. 1) cycle
         n1  =0  
         do i = 1, nbin
            if(xin(i) .gt. x2(region) ) exit
            if(xin(i) .lt. x1(region) ) cycle
            if(code .eq. 2 .and. xin(i) .lt. x1h ) cycle
            n1 = n1 +1
            xuse(n1) = xin(i)
            yuse(n1) = yin(i)
         enddo
c////////////
c         write(0,*) ' region=',region, ' points=', n1
c         write(0,*) ' param(1,region)=', param(1, region)
c         write(0,*) ' param(2,region)=', param(2, region)
c         write(0,*) ' param(3,region)=', param(3, region)
c//////////
c                       region,  x,   y,   #, output param
         call fittime1(region, code,  xuse, yuse, n1, param(1, region),
     *          prmout(1, region))

         if(cond .eq. 1) then
c                only coeff. is put on stdout
            write(*,'(3g12.4)')  prmout(1, region), prmout(2, region),
     *      prmout(3, region)
            if(code .eq. 2) then
c                for hadron, region 1 is mssing so we repeat region 2 data
               write(*,'(3g12.4)')
     *          prmout(1, region), prmout(2, region), prmout(3, region)
            endif
         else
c             coeff is put on stderr
            write(0,'(3g12.4)')  prmout(1, region), prmout(2, region),
     *      prmout(3, region)
         endif
         if(cond .eq. 0) then
c               to see fitted result (r, t) is put on stdout  
            xx = drx1( region)
            do while ( xx .le. drx2(region) )
               f=prmout(1,region)*xx**(prmout(2,region) +
     *              prmout(3,region)*log(xx))
               write(*,*) xx, f
               xx = xx*10.0**0.02
            enddo
         endif
      enddo
      end
c     ***************************************************
      subroutine fittime1(region,  code, xin, yin, n, prmin, prmout )
      implicit none
      include "Zfit.h"
      integer region, code
      real*8  prmin(nparam), prmout(nparam)
      integer  n
      real*8 xin(n), yin(n)

      integer nlabel(nparam)
      character*10  pname(nparam)
      real*8 initval(nparam)
      real*8 step(nparam)

      data nlabel/ 1,  2, 3/
      data pname/ 'p',  'q', 'r'/
      data step/  1.,  0.001d0, 0.0001d0/
      real*8 zero, one, three, four, five
      data zero,one,three,four, five / 0., 1., 3.,4., 5. /
      real*8 fval, xx
      integer i, ierflg

      external timefnc
      
c
c           in fortran mode, this must be called for a new fnc
c
      npoint = n
      do i = 1, npoint
         x(i) = xin(i)
         y(i) = yin(i)
      enddo

      do i = 1, nparam
         initval(i) = prmin(i)
      enddo

      call mninit( 5, minout, minsave)

      do  i= 1, nparam
c        nprm: a number given to a parameter: (label)
c        pnam: name of the parameer
c        vstrt: initial value of the parameter
c        stp:   initial step size of the //
c        next two: zero-->the parameter is not bounded (lower or upper)
c        ierflg: retrun value; cond code. 0--> ok

         call mnparm(nlabel(i), pname(i), initval(i), step(i),
     *     low(i, code), up(i, code), ierflg)

         if (ierflg .ne. 0)  then
            write (0,'(a,i3)')  ' unable to define parameter no.',i
            stop
         endif
      enddo
c
      call mnseti('tf as a function of core distance')
c       request fcn to read in (or generate random) data (iflag=1)
c            fcnk0: function to be minimuzed is calculated. also 
c              there are other funcitons
c            one is the  argument to fcnk0.  seems to be converted to
c            integer inside.
c            1 number of argument in one  (one could be array)
c           ierflf: ouptut. 0-->ok
c            0: no external function is used in fcnk0

      call mnexcm(timefnc, 'call fcn', one ,1,ierflg, 0)
c        fix the  3,4,5-th parameters,  
c      call mnexcm(timefnc,'fix', fixlist ,3, ierflg,0)
c       print minumum things   
      call mnexcm(timefnc,'set print', zero ,1,ierflg,0)
c                use migrad method for minimization
c                with default condtions
      call mnexcm(timefnc,'migrad', zero ,0,ierflg,0)
c                analysis of errors for all parameters
      call mnexcm(timefnc,'minos', zero ,0,ierflg,0)
c                5th parameter is now set to a variable parameter
c         call mnexcm(timefnc,'release', five ,1,ierflg,0)
c            release 3,4
c         call mnexcm(timefnc,'release', fixlist ,2,ierflg,0)
c
c                  and use migrad again
c         call mnexcm(timefnc,'migrad', zero ,0,ierflg,0)
c                  error analysis
c         call mnexcm(timefnc,'minos',  zero ,0,ierflg,0)
c
c
c      call mnexcm(timefnc,'migrad', zero ,0,ierflg,0)
c                  error analysis
c      call mnexcm(timefnc,'minos',  zero ,0,ierflg,0)
c
c                 call fcn with 3. i.e, ouput etc.  
      call mnexcm(timefnc,'call fcn', three , 1,ierflg, 0)

      do i = 1, nparam
         prmout(i) = oparam(i)
      enddo

      end
