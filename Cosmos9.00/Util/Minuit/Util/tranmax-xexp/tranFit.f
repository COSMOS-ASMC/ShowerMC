#include "BlockData/cblkGene.h"
      include "ZfitBD.h"
      program tranFit
      implicit none
      include "Zfit.h"

      real*8 xin(maxbin), yin(maxbin)
      character*128 buf 
      integer nbin0, count, icon, status

#if defined (MACOSX)
      integer iargc
      count = iargc() +1
#else
      count = NARGS()
#endif
      if(count .ne. 4) then
         write(0,*) ' command line arg is ', count
         write(0,*) 'Usage: ./tranaFit  op  < inputfile > outfile'
         write(0,*) '  op 1: # of layers to used before max '
         write(0,*) '  op 2: # of layers to used after  max '
         write(0,*) '  third option is stdout contnent: if it is'
         write(0,*) " 1 stdout is number of x,y's"
         write(0,*) " 2 stdout is fitting coeff. a,b,c,d,x0"
         write(0,*) " 3 stdout is xmax, apparent_xmax"
         write(0,*) " 4 stdout is xmax, apparent_xmax, a,b,c,d,x0"
         write(0,*) " 5 stdout is xmax, apparent_xmax, a,b,c,d,x0"
         write(0,*) "             x1 y1"
         write(0,*) "             x2 y2"
         write(0,*) "             ...  "
         write(0,*) ' inputfile: depth vs flux'
         stop
      endif
#if defined (MACOSX)
      call getarg(1, buf)
      read(buf, *) maxleft
      call getarg(2, buf)
      read(buf, *) maxright
      call getarg(3, buf)
      read(buf, *) op
#else
      call getarg(1, buf, status)
      read(buf, *) maxleft
      call getarg(2, buf, status)
      read(buf, *) maxright
      call getarg(3, buf, status)
      read(buf, *) op
#endif
      nbin0 = 0
      do while(.true.)
         read(*,*,end=100) xin(nbin0+1), yin(nbin0+1)
         nbin0 = nbin0 +1
      enddo
100   continue
      if(nbin0 .gt. maxbin) then
         write(0,*) ' too many data> ', maxbin
         stop
      endif


      call copenfw2(minout, "/dev/null",   1,  icon)


      call fitTran0(xin, yin, nbin0)
      end

      subroutine fitTran0(xin, yin, nbin)
      implicit none
      include "Zfit.h"
      integer nbin
      real*8 xin(nbin), yin(nbin)
      real*8 xuse(nbin), yuse(nbin)
      real*8 prmout(nparam)

      integer i 
      integer n1
      real*8 xx, f, xb, xmax
      real*8 a, b, c, d, x0




c          find max pos and use max of 
c          maleft+maxright+1  points
c          around max.
c     
      maxval = yin(1)-1.0e10
      do i = 1, nbin
         if(maxval .lt. yin(i))  then
            maxpos = i
            maxval = yin(maxpos)
         endif
      enddo

      maxdep = xin(maxpos)
      from = max(maxpos-maxleft, 1)
      to = min(maxpos+maxright, nbin)
c   
      write(0,*) ' max depth index and depth=', maxpos, maxdep
      write(0,*) ' layers for fitting are from', from, ' to ', to
c            a*z**b * exp(-c z **d); z = x/x0
c  
      param(1) = yin(maxpos)   ! a
      param(2) = maxdep       !  x0  
      param(3) = 10.          ! b
      param(4) = 10.          ! c
      param(5) = 1.0          ! d 
      n1 = 0
      do i = from, to
         n1 = n1 + 1
         xuse(n1) = xin(i)
         yuse(n1) = yin(i)
      enddo
c 
      call fitTran( xuse, yuse, n1, param, prmout)
c

c               to see fitted result

      a = prmout(1)
      x0 = prmout(2)
      b = prmout(3)
      c = prmout(4)
      d = prmout(5)

c             coeff is put on stderr
      write(0,'(5g12.4)')  prmout
      xmax = x0* (b/c/d)**(1./d)
      write(0,'("max dep=",g12.3)') xmax

      if(op .eq. 1) then
         xx = drx1
         do while ( xx .le. drx2 )
            f= a*xx**b *exp(-c*xx**d)
            write(*,*) xx*x0, f
            xx = xx*10.0**0.002
         enddo
      elseif(op .eq. 2) then
         write(*,'(5G14.4)')  a, b, c, d, x0
      elseif(op .eq. 3) then 
         write(*,*)  xmax, maxdep
      elseif(op .eq. 4) then
         write(*,'(7g14.4)') xmax, maxdep, a,b,c,d,x0
      elseif(op .eq. 5) then
         write(*,'(7g14.4)') xmax, maxdep, a,b,c,d,x0
         xx = drx1
         do while ( xx .le. drx2 )
            f= a*xx**b *exp(-c*xx**d)
            write(*,*) xx*x0, f
            xx = xx*10.0**0.002
         enddo
      endif
      end
c     ***************************************************
      subroutine fitTran(xin, yin, n, prmin, prmout )
      implicit none
      include "Zfit.h"
      real*8  prmin(nparam), prmout(nparam)
      integer  n
      real*8 xin(n), yin(n)

      integer nlabel(nparam)
      character*10  pname(nparam)
      real*8 initval(nparam)
      real*8 step(nparam)

      data nlabel/ 1,  2, 3, 4, 5/
      data pname/ 'a',  'x0', 'b', 'c', 'd'/
      data step/  1.,  0.01d0, 0.1d0, 0.1d0, 0.1d0/
      real*8 zero, one, three, four, five
      data zero,one,three,four, five / 0., 1., 3.,4., 5. /
      real*8 low(nparam), high(nparam)
      real*8 fval, xx
      integer i, ierflg

      external tranfnc
      
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
      low(1) = prmin(1)
      high(1) = prmin(1)*1.d7
      low(2) = prmin(2)*0.9
      high(2) = prmin(2)*1.1
      low(3) = prmin(3)*0.5
      high(3) = prmin(3)*2.
      low(4) = prmin(4)*0.5
      high(4) = prmin(4)*2.
      low(5) = prmin(5)*0.5
      high(5) = prmin(5)*2.

      call mninit( 5, minout, minsave)

      do  i= 1, nparam
c        nprm: a number given to a parameter: (label)
c        pnam: name of the parameer
c        vstrt: initial value of the parameter
c        stp:   initial step size of the //
c        next two: zero-->the parameter is not bounded (lower or upper)
c        ierflg: retrun value; cond code. 0--> ok

         call mnparm(nlabel(i), pname(i), initval(i), step(i),
c     *    low(i, code), up(i, code), ierflg)
     *     low(i), high(i), ierflg)

         if (ierflg .ne. 0)  then
            write (0,'(a,i3)')  ' unable to define parameter no.',i
            stop
         endif
      enddo
c
      call mnseti('tranfit')
c       request fcn to read in (or generate random) data (iflag=1)
c            fcnk0: function to be minimuzed is calculated. also 
c              there are other funcitons
c            one is the  argument to fcnk0.  seems to be converted to
c            integer inside.
c            1 number of argument in one  (one could be array)
c           ierflf: ouptut. 0-->ok
c            0: no external function is used in fcnk0

      call mnexcm(tranfnc, 'call fcn', one ,1,ierflg, 0)
c        fix the  3,4,5-th parameters,  
c      call mnexcm(timefnc,'fix', fixlist ,3, ierflg,0)
c       print minumum things   
      call mnexcm(tranfnc,'set print', zero ,1,ierflg,0)
c                use migrad method for minimization
c                with default condtions
      call mnexcm(tranfnc,'migrad', zero ,0,ierflg,0)
c                analysis of errors for all parameters
      call mnexcm(tranfnc,'minos', zero ,0,ierflg,0)
c                 call fcn with 3. i.e, ouput etc.  
      call mnexcm(tranfnc,'call fcn', three , 1,ierflg, 0)

      do i = 1, nparam
         prmout(i) = oparam(i)
      enddo

      end

