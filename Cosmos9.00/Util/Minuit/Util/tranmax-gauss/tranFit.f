#include "BlockData/cblkGene.h"
      include "ZfitBD.h"
      program tranFit
      implicit none
      include "Zfit.h"

      real*8 xin(maxbin), yin(maxbin)
      character*128 buf 
      integer nbin0, count, icon, status


      count = NARGS()
      if(count .ne. 2) then
         write(0,*) 'Usage: tranaFit... op  < inputfile > outfile'
         write(0,*) "op: 1 stdout:  fitted x,y line by line "
         write(0,*) "op: 2 stdout: fitting coeff. a,s,xmax"
         write(0,*) "op: 3 stdout: xmax, apparent_xmax"
         write(0,*) 
     *       "op: 4 stdout: xmax, apparent_xmax, a,s,xmax"
         write(0,*)
     *       "op: 5 stdout: xmax, apparent_xmax, a,s,xmax"
         write(0,*) "    fitted   x y line by line  "
         write(0,*) ' inputfile: depth vs flux'
         stop
      endif
      call getarg(1, buf, status)
      read(buf, *) op
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
      real*8 a, s, x0




!          find max pos and use max of 9 points
!          around max.
!     
      maxval = yin(1)-1.0e10
      do i = 1, nbin
         if(maxval .lt. yin(i))  then
            maxpos = i
            maxval = yin(maxpos) 
         endif
      enddo

      maxdep = xin(maxpos)
      from = max(maxpos-2, 1)
      to = min(maxpos+2, nbin)
!   
      write(0,*) ' max depth index and depth=', maxpos, maxdep
      write(0,*) ' layers for fitting are from', from, ' to ', to
!            a* exp(-((x-x0)/s)**2/2)
!  
      param(1) = yin(maxpos)   ! a
      param(2) =  50.          ! s
      param(3) = maxdep       !  x0 
      n1 = 0
      do i = from, to
         n1 = n1 + 1
         xuse(n1) = xin(i)
         yuse(n1) = yin(i)
      enddo
! 
      call fitTran( xuse, yuse, n1, param, prmout)
!

!               to see fitted result

      a = prmout(1)
      s = prmout(2)
      x0 = prmout(3)
!             coeff is put on stderr
      write(0,'(3g12.4)')  prmout
      xmax = x0
      write(0,'("max dep=",g12.3)') x0

      if(op .eq. 1) then
         xx = drx1*xmax
         do while ( xx .le. drx2*xmax )
            f= a*exp( - ((xx-x0)/s)**2/2.)
            write(*,*) xx, f
            xx = xx*10.0**0.002
         enddo
      elseif(op .eq. 2) then
         write(*,'(5G14.4)')  a, s, x0
      elseif(op .eq. 3) then 
         write(*,*)  xmax, maxdep
      elseif(op .eq. 4) then
         write(*,'(7g14.4)') xmax, maxdep, a,s,x0
      elseif(op .eq. 5) then
         write(*,'(7g14.4)') xmax, maxdep, a,s,x0
         xx = drx1*xmax
         do while ( xx .le. drx2*xmax )
            f= a*exp( - ((xx-x0)/s)**2/2.)
            write(*,*) xx, f
            xx = xx*10.0**0.002
         enddo
      endif
      end
!     ***************************************************
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

      data nlabel/ 1,  2, 3/
      data pname/ 'a',  's', 'x0'/
      data step/  10.,  1.0, 1.0/
      real*8 zero, one, three, four, five
      data zero,one,three,four, five / 0., 1., 3.,4., 5. /
      real*8 low(nparam), high(nparam)
      real*8 fval, xx
      integer i, ierflg

      external tranfnc
      
!
!           in fortran mode, this must be called for a new fnc
!
      npoint = n
      do i = 1, npoint
         x(i) = xin(i)
         y(i) = yin(i)
      enddo

      do i = 1, nparam
         initval(i) = prmin(i)
      enddo
      low(1) = prmin(1)*0.1
      high(1) = prmin(1)*4
      low(2) = prmin(2)*0.1
      high(2) = prmin(2)*10.
      low(3) = prmin(3)-12.5
      high(3) = prmin(3)+12.5

      call mninit( 5, minout, minsave)

      do  i= 1, nparam
!        nprm: a number given to a parameter: (label)
!        pnam: name of the parameer
!        vstrt: initial value of the parameter
!        stp:   initial step size of the //
!        next two: zero-->the parameter is not bounded (lower or upper)
!        ierflg: retrun value; cond code. 0--> ok

         call mnparm(nlabel(i), pname(i), initval(i), step(i),
!     *    low(i, code), up(i, code), ierflg)
     *     low(i), high(i), ierflg)

         if (ierflg .ne. 0)  then
            write (0,'(a,i3)')  ' unable to define parameter no.',i
            stop
         endif
      enddo
!
      call mnseti('tranfit')
!       request fcn to read in (or generate random) data (iflag=1)
!            fcnk0: function to be minimuzed is calculated. also 
!              there are other funcitons
!            one is the  argument to fcnk0.  seems to be converted to
!            integer inside.
!            1 number of argument in one  (one could be array)
!           ierflf: ouptut. 0-->ok
!            0: no external function is used in fcnk0

      call mnexcm(tranfnc, 'call fcn', one ,1,ierflg, 0)
!        fix the  3,4,5-th parameters,  
!      call mnexcm(timefnc,'fix', fixlist ,3, ierflg,0)
!       print minumum things   
      call mnexcm(tranfnc,'set print', zero ,1,ierflg,0)
!                use migrad method for minimization
!                with default condtions
      call mnexcm(tranfnc,'migrad', zero ,0,ierflg,0)
!                analysis of errors for all parameters
      call mnexcm(tranfnc,'minos', zero ,0,ierflg,0)
!                 call fcn with 3. i.e, ouput etc.  
      call mnexcm(tranfnc,'call fcn', three , 1,ierflg, 0)

      do i = 1, nparam
         prmout(i) = oparam(i)
      enddo

      end

