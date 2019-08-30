!c  include "BlockData/cblkGene.h"
    
      module modxsecFit
      implicit none

      integer,parameter:: maxbin=100
      integer,parameter:: nparam=3
      integer,parameter::minout=56
      integer,parameter:: minsave=7
      real(8):: oparam(nparam)	
      real(8):: param(nparam)
      real(8):: x(maxbin), y(maxbin)
      integer npoint
      real(8):: chisq
      character(2):: id
      integer:: TargA, Eregion
      end module modxsecFit

      program xsecFit
      use modxsecFit
      implicit none
      real*8 xin(maxbin), yin(maxbin)
      character*128 buf 
      integer nbin0, count, icon, status
      integer:: n1, n2 
#if defined (MACOSX)
      integer iargc
      count = iargc() +1
#else
      count = NARGS()
#endif
      if(count .ne. 3) then
         write(0,*) ' command line arg is ', count
         write(0,*)
     *      'Usage: ./xsecFit  id A  < inputfile > outfile'
         write(0,*)    ' Id:  p or pi or K'
         write(0,*)    ' A   target mass #'
         write(0,*) ' inputfile: E(GeV)  vs Cross-sec(mb)'
         stop
      endif
#if defined (MACOSX)
      call getarg(1, buf)
      read(buf, '(a)') id
      call getarg(2, buf)
      read(buf, *) TargA
#else
      call getarg(1, buf, status)
      read(buf, '(a)') id
      call getarg(2, buf, status)
      read(buf, *) TargA
#endif
!/////////////
      write(0,*) ' id, TargA=',id, TargA
!////////////
      nbin0 = 0
      do while(.true.)
         read(*,*,end=100) xin(nbin0+1), yin(nbin0+1)
!         write(0,*) xin(nbin0+1), yin(nbin0+1)
         nbin0 = nbin0 +1
      enddo
100   continue
      if(nbin0 .gt. maxbin) then
         write(0,*) ' too many data> ', maxbin
         stop
      endif


      call copenfw2(minout, "/dev/null",   1,  icon)
          ! nbin0=86
      n1 = 1
      n2 = nbin0/3 + 5   ! 28+5=33
      Eregion = 1
      call fitXsec0(xin(n1), yin(n1),  n2)

      n1 = n2 -10      !    23
      n2 = n1 + nbin0/3+10   ! 23 +28+10 = 61
      Eregion = 2
      call fitXsec0(xin(n1), yin(n1), n2-n1+1)
      n1 = n2 - 10   ! 51
      n2 = nbin0     ! 86
      Eregion = 3
      call fitXsec0(xin(n1), yin(n1), n2-n1+1)
      end

      subroutine fitXsec0(xin, yin, nbin)
      use modxsecFit
      implicit none
      integer nbin
      real*8 xin(nbin), yin(nbin)
      real*8 prmout(nparam)

      integer i 

      real*8 temp, f
      real*8 a, b, c
!  
      param(1) = yin(1)/10.  ! a
      param(2) = 0.1     ! b 
      param(3) = 0.1     ! c
! 
      call fitXsec( xin, yin, nbin, param, prmout)
!
!               to see fitted result

      a = prmout(1)
      b = prmout(2)
      c = prmout(3)
!             coeff is put on stderr
      write(*,'(a, 2i4, 1p, 3g14.4)') id, TargA, Eregion, prmout(:)
      do i = 1, nbin
         temp = log(x(i))
         f = a+ (b + c*temp) * temp
!         write(*,'(1p, 2g14.4)') x(i), f
      enddo
      end
!     ***************************************************
      subroutine fitXsec(xin, yin, n, prmin, prmout )
      use modxsecFit
      implicit none
      real*8  prmin(nparam), prmout(nparam)
      integer  n
      real*8 xin(n), yin(n)

      integer nlabel(nparam)
      character*10  pname(nparam)
      real*8 initval(nparam)
      real*8 step(nparam)

      data nlabel/ 1,  2, 3/
      data pname/ 'a',  'b', 'c'/
      data step/  1.,  0.1d0, 0.1d0/
      real*8 zero, one, three
      data zero,one,three / 0., 1., 3./
      real*8 low(nparam), high(nparam)
      real*8 fval, xx
      integer i, ierflg

      external xsecfnc
      
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
!      low(1) = prmin(1)*0.1
      low(1) = 0.
!      high(1) = prmin(1)*1d4
      high(1) = 0.
!      low(2) = prmin(2)*0.1
      low(2) = 0.
!      high(2) = prmin(2)*10
      high(2) = 0.
!      low(3) = prmin(3)*0.1
      low(3) = 0.
!      high(3) = prmin(3)*10.
      high(3) =0.
      call mninit(5, minout, minsave)

      do  i= 1, nparam
!        nprm: a number given to a parameter: (label)
!        pnam: name of the parameer
!        vstrt: initial value of the parameter
!        stp:   initial step size of the //
!        next two: zero-->the parameter is not bounded (lower or upper)
!        ierflg: retrun value; cond code. 0--> ok
         call mnparm(nlabel(i), pname(i), initval(i), step(i),
     *     low(i), high(i), ierflg)
         if (ierflg .ne. 0)  then
            write (0,'(a,i3)')  ' unable to define parameter no.',i
            stop
         endif
      enddo
!
      call mnseti('xsecfit')
!       request fcn to read in (or generate random) data (iflag=1)
!            fcnk0: function to be minimuzed is calculated. also 
!              there are other funcitons
!            one is the  argument to fcnk0.  seems to be converted to
!            integer inside.
!            1 number of argument in one  (one could be array)
!           ierflf: ouptut. 0-->ok
!            0: no external function is used in fcnk0

      call mnexcm(xsecfnc, 'call fcn', one ,1,ierflg, 0)
!        fix the  3,4,5-th parameters,  
!      call mnexcm(timefnc,'fix', fixlist ,3, ierflg,0)
!       print minumum things   
      call mnexcm(xsecfnc,'set print', zero ,1,ierflg,0)
!                use migrad method for minimization
!                with default condtions
      call mnexcm(xsecfnc,'migrad', zero ,0,ierflg,0)
!                analysis of errors for all parameters
      call mnexcm(xsecfnc,'minos', zero ,0,ierflg,0)
!                 call fcn with 3. i.e, ouput etc.  
      call mnexcm(xsecfnc,'call fcn', three , 1,ierflg, 0)

      do i = 1, nparam
         prmout(i) = oparam(i)
      enddo
      end
      subroutine xsecfnc(npar,gin,f,paramx,iflag)
      use modxsecFit 
      implicit none
! npar: input. number of current variable  parameters
! gin:  optional output of gradient
!  f:  output. function value to be minimized 
! paramx: input. vector of const and variable parameters
!  iflag: input. depending on this value, what to do 
!        is determined.
!       1--> some preparation for computing f
!       2--> compute gin
!       3--> fittng finished
!       for all other cases: we  must compute f
!
      integer  npar
      real*8   gin(*)
      real*8   f
      real*8   paramx(*)
      integer iflag
!
      integer i
      real*8 fval,  xx, temp
      real*8  a, b, c
      save
      if (iflag .eq. 1) then
!         write(0,*) 'current npar=',npar,  ' npoint=',npoint
!         write(0,*) ' param=',(paramx(i), i=1,npar)
      endif
      if(iflag .eq. 2 ) then
         write(0,*) ' no grad computed'
         stop
      endif
!        compute f

      chisq = 0.
      do  i= 1, npoint
         a = paramx(1)
         b =  paramx(2)
         c = paramx(3)
         temp = log(x(i))
         fval =  a + (b+c*temp)* temp
         chisq= chisq + (fval-y(i))**2/y(i)
      enddo
      f = chisq

      if(iflag .eq. 3) then
!         do i = 1, npoint
!            fval =   paramx(1)*x(i)**(paramx(2) +
!            write(*,*) x(i), y(i), fval
!         enddo
         do i = 1, npar
            oparam(i) = paramx(i)
         enddo
      endif
      end

