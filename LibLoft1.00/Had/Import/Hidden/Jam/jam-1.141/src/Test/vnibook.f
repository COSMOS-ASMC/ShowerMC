C $Id: vnibook.f,v 3.2.5.0 1997/08/30  2:13:42 klaus Exp klaus $
C*********************************************************************
C*********************************************************************
C
C
C                   VNIBOOK - HISTOGRAM PACKAGE
C
C                                                                
C                       Klaus Kinder-Geiger                        
C                                                                   
C                 Brookhaven National Laboratory                      
C                 Physics Department, Bldg 510-A                      
C                    Upton, N.Y. 11973, U.S.A.                        
C                     phone: (516) 344-3791                           
C                     e-mail: klaus@bnl.gov                          
C               http:/penguin.phy.bnl.gov/~klaus                     
C                                                                     
C                                                                     
C********************************************************************* 
C********************************************************************* 
C                                                                      
C
C  VNIBOOK is a simple histogram package consisting of a few fortran
C  subroutines. It provides a `portable' tool of filling and printing
C  out 1- or 2-dimensional histograms. The VNIBOOK package is based
C  on a subroutine package written by T. Sjoestrand, and is similar 
C  to the CERN library HBOOK package. The use of the subroutines may
C  be divided into 4 steps: booking, filling editing and printing.
C  All subroutines are written in fortran 77 standard. They are
C  briefly discussed below.
C
C  The common block   COMMON /VBOOK/ A(10000) is used to store histo-
C  gram information. The histogram index takes 102 positions. Each 
C  booked 1- (2-) dimensional histogram takes an additional 38+NX 
C  (38+NX+NY) positions. If more space is required, the program has 
C  to be recompiled with extended common block.
C
C***  booking  *******************************************************
C
C  There are 100 histograms at the disposal of the user, each given by
C  a number between 1 and 100. Before a histogram can be used, space 
C  must be reserved for it.
C
C  CALL VBOOK1(ID,TITLE,NX,XL,XU)
C  Purpose: to book a one-dimensional histogram.
C
C  ID: histogram number, integer between 1 and 100.
C  TITLE: histogram title, at most 60 characters.
C  NX: number of bins (in X direction) in histogram, integer between
C     1 and 100.
C  XL, XU: lower and upper bound, respectively, on the (X) range
C     covered by the histogram.
C
C  CALL VBOOK2(ID,TITLE,NX,XL,XU,NY,YL,YU)
C  Purpose: to book a two-dimensional histogram.
C
C  ID: histogram number, integer between 1 and 100.
C  TITLE: histogram title, at most 60 characters.
C  NX: number of bins in X direction, integer between 1 and 50.
C  XL, XU: lower and upper bound on the X range of histogram.
C  NY: number of bins in Y direction, arbitrary positive integer.
C  YL, YU: lower and upper bound on the Y range of histogram.
C
C***  filling  *******************************************************
C
C  For booked histograms weights can be filled at given coordinates.
C
C  CALL VFILL1(ID,X,W)
C  Purpose: to fill in a one-dimensional histogram.
C  ID: histogram number.
C  X:  X coordinate of point.
C  W:  weight to be added in this point.
C
C  CALL VFILL2(ID,X,Y,W)
C  Purpose: to fill in a two-dimensional histogram.
C  ID: histogram number.
C  X:  X coordinate of point.
C  Y:  Y coordinate of point.
C  W:  weight to be added in this point.
C
C***  editing  *******************************************************
C
C  For editing of histograms before printout, 2 routines are available.
C
C  CALL VSCALE(ID,F)
C  Purpose: to rescale the contents of a histogram.
C  ID: histogram number.
C  F:  rescaling factor, i.e. factor that all bin contents (including
C      overflow etc.) are multiplied by.
C  Remark: A typical rescaling factor for a one-dimensional histogram
C      could be 
C       F = 1/(bin size * no. events) = NX/(XU-XL) * 1/(no. events).
C
C  CALL VOPERA(ID1,OPER,ID2,ID3,F1,F2)
C  Purpose: to provide general purpose editing of one or several 
C     histograms, which all are assumed to have the same number of 
C     bins. Operations are carried out bin by bin, including overflow
C     bins etc.
C  OPER: gives the type of operation to be carried out, a one-character
C     string or a character*1 variable.
C  OPER= '+', '-', '*', '/': add, subtract, multiply, or divide the
C     contents in ID1 and ID2 and put the result in ID3. F1 and F2, if
C     not 1., give factors by which the ID1 and ID2 bin contents are
C     multiplied brfore the indicated operation. (Division by a
C     vanishing bin content will give 0).
C  OPER= 'A', 'S', 'L': for 'S' the square root of the content in ID1
C     is taken (result 0 for negative bin contents) and for 'L' the
C     10-logarithm is taken (a non-positive bin content is before
C     replaced by 0.8 times the smalles positive bin content).
C     Thereafter, in all three cases, the content is multiplied by F1
C     and added with F2, and the result is placed in ID3. Thus ID2
C     is dummy in this case.
C  OPER= 'M': intended for statistical analysis, bin-by-bin mean and
C     standard deviation of a variable, assuming that ID1 contains
C     accumulated weights, ID2 accumulated weights*variable and ID3
C     accumulated weight*variable-squared. Afterwards ID2 will
C     contain the mean values (=ID2/ID1) and ID3 the standard
C     deviations (=SQRT(ID3/ID1-(ID2/ID1)**2)). In the end, F1
C     multiplies ID1 (for normalization purposes), while F2 is dummy.
C
C***  printing  ******************************************************
C
C  At printing, the bin contents are listed in tabular form, and the 
C  mean-values, standard deviations and overflow are given at the end.
C  The printout is directed to be written to disc, into files named 
C   "vbook.???", C  where ???=ID=histogram number (e.g. vbook.001).
C
C  CALL VCLEAR
C  Purpose: to print out all histograms that have been filled, and
C     reset them thereafter to 0.
C
C  CALL VPRINT(ID)
C  Purpose: to print out a single histogram.
C  ID: histogram to be printed.
C
C  CALL VRESET(ID)
C  Purpose: to reset all bin contents, including overflow etc. to 0.
C  ID: histogram to be reset.
C
C*********************************************************************
C*********************************************************************
 
      subroutine vbook1(id,title,nx,xl,xu)
 
C...Double precision declaration.
      implicit double precision(a-h, o-z)
C...Commonblock.
      common/vnibins/ihist(4),indx(1000),bin(20000)
      save /vnibins/
C...Local character variables.
      character title*(*), titfx*60,check*70
 
C...Check that input is sensible. Find initial address in memory.
      if(id.le.0.or.id.gt.ihist(1)) then
         write(check,*)id,ihist(1)
          call pjerrm(28,
     &       '(VBOOK1:) not allowed histogram number'//check)
      endif

      if(nx.le.0.or.nx.gt.100) then
         write(check,*)' id nx=',id,nx
         call pjerrm(28,'(VBOOK1:) not allowed number of bins'//check)
      endif

      if(xl.ge.xu) then
         write(check,*)id,nx,xl,xu
         call pjerrm(28,'(VBOOK1:) x limits in wrong order'//check)
      endif

      indx(id)=ihist(4)
      ihist(4)=ihist(4)+38+nx
      if(ihist(4).gt.ihist(2)) call pjerrm(28,
     &'(VBOOK1:) out of histogram space')
      is=indx(id)
 
C...Store histogram size and reset contents.
      bin(is+1)=nx
      bin(is+2)=xl
      bin(is+3)=xu
      bin(is+4)=(xu-xl)/nx
      bin(is+5)=1.0d0
      call vreset(id)
 
C...Store title by conversion to integer to double precision.
      titfx=title//' '
      do it=1,20
        bin(is+18+nx+it)=256**2*ichar(titfx(3*it-2:3*it-2))+
     &  256*ichar(titfx(3*it-1:3*it-1))+ichar(titfx(3*it:3*it))
      end do
 
      return
      end

 
C*********************************************************************
 
      subroutine vbook2(id,title,nx,xl,xu,ny,yl,yu)

C...Purpose: to reset and book 2-dim histogram.
C...Commonblock.
      implicit double precision(a-h, o-z)
      common/vnibins/ihist(4),indx(1000),bin(20000)
      save /vnibins/
C...Local character variables.
      character title*(*), titfx*60,check*80

C...Identification; store histogram specifics.
  
      indx(id)=ihist(4)
      ihist(4)=ihist(4)+38+nx*ny
      if(ihist(4).gt.ihist(2)) then
        write(check,*)ihist(4),ihist(2)
        call pjerrm(28,'(VBOOK2:) out of histogram space'//check)
      endif
      is=indx(id)
      bin(is+1)=nx
      bin(is+2)=xl
      bin(is+3)=xu
      bin(is+4)=(xu-xl)/nx
      bin(is+5)=ny
      bin(is+6)=yl
      bin(is+7)=yu
      bin(is+8)=(yu-yl)/ny

C...Reset.
      call vreset(id)

C...Store title.
      titfx=title//' '
      do it=1,20
        bin(is+18+nx*ny+it)=256**2*ichar(titfx(3*it-2:3*it-2))+
     &  256*ichar(titfx(3*it-1:3*it-1))+ichar(titfx(3*it:3*it))
      end do

      return
      end
 
C*********************************************************************
 
      subroutine vfill1(id,x,w)

C...Purpose: to accumulate statistics and fill 1-dim. histogram.
      implicit double precision(a-h, o-z)
      common/vnibins/ihist(4),indx(1000),bin(20000)
      common /vcall/ nevt(1000),nhit(1000)
      save /vnibins/,/vcall/
      character check*80

C...Counter.
      nevt(id)=nevt(id)+1
      
C...Identification.
      is=indx(id)
      if(is.eq.0) then
        write(check,*)id
        call pjerrm(28,'(VFILL1:) filling unbooked histogram id='//
     $         check)
      endif
      bin(is+9)=bin(is+9)+1.d0

C...Find bin in x, including under/overflow, and fill.
      if(x.lt.bin(is+2)) then
        bin(is+13)=bin(is+13)+w
      elseif(x.ge.bin(is+3)) then
        bin(is+15)=bin(is+15)+w
      else
        bin(is+14)=bin(is+14)+w
        ix=(x-bin(is+2))/bin(is+4)
        ix=max(0,min(nint(bin(is+1))-1,ix))
        bin(is+19+ix)=bin(is+19+ix)+w
        nhit(id)=nhit(id)+1
      endif


      return
      end
 
C*********************************************************************
 
      subroutine vfill2(id,x,y,w)

C...Purpose: to accumulate statistics and fill 2-dim. histogram.
      implicit double precision(a-h, o-z)
      common/vnibins/ihist(4),indx(1000),bin(20000)
      common /vcall/ nevt(1000),nhit(1000)
      save /vnibins/,/vcall/
      character check*80

C...Counter.
      nevt(id)=nevt(id)+1

C...Identification.
      is=indx(id)
      if(is.eq.0) then
        write(check,*)id
        call pjerrm(28,'(VFILL2:) filling unbooked histogram id='//
     $         check)
      endif

C...Check over/underflow.
      bin(is+9)=bin(is+9)+1.d0
      iox=2
      if(x.lt.bin(is+2)) iox=1
      if(x.ge.bin(is+3)) iox=3
      ioy=2
      if(y.lt.bin(is+6)) ioy=1
      if(y.ge.bin(is+7)) ioy=3
      bin(is+6+3*ioy+iox)=bin(is+6+3*ioy+iox)+w
      if(iox.ne.2.or.ioy.ne.2) return

C...Fill bins.
      nhit(id)=nhit(id)+1
      ix=(x-bin(is+2))/bin(is+4)
      iy=(y-bin(is+6))/bin(is+8)
      ic=int(bin(is+1)+0.5d0)*iy+ix
      bin(is+19+ic)=bin(is+19+ic)+w

      return
      end
 
C*********************************************************************
 
      subroutine vscale(id,f)

C...Purpose: to rescale histogram.
      implicit double precision(a-h, o-z)
      common/vnibins/ihist(4),indx(1000),bin(20000)
      save /vnibins/
      character check*80

C...Multiply entries by factor F.
      is=indx(id)
      if(is.eq.0) then
        write(check,*)id
        call pjerrm(28,'(VSCALE:) filling unbooked histogram id='//
     $         check)
      endif
      do ix=is+10,is+18+nint(bin(is+1))*nint(bin(is+5))
       bin(ix)=f*bin(ix)
      end do

      return
      end
 
C*********************************************************************
 
      subroutine vopera(id1,oper,id2,id3,f1,f2)

C...Purpose: to perform various editing operations OPER on  histogram.
      implicit double precision(a-h, o-z)
      common/vnibins/ihist(4),indx(1000),bin(20000)
      save /vnibins/
      character oper*(*)

C...Identification.
      is1=indx(id1)
      is2=indx(min(ihist(1),max(1,id2)))
      is3=indx(min(ihist(1),max(1,id3)))
      nc=nint(bin(is3+1))*nint(bin(is3+5))

C...Type of operation.
      if(oper.eq.'+'.or.oper.eq.'-'.or.oper.eq.'*'.or.oper.eq.'/') then
        bin(is3+9)=bin(is1+9)+bin(is2+9)
      else if(oper.eq.'A'.or.oper.eq.'S'.or.oper.eq.'L') then
        bin(is3+9)=bin(is1+9)
      endif

C...Addition, subtraction, multiplication, division.
      if(oper.eq.'+') then
        do ic=10,18+nc
        bin(is3+ic)=f1*bin(is1+ic)+f2*bin(is2+ic)
        end do
      elseif(oper.eq.'-') then
        do ic=10,18+nc
        bin(is3+ic)=f1*bin(is1+ic)-f2*bin(is2+ic)
        end do
      elseif(oper.eq.'*') then
        do ic=10,18+nc
        bin(is3+ic)=f1*bin(is1+ic)*f2*bin(is2+ic)
        end do
      elseif(oper.eq.'/') then
        do ic=10,18+nc
          fa2=f2*bin(is2+ic)
          if(abs(fa2).le.1d-10) bin(is3+ic)=0.d0
          if(abs(fa2).gt.1d-10) bin(is3+ic)=f1*bin(is1+ic)/fa2
        end do

C...Statistical analysis.
      elseif(oper.eq.'A') then
        do ic=10,18+nc
        bin(is3+ic)=f1*bin(is1+ic)+f2
        end do
      elseif(oper.eq.'S') then
        do ic=10,18+nc
        bin(is3+ic)=f1*sqrt(max(0.d0,bin(is1+ic)))+f2
        end do
      elseif(oper.eq.'L') then
        zmin=1d30
        do ic=19,18+nc
        if(bin(is1+ic).lt.zmin.and.bin(is1+ic).gt.1d-20) 
     &  zmin=0.8d0*bin(is1+ic)
        end do
        do ic=10,18+nc
        bin(is3+ic)=f1*dlog10(max(bin(is1+ic),zmin))+f2
        end do
      elseif(oper.eq.'M') then
        do ic=10,18+nc
          if(abs(bin(is1+ic)).le.1d-10) then
            bin(is2+ic)=0.d0
          else if(abs(bin(is1+ic)).gt.1d-10) then
            bin(is2+ic)=bin(is2+ic)/bin(is1+ic)
          endif
          if(id3.ne.0) then
            if(abs(bin(is1+ic)).le.1d-10) then
              bin(is3+ic)=0.d0
            else if(abs(bin(is1+ic)).gt.1d-10) then
              bin(is3+ic)=
     &        sqrt(max(bin(is3+ic)/bin(is1+ic)-bin(is2+ic)**2,0.d0))
            endif
          endif
          bin(is1+ic)=f1*bin(is1+ic)
        end do
      endif

      return
      end
 
C*********************************************************************
 
      subroutine vclear(mnorm,mform)

C...Purpose: to print out and then reset all bin contents to 0.
      implicit double precision(a-h, o-z)
      common/vnibins/ihist(4),indx(1000),bin(20000)
      save /vnibins/

      do 100 id=1,1000
        is=indx(id)
        if(is.eq.0.or.bin(is+9).lt.0.5d0) goto 100
        call vprint(id,mnorm,mform)
        call vreset(id)
  100 continue

      return
      end
 
C*********************************************************************
 
      subroutine vprint(id,mnorm,mform)

C...Purpose: to print out accumulated statistics in table format.
C...MNORM = 0 ( 1 ) gives unnormalized (normalized) distribution,
C...MFORM = 0 ( 1 ) prints F ( E ) format of values.
      implicit double precision(a-h, o-z)
      common/vnibins/ihist(4),indx(1000),bin(20000)
      common /vcall/ nevt(1000),nhit(1000)
      save /vnibins/,/vcall/
      character title*60,cha(40)*1,chfile*9
C...Internal arrays.
      parameter (nz=1000)
      dimension xz(nz),yz(nz),zz(nz,nz)

C...Variables for time: NORSK DATA, VAX, IBM, CDC, SUN
C     DIMENSION IDAT(7)
      character ctime*8
C     CHARACTER DATE*8,TIME*8
C     CHARACTER DATE*10,TIME*10,CDATE*10,CTIME*10
C     DIMENSION IDAT(3),ITIM(3)
      dimension dyac(10),ev(20)
      character check*80
C...LIN gives maximum length for histogram proper.
      data lin/37/, dyac/.04d0,.05d0,.06d0,.08d0,.10d0,.12d0,.15d0, 
     & .20d0,.25d0,.30d0/
      data cha/' ','0','1','2','3','4','5','6','7','8','9','A','B',
     &'C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q',
     &'R','S','T','U','V','W','X','Y','Z','*','-','!'/

C...Indentification, reset internal arrays.
      is=indx(id)
      if(is.eq.0) then
        write(check,*)id
        call pjerrm(28,'(VPRINT:) filling unbooked histogram id='//
     $         check)
      endif
      nx=nint(bin(is+1))
      ny=nint(bin(is+5))
      if((nx.le.0.or.nx.gt.1000).or.(ny.le.0.or.ny.gt.1000)) then
        write(check,*)id,nx,ny
        call pjerrm(28,'(VPRINT:) invalid nx ny, id nx ny='//
     $         check)
      endif
      do ix=1,nx
      xz(ix)=0.d0
      end do
      do iy=1,ny
      yz(ix)=0.d0
      end do
      do ix=1,nx
      do iy=1,ny
      zz(ix,iy)=0.d0
      end do
      end do
      do it=1,20
      ev(it)=0.d0
      ieq=nint(bin(is+18+nx*ny+it))
      title(3*it-2:3*it)=char(ieq/256**2)
     &//char(mod(ieq,256**2)/256)//char(mod(ieq,256))
      end do

C...Open output file.
      iunit=10+id
      if(iunit.gt.99) iunit=iunit-89
      if(id.lt.10) write(chfile,400) id
      if(id.ge.10.and.id.lt.100) write(chfile,410) id
      if(id.ge.100.and.id.lt.1000) write(chfile,420) id
      open(unit=iunit,file=chfile,status='UNKNOWN')

C...Exit, if not data have been filled.
      if(bin(is+9).lt.0.5d0) then
        write(iunit,500) title,id
        return
      endif

C...Writing of title.
      write(iunit,1000) title,id,nevt(id),nhit(id)

C...Writing of time: NORSK DATA, VAX, IBM, CDC, SUN
C     CALL CLOCK(IDAT)
C     WRITE(IUNIT,1100) (IDAT(8-I),I=1,5)
C      CALL IDATE(IMON,IDAY,IYEAR)
C      CALL TIME(CTIME)
c     write(iunit,1100) iyear, imon, iday, ctime(1:5)
C     CALL DAY(DATE,TIME)
C     WRITE(IUNIT,1100) DATE(7:8), DATE(4:5), DATE(1:2),
C    &TIME(1:2), TIME(4:5)
C     CALL DATIME(IDAT,ITIM)
C     WRITE(IUNIT,1100) 1900+IDAT/10000,MOD(IDAT/100,100),
C    &MOD(IDAT,100),ITIM/100,MOD(ITIM,100)
C     CDATE=DATE()
C     CTIME=TIME()
C     WRITE(IUNIT,1100) CDATE(2:3),CDATE(5:6),CDATE(8:9),
C    &CTIME(2:3), CTIME(5:6)
C     CALL IDATE(IDAT)
C     CALL ITIME(ITIM)
C     WRITE(IUNIT,1100) IDAT(3),IDAT(2),IDAT(1),ITIM(1), 
C    &ITIM(2)

C......Normalization.
        fev=1.d0
        fac=1.d0
        if(mnorm.eq.1) then
          fev=1.d0/max(1,nevt(id))
          fac=1.d0/max(1,nhit(id))
        endif

C...1-dimensional histogram:
      if(ny.eq.1) then

C......Print out accumulated statistics.
        write(iunit,1200) 
        is=indx(id)
        nx=int(bin(is+1)+0.5d0)
        ys=0.d0
        ya=0.d0
        do 200 ix=1,nx
        xx=bin(is+2)+(ix-0.5d0)*bin(is+4)
        yx=bin(is+18+ix)
        ys=ys+yx
        ya=ya+yx*bin(is+4)
        if(mform.eq.0) write(iunit,1300) xx,fev*yx
        if(mform.eq.1) write(iunit,1400) xx,fev*yx
  200   continue
        if(mform.eq.0) write(iunit,1600) fac*ys
        if(mform.eq.1) write(iunit,1800) fac*ys

C...Statistical analysis.
        do 210 iy=1,ny
        y=bin(is+6)+(iy-0.5d0)*bin(is+8)
        do 210 ix=1,nx
        x=bin(is+2)+(ix-0.5d0)*bin(is+4)
        cta=abs(bin(is+18+nx*(iy-1)+ix))
        ev(1)=ev(1)+cta
        ev(2)=ev(2)+cta*x
        ev(3)=ev(3)+cta*x**2
        ev(4)=ev(4)+cta*y
        ev(5)=ev(5)+cta*y**2
        ev(6)=ev(6)+cta*x*y
  210   continue
        xmean=ev(2)/max(ev(1),1d-20)
        xrms=sqrt(max(0.d0,ev(3)/max(ev(1),1d-20)-xmean**2))
        ymean=ev(4)/max(ev(1),1d-20)
        yrms=sqrt(max(0.d0,ev(5)/max(ev(1),1d-20)-ymean**2))
        xycor=(ev(6)/max(ev(1),1d-20)-xmean*ymean)/max(1d-20,xrms*yrms)
        write(iunit,1900) int(bin(is+9)+0.5d0),xmean,
     $   bin(is+2),xrms,bin(is+3),
     &  ymean,bin(is+6),yrms,bin(is+7),xycor
        write(iunit,2000) (bin(is+j),j=16,18),(bin(is+j), j=13,15),
     &  (bin(is+j),j=10,12)


C...2-dimensional histogram:
      else

C......Print out accumulated statistics.
c       is=bin(id+2)+0.5
        is=indx(id)
        nx=int(bin(is+1)+0.5d0)
        ny=int(bin(is+5)+0.5d0)
        do ix=1,nx
        xz(ix)=bin(is+2)+(ix-0.5d0)*bin(is+4)
        end do
        do iy=1,ny
        yz(iy)=bin(is+6)+(iy-0.5d0)*bin(is+8)
        end do
        do ix=1,nx
        do iy=1,ny
        zz(ix,iy)=bin(is+18+nx*(iy-1)+ix)
        end do
        end do
        nys=ny
        if(nys.gt.15) ny=20
        if(nys.le.15) ny=15
        if(nys.le.10) ny=10
        if(mform.eq.0) then
          if(ny.eq.10) then
            write(iunit,2300) (yz(iy),iy=1,ny)
            write(iunit,2400) (xz(ix),(fev*zz(ix,iy),iy=1,ny),
     &      ix=1,nx)
          elseif(ny.eq.15) then
            write(iunit,2500) (yz(iy),iy=1,ny)
            write(iunit,2600) (xz(ix),(fev*zz(ix,iy),iy=1,ny),
     &      ix=1,nx)
          elseif(ny.eq.20) then
            write(iunit,2700) (yz(iy),iy=1,ny)
            write(iunit,2800) (xz(ix),(fev*zz(ix,iy),iy=1,ny),
     &      ix=1,nx)
          endif
        elseif(mform.eq.1) then
          if(ny.eq.10) then
            write(iunit,2900) (yz(iy),iy=1,ny)
            write(iunit,3000) (xz(ix),(fev*zz(ix,iy),iy=1,ny),
     &      ix=1,nx)
          elseif(ny.eq.15) then
            write(iunit,3100) (yz(iy),iy=1,ny)
            write(iunit,3200) (xz(ix),(fev*zz(ix,iy),iy=1,ny),
     &      ix=1,nx)
          elseif(ny.eq.20) then
            write(iunit,3300) (yz(iy),iy=1,ny)
            write(iunit,3400) (xz(ix),(fev*zz(ix,iy),iy=1,ny),
     &      ix=1,nx)
          endif
        endif

C...Statistical analysis.
        do 250 iy=1,ny
        y=bin(is+6)+(iy-0.5d0)*bin(is+8)
        do 250 ix=1,nx
        x=bin(is+2)+(ix-0.5d0)*bin(is+4)
        cta=abs(bin(is+18+nx*(iy-1)+ix))
        ev(1)=ev(1)+cta
        ev(2)=ev(2)+cta*x
        ev(3)=ev(3)+cta*x**2
        ev(4)=ev(4)+cta*y
        ev(5)=ev(5)+cta*y**2
        ev(6)=ev(6)+cta*x*y
  250   continue
        xmean=ev(2)/max(ev(1),1d-20)
        xrms=sqrt(max(0.d0,ev(3)/max(ev(1),1d-20)-xmean**2))
        ymean=ev(4)/max(ev(1),1d-20)
        yrms=sqrt(max(0.d0,ev(5)/max(ev(1),1d-20)-ymean**2))
        xycor=(ev(6)/max(ev(1),1d-20)-xmean*ymean)/max(1d-20,xrms*yrms)
        write(iunit,1900) int(bin(is+9)+0.5d0),xmean,bin(is+2)
     $                    ,xrms,bin(is+3),
     &  ymean,bin(is+6),yrms,bin(is+7),xycor
        write(iunit,2000) (bin(is+j),j=16,18),(bin(is+j), j=13,15),
     &  (bin(is+j),j=10,12)

      endif


C...Close output file.
      close(iunit)


C...Format statements for output file.
  400 format('vbook.00',i1)
  410 format('vbook.0',i2)
  420 format('vbook.',i3)

C...Format statements for output file. information on results
c...and errors.
  500 format('#'/'#',5x,'VBOOK: ',a60/'#'/'#',5x,'histogram no.: ',
     &i8,' - no entries!')

C...Format statement for title.
 1000 format('#'/'#',5x,'VBOOK: ',a60/'#'/'#',5x,'histogram no.: ',
     &i9/'#',5x,'no. of events: ',i9/'#',5x,'no. of hits:   ',i9/'#')

C...Format statement for time: NORSK DATA, IBM VS, VAX, IBM SIE, CDC
C 1100 FORMAT('#'/'#',5X,'date: ',I4,'-',I2,'-',I2,I3,':',I2/'#')
 1100 format('#'/'#',5x,'date: ','19',i2,'-',i2,'-',i2,1x,a5/'#')
C 1100 FORMAT('#'/'#',5X,'date: ','19',A2,'-',A2,'-',A2,1X,A2,':',
C      A2/'#')
C 1100 FORMAT('#'/'#',5X,'date: ','19',A2,'-',A2,'-',A2,1X,A2,':',
C      A2/'#')

C...Format statements for data output.
 1200 format('#',10x,'   x   ',10x,'    y(x)    '/'#')
 1300 format(1x,f15.5,5x,f15.5)
 1400 format(1x,e15.5,5x,e15.5)
 1600 format('#',9x,'sum y:',5x,f15.5)
 1800 format('#',9x,'sum y:',5x,e15.5)
 1900 format('#'/'#'/'#',5x,'entries =',i10,1p,5x,'xmean =',e10.3/
     &'#',5x,'xmin    =',e10.3,5x,'xrms  =',e10.3/'#',5x,'xmax    =',
     &e10.3,5x,'ymean =',e10.3/'#',5x,'ymin    =',e10.3,5x,'yrms  =',
     &e10.3/'#',5x,'ymax    =',e10.3,5x,'xycor =',e10.3)
 2000 format('#'/'#'/'#',10x,e10.3,' ! ',e10.3,' ! ',e10.3,/'#',10x,
     &36('-')/'#',1x,'overflow',1x,e10.3,' ! ',e10.3,' ! ',e10.3,7x/
     &'#',10x,36('-')/'#',10x,e10.3,' ! ',e10.3,' ! ',e10.3)
 2300 format('#','       y  ',10f9.4/'#','  x   '/'#')
 2400 format(f9.4,2x,10f9.4)
 2500 format('#','     y  ',15f7.3/'#','  x   '/'#')
 2600 format(f7.3,2x,15f7.3)
 2700 format('#','   y  ',20f5.1/'#','  x   '/'#')
 2800 format(f5.1,2x,20f5.1)
 2900 format('#','       y     ',10e10.3/'#','  x   '/'#')
 3000 format(e10.3,4x,10e10.3)
 3100 format('#','     y     ',15e8.1/'#','  x   '/'#')
 3200 format(e8.1,4x,15e8.1)
 3300 format('#','   y     ',20e8.1/'#','  x   '/'#')
 3400 format(e8.1,4x,20e8.1)


      return
      end

C*********************************************************************
 
      subroutine vreset(id)

C...Purpose: to reset histogram contents to 0.
      implicit double precision(a-h, o-z)
      common/vnibins/ihist(4),indx(1000),bin(20000)
      save /vnibins/

      is=indx(id)
      do ic=is+9,is+18+nint(bin(is+1))*nint(bin(is+5))
      bin(ic)=0.d0
      end do

      return
      end
 
C*********************************************************************
 
      subroutine vdump(file)

C...Purpose: to dump histograms into file.
      implicit double precision(a-h, o-z)
      common/vnibins/ihist(4),indx(1000),bin(20000)
      save /vnibins/
      character*(*) file
      character title*60,filnam*60
      equivalence (req,ieq)
      dimension val(5)

C...Find and open file for dump.
c     filnam='/home/surya11/klaus'//file//' '
      filnam=file
C***STATUS for SUN and DEC.
C     open(23,file=filnam,status='replace',form='formatted') 
      open(23,file=filnam,status='unknown',form='formatted')

C...Loop over histograms and find which booked.
      do 130 id=1,1000
      is=indx(id)
      if(is.eq.0) goto 130

C...Write title and histogram size.
      nx=int(bin(is+1)+0.5d0)
      do 100 it=1,20
      req=bin(is+18+nx+it)
  100 title(3*it-2:3*it)=char(ieq/256**2)//char(mod(ieq,256**2)/256)
     &//char(mod(ieq,256))
      write(23,'(I5,5X,A60)') id,title
      write(23,'(I5,2D14.6)') nx,bin(is+2),bin(is+3)

C...Write histogram contents, in groups of five.
      do 120 ixg=1,(nx+4)/5
      do 110 ixv=1,5
      ix=5*ixg+ixv-5
      if(ix.le.nx) then
        val(ixv)=bin(is+18+ix)
      else
        val(ixv)=0.d0
      endif
  110 continue
  120 write(23,'(5e14.6)') (val(ixv),ixv=1,5)

C...Go to next histogram; finish.
  130 continue


      return
      end

C*********************************************************************
 
      block data vdata

C...Purpose: to set array A equal to 0 initially.
      implicit double precision(a-h, o-z)
      common/vnibins/ihist(4),indx(1000),bin(20000)
      save /vnibins/

      data indx/1000*0/
      data ihist/1000,20000,55,1/

      end
 
C*********************************************************************
