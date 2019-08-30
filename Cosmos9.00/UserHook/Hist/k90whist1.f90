module modHistogram
  implicit none

  integer:: BinWrite

  type histogram
#if defined  (KEKA) || defined (KEKB)
     sequence
#endif
     character(4):: label
     character(8):: unit
     real::  inc
     real:: bin
     real:: xmin
     real:: xm
     real:: sumw 
     integer(2):: nhist
     integer(2):: cent
     integer(2):: imin
     integer(2):: imax
     integer(2):: step
     logical(1):: tklg
     logical(1):: ufl
     logical(1):: ofl
  end type histogram

  type histoc
#if defined  (KEKA) || defined (KEKB)
     sequence
#endif
     real norm
     real pw
     integer*2 eventno
     logical*1 logv
     character*8  init
     character*128 title
     character*8  categ
     character*96 id 
     character*32 dNunit
     character*128 dir
  end type histoc
end module modHistogram

module modHistogram1
  use modHistogram
  implicit none
  
  type histogram1
#if defined  (KEKA) || defined (KEKB)
     sequence  
#endif
     type(histoc) c
     type(histogram) x
     real, allocatable ::  xw(:)
     real, allocatable ::  dnw(:)
     real, allocatable ::  dndx(:)
     real, allocatable ::  mean(:)
     real, allocatable ::  xforI(:)
     real, allocatable ::  integ(:)
  end type histogram1
  contains
!     
!            weighted histograming  fortan 90 version
!            (Not work under Absoft fortran 90) 
!      Usage:  kwhisti:   insatnciate one histogram
!              kwhistc:   clear histogram area
!              kwhist:    take histogram
!              kwhists0:  specify integral to be from small or large
!              kwhists:   compute statistical result.
!                        This can be used more than once 
!                        with a differennt normalization factor
!                        for the same histogram.   
!                 Therefore, in bin2ascii you may use inorm=-1 to 
!                    keep the previous value.

!              kwhistp:   print statistical result or binary write histogram result
!                         by calling kwhistpr or kwhistw
!              kwhistw:  write histogram with binary format
!                        for later use.
!              kwhistr:  read histogram written by kwhistw
!              kwhista:  add two histograms with identical
!                        structure.
!              kwhistpr: print histogram 
!              print format
!         #hist1 

!         ...
!         ...
!         0 0 0 0 0 0
!
!      MinIndex: min. bin index where non zero data is stored
!      MaxIndex: max. bin index where non zero data is stored
!
    subroutine kwhisti(h, ixmin, ibinORxmax, inbin, itklg )
      use modHistogram
      implicit none
!         instanciate
      integer inbin  ! input. request inbin histogram area ( actual bin number
                     ! will be <= inbin-1. ).
      real ixmin     ! input. xmin. not in log even if log10(variable) is taken
                     !         see itklg 
      real ibinORxmax  ! input. bin or ixmax. depends on itklg.
                     !  If bin and log10 is taken, bin is for log10 

      integer itklg  ! input.  bit pattern. give it like b'10001'
                     !         bit 1 is LSB.
                     !         bit 1: 0--> not take log10 of variable
                     !                1--> take log10   //
                     !             2: 0--> ixmin is the min of lowest bin
                     !                    |---|---|---|....     |...|
                     !                    |                         |
                     !                    ixmin                     ixmax
                     !                1--> ixmin is the center of the lowest bin
                     !                  |--*--|-----|-----|....    |--*--|
                     !                     |                          |
                     !                     ixmin                      ixmax
                     !            max follows the same rule.
   
                     !             3: 0--> neglect underflow
                     !                1--> underflow is put in lowest bin
                     !                     mean bin value is affected by
                     !                     those with underflowed values
                     !             4: 0--> neglect overflow
                     !                1--> overflow is put in the highest bin
                     !                     mean bin value is affected by
                     !                     overflowed ones    
                     !             5: 0-->ibinORxmax  is the bin
                     !                        xmax is determined by bin,
                     !                        xmin and inbin
                     !                1-->ibinORxmax  is ixmax. 
                     !                        bin is determined by xmax xmin
                     !                        and inbin.


      type(histogram1) h, h1, h2
!     ====================      
      integer fno  !  if < 0, standard output is used else fno is used for histogram output
                   !  fno must be opened by the user beforehand.
      real inorm   !  input. used in the normalization as dN/dx/inorm
                   !  if 0, area normalization is tried.
!
      real  x, w 
      real  xx
      integer*2  nbin
      integer   i, ndiv, dealloc
      logical  asmax
      real*8 isumw
      real dx, normf
      real xmin, xmax
      logical oldfmt/.false./, old
      integer iofl, iufl
      integer bfnow  !  binary  write file no.
      integer bfnor  !  binary  read file no.
      integer icon   !  0; binary read was successful
                     !  1; unexpected EOF
!
!       next one will be used as key (gnuplot); it can be different
!       from graph to graph
      character *(*) id  ! input (max 96 char)
!       
!        all the rest  input below will be used only once 
!        (first one will be used).  if graphs are overlapped.
! 
      integer evno   !  input event no. 
      character*(*)  title  !  input. graph title.   max 128 char.
      character*(*)  categ   !  input. max 8 char. representing histogram.
                            !         say, lat, lateral, energy, etc
                            ! which can be a directory name.
      character*(*)  dNunit !  input. max 32 char. what is dN. May have  () etc
                            !   say, GeV/(kg/m2)
      real pw        !  input. vertical scale is displayed by x**pw dN/dx in default
      logical logv   !  input. if T, vertical scale is displayed in log in default
      character*(*)  label  !  input. a few char string for x-axis  max 4 char.
                   !       r, E, Erg, L, T, time,  l, a-b  are valid.
                   !       (r),  E! E*  etc are not valid
      character*(*)  unit   !  input. a few char string for x-axis unit
            !  max of 8 char.  m2, m^2  cm, s, GeV etc are valid
      character*(*)  dir   ! directory where split ascii histogram data is
                           ! is stored.  maindir/categ/dir/file.dat is the 
                           ! actual place. maindir is determined when
                           ! ascii file is split. categ is as above.
!   -------------------------
      integer itempx, itempv, nstr
      character*96  dirstr
      real*8 integral, xl
      integer klena
      integer from, fromwhich/0/   ! spcify integration direction. 
                   !  0--> from small.  !=0 --> from  large.
!      save normf, fromwhich, oldfmt
      save

      if( h%c%init .eq. 'initend') then
         write(0, *) '1D hist already instanciated'
         write(0, *) ' title=',h%c%title
         write(0, *) ' category=',h%c%categ
         write(0, *) ' id=',h%c%id
         stop 9999
      else
         h%c%init = 'initend'
      endif

      h%x%nhist = inbin 
      allocate( h%xw(inbin) )
      allocate( h%dnw(inbin) )
      allocate( h%mean(inbin) )
      allocate( h%dndx(inbin) )
      allocate( h%xforI(inbin) )
      allocate( h%integ(inbin) )


      h%x%tklg  =( itklg - (itklg/2)*2  ) .ne. 0
!     h%x%tklg  = bit(0, itklg)   ! etc may be used; logical bit needed
      h%x%cent  =( (itklg/2)*2 - (itklg/4)*4 ) /2    ! integer
      h%x%ufl  = ( (itklg/4)*4 - (itklg/8)*8 ) .ne. 0
      h%x%ofl  = ( (itklg/8)*8 - (itklg/16)*16 ) .ne. 0
      asmax = ( (itklg/16)*16 - (itklg/32)*32 ) .ne. 0
      

      h%x%xmin = ixmin    !  not used at present
      if(asmax) then
         if(ixmin .ge. ibinORxmax ) then
            write(0,*) ' ibinORxmax is regarded as ixmax but <= ixmin'
            stop 99999
         else
            if( h%x%cent .eq. 1 ) then
               ndiv= inbin - 2
            else
               ndiv = inbin-1
            endif
            if(h%x%tklg) then
               h%x%bin = log10(ibinORxmax/ixmin)/ndiv
            else
               h%x%bin = (ibinORxmax - ixmin )/ndiv
            endif
         endif
      else
         h%x%bin = ibinORxmax
      endif

      if( h%x%tklg  ) then
         if( h%x%xmin <= 0.0 )  then
            write(0, '("min must be > 0 for log option")')
            stop
         endif
         h%x%xm = log10(h%x%xmin) - h%x%cent * h%x%bin/2
         h%x%inc = 10.**h%x%bin
      else 
         h%x%xm = h%x%xmin  -  h%x%cent * h%x%bin/2
         h%x%inc  = h%x%bin
      endif

      h%c%id = ' '
!
      h%c%eventno = 1
      h%x%label =' '
      h%x%unit = ' '
      h%c%title = ' '
      h%c%dNunit=' '
      h%c%categ = ' '
      h%c%dir = ' '
      h%c%pw = 0
      h%c%logv = .true.
      h%c%norm = -1.0
      return
!    ************************
      entry kwhistc(h)
!    ************************

      do i = 1, h%x%nhist
         h%xw(i) = 0.
         h%dnw(i) = 0.
         h%mean(i) =0. 
      enddo
      return

      entry kwhistfmt(old)
      oldfmt=old
      return
!    *************************
      entry kwhist( h, x, w )
!    *************************
      if( h%x%tklg  .and. x .le. 0.) then
!         neglect this data
      else
         if( h%x%tklg  ) then
            xx = log10(x)
         else
            xx = x
         endif
         i = ( xx-h%x%xm ) / h%x%bin  + 1

         if(i .le. 0 .and. h%x%ufl ) then
            i = 1
         elseif(i .ge. h%x%nhist .and. h%x%ofl ) then
            i = h%x%nhist-1
         endif
         if(i .ge. 1 .and.  i  .lt. h%x%nhist )  then
            h%xw(i) = h%xw(i)  +  x*w
            h%dnw(i) = h%dnw(i) + w
            h%mean(i) = h%mean(i) + 1.0
         endif

      endif
      return
!     ******************
      entry kwhists0(from)
      fromwhich = from
      return
!     ***********************
      entry kwhists( h, inorm )
!     ************* take statistics
!         if inorm = -1.0, use alredy fixed one.
      if( inorm .ne. -1.0) then
         h%c%norm = inorm
      endif

      h%x%imin = 1
      do while( h%x%imin .lt. h%x%nhist .and.  h%dnw(h%x%imin) .eq. 0.) 
         h%x%imin = h%x%imin + 1
      enddo
      h%x%imax = h%x%nhist-1
      do while (h%x%imax .gt. 1 .and.  h%dnw(h%x%imax) .eq.  0.)  
         h%x%imax = h%x%imax -1
      enddo


      isumw = 0.
      do i = h%x%imin, h%x%imax
         isumw = isumw +  h%dnw(i)
      enddo
      h%x%sumw =  isumw

      if(h%c%norm .eq. 0. .and.  h%x%sumw  .gt. 0.) then
        normf = h%x%sumw 
      elseif(h%c%norm .le. 0. ) then
         normf = 1.0
      else
         normf = h%c%norm
      endif

!        bin center value and left boundary
      if( h%x%tklg ) then
         xx =10**(h%x%xm + h%x%bin/2.) * h%x%inc**(h%x%imin-1)
         xl =10**h%x%xm * h%x%inc**(h%x%imin-1)
      else
         xx = h%x%xm +   h%x%bin/2 + h%x%inc*(h%x%imin-1)
         xl = h%x%xm + h%x%inc*(h%x%imin-1)
      endif

         
      if(fromwhich .eq. 0 ) then
         integral = 0.
      else
         integral = isumw
      endif
      

      dx = h%x%bin      
      do i = h%x%imin, h%x%imax
         if( h%x%tklg ) then
            dx  = 10.0**(h%x%xm + i * h%x%bin) -  &
                 10.0**(h%x%xm + (i-1)*h%x%bin)
         endif
!         if( inorm .eq. -1.0) then
!            data is probbly read from file so mean exists
!cc            if( h%dnw(i) .eq. 0) then
!cc               h%xw(i) = xx
!cc            else
!cc               h%xw(i) = h%mean(i)* h%dnw(i)
!cc            endif
!cc         else
!c            if(h%dnw(i) .eq. 0) then
!c               h%mean(i) = xx
!c            else
!c               h%mean(i) = h%xw(i)/h%dnw(i)
!c            endif
!c         endif
!
         h%xforI(i) = xl
         h%integ(i) = integral/normf
         h%dndx(i) = h%dnw(i)/dx/normf
         if( h%x%tklg ) then
            xx = xx * h%x%inc
            xl = xl * h%x%inc
         else
            xx = xx + h%x%inc
            xl = xl + h%x%inc
         endif
         if(fromwhich .eq. 0) then
            integral = integral + h%dnw(i)
         else
            integral = integral - h%dnw(i)
         endif
      enddo
!      if(.not. oldfmt) then
         h%xforI(h%x%imax+1) =xl
         h%integ(h%x%imax+1) =integral/normf
!      endif
      h%dndx(h%x%imax+1) =0.
      h%dnw(h%x%imax+1) =0.
!c      h%mean(h%x%imax+1) =xx
      h%mean(h%x%imax+1) =0.
      return
!     ********************
      entry kwhistev(h, evno)
!     ********************
      h%c%eventno = evno
      return
!     ********************
      entry kwhistid(h,  id )
!     *******************
      h%c%id = id
      return
!     ********************
      entry kwhistai(h,  title, categ, dNunit, logv, pw, label, unit)
!     *******************
!       additional info.

      h%c%title = title
      h%c%categ =  categ
      h%c%dNunit = dNunit
      h%c%pw = pw
      h%c%logv = logv
      h%x%label = label
      h%x%unit =  unit
      return
!     *******************
      entry kwhistdir(h, dir)
!     ******************must be called after kwhistai is called
      dirstr = dir
      call kseblk(dirstr, "|", nstr)
      if( klena(h%c%categ) .gt. 0 )  then
         h%c%dir=h%c%categ(1:klena(h%c%categ))//"/" &
         //dirstr(1:nstr)
      else
         h%c%dir=dirstr(1:nstr)
      endif
      return
!
!     *********************
      entry kwhistpr( h, fno )
!     ****************print  hist

!       next lines for getting sumw is needed
!       sumv read from ascii file has some error
!       and results in 10^-5 diff. eventually
!
      isumw = 0.
      do i = h%x%imin, h%x%imax
         isumw =  isumw + h%dnw(i)
      enddo
      h%x%sumw =  isumw 
    
!        next must be given ; normf is undef. when
!        ascii file is read.
      if(h%c%norm .eq. 0. .and.  h%x%sumw  .gt. 0.) then
         normf = isumw
      elseif(h%c%norm .le. 0. ) then
         normf = 1.0
      else
         normf = h%c%norm
      endif
!      ......................
      if( h%x%tklg ) then
         xx = 10.0**(h%x%xm + h%x%bin/2.0) * h%x%inc**(h%x%imin-1)
      else
         xx = h%x%xm + h%x%bin/2. + h%x%inc*(h%x%imin-1)
      endif
!        header
      iufl = 0
      if(h%x%ufl) iufl = 1
      iofl = 0
      if(h%x%ofl) iofl = 1
 
      if(fno .lt. 0) then
         write(*, '(a, i5, i7, 3i3, g15.6)')  &
         '#hist1 ', h%c%eventno, h%x%nhist, &
         h%x%cent, iufl, iofl, h%x%bin 
         write(*, '(a,a)') '#t ', h%c%title(1:klena(h%c%title))
         write(*, '(a,a)') '#c ', h%c%categ(1:klena(h%c%categ))

         write(*, '(a, a,1x, a)')  '#x ', h%x%label, h%x%unit
         write(*, '(a,f10.2)') '#pw ', h%c%pw

         write(*, '(a,a)') '#dN ', h%c%dNunit(1:klena(h%c%dNunit))
         write(*, '(a,a)') '#k ', h%c%id(1:klena(h%c%id))
         itempx  = 0
         if( h%x%tklg )  itempx = 1
         itempv = 0
         if( h%c%logv) itempv = 1
         write(*,'(a, 2i3)')'#l ',  itempx, itempv
         write(*, '(a, 1pE15.8, g15.6)')'#n ', isumw, normf
         write(*, '(a,2i6,2g15.6)') &
        '#o ', h%x%imin, h%x%imax, &
               h%x%xm, h%x%inc
      else
         write(fno, '(a, i5, i7, 3i3, g15.6)')  &
          '#hist1 ', h%c%eventno, h%x%nhist,  &
          h%x%cent, iufl, iofl, h%x%bin 
         write(fno, '(a,a)') '#t ', h%c%title(1:klena(h%c%title))
         write(fno, '(a,a)') '#c ', h%c%categ(1:klena(h%c%categ))

         write(fno, '(a, a,1x, a)')  '#x ', h%x%label, h%x%unit
         write(fno, '(a,f10.2)') '#pw ', h%c%pw

         write(fno, '(a,a)') '#dN ', h%c%dNunit(1:klena(h%c%dNunit))
         write(fno, '(a,a)') '#k ', h%c%id(1:klena(h%c%id))
         itempx = 0
         if( h%x%tklg )  itempx = 1
         itempv = 0
         if( h%c%logv) itempv = 1
         write(fno,'(a,2i3)') '#l ',  itempx, itempv
         write(fno, '(a, 1pE15.8, g15.6)') '#n ', isumw, normf
         write(fno, '(a,2i6,2g15.6)') &
        '#o ', h%x%imin, h%x%imax, &
               h%x%xm,  h%x%inc
      endif

!          directory where data is saved
      dirstr= h%c%dir(1:klena(h%c%dir))
      call kseblk(dirstr, "|", nstr)
      if(fno .lt. 0 ) then
         write(*,'(a, a)') '#d ', dirstr(1:nstr)
      else
         write(fno,'(a, a)') '#d ', dirstr(1:nstr)
      endif

         

      do i = h%x%imin, h%x%imax
         if(fno .lt. 0) then
            write(*, '(i5, 1p6E13.5)')   i,  &
           xx,  h%dndx(i), h%dnw(i), h%mean(i),  h%xforI(i), h%integ(i)
         else
            write(fno, '(i5, 1p6E13.5)')  i,  &
           xx,  h%dndx(i), h%dnw(i), h%mean(i), &
           h%xforI(i), h%integ(i)
         endif            

         if( h%x%tklg ) then
            xx =  xx * h%x%inc 
         else
            xx = xx + h%x%inc
         endif
      enddo


      if(fno .lt. 0) then
         write(*, '(i5, 1p6E13.5)')    h%x%imax+1,  &
          xx,  0.,  0., 0.,  &
           h%xforI( h%x%imax+1), h%integ( h%x%imax+1)
      else
         write(fno, '(i5, 1p6E13.5)')  h%x%imax+1, &
             xx,  0.,  0., 0.,  &
           h%xforI( h%x%imax+1), h%integ( h%x%imax+1)
      endif            
      
!       trailer
      if(fno .lt. 0) then
         write(*,'(7i3)')  0,0,0,0,0,0, 0
      else
         write(fno,'(7i3)')  0,0,0,0,0,0, 0
      endif
      return
!     *********************
      entry kwhistw(h, bfnow)
!     ********************      
!       binary write of h to bfnow
      write(bfnow) '#hist1'
      write(bfnow) h%x%nhist
      write(bfnow) h%x, h%c
      write(bfnow) h%xw, h%dnw, h%mean, h%dndx, &
                  h%xforI, h%integ
      return
!     *********************
      entry kwhistr(h, bfnor, icon)
!     ********************
!        #hist1 must be read outside
!
      read(bfnor, end =222)  nbin 
!/////////////////////
!      write(0,*) ' nbin=',nbin
!///////////////////////
      allocate( h%xw(nbin) )
      allocate( h%dnw(nbin) )
      allocate( h%mean(nbin) )
      allocate( h%dndx(nbin) )
!      if(.not. oldfmt) then
         allocate( h%xforI(nbin) )
         allocate( h%integ(nbin) )
!      endif

      read(bfnor, end=222)  h%x, h%c
!////////////////////////
!      write(0,*)' h%x, h%c read'
!/////////////////
      if(oldfmt) then
         read(bfnor, end= 222)  h%xw, h%dnw, h%mean, h%dndx
      else
        read(bfnor, end= 222)  h%xw, h%dnw, h%mean, &
        h%dndx, h%xforI, h%integ
!///////////////////////
!        write(0,*) ' xw.. read'
!//////////////
      endif
      icon = 0
      return
 222  continue
      write(0,*) ' kwhistr reached EOF unexpectedly'
      icon = 1
      return
!     ****************
      entry kwhistd(h)
!     ***************
!      deallocate histogram area
!
      h%c%init = ' '
!      if(.not. oldfmt) then
         deallocate( h%xw, h%dnw,  h%mean,  h%dndx,  &
                 h%xforI, h%integ, stat=dealloc)
!      else
!         deallocate( h%xw, h%dnw,  h%mean,  h%dndx, 
!     *     stat=dealloc)
!      endif
      if(dealloc .ne. 0) then
         write(0,*) ' in kwhist1: dealloc failed =',dealloc
         stop 12345
      endif
      return
!     ********************
      entry kwhista(h1, h2, h)
!     ******************
!      h = h1 + h2  of bin area. For others, h1 is inherited
!      h,  h1 and h2 must be the same size  histogram  of same 
!     type.  h can be h1
!
      if( h1%x%nhist .ne. h2%x%nhist) then
         write(0, *) ' h1 and h2 diff. size histogram in kwhista'
         stop 9876
      endif
      if( h%c%init .ne. 'initend') then
!           not yet initialized.
         nbin = h1%x%nhist
         allocate( h%xw(nbin) )
         allocate( h%dnw(nbin) )
         allocate( h%mean(nbin) )
         allocate( h%dndx(nbin) )
!         if(.not. oldfmt) then
            allocate( h%xforI(nbin) )
            allocate( h%integ(nbin) )
!         endif
         h%c%init = 'initend'
      endif
      h%x = h1%x
      h%c = h1%c
      do i = 1, h%x%nhist
         h%xw(i) = h1%xw(i) + h2%xw(i)
         h%dnw(i) = h1%dnw(i) + h2%dnw(i)
      enddo
                                          
    end subroutine kwhisti
    subroutine kwhistso( binw )
        use modHistogram
!        specify output method
      implicit none
      integer binw  ! input.  1--> ascii write
                    !         2--> binary write
      BinWrite = binw
      if(binw .ne. 1 .and. binw .ne. 2) then
         write(0,*) 'binw=',binw,' for kwhistso is invalid'
         stop
      endif
    end subroutine kwhistso

    subroutine kwhistp( h, fno )
      use modHistogram
!
!         print or binary write histogram
!         kwhistso must be called to
!         fix binary write or print
!
      implicit none

      type(histogram1) h
      integer fno

      if( BinWrite .eq. 2 ) then
         call kwhistw(h, fno)
      elseif(BinWrite .eq. 1 ) then
         call kwhistpr(h, fno)
      endif
    end subroutine kwhistp
    
    integer function kwhistIxy(h, x, y, n)
      use modHistogram
      implicit none
      type(histogram1) h
!     ==================
      integer n  ! input size of x,y
      real*8 x(n), y(n) ! output, (normalized) 
                      ! integral value y at x (i=1, m)
      integer m ! function value. number of data points in x,y
                ! if number of x,y must be > n, m=-1 will be returned.
                ! (error)
      integer i

      m = h%x%imax - h%x%imin +2
      if(m .gt. n) then
         m = -1
      else
         m = 0
         do i = h%x%imin,  h%x%imax + 1
            m = m + 1
            x(m) = h%xforI(i)
            y(m) = h%integ(i)
         enddo
      endif
      kwhistIxy = m
    end function kwhistIxy


    integer function kwhistxy(h, x, y, n)
      use modHistogram
      implicit none

      type(histogram1) h
!     ==================
      integer n  ! input size of x,y
      real*8 x(n), y(n) ! output, normalzied dn/dx value y at x (i=1, m)
                      ! x is bin center.
      integer m ! function value. number of data points in x,y
                ! if number of x,y must be > n, m=-1 will be returned.
                ! (error)
      integer i
      real*8 xx

      if( h%x%tklg ) then
         xx = 10.0**(h%x%xm + h%x%bin/2.0) * h%x%inc**(h%x%imin-1)
      else
         xx = h%x%xm + h%x%bin/2. + h%x%inc*(h%x%imin-1)
      endif


      m = h%x%imax - h%x%imin + 1
      if(m .gt. n) then
         m = -1
      else
         m = 0
         do i = h%x%imin,  h%x%imax
            m = m + 1
            x(m) = xx
            y(m) = h%dndx(i)
            if( h%x%tklg ) then
               xx =  xx * h%x%inc 
            else
               xx = xx + h%x%inc
            endif
         enddo
      endif
      kwhistxy = m
    end function kwhistxy
  end module modHistogram1
