!     
!            weihted histograming  fortan 90 version. 2D case 
!            (Not work with Absoft) 
!      Usage:  kwhisti2:   instanciate one histogram
!              kwhistc2:   clear histogram area
!              kwhist2:    take histogram
!              kwhists2:   compute statistical result 
!              kwhistp2:   print statistical resutl
!
  !
module modHistogram2
  use modHistogram
  implicit none

  type histogram2
#if defined  (KEKA) || defined (KEKB)
     sequence
#endif
     type(histoc)  c
     type(histogram) x
     type(histogram) y
     real, allocatable ::  dnw(:,:)
     real, allocatable ::  dndxdy(:,:)
  end type histogram2

contains

subroutine kwhisti2( h, &
          ixmin, ixbinORxmax, ixnbin, ixtklg,  &
          iymin, iybinORymax, iynbin, iytklg )
!  use modHistogram
  implicit none
!         initialize 
  integer,intent(in):: ixnbin  !  request inbin histogram area
  real,intent(in):: ixmin     ! . xmin. not in log even if log10(variable) is taken
                     !         see itklg 
  real,intent(in):: ixbinORxmax  !  bin or ixmax. depends on itklg.
                     !  If bin and log10 is taken, bin is for log10 

  integer,intent(in):: ixtklg  !  bit pattern. give it like b'10001'
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



  integer,intent(in):: iynbin  !   request inbin histogram size
  real,intent(in):: iymin     !  xmin. not in log even itklg ==1
  real,intent(in):: iybinORymax   ! . bin.  If log, bin is for log10
  integer,intent(in):: iytklg    !  same as ixtklg.

  integer dealloc

  real sum
  logical asmax
  real  normf

  character*(*)  id
  integer  evno  !  input event no.
  character*(*)  title      !  input. graph title
  character*(*)  categ      !  input. category of the histogram. such as lat re
  real pw                   !  input. vertical scale is displayed by x**pw dN/dx
                                !  in default
  logical logv              !  input. if T, vertical scale is displayed in
                                !  log in default
  character*(*)  labelx      !  input. a few char string for x-axis
                                !   (must not include special char).
  character*(*)  unitx       !  input. a few char string for x-axis unit
  character*(*)  labely      !  input. same as
  character*(*)  unity       !  input.    above.
  character*(*)  dir         ! input
  character*(*)  dNunit      ! input
  integer  xstep   ! input
  integer nstr
  character*96 dirstr
  integer klena

  type(histogram2) h, h1, h2
!     ====================      
  real inorm   !  input. used in the normalization as dN/dxdy/inorm
                   ! if 0, area normalization is tried.  Here area is not
                   !  whole (x,y) region but, x is regarded as a fixed paramter
                   !  and area is computed  along y at a given x.
  integer fno  !  if < 0, standard output is used else fno is used for histogram output
                   !  fno must be opened by the user beforehand.
  integer bfnow  ! input. binary  write file no.
  integer bfnor  ! input. binary  read file no.
  integer icon   ! output. 0; binary read was successful
                     !  1;unexpected EOF
  integer itempy, itempv
  
  character*32 filen
  integer nchar


  integer*2 xnbin, ynbin
  real  x, y,  w 
  real  xx, yy
  integer  i, j, ndiv
  integer istep
  real dx,  dy

  if( h%c%init == 'initend') then
     write(0, *) '2D hist  already instanciated; '
     write(0, *) ' title=',h%c%title
     write(0, *) ' category=',h%c%categ
     write(0, *) ' id=',h%c%id
     stop 9999
  else
     h%c%init = 'initend'
  endif
      
  h%x%nhist = ixnbin
  h%y%nhist = iynbin

  allocate( h%dnw(ixnbin, iynbin) )
  allocate( h%dndxdy(ixnbin, iynbin) )
      
  h%x%tklg  = ( ixtklg - (ixtklg/2)*2 ) .ne. 0
  h%x%cent  = ( (ixtklg/2)*2 - (ixtklg/4)*4 )/2
  h%x%ufl  = ( (ixtklg/4)*4 - (ixtklg/8)*8) .ne. 0
  h%x%ofl  = ( (ixtklg/8)*8 - (ixtklg/16)*16 ) .ne. 0
                        
  h%x%xmin = ixmin
  asmax =  ( (ixtklg/16)*16 - (ixtklg/32)*32 ) .ne. 0
  
  if(asmax) then
     if(ixmin .ge. ixbinORxmax) then
        write(0,*) ' ixbinORxmax is regarded as ixmax but <= ixmin'
        stop 99999
     else
        if( h%x%cent .eq. 1 ) then
           ndiv= ixnbin - 1
        else
           ndiv = ixnbin
        endif
        if(h%x%tklg) then
           h%x%bin = log10(ixbinORxmax/ixmin)/ndiv
        else
           h%x%bin = (ixbinORxmax - ixmin )/ndiv
        endif
     endif
  else
     h%x%bin = ixbinORxmax
  endif


  h%y%tklg  = ( iytklg - (iytklg/2)*2 ) .ne. 0
  h%y%cent  = ( (iytklg/2)*2 - (iytklg/4)*4 )/2
  h%y%ufl  = ( (iytklg/4)*4 - (iytklg/8)*8 ) .ne. 0
  h%y%ofl  = ( (iytklg/8)*8 - (iytklg/16)*16 ) .ne. 0
  
  h%y%xmin = iymin

  asmax =  ( (iytklg/16)*16 - (iytklg/32)*32 ) .ne. 0
  if(asmax) then
     if(iymin .ge. iybinORymax) then
        write(0,*) ' iybinORymax is regarded as iymax but <= iymin'
        stop 99999
     else
        if( h%y%cent .eq. 1 ) then
           ndiv= iynbin - 1
        else
           ndiv = iynbin
        endif
        if(h%y%tklg) then
           h%y%bin = log10(iybinORymax/iymin)/ndiv
        else
           h%y%bin = (iybinORymax - iymin )/ndiv
        endif
     endif
  else
     h%y%bin = iybinORymax
  endif

!-------------------------------------

  if( h%x%tklg  ) then
     if( h%x%xmin <= 0.0 )  then
        write(0,  '("min must be > 0 for log option")')
        stop
     endif
     h%x%xm = log10(h%x%xmin)-h%x%cent * h%x%bin/2.
     h%x%inc = 10.**h%x%bin
  else
     h%x%xm = h%x%xmin - h%x%cent * h%x%bin/2.
     h%x%inc = h%x%bin
  endif

  if( h%y%tklg  ) then
     if( h%y%xmin <= 0.0 )  then
        write(0, '("min must be > 0 for log option")')
        stop
     endif
     h%y%xm = log10(h%y%xmin) - h%y%cent*h%y%bin/2.
     h%y%inc = 10.**h%y%bin
  else 
     h%y%xm = h%y%xmin - h%y%cent * h%y%bin/2.
     h%y%inc  = h%y%bin
  endif
  h%c%id = ' '
  h%c%eventno = 1
  h%x%label =' '
  h%x%unit = ' '
  h%y%label =' '
  h%y%unit = ' '
  h%c%title = ' '
  h%x%step = 1
  h%c%categ = ' '
  h%c%pw = 0
  h%c%logv = .true.

  return
!     ******************
  entry kwhistc2(h)
  !     ******************
  do i = 1, h%x%nhist
     do j = 1, h%y%nhist
        h%dnw(i,j) = 0.
     enddo
  enddo
  
  return
!    *************************
  entry kwhist2( h, x, y,  w )
!    *************************
  if( h%x%tklg  .and. x .le. 0.) then
!         neglect this data
  elseif( h%y%tklg  .and. y .le. 0.) then
!         neglect this data
  else
     if( h%x%tklg  ) then
        xx = log10(x)
     else
        xx = x
     endif
     if( h%y%tklg  ) then
        yy = log10(y)
     else
        yy = y
     endif
     
     i = ( xx-h%x%xm ) / h%x%bin  + 1
     j = ( yy-h%y%xm ) / h%y%bin  + 1

     if(  i .lt. 1 .and. h%x%ufl ) i = 1
     if(  i .gt. h%x%nhist .and. h%x%ofl ) i = h%x%nhist
     if(  j .lt. 1 .and. h%y%ufl ) j = 1
     if(  j .gt. h%y%nhist .and. h%y%ofl ) i = h%y%nhist

     if(  i .ge. 1 .and.  i  .le. h%x%nhist   .and. &
          j .ge. 1 .and.  j  .le. h%y%nhist )  then
        h%dnw(i,j) = h%dnw(i,j) + w
     endif
  endif
  return

!     ***********************
  entry kwhists2( h, inorm )
!     ************* take statistics
  if( inorm .ne. -1.0) then
     h%c%norm = inorm
  endif

  h%x%imin = h%x%nhist
  h%x%imax = 1

  h%y%imin = h%y%nhist
  h%y%imax = 1


  do i = 1, h%x%nhist
     do j = 1, h%y%nhist
        if( h%dnw(i,j) .ne. 0) then
!                 IBM xl fortran will detect error if int() below is 
!                 omitted.  
           h%x%imin = min(int(h%x%imin), i)
           h%y%imin = min(int(h%y%imin), j)
           h%x%imax = max(int(h%x%imax), i)
           h%y%imax = max(int(h%y%imax), j)
        endif
     enddo
  enddo

!      write(0,*) ' min,max=', h%x%imin,h%x%imax, h%y%imin,h%y%imax
  dy = h%y%bin
  do i = h%x%imin, h%x%imax
     call kwhistgetnorm2(h, i, dx, sum, normf)
     do j = h%y%imin, h%y%imax
        if(h%y%tklg ) then
           dy  = 10.0**(h%y%xm + j * h%y%bin) - &
                10.0**(h%y%xm + (j-1)*h%y%bin)
        endif
        h%dndxdy(i,j) = h%dnw(i,j)/dx/dy/normf
     enddo
  enddo
  return
!     *******************
  entry kwhistev2(h, evno)
!     ****************
  h%c%eventno = evno
  return
!     **********************
  entry kwhistid2(h, id)
!     ********************
  h%c%id = id
  return
!     ********************
  entry kwhistai2(h,  title, categ, dNunit, logv, pw,  &
         labelx, unitx,   labely, unity)
!     *******************
!      additional info

  h%c%title = title
  h%c%categ =  categ
  h%c%logv = logv
  h%c%pw = pw
  h%c%dNunit = dNunit
  h%x%label = labelx
  h%x%unit =  unitx
  h%y%label = labely
  h%y%unit =  unity
  return
!     *******************
  entry kwhistdir2(h, dir)
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
!     ********************
  entry kwhiststep2(h, xstep)
!     *******************
  if(xstep .gt. 0) then
     h%x%step = xstep 
  endif
  return
!     *********************
  entry kwhistpr2( h, fno )
!     ****************print  hist


  if(h%x%tklg ) then
     xx = 10.0**( h%x%xm + h%x%bin/2) * h%x%inc**(h%x%imin-1)
  else
     xx = h%x%xm + h%x%bin/2 + h%x%inc*(h%x%imin-1)
  endif

  do i = h%x%imin, h%x%imax, h%x%step
     call kwhistgetnorm2(h, i, dx, sum,  normf )
!         header
     dirstr =  h%c%dir(1:klena(h%c%dir))
     call kseblk(dirstr,"|", nstr)
     if(fno .lt. 0) then
        write(*,'(a, i3)') '#hist2 ', h%c%eventno
        write(*, '(a,a)')  '#t ', h%c%title(1:klena(h%c%title))
        write(*, '(a,a)')  '#c ', h%c%categ(1:klena(h%c%categ))
        write(*, '(a, a,1x, a)')  '#x ', h%x%label, h%x%unit
        write(*, '(a, a,1x, a)')  '#y ', h%y%label, h%y%unit
        write(*, '(a,f10.2)') '#pw ', h%c%pw
        write(*, '(a,a)') '#dN ', h%c%dNunit(1:klena(h%c%dNunit))
        write(*, '(a,a,1pE11.3)')  &
             '#k ', h%c%id(1:klena(h%c%id)), xx
        itempy  = 0
        if( h%y%tklg )  itempy = 1
        itempv = 0
        if( h%c%logv) itempv = 1
        write(*,'(a, 2i3)') '#l ',  itempy, itempv
        write(*, '(a, 1p3E11.3)') '#n ', sum, normf, dx
        write(*, '(a,a)')  "#d ",dirstr(1:nstr)
        write(filen, '(a,i2,a)') h%x%label, i, ".dat"
        call kseblk(filen, "|", nchar)
        write(*, '(a,a)') "#f ",filen(1:klena(filen))
     else
        write(fno,'(a,i3)')  '#hist2 ', h%c%eventno
        write(fno, '(a,a)') '#t ', h%c%title(1:klena(h%c%title))
        write(fno, '(a,a)') '#c ', h%c%categ(1:klena(h%c%categ))
        write(fno, '(a, a,1x, a)')  '#x ', h%x%label, h%x%unit
        write(fno, '(a, a,1x, a)')  '#y ', h%y%label, h%y%unit
        write(fno, '(a,f10.2)') '#pw ', h%c%pw
        write(fno, '(a,a)') '#dN ', h%c%dNunit(1:klena(h%c%dNunit))
        write(fno, '(a,a,1pE11.3)') &
             '#k ', h%c%id(1:klena(h%c%id)), xx
        itempy  = 0
        if( h%y%tklg )  itempy = 1
        itempv = 0
        if( h%c%logv) itempv = 1
        write(fno,'(a, 2i3)')'#l ',  itempy, itempv
        write(fno, '(a, 1p3E11.3)')'#n ', sum, normf, dx
        write(fno, '(a,a)')  "#d ", dirstr(1:nstr)
        write(filen, '(a,i2,a)') h%x%label, i, ".dat"
        call kseblk(filen, "|", nchar)
        write(fno, '(a,a)') "#f ",filen(1:klena(filen))
     endif

     if(h%y%tklg ) then
        yy =10.**( h%y%xm + h%y%bin/2) * h%y%inc**(h%y%imin-1)
     else
        yy = h%y%xm + h%y%bin/2 + h%y%inc*(h%y%imin-1)
     endif
     dy = h%y%bin
     do j = h%y%imin, h%y%imax
        if(h%y%tklg ) then
           dy  = 10.0**(h%y%xm + j * h%y%bin) - &
                 10.0**(h%y%xm + (j-1)*h%y%bin)
        endif

        if( fno .lt. 0 ) then
           write(*, '( 2i4,1p5E11.3)')  i, j, &
                xx, yy, h%dndxdy(i,j), h%dnw(i,j), dy
        else
           write(fno, '( 2i4, 1p5E11.3)') i, j, &
                xx, yy, h%dndxdy(i,j), h%dnw(i,j), dy
        endif
 
        if( h%y%tklg ) then
           yy =  yy * h%y%inc 
        else
           yy =  yy + h%y%inc
        endif
     enddo
     do istep = 1, h%x%step
        if( h%x%tklg ) then
           xx =  xx * h%x%inc 
        else
           xx = xx + h%x%inc
        endif
     enddo
     if(fno .lt. 0) then
        write(*,'(7i3)')  0, 0, 0, 0, 0, 0, 0
     else
        write(fno,'(7i3)') 0, 0, 0, 0, 0, 0, 0
     endif
  enddo
!        trailer
  if(fno .lt. 0) then
     write(*,'(7i3)')  0, 0, 0, 0, 0, 0, 0
  else
     write(fno,'(7i3)')  0, 0, 0, 0, 0, 0, 0
  endif
  return
!     *********************
  entry kwhistw2(h, bfnow)
!     *****************
!       binary write of h to bfnow
  write(bfnow) '#hist2'
  write(bfnow) h%x%nhist, h%y%nhist
  write(bfnow) h%x, h%y, h%c
  write(bfnow) h%dnw,  h%dndxdy
  return
!     *********************
  entry kwhistr2(h, bfnor, icon)
!     *********************
!        #hist2 must be read outside
      
  read(bfnor, end =222) xnbin, ynbin

  allocate( h%dnw(xnbin, ynbin) )
  allocate( h%dndxdy(xnbin, ynbin) )

  read(bfnor, end=222)  h%x, h%y,  h%c
  read(bfnor, end= 222) h%dnw, h%dndxdy
  icon = 0
  return
 222  continue
  write(0,*) ' kwhistr2 reached EOF unexpectedly'
  icon = 1
  return
!     ******************
  entry kwhistd2(h)
!     ****************
!        deallocate histo area
  h%c%init = ' '
  deallocate(h%dnw, h%dndxdy, stat=dealloc)
  if(dealloc .ne. 0 ) then
     write(0,*) ' failed deallocation=' , dealloc
     stop 9999
  endif
  return
!     ********************
  entry kwhista2(h1, h2, h)
!     ******************
!      h = h1 + h2  of bin area. For others, h1 is inherited
!      h,  h1 and h2 must have the identical structure
!      h can be h1
!
  if( h1%x%nhist .ne. h2%x%nhist .or. &
       h1%y%nhist .ne. h2%y%nhist ) then
     write(0, *) &
         ' h1 and h2 diff. size histogram in kwhista2'
     write(0,*) " h1x, h2x=", h1%x%nhist, h2%x%nhist 
     write(0,*) " h1y, h2y=", h1%y%nhist, h2%y%nhist 
     stop 9876
  endif
  if( h%c%init .ne. 'initend') then
!           not yet initialized.
     xnbin = h1%x%nhist
     ynbin = h1%y%nhist
     allocate( h%dnw(xnbin, ynbin) )
     allocate( h%dndxdy(xnbin, ynbin) )
     h%c%init = 'initend'
  endif
  h%x = h1%x
  h%y = h1%y
  h%c = h1%c

  do i = 1, h%x%nhist
     do j = 1, h%y%nhist
        h%dnw(i,j) = h1%dnw(i,j) + h2%dnw(i,j)
     enddo
  enddo
end subroutine kwhisti2
end module modHistogram2

subroutine kwhistp2( h, fno )
!  use modHistogram
  use modHistogram2

  implicit none

  type(histogram2) h
  integer fno
      
  if( BinWrite .eq. 1) then
     call kwhistpr2( h, fno )
  else
     call kwhistw2( h, fno )
  endif
end subroutine kwhistp2


subroutine kwhistgetnorm2(h, i, dx,  sum, normf)
!  use modHistogram
  use modHistogram2
  implicit none

  type(histogram2) h
!     ====================      
  integer i   ! input. index of x for histogram h(x,y)
  real normf   ! output. normalizaton factor to be used.
  real dx      !  //     bin of x  
  real  sum    !  //     sum to be used for area normalization
  integer j

  dx = h%x%bin      

  if(h%x%tklg ) then
     dx  = 10.0**(h%x%xm + i * h%x%bin) - &
          10.0**(h%x%xm + (i-1)*h%x%bin)
  endif


  sum = 0.
  do j = h%y%imin, h%y%imax
     sum = sum +  h%dnw(i,j)
  enddo
  if( sum .eq. 0. ) then
     normf= 1.0
  else 
     sum = sum / dx
  endif
  if(h%c%norm .eq. 0.) then
     normf = sum
  elseif(h%c%norm .lt. 0.) then
!          not initilized. strange
     normf = 1.
  else
     normf = h%c%norm
  endif
end subroutine kwhistgetnorm2


