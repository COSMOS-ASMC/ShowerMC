!     
!            weihted histograming  fortan 90 version. 3D case 
!            (Not work with Absoft) 
!
  !
module modHistogram3
  use modHistogram
  implicit none
  type histogram3
#if defined  (KEKA) || defined (KEKB)
     sequence
#endif
     type(histoc) c
     type(histogram) x
     type(histogram) y
     type(histogram) z
     real, allocatable ::  dnw(:,:,:)
     real, allocatable ::  dndxdydz(:,:,:)
  end type histogram3

contains
  
subroutine kwhisti3( h, &
     ixmin, ixbinORxmax, ixnbin, ixtklg,   &
     iymin, iybinORymax, iynbin, iytklg,   &
     izmin, izbinORzmax, iznbin, iztklg )

  use modHistogram
  
  implicit none
!         initialize 
  integer ixnbin  ! input.  request ixnbin histogram size
  real ixmin     ! input. xmin. not in log even log10
  real ixbinORxmax    ! input. bin or xmax.  see k90whist2.f
  integer ixtklg    ! input. same as k90wist or k90whist2

  integer iynbin  ! input.  request iynbin histogram size
  real iymin     ! input. xmin. not in log even log10
  real iybinORymax      ! input. bin.  If log, bin is for log10
  integer iytklg    ! input.  same as ixtklg 
                       !         
  integer iznbin  ! input.  request iznbin histogram size
  real izmin     ! input. xmin. not in log even log10
  real izbinORzmax      ! input. bin.  If log, bin is for log10
  integer iztklg    ! input.   same as ixtklg 
  integer xstep, ystep ! input

  integer dealloc
  real sum 
  character*64 item(4)
  integer nitem
  integer nc 
  real normf
  integer*2   xnbin, ynbin, znbin
                       !         
  character*(*)  id
  character*(*)  title      !  input. graph title
  character*(*)  categ      !  input. max 8 char. category of histogram suc as lateral
  real pw                   !  input. vertical scale is displayed by x**pw dN/dx
                                !  in default
  logical logv              !  input. if T, vertical scale is displayed in
                                !  log in default
  character*(*)  labelx     !  input. a few char string for x-axis
                                !   (must not include special char).
  character*(*)  unitx      !  input. a few char string for x-axis unit
  character*(*)  labely     !  input. same as
  character*(*)  unity      !  input.    above.
  character*(*)  labelz     !  input. same as
  character*(*)  unitz      !  input.    above.
  character*(*)  dir        ! input.  
  character*(*)  dNunit ! input.

  integer klena
!     ******************

  type(histogram3) h, h1, h2
!     ====================      
  real inorm   !  input. used in the normalization as dN/dx/inorm
                   ! if 0, area normalization is tried. Hear the area is
                   ! not whole (x,y,z) region, but along  z at given (x,y)
                   ! ((x,y) are regarded as fixed parameter.
  integer fno  !  if < 0, standard output is used else fno is used for histogram output
                   !  fno must be opened by the user beforehand.
  integer bfnow  ! input. binary  write file no.
  integer bfnor  ! input. binary  read file no.
  integer icon   ! output. 0; binary read was successful
                     !  1;unexpected EOF
  integer evno   ! input
  !
  integer  itempz, itempv
  integer istep
  character*96 dirstr
  character*32 filen
  integer nchar
  real  x, y, z,  w 
  real  xx, yy, zz
  integer  i, j, k, ndiv
  logical asmax
  real dx,  dy, dz



  if( h%c%init .eq. 'initend') then
     write(0, *) '3D hist.  already instanciated;'
     write(0, *) ' title=',h%c%title
     write(0, *) ' category=',h%c%categ
     write(0, *) ' id=',h%c%id
     stop 9999
  else
     h%c%init = 'initend'
  endif
                                                                              
  h%x%nhist = ixnbin
  h%y%nhist = iynbin
  h%z%nhist = iznbin
  h%x%step = 1
  h%y%step = 1
  allocate( h%dnw(ixnbin, iynbin, iznbin) )
  allocate( h%dndxdydz(ixnbin, iynbin, iznbin) )

  h%x%tklg  = ( ixtklg - (ixtklg/2)*2 ) .ne. 0
  h%x%cent  = ( ( ixtklg/2)*2 - (ixtklg/4)*4 )/2
  h%x%ufl  = ( (ixtklg/4 )*4 - (ixtklg/8)*8) .ne. 0
  h%x%ofl  = ( (ixtklg/8 )*8 - (ixtklg/16)*16 ) .ne. 0
                        
  h%x%xmin = ixmin

  asmax  = ( (ixtklg/16 )*16 - (ixtklg/32)*32 ) .ne. 0
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


  h%y%tklg  = (iytklg - (iytklg/2)*2) .ne. 0
  h%y%cent  = ( (iytklg/2)*2 - (iytklg/4)*4) /2
  h%y%ufl  = ( (iytklg/4)*4 - (iytklg/8)*8 ) .ne. 0
  h%y%ofl  = ( (iytklg/8)*8 - (iytklg/16)*16 ) .ne. 0

  h%y%xmin = iymin
  
  asmax  = ( (iytklg/16 )*16 - (iytklg/32)*32 ) .ne. 0
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



  h%z%tklg  =( iztklg - (iztklg/2)*2 ) .ne. 0
  h%z%cent  =( (iztklg/2)*2 - (iztklg/4)*4 ) /2
  h%z%ufl  = ( (iztklg/4)*4 - (iztklg/8)*8 ) .ne. 0
  h%z%ofl  = ( (iztklg/8)*8 - (iztklg/16)*16 ) .ne. 0
  
  h%z%xmin = izmin

  asmax  = ( (iztklg/16 )*16 - (iztklg/32)*32 ) .ne. 0
  if(asmax) then
     if(izmin .ge. izbinORzmax) then
        write(0,*) ' izbinORzmax is regarded as izmax but <= izmin'
        stop 99999
     else
        if( h%z%cent .eq. 1 ) then
           ndiv= iznbin - 1
        else
           ndiv = iznbin
        endif
        if(h%z%tklg) then
           h%z%bin = log10(izbinORzmax/izmin)/ndiv
        else
           h%z%bin = (izbinORzmax - izmin )/ndiv
        endif
     endif
  else
     h%z%bin = izbinORzmax
  endif


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
     h%y%xm = log10(h%y%xmin) - h%y%cent * h%y%bin/2.
     h%y%inc = 10.**h%y%bin
  else 
     h%y%xm = h%y%xmin - h%y%cent * h%y%bin/2.
     h%y%inc  = h%y%bin
  endif

  if( h%z%tklg ) then
     if( h%z%xmin <= 0.0 )  then
        write(0, '("min must be > 0 for log option")')
        stop
     endif
     h%z%xm = log10(h%z%xmin) - h%z%cent * h%z%bin/2.
     h%z%inc = 10.**h%z%bin
  else 
     h%z%xm = h%z%xmin - h%z%cent * h%z%bin/2.
     h%z%inc  = h%z%bin
  endif

  h%c%eventno = 1
  h%x%label =' '
  h%x%unit = ' '
  h%y%label =' '
  h%y%unit = ' '
  h%z%label =' '
  h%z%unit = ' '
  h%c%title = ' '
  h%c%categ = ' '
  h%c%pw = 0
  h%c%logv = .true.
                                          
  return

!     ***************************
  entry kwhistc3(h)
!     ***************************
  do i = 1, h%x%nhist
     do j = 1, h%y%nhist
        do k = 1, h%z%nhist
           h%dnw(i,j,k) = 0.
        enddo
     enddo
  enddo
  return

!    *************************
  entry kwhist3( h, x, y, z,  w )
!    *************************
  if( h%x%tklg  .and. x .le. 0.) then
!         neglect this data
  elseif( h%y%tklg  .and. y .le. 0.) then
!         neglect this data
  elseif( h%z%tklg  .and. z .le. 0.) then
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
     if( h%z%tklg  ) then
        zz = log10(z)
     else
        zz = z
     endif

     i = ( xx-h%x%xm ) / h%x%bin  + 1
     if( i .le.  0 .and. h%x%ufl ) i = 1
     if( i .gt.  h%x%nhist  .and. h%x%ofl ) i = h%x%nhist
     
     j = ( yy-h%y%xm ) / h%y%bin  + 1
     if( j .le.  0 .and. h%y%ufl ) j = 1
     if( j .gt.  h%y%nhist  .and. h%y%ofl ) j = h%y%nhist

     k = ( zz-h%z%xm ) / h%z%bin  + 1
     if( k .le.  0 .and. h%z%ufl ) k = 1
     if( k .gt.  h%z%nhist  .and. h%z%ofl ) k = h%z%nhist

     if(  i .ge. 1 .and.  i  .le. h%x%nhist .and. &
          j .ge. 1 .and.  j  .le. h%y%nhist   .and. &          
          k .ge. 1 .and.  k  .le. h%z%nhist )  then
!            h%dn(i,j, k) = h%dn(i,j, k)  + 1
        h%dnw(i,j,k) = h%dnw(i,j,k) + w
     endif
  endif
  return

!     ***********************
  entry kwhists3( h, inorm )
!     ************* take statistics
  if( inorm .ne. -1.0) then
     h%c%norm = inorm
  endif

  h%x%imin = h%x%nhist
  h%x%imax = 1

  h%y%imin = h%y%nhist
  h%y%imax = 1

  h%z%imin = h%z%nhist
  h%z%imax = 1

  do i = 1, h%x%nhist
     do j = 1, h%y%nhist
        do k = 1, h%z%nhist
           if( h%dnw(i,j,k) .ne. 0) then
!                    IBM XL fortran complains without int below
              h%x%imin = min(int(h%x%imin), i)
              h%y%imin = min(int(h%y%imin), j)
              h%z%imin = min(int(h%z%imin), k)
              h%x%imax = max(int(h%x%imax), i)
              h%y%imax = max(int(h%y%imax), j)
              h%z%imax = max(int(h%z%imax), k)
           endif
        enddo
     enddo
  enddo

  dz = h%z%bin
  do i = h%x%imin, h%x%imax
     do j = h%y%imin, h%y%imax
        call kwhistgetnorm3(h, i, j, dx, dy, sum,  normf)
        do k = h%z%imin, h%z%imax
           if(h%z%tklg ) then
              dz  = 10.0**(h%z%xm + k * h%z%bin) - &
                                10.0**(h%z%xm + (k-1)*h%z%bin)
           endif
           h%dndxdydz(i,j,k) = h%dnw(i,j,k)/dx/dy/dz/normf
        enddo
     enddo
  enddo
  return
!     ****************
  entry kwhistev3(h, evno)
!     ********************
  h%c%eventno = evno
  return
!                        
!     *********************
  entry kwhistid3( h, id)
  !     ********************
  h%c%id = id
  return
!     ********************
  entry kwhistai3(h,  title, categ, dNunit, &
       logv, pw, labelx, unitx, &
       labely, unity, labelz, unitz )
!     *******************
!       additional info.

  h%c%title = title
  h%c%categ =  categ
  h%c%dNunit = dNunit
  h%x%label = labelx
  h%x%unit =  unitx
  h%y%label = labely
  h%y%unit =  unity
  h%z%label = labelz
  h%z%unit =  unitz
  h%c%pw = pw
  h%c%logv = logv
  return
!     *******************
  entry kwhistdir3(h, dir)
!     ******************must be called after kwhistai is called
  if( klena(h%c%categ) .gt. 0 )  then
     h%c%dir=h%c%categ(1:klena(h%c%categ))//"/" &
         //dir(1:klena(dir))
  else
     h%c%dir=dir(1:klena(dir))
  endif
  return
!     ********************
  entry kwhiststep3(h, xstep, ystep)
!     ******************
!        give step size of indexes for parameter X and Y
!
  if(xstep .gt. 0) then 
     h%x%step = xstep
  endif
  if(ystep .gt. 0) then 
     h%y%step = ystep
  endif
  return
                                                            
!     *********************
  entry kwhistpr3( h, fno )
!     ****************print  hist
!
!
  if(h%x%tklg ) then
     xx = 10.0**(h%x%xm + h%x%bin/2) * h%x%inc**(h%x%imin-1)
  else
     xx = h%x%xm + h%x%bin/2  + h%x%inc*(h%x%imin-1)
  endif
  do i = h%x%imin, h%x%imax, h%x%step
     if(h%y%tklg ) then
        yy = 10.**(h%y%xm + h%y%bin/2) * h%y%inc**(h%y%imin-1)
     else
        yy = h%y%xm + h%y%bin/2 + h%y%inc*(h%y%imin-1)
     endif
     do j = h%y%imin, h%y%imax, h%y%step
        if(h%z%tklg ) then
           zz =10.0**( h%z%xm + h%z%bin/2) * &
                h%z%inc**(h%z%imin-1)
        else
           zz = h%z%xmin + h%z%bin/2 + h%z%inc*(h%z%imin-1)
        endif
!              header
        call kwhistgetnorm3(h, i, j, dx, dy, sum, normf)
        if( fno .lt. 0 ) then
           write(*, '(a,i3)') '#hist3 ', h%c%eventno
           write(*, '(a,a)') '#t ', h%c%title(1:klena(h%c%title))
           write(*, '(a,a)') '#c ', h%c%categ(1:klena(h%c%categ))
           write(*, '(a, a,1x, a)')  '#x ', h%x%label, h%x%unit
           write(*, '(a, a,1x, a)')  '#y ', h%y%label, h%y%unit
           write(*, '(a, a,1x, a)')  '#z ', h%z%label, h%z%unit
           write(*, '(a,f10.2)') '#pw ', h%c%pw
           write(*, '(a,a)') '#dN ', &
                 h%c%dNunit(1:klena(h%c%dNunit))
           write(*, '(a,a,1p2E11.3)') '#k ', &
                  h%c%id(1:klena(h%c%id)), xx, yy
           itempz  = 0
           if( h%z%tklg )  itempz = 1
           itempv = 0
           if( h%c%logv) itempv = 1
           write(*,'(a,2i3)')'#l ',  itempz, itempv
           write(*,'(a, 1p4E11.3)') "#n ", sum,  normf, dx, dy
        else
           write(fno, '(a,i3)') '#hist3 ',  h%c%eventno
           write(fno, '(a,a)') '#t ', h%c%title(1:klena(h%c%title))
           write(fno, '(a,a)') '#c ', h%c%categ(1:klena(h%c%categ))
           write(fno, '(a, a,1x, a)')  '#x ', h%x%label, h%x%unit
           write(fno, '(a, a,1x, a)')  '#y ', h%y%label, h%y%unit
           write(fno, '(a, a,1x, a)')  '#z ', h%z%label, h%z%unit
           write(fno, '(a,f10.2)') '#pw ', h%c%pw
           write(fno, '(a,a)') '#dN ', &
                 h%c%dNunit(1:klena(h%c%dNunit))
           write(fno, '(a,a,1p2E11.3)') '#k ',  &
                h%c%id(1:klena(h%c%id)), xx, yy
           itempz  = 0
           if( h%z%tklg )  itempz = 1
           itempv = 0
           if( h%c%logv) itempv = 1
           write(fno,'(a,2i3)')'#l ',  itempz, itempv
           write(fno,'(a, 1p4E11.3)') "#n ", sum,  normf, dx, dy
        endif

        if(fno .lt. 0 ) then
           write(dirstr,'(a,i2)') &
                h%c%dir(1:klena(h%c%dir))//h%x%label,i
           call kseblk(dirstr,"|", nchar)
           write(*,'(a,a)') "#d ",dirstr(1:klena(dirstr))

           write(filen, '(a, i3, a)') h%y%label, j, ".dat"
           call kseblk(filen, "|", nchar)
           write(*,'(a, a)')   '#f ',  filen(1:klena(filen))
        else
           write(dirstr,'(a,i2)') &
              h%c%dir(1:klena(h%c%dir))//h%x%label,i
           call kseblk(dirstr,"|", nchar)
           write(fno,'(a,a)') "#d ",dirstr(1:klena(dirstr))

           write(filen, '(a, i3, a)') h%y%label, j, ".dat"
           call kseblk(filen, "|", nchar)
           write(fno,'(a, a)')   '#f ',  filen(1:klena(filen))
        endif

        dz = h%z%bin
        do k = h%z%imin, h%z%imax
           if(h%z%tklg ) then
              dz  = 10.0**(h%z%xm + k * h%z%bin) - &
                   10.0**(h%z%xm + (k-1)*h%z%bin)
           endif
           if(fno .lt. 0 ) then
              write(*, '(2i4,i5, 1p6E11.3)') &
                   i, j, k, &
                   xx, yy, zz,  h%dndxdydz(i,j,k), h%dnw(i,j,k), dz
           else
              write(fno, '(2i4, i5, 1p6E11.3)') &
                 i, j, k, &
                 xx, yy, zz,  h%dndxdydz(i,j,k), h%dnw(i,j,k), dz
           endif
           if( h%z%tklg ) then
              zz =  zz * h%z%inc 
           else
              zz =  zz + h%z%inc
           endif
        enddo
        if(fno .lt. 0) then
           write(*, '(9i3)')  0, 0, 0, 0, 0, 0, 0, 0, 0
        else
           write(fno, '(9i3)')  0, 0, 0, 0, 0, 0, 0, 0, 0  
        endif
        do istep = 1, h%y%step
           if( h%y%tklg ) then
              yy =  yy * h%y%inc 
           else
              yy =  yy + h%y%inc
           endif
        enddo
     enddo

     if( fno .lt. 0 ) then
        write(*, '(9i3)')  0,  0, 0, 0, 0, 0, 0, 0, 0
     else
        write(fno,'(9i3)')  0,  0, 0, 0, 0, 0, 0, 0, 0 
     endif
     do istep = 1, h%x%step
        if( h%x%tklg ) then
           xx =  xx * h%x%inc 
        else
           xx = xx + h%x%inc
        endif
     enddo
  enddo
  if( fno .lt. 0 ) then
     write(*, '(9i3)') 0,  0, 0, 0, 0, 0, 0, 0, 0 
  else
     write(fno,'(9i3)') 0,  0, 0, 0, 0, 0, 0, 0, 0 
  endif
  return
!     *********************
  entry kwhistw3(h, bfnow)
!     *****************
!       binary write of h to bfnow
  write(bfnow) '#hist3'
  write(bfnow) h%x%nhist, h%y%nhist, h%z%nhist 
  write(bfnow) h%x, h%y, h%z, h%c
  write(bfnow) h%dnw,  h%dndxdydz
  return

!     *********************
  entry kwhistr3(h, bfnor, icon)
!     *********************
!        #hist3 must be read outside

  read(bfnor, end =222) xnbin, ynbin, znbin
  allocate( h%dnw(xnbin, ynbin, znbin) )
  allocate( h%dndxdydz(xnbin, ynbin, znbin) )
  read(bfnor, end=222)  h%x, h%y, h%z, h%c
  read(bfnor, end= 222) h%dnw, h%dndxdydz
  icon = 0
  return
222 continue
  write(0,*) ' kwhistr3 reached EOF unexpectedly'
  icon = 1
  return
!     *****************
  entry kwhistd3(h)
!     ***************
!       deallocate histo area
  h%c%init = ' '
  deallocate(h%dnw,h%dndxdydz, stat=dealloc)
  if(dealloc .ne. 0) then
     write(0,*) ' failed to deallocated hist 3=',dealloc
     stop 3333
  endif
  return
!     ********************
  entry kwhista3(h1, h2, h)
!     ******************
!      h = h1 + h2  of bin area. For others, h1 is inherited
!      h,  h1 and h2 must have the identical structure
!      h can be h1
!
  if( h1%x%nhist .ne. h2%x%nhist .or. &
       h1%y%nhist .ne. h2%y%nhist .or. &
       h1%z%nhist .ne. h2%z%nhist ) then
     write(0, *) &
          ' h1 and h2 diff. size histogram in kwhista3'
     stop 9876
  endif
  if( h%c%init .ne. 'initend') then
!           not yet initialized.
     xnbin = h1%x%nhist
     ynbin = h1%y%nhist
     znbin = h1%z%nhist
     allocate( h%dnw(xnbin, ynbin, znbin) )
     allocate( h%dndxdydz(xnbin, ynbin, znbin) )
     h%c%init = 'initend'
  endif
  h%x = h1%x
  h%y = h1%y
  h%z = h1%z
  h%c = h1%c
  do i = 1, h%x%nhist
     do j = 1, h%y%nhist
        do k = 1, h%z%nhist
           h%dnw(i,j,k) = h1%dnw(i,j,k) + h2%dnw(i,j,k)
        enddo
     enddo
  enddo

end subroutine kwhisti3
end module modHistogram3

subroutine kwhistp3( h, fno )
!  use modHistogram
  use modHistogram3

  implicit none

  type(histogram3) h

  integer fno
      
  if( BinWrite .eq. 1) then
     call kwhistpr3( h, fno )
  else
     call kwhistw3( h, fno )
  endif
end subroutine kwhistp3

subroutine kwhistgetnorm3(h, i, j,  dx, dy,  sum, normf)
!  use modHistoram
  use modHistogram3

  implicit none
  type(histogram3) h
!     ====================
  integer i, j                 ! input. index of x for histogram h(x,y)
  real normf                ! output. normalizaton factor
  real dx, dy               !  //     bin of x
  
  integer k
  real  sum

  dx = h%x%bin
  dy = h%y%bin
  if(h%x%tklg ) then
     dx  = 10.0**(h%x%xm + i * h%x%bin) -  &
          10.0**(h%x%xm + (i-1)*h%x%bin)
  endif
         
      
  if(h%y%tklg ) then
     dy  = 10.0**(h%y%xm + j * h%y%bin) - &
             10.0**(h%y%xm + (j-1)*h%y%bin)
  endif
                                                                                                                            

  sum = 0.
  do k = h%z%imin, h%z%imax
     sum = sum + h%dnw(i,j,k)
  enddo
  if( sum .eq. 0. ) then
     sum = 1.
  else 
     sum = sum/dx/dy
  endif

  if(h%c%norm .eq. 0.) then
     normf = sum
  elseif(h%c%norm .lt. 0.) then
     normf = 1.0  ! strange. not initialize
  else
     normf = h%c%norm
  endif
end subroutine kwhistgetnorm3


                                                            
                                                
