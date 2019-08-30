integer function  kwhistReadAscii(h, fno)
!  use modHistogram
  use modHistogram1
!           this is only for 1D histogram
!        data is read 1 histogram by 1 in the file
!        specified by fno.  That is, when one histogram
!        is read, return is made so that the user
!        must use this one again if more histogram is
!        contained in the same file.  
!        instanciation of h is done in this prog. 
!        deallocaion must be done by the user
  implicit none
!     ******************        
  type(histogram1) h
  integer fno
  integer icon
!      -1: requested data not exist
!      >=0: some data read. the number of bins
!      with non zero data is given
  integer MAXLENG
  parameter (MAXLENG=121)

  character*(MAXLENG) string
  integer nch, ndim

  integer MAXL
  parameter( MAXL = MAXLENG-1 )

  parameter (nch=MAXL, ndim=15)

  character*(nch)  head


  integer nbin, i, iofl, iufl
  real xx, xx2, xx3
  real tt, tt1, tt2
  character*(nch) term(ndim)
  integer nr
  integer klena
  save


  icon = -1
  do while (.true.)
     string=" "
     read(fno, '(a)', end=100) string
     do i = 1, ndim
        term(i) = " "
     enddo
     call ksplit(string, nch, ndim, term, nr)
     head=" "
     head=term(1)
!
!         write(0,*) ' head=', head
!
     if( head .eq. "#hist1" ) then
        icon = 0
        string(1:6)="      "
        
        read(string, *)  h%c%eventno, h%x%nhist, &
             h%x%cent, iufl, iofl,   h%x%bin
        
        h%x%ufl = iufl == 1
        h%x%ofl = iofl == 1
!
!           write(0,*) ' header read;  #of bin=',h%x%nhist
!           write(0,*) ' init cond=', h%c%init

        if( h%c%init .ne. "initend") then
           nbin = h%x%nhist
           allocate( h%xw(nbin) )
           allocate( h%dnw(nbin) )
           allocate( h%mean(nbin) )
           allocate( h%dndx(nbin) )
           allocate( h%xforI(nbin) )
           allocate( h%integ(nbin) )
           !               write(0,*) ' histo region obained'
!               write(0, *) ' ev#=', h%c%eventno,  ' cent=',
!     *          h%x%cent, ' uf=',  h%x%ufl, ' of=',h%x%ofl,
!     *          ' bin=', h%x%bin
!               write(0,*) ' going to clear ', h%x%nhist,' bins'
           call kwhistc( h )

           h%c%init = "initend"
!
!               write(0,*) ' hist init end'
! 
        endif
     elseif( head .eq. "#t") then
        h%c%title=string(klena(head)+2:klena(string))
     elseif( head .eq. "#c" )   then
        h%c%categ = string(klena(head)+2:klena(string))
     elseif( head .eq. "#x")  then
        h%x%label=term(2)
        h%x%unit =term(3)
     elseif( head .eq. "#pw")  then
        read(term(2), *)   h%c%pw
     elseif( head .eq. "#dN" ) then
        h%c%dNunit =string(klena(head)+2:klena(string))
     elseif( head .eq. "#k") then
        h%c%id = string(klena(head)+2:klena(string))
     elseif( head .eq. "#l" ) then
        read(term(2), *) i
        h%x%tklg = i .eq. 1
        read(term(3), *) i
        h%c%logv = i .eq. 1
     elseif( head .eq. "#n") then
        read(term(2), *) h%x%sumw
        read(term(3), *) h%c%norm
     elseif( head .eq. "#o" ) then
        read(term(2), *) h%x%imin
        read(term(3), *) h%x%imax
        read(term(4), *) h%x%xm
        read(term(5), *) h%x%inc
     elseif( head .eq. "#d") then
        h%c%dir =  string(klena(head)+2:klena(string))
     elseif( head .eq. "0" ) then
        exit   ! **************************
     else
        read(string,*) i, xx, tt, tt1, tt2, xx2,  xx3
!               correct data contains max + 1 
!               but need not store last one
        if(i .le. h%x%imax+1) then
           h%dndx(i) = tt
           h%dnw(i) = tt1
           h%mean(i) = tt2
           h%xforI(i) = xx2
           h%integ(i) = xx3
           icon= icon +1
        endif
     endif
  enddo
100 continue
  kwhistReadAscii=icon
end function kwhistReadAscii

