      module nrl_atmos
        implicit none
        integer,save::nsize, nsizem
        logical,save::first=.true.

        real(4),save::glat, glong
        integer,save:: gday1,  gday2
        integer,save:: ghour1, ghour2

         !    d(drho/dh)/dh
        real(8)::drhodh2
        real(4),save::heightx=499d3
         !    above h> hightx rho=a exp(-x/b) is used
         !    the values are typical ones and will be changed  by actual data.
        real(4),save::thickx 
        real(4),save::a=3.6399d-09
        real(4),save::b=63000.0
        real(4),save::al, bl  ! same for h< h(1)
        real(8),parameter:: heightSpace=1000.d3 ! abvoe this, space




        real(4),allocatable,save::h(:), den(:), thick(:) !!, denp(:)
        real(4),allocatable,save::temp(:)  ! temperature  K
        real(4),allocatable,save::h2denCoef(:,:)
        real(4),allocatable,save::h2thickCoef(:,:)
        real(4),allocatable,save::h2tempCoef(:,:)
        real(4),allocatable,save::thick2hCoef(:,:)
        real(4),allocatable,save::thick2denCoef(:,:)
!!        real(4),allocatable,save::h2denpCoef(:,:)
      end module nrl_atmos
      
      function cvh2den(vh) result(ans)
      use nrl_atmos
      implicit none
      real(8),intent(in)::vh    !  m
      real(8)::ans
      
      real(4):: svh, sans
      
      svh = vh
      if( vh < heightx ) then
         if( vh >= h(1) ) then
            call kScsplIntp(h, den, nsize, h2denCoef, nsizem, svh, sans)
            ans = sans
         else
            ans = al*exp(-vh/bl)
         endif
      else
         ans = a*exp(-vh/b)
      endif

      end function cvh2den

      function cvh2temp(vh) result(ans)
      use nrl_atmos
      implicit none
      real(8),intent(in):: vh  ! m
      real(8):: ans !  temperature Kelvin

      real(4):: svh, sans

      svh = vh
      if( vh < heightx ) then
         if( vh >= h(1) ) then
            call kScsplIntp(h, temp, nsize, h2tempCoef, nsizem,
     *       svh, sans)
            ans = sans
         else
            ans = temp(1)-(vh-h(1))*0.3/100.
         endif
      else
         ans =  temp(nsize)
      endif
      end function cvh2temp


      function cvh2denp(vh) result(ans)
!               get drho/dh and its derivative 
      use nrl_atmos
      implicit none
      real(8),intent(in)::vh    !  m
      real(8)::ans
      
      real(4)::svh, sans, dum
      real(8):: rho 
      svh = vh
      if( vh < heightx ) then
         if( vh >= h(1) ) then
            call kScsplDif(h, nsize, h2denCoef, nsizem, svh, sans, dum)
            drhodh2 = dum
            ans = sans
            if( vh > 65.0d3) then
!
!           2nd derivative (dum) is not stable at ~ 72 km or so.
!           and h> 200km 
!           This is rather due to original data. 
!!
!           It's very small (|d(drho/dh)/dh| ~ 10^-(3~4)|(drho/dh)|)
!           so we put it 0 over 65 km.  Also if we take small
!           h step (100 m ), result becomes worse than   500 m step
!           so we use 500 m step.  single precision and double 
!           precision does  not change the situaiton. So we use
!           single precision. Also taking log of density cannot
!           improve the situaiton.

!!!!             drhodh2 = 0.    ! at prsent don't do this
            endif
         else
            call kScsplDif(h, nsize, h2denCoef, nsizem, h(1), sans, dum)
              ! not so good but don't care 
            ans = -al/bl*exp(-vh/bl) + al/bl*exp(-h(1)/bl) + sans
            drhodh2 = al/bl/bl*(exp(-vh/bl) - exp(-h(1)/bl))
     *             + dum
         endif
      else
         ans = -a/b*exp(-vh/b)
         drhodh2 = -ans/b
      endif

      end function cvh2denp
      
      function cvh2den2p(vh) result(ans)  
!                this must be called after cvh2denp
      use nrl_atmos
      implicit none
      real(8),intent(in)::vh    !  m
      real(8)::ans

      ans = drhodh2
      end function cvh2den2p
      
            
      subroutine cNRLdataRead(devno,filename)
      use nrl_atmos
      implicit none
      integer,intent(in):: devno ! temporary  logical dev. # of NRL data 
                          ! container
      character(*),intent(in)::filename ! file name of the container
                 ! file should contain height (m) vs density (kg/m3) table
                 ! first  h is -400 m.  step h should be ~ 500 m.  (100 m is
                 ! too small and spline interpolation becomes NG for the 
                 ! 2nd derivative 500 m gives as good as 100 m in the stable
                 ! region. 
      character(len=64):: line  ! 1 line should be <= 64 cha.

      integer:: i, icon
      call copenf(devno, filename, icon)
      if(icon /= 0 ) then
         write(0,*) ' file=',filename
         write(0,*) ' cannot be opened for cNRLdataRead'
         stop
      endif
      call cNRLHeaderRead(devno)

      call cskipComment(devno, icon)
!             find # of lines
      if(icon /= 0 ) then
         write(0,*) ' filename=',trim(filename)
         write(0,*) ' has no #------------- line'
         stop
      endif
      nsize=0
      do while ( .true. )
         read(devno,'(a)', end=100)  line
         nsize = nsize + 1
      enddo
 100  continue

      rewind devno

      nsizem = nsize - 1
      write(0,*) ' nsize=',nsize
      call cNRLalloc            ! allocate arrays

      call cskipComment(devno, icon)
      do i = 1, nsize
         read(devno, *)  h(i), den(i), temp(i)
      enddo
      close(devno)
      write(0,*) ' NRL atmospheric data has been read from'
      write(0,*)  filename
      end   subroutine cNRLdataRead

      subroutine cNRLHeaderRead(io)
      use nrl_atmos
      implicit none
      integer,intent(in):: io  ! read dev #

      character(len=10)::term1, term2
      character(len=20)::term3

      read(io,*) term1, term2, term3
      if( term2 /= "NRL") then
         write(0,*)
     *    'NRL data is needed since ATMOSPHERE 3 is specified'
         write(0,*)
     *    'in Zcondc.h '    
         write(0,*)
     *    'For the NRL atmosphere data: 1st line must be like'
         write(0,*) '# NRL  atmos...  but we have'
         write(0,'(a," ", a," ", a)') 
     *     trim(term1), trim(term2), trim(term3)
         stop
      endif
      read(io,*) term1, term2, glat
      if( term2 /= "lat") then
         write(0,*)
     *   'Atmosphere data for NRL: 2nd line must be like'
         write(0,*) '# lat  32.0'
         stop
      endif
      read(io,*) term1, term2, glong
      if( term2 /= "long") then
         write(0,*)
     *   'Atmosphere data for NRL: 3rd line must be like'
         write(0,*) '# long  132.0'
         stop
      endif
      read(io,*) term1, term2, gday1
      if( term2 /= "day1") then
         write(0,*)
     *   'Atmosphere data for NRL: 4th line must be like'
         write(0,*) '# day1  125'
         stop
      endif

      read(io,*) term1, term2, gday2
      if( term2 /= "day2") then
         write(0,*)
     *   'Atmosphere data for NRL: 6th line must be like'
         write(0,*) '# day2  135'
         stop
      endif

      read(io,*) term1, term2, ghour1
      if( term2 /= "hour1") then
         write(0,*)
     *   'Atmosphere data for NRL: 5th line must be like'
         write(0,*) '# hour1  3'
         stop
      endif

      read(io,*) term1, term2, ghour2
      if( term2 /= "hour2") then
         write(0,*)
     *   'Atmosphere data for NRL: 7th line must be like'
         write(0,*) '# hour2  1.0'
         stop
      endif
      end       subroutine cNRLHeaderRead

      subroutine cNRLHeaderWrite(io)
      use nrl_atmos
      implicit none
      
      integer,intent(in):: io    ! write dev #

      call cNRLHeaderW0(io)

      write(io,'(a)') 
     * "# 3 terms of the 1st 7 lines must be as above (except numbers)"
      write(io,'(a)')
     * "# any comment may be put after 3 terms. and lines below"
      write(io,'(a)')
     * "# until #-----------------  line "
      write(io,'(a)') "# h(m) den(kg/m3) T(K)"
      write(io,'(a)') "#-------------------------------"
      end  subroutine cNRLHeaderWrite

      subroutine cNRLHeaderW0(io)
      use nrl_atmos
      implicit none
      
      integer,intent(in):: io    ! write dev #

      write(io,'(a)') "# NRL atmosphere: # 3"
      write(io,'(a,f10.1)')  "# lat ", glat
      write(io,'(a, f10.1)') "# long ", glong
      write(io,'(a, i3)') "# day1 ", gday1
      write(io,'(a, i3)') "# day2 ", gday2
      write(io,'(a, i3)') "# hour1 ", ghour1
      write(io,'(a, i3)') "# hour2 ", ghour2
      end  subroutine cNRLHeaderW0

      subroutine cNRLLatLongCheck(lat, long)
!        should be called  when AtmosFile /= " "
!     for ATMOSPHERE 3 (NRL atmosphere)
!     if lat, long in AtmosFile differ from
!     this input, stop is made.
      use nrl_atmos
      implicit none
      real(8),intent(in):: lat   ! LatitOfSite in deg.
      real(8),intent(in):: long ! LongitOfSite in deg
! these two should be the same those in the Namlist input.
!     and if differ from glat glong  stop is made.
!     We may replace lat and long by glat glong to keep
!     consisency with Namelist input. But we worry about
!     doing so. (future program update may use LatitOfSite
!     and LongitOfSite before this rouine is called. etc).
!
      logical::warning
      warning =  abs(lat -glat) > 1.0
      warning= warning .or. abs(long -glong) > 2.0
      if( warning ) then
         write(0,*) 'ATMOSPHERE 3 (NRL atmosphere) is being used'
         write(0,*) 'and AtmosFile contains filename in which '
         write(0,*) 'latitude and longitude are gvien.'
         write(0,*) 'They differ from LatitOfSite and LongitOfSite'
         write(0,*) 'more than 1 or 2 deg.;diff must be < 1 or 2 deg'
         write(0,*) '              latitude         longitude '
         write(0,*) 'AtmosFile  ', glat, glong
         write(0,*) 'param      ', sngl(lat), sngl(long)
         write(0,*)
     *        'Give LatitOfSite and LongitOfSite values close to'
         write(0,*) 'those in AtmosFile; Or '
         write(0,*)
     *   'Consider using NRL_period without using AtmosFile' 
         stop
      endif
      end  subroutine cNRLLatLongCheck
      
      subroutine cNRLdataWrite(io, filename)
      use nrl_atmos
      implicit none
      integer,intent(in):: io  ! write dev #. If this is 6,
             ! filename is not referred.
      character(*),intent(in):: filename  ! filename path

      integer:: i, icon
      
      if( io /=  6) then 
         call copenfw(io, filename, icon)
         if( icon /= 0 ) then
            write(0,*) ' file=',trim(filename)
            write(0,*) 'cannot be opened for cNRLdataWrite'
            stop
         endif
      endif
      call cNRLHeaderWrite(io)
      do i = 1, nsize
         write(io, '(i7, 1p, E12.5, 0p, f8.1)') 
     *     int( h(i) ), den(i), temp(i)
      enddo
      if( io /= 6 ) then
         close(io)
      endif
      end      subroutine cNRLdataWrite
      
      subroutine cNRLGenData(lat, long, period)
      use nrl_atmos
!        get the  average of height vs air density table of the atmosphere 
!        at a given place during the given period specified by
!        starting time and ending time.  
      implicit none
      real(8),intent(in):: lat  ! latitute in deg. of the place
      real(8),intent(in):: long ! longitude in deg. of the place
!                  time info below is for the local time of the place
      integer,intent(in):: period(4) !  period(1):  starting day (Jan.1 is 1)  
                                   !                Dec. 31 is 365 (366 for leap
                                   !                year. but at present, year
                                   !                is not considered. so 365 is
                                   !                better?)    
                                   !    period(2):  endign day 1 ~ 365
                                   !    period(3):  starting hours
                                   !                0,24 is midnight. 12 is noon.
                                   !    period(4)   same for   ending hour
                                   !   
!          Let dayi=period(i) (i=1,2), hourj-2=period(j) (j=3,4)
!          If day1 > day2, the period is understood as straddling Dec. to Jan.
!          If hour1 > hour2, the period is understood as straddling midnight
!           
!          The data is generated by taking the average of data for sampled
!          times in the period:
!          a)  Samples days from day1 to day2 (7 day step; day1 and day2 is always
!              included. If the last but one day is close to day2 (< nearday) 
!              such day is skipped. 
!          b)  For each sampled day, get data for sampled hours between 
!          　hour1 and hour2. (4 hour step; hour1 and hour2 are always included.
!            If the last but one hour is close to hour2 (< nearhour) such  hour
!            is skipped.
!          The user can generate data, e.g, 
!          the average during day time of  the winter season (say, from Dec to Feb),
!          or night time of the same period.  Also it is possible to take the 
!          average of  whole days during the summer season etc.
!     As to the height, data is generated  from lowh to heightx, 
!         above heightx extrapolation by exp formula  is used.
!         Also below lowh, extrapolation is made but it must not be large.
 
      real(8),parameter:: lowh=-400d0  ! generate data from this height (m a.s.l)
      real(8),parameter:: step=500d0   ! height step m. this is enough      
      real(8),parameter:: dayinsec = 24*3600
      integer,parameter:: daystep  = 7  !7 day step in whole days
      integer,parameter:: nearday  = 3   !
      integer,parameter:: hourstep = 4   ! 4 hour step in each  day
      integer,parameter:: nearhour = 2
      integer,save::iday, day, tempday
      integer,save::hour, temphour
      real(8):: dalti

      real(4):: totalden(9), tempera(2)

      integer::nsum
!                usec: UT in sec in the day
!                stl:  apparent local soloar time in the day in hours
!               alti: atltitude in km
      real(4),save::usec, stl, alti
      integer::i

      nsize = (heightx-lowh)/step+1   
      nsizem = nsize - 1

      call cNRLalloc  ! allocate arrays

      den(:) = 0.     ! clear rho and temperature area
      temp(:) = 0.

      glat = lat       ! save lat, long
      glong = long

      gday1 = period(1)
      gday2 = period(2)
      ghour1 = period(3)
      ghour2 = period(4)

      nsum = 0
      day=gday1

      if(day > gday2) then
         day = day-366
      endif

      do while ( day <=   gday2 )
         hour = ghour1
         if( hour > ghour2) then
            hour = hour - 24
         endif

         do while ( hour <= ghour2 )
            if( day > 0 ) then
               tempday = day
            elseif( day < 0 ) then
               tempday = day+366
            elseif( day == 0 ) then
               tempday = 1
            endif
            
            if( hour > 0 ) then
               temphour = hour
            else
               temphour = hour + 24
            endif
 !             local apparent solar time in hours
            stl = temphour
             !  STL=UT(sec)/3600+GLONG/15   hours
            usec = (3600*stl - long*240)  ! UT in sec
            iday = tempday      ! day; year may be added; if y2015, 15*1000 may be
                           ! added
            dalti = lowh      ! lowest height
            nsum = nsum + 1
                ! for given place and time, compute 
                ! density at all heights
            do  i=1,  nsize
               alti = dalti/1000.  ! to km
!                         day  UT(sec) km   deg   deg   hours
               call GTD7(iday, usec, alti, glat, glong, stl,
     &              150., 150., 4., 48, totalden, tempera)
               den(i) = den(i) + totalden(6)   ! g/cm3
               temp(i) = temp(i) +  tempera(2) ! K
               dalti = dalti + step
            enddo

            if( hour == ghour2 ) exit

            hour = hour + hourstep
            if( ghour2 - hour <= nearhour) then
               hour = ghour2
            endif
         enddo

         if( day == gday2 ) exit

         day  = day + daystep
         if( gday2 -  day <= nearday ) then
            day = gday2
         endif
      enddo
!            get average
      den(:) = den(:)*1000./nsum  ! kg/m3
      temp(:) = temp(:) / nsum  ! K

!        make height table
      dalti = lowh
      do i=1, nsize
         h(i) = dalti      ! m
         dalti = dalti + step
      enddo
!           generate other data.
      call cNRLdataManip
      write(0,*) 'NRL atmospheric data has been made for'
      write(0,*) 'Latitute=',glat, ' Londitude=',glong, ' deg'
      write(0,*) 'for period=',period(:)
      end        subroutine cNRLGenData


      subroutine cNRLalloc
      use nrl_atmos
      implicit none
      if( allocated( h ) ) return !   !!!!
      allocate( h(nsize) )
      allocate( den(nsize) )
      allocate( temp(nsize) )
!!        allocate( denp(nsize) )
      allocate( thick(nsize) )
      allocate( h2denCoef(nsizem,3) )
      allocate( h2tempCoef(nsizem,3) )
!!        allocate( h2denpCoef(nsizem,3) )
      allocate( h2thickCoef(nsizem,3) )
      allocate( thick2denCoef(nsizem,3) )
      allocate( thick2hCoef(nsizem,3) )
      end  subroutine cNRLalloc
        
      subroutine cNRLdataManip
      use nrl_atmos
      implicit none
      integer:: i
      
      real(4)::dum

      call kScsplCoef(h, den, nsize, h2denCoef, nsizem)


      do i = nsize/2, nsize
         if( h(i) > heightx-50.0e3 ) exit
      enddo
      b = (h(nsize)-h(i)) /log(den(i)/den(nsize))
      a = den(i)*exp(h(i)/b)
      write(0,*) ' a,b=',a,b

      bl = (h(2)-h(1))/log(den(1)/den(2)) 
      al = den(1)*exp(h(1)/bl)
!         integral of den at each segment
      thick(:) = 0.
      do i =  1, nsizem
         call kScsplInteg(h, den, nsize, h2denCoef, nsizem,  
     *        h(i), h(i+1), thick(i) )
      enddo
           !  thickness above h> h(nsize) 
           !  ~  3.e-3  thickness at h(nsize) kg/m^2.  very rough
      thickx =  a*b*exp(-h(nsize)/b)
      thick(nsize) = thickx
!       make  cummulative 
      do i = nsizem, 1, -1
         thick(i) = thick(i) + thick(i+1)
      enddo
      
      call kScsplCoef(h, thick, nsize, h2thickCoef, nsizem)
      call kScsplCoef(thick, den, nsize, thick2denCoef, nsizem)
      call kScsplCoef(thick, h, nsize, thick2hCoef, nsizem)
      call kScsplCoef(h, temp, nsize, h2tempCoef, nsizem)
      end       subroutine cNRLdataManip
      
      function cthick2h(t) result(ans)
      use nrl_atmos
      implicit none
      real(8),intent(in)::t     ! kg/m2
      real(8)::ans
      
      real(4):: sans, st
      
      if( t > thickx ) then
         if( t <= thick(1) ) then
            st = t 
            call kScsplIntp(thick, h, nsize, thick2hCoef, nsizem, st, 
     *         sans)
            ans = sans
         else
            ans = -bl* log( (t-thick(1))/al/bl + exp(-h(1)/bl ))
         endif
      else
         ans = -b*log(t/a/b)
      endif
      end function cthick2h
      
      function cthick2den(t) result(ans)
      use nrl_atmos
      implicit none
      real(8),intent(in):: t    ! kg/m2
      
      real(8)::ans              !  kg/m3
      
      real(4):: sans, st
      if( t > thickx ) then
         if( t<= thick(1) ) then
            st = t
            call kScsplIntp(thick, den, nsize, thick2denCoef, nsizem,
     *           st, sans)
            ans = sans
         else
            ans= (t-thick(1))/bl + al*exp(-h(1)/bl)
         endif
      else
       !   ans = a*exp( b*log(t/a/b)/b) =a t/(ab) = t/b
         ans = t/b
      endif
      end function cthick2den
      
      function cvh2thick(vh) result(ans)
      use nrl_atmos
      implicit none
      real(8),intent(in):: vh   ! m
      real(8)::ans              !  kg/m2
      
      real(4):: sans, svh
      if( vh < heightx ) then
         if( vh >= h(1) ) then
            svh = vh
            call kScsplIntp(h, thick, nsize, h2thickCoef, nsizem, 
     *           svh, sans)
            ans = sans
         else
            ans = al*bl*(exp(-vh/bl) -exp(-h(1)/bl)) + thick(1)
         endif
      else
         ans =  a*b*exp(-vh/b)
      endif
      end function cvh2thick

      function cvh2scaleh(vh) result(ans)
      use nrl_atmos
      implicit none
      real(8),intent(in):: vh   ! m
      real(8)::ans              ! m 

      real(8),external:: cvh2den, cvh2denp
!
      ans = - cvh2den(vh)/ cvh2denp(vh)
      end    function cvh2scaleh
      
      include "./nrlmsise00_sub.f"

