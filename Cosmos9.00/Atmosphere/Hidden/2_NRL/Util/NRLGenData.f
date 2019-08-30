!  useage
!  0) prepare a file.  For that copy template.d to some file and
!     fill the data in that file.
!  1)  make clean; make
!  2)  ./a.out < that input_file  > output_datafile 
!     this  output contains   
!  
!        height, rho, thickness, hh, density, drhodh, drhodh2,
!
!      hh should be the same as height
!      density should be the same as rho
!      drhodz is drho/dh. 
!      drhodz2 is d(drhodz)/dh. This will show some irregular
!              behaviour but don't worry.
!     If the basic atmospheric  data is  to be read from a file
!        no other output is made.
!     else the basic atmosheric data is created as ./NRLAtmos.d
!     which can be used as the input in another (Cosmos) run.
!
      program main
        use nrl_atmos
        implicit none
        real(8):: height, rho, thickness, hh, density, drhodz, drhodz2
        real(8):: tk
        integer::i
        real(4):: dh=10.
        real(8):: lat, long
        real(8),external:: cvh2den, cthick2h, cvh2denp, cvh2den2p,
     *                    cvh2thick, cthick2den, cvh2temp 

        logical:: fromfile
        character(len=128) filepath
        integer::period(4)
        integer::icon

        call cskipComment(5,icon)
        read(*, *) fromfile 
        if( fromfile ) then
           read(*,*) filepath
        else
           read(*,*)   ! skip a line
           read(*,*) lat, long
           read(*,*) period
        endif
        if( fromfile ) then 
           call cNRLdataRead(11, filepath)
           call cNRLdataManip
        else
           call cNRLGenData(lat, long, period)
        endif
        write(*,'(a)')
     *  "H(m) rho(kg/m3) depth(kg/m2) H(m)  rho drhodz drhodz2 T(K)"
        height = -400. -1000. !  m
        do while (height < 1000.d3) 
           rho= cvh2den(height)
           thickness = cvh2thick(height)
           hh = cthick2h(thickness)
           density = cthick2den(thickness)
           drhodz = cvh2denp(height)
           drhodz2 = cvh2den2p(height)
           tk = cvh2temp(height)
           write(*,'(f10.2,  1p,6g15.6,0p, f10.2)' ) 
     *      height, rho, thickness, hh, density, drhodz, drhodz2,
     *      tk 
           height = height + dh
        enddo
        if(.not. fromfile)  then
           call cNRLdataWrite(20,"NRLAtmos.d")
        endif
      end program main
      
