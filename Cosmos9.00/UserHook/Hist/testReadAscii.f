!  compile: 
!  make -f testReadAscii.mk
!  exectution:
!  ./testReadAsciiPcLinuIFC < inputAscii.hist 
!  output is temp.hist (binary). To change the file
!  change the program.  temp.hist must not be existent.
!  
!
!  Make some ascii histogram data file: 
!  For that test purpose, you can use
!  test1.c which generate 3 histograms of power spectra
!  and 10 histograms by  gaussian distributions.
!  you can create it as ascii or binary file.
!  In case of binary file, you can use bin2ascii.c
!  to convert it to ascii file.
!       For usage see Readme.
	program main
	use modHistogram1
	implicit none

	integer kwhistReadAscii

	type(histogram1) h
	integer fni, fno
	integer icon

	fni=5
	fno=21
	call copenfw2(fno, "temp.hist", 2, icon)
	if(icon .ne. 0) then
	   write(0,*)  ' output file err'
	   stop
	endif
	call kwhistso(1)

	icon = 0
	do while ( icon .ge. 0  )
	   icon = kwhistReadAscii( h, fni) 
	   if(icon .ge. 0 ) then
	      write(0,*)
     *	      "None zero hist rec=",icon
!	      call kwhistp( h, -6) ! print ascii
	      call kwhistw( h, fno)  ! bin out
!                free h
	      call kwhistd( h)
	   endif
	enddo
	write(0,*) "last icon=", icon
	end
