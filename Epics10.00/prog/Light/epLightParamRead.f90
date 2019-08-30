!     *****************************************************************
!     * epLightParmRead; Read parameters for light ray tracing 
!    *****************************************************************
!


subroutine epLightParamRead0(dir, mn, parm)

  use modepLightPty
  implicit none
  
  character(len=*),intent(in)::dir   ! directory where Light10.dat etc exists
  integer,intent(in)::mn       ! mn for Lightmn.dat or Sensormn.dat
  type(property),intent(out)::parm

  character(len=256)::path

  
  if(mn < SensorNoMin ) then
     ! for  Lightmn.dat.   i2.2 means non zero suppress in head part
     write(path,'(a,"/Light",i2.2,".dat")') trim(dir), mn
  elseif(mn < SensorNoMax) then
     ! for Sensormn.dat
     write(path,'(a,"/Sensor",i2.2,".dat")') trim(dir), mn
  else
     write(0,*) ' mn=',mn,'  for sensor in CountDE is wrong'
     write(0,*) ' detected in epLightParamRead0'
     stop 1234
  endif
  call epLightParamRead(path, parm)
end subroutine epLightParamRead0



subroutine epLightParamRead(file, parm)
  use modepLight
  use modepLightPty
  implicit none
  type(property),intent(out)::parm
  character(len=*) file


  integer::icon
  character(24)::vname
  character(100):: vvalue
  logical::epgetParmN
  call copenf(iowk, file, icon)
  if(icon .ne. 0) then
     write(0,*) ' file:',trim(file),' could not be opend'
     stop 7979
  endif
! skip separator
  call afsep(iowk)
!             read each item
  do while ( epgetParmN(iowk, vname, vvalue) )
     if( vname == 'refracIndex' ) then
         parm%refracIndex = 0  ! will be replaced by media value if  
                               ! vvalue is "   /".
        read(vvalue, *) parm%refracIndex
     elseif( vname == 'refracFile' ) then
        read(vvalue, *) parm%refracFile 
     elseif( vname == 'refMode') then
        read(vvalue, *)  parm%refMode
     elseif( vname == 'wrapper') then
        read(vvalue, *)  parm%wrapper
     elseif( vname == 'wrapperReflecFile') then
        read(vvalue, *) parm%wrapperReflecFile 
     elseif( vname == 'wrapperRefracFile') then
        read(vvalue, *) parm%wrapperRefracFile 
     elseif( vname  == 'mirror') then
        read(vvalue, *)  parm%mirror
     elseif( vname == 'fuzzy') then
        read(vvalue, *)  parm%fuzzy
     elseif( vname ==  'maxPathFac') then
        read(vvalue, *)  parm%maxPathFac
     elseif( vname ==  'Rayleigh') then
        read(vvalue, *)  parm%Rayleigh
     elseif( vname ==  'Mie') then
        read(vvalue, *)  parm%Mie
     elseif( vname .eq. 'attenL' ) then
        read(vvalue, *)  parm%attenL
     elseif( vname .eq. 'attenFile') then
        read(vvalue, *)  parm%attenFile
     elseif( vname  == 'WLS' ) then
        read(vvalue, *)  parm%WLS
     elseif( vname .eq. 'NpPerMeV') then
        read(vvalue,*)  parm%NpPerMeV
     elseif( vname .eq. 'waveLenFile') then
        read(vvalue,*)  parm%waveLenFile
     elseif( vname .eq. 'peakWL') then
        read(vvalue,*)  parm%peakWL
     elseif( vname .eq. 'quench') then
        write(0,*) '***********************************' 
        write(0,*) '***********************************' 
        write(0,*) ' In ', trim(file) 
        write(0,*) ' "quench" is obsolete: use "Modify" file'
        write(0,*) ' to give comp. specific quenching factor'
!        read(vvalue,*)  parm%quench
        write(0,*) '***********************************' 
        write(0,*) '*************sorry*****************' 
        stop
     elseif( vname .eq. 'minmaxWL' ) then
        read(vvalue,*)  parm%minmaxWL
     elseif( vname .eq. 'NpSample') then
        read(vvalue,*)  parm%NpSample
     elseif( vname .eq. 'CellSize') then
        read(vvalue,*)  parm%CellSize
     elseif( vname .eq. 'Qeff') then
        read(vvalue,*)  parm%Qeff
     elseif( vname .eq. 'QeffFile') then
        read(vvalue,*)  parm%QeffFile
     elseif( vname .eq. 'mnOfScinti') then
        read(vvalue,*)  parm%mnOfScinti
     elseif( vname .eq. 'cf') then
        read(vvalue,*)  parm%cf
     else
        write(0, *) 'error in file=', file
        write(0, *) vname, ' is not defined '
        stop
     endif
  end do
  close(iowk)
end subroutine epLightParamRead

subroutine epLightParamWrite(io, parm)
  use modepLight
  use modepLightPty
  implicit none
  integer, intent(in) :: io  !   output logical dev. # . must have been
                             ! opened. 
  type(property),intent(in)::parm  !  comp & media property parameter

  write(io, '(a)') '--------------------------'

  call awprmr(io,'refracIndex', parm%refracIndex)
  call awprmc(io, 'refracFile', parm%refracFile )
  call awprmra(io, 'refMode', parm%refMode, maxSurfaces)
  call awprmra(io, 'wrapper', parm%wrapper, maxSurfaces)
  call awprmc(io, 'wrapperReflecFile', parm%wrapperReflecFile )
  call awprmc(io, 'wrapperRefracFile', parm%wrapperRefracFile )
  call awprmra(io, 'mirror',   parm%mirror, maxSurfaces)
  call awprmra(io,'fuzzy',  parm%fuzzy, maxSurfaces)
  call awprmr(io,'maxPathFac',  parm%maxPathFac)
  call awprmr(io,'Rayleigh',  parm%Rayleigh)
  call awprmr(io,'Mie',  parm%Mie)
  call awprmr(io,'attenL', parm%attenL)
  call awprmc(io,'attenFile', parm%attenFile)
  call awprmra(io,'WLS', parm%WLS, 4)
  call awprmr(io,'NpPerMeV',  parm%NpPerMeV)
  call awprmc(io, 'waveLenFile',parm%waveLenFile)
  call awprmr(io, 'peakWL', parm%peakWL)
  call awprmr(io, 'quench', parm%quench )
  call awprmra(io, 'minmaxWL',  parm%minmaxWL, 2)
  call awprmr(io, 'NpSample', parm%NpSample)
  call awprmra(io, 'CellSize',   parm%CellSize, 3)
  call awprmr(io, 'Qeff', parm%Qeff)
  call awprmc(io,  'QeffFile', parm%QeffFile)
  call awprmi(io,  'mnOfScinti', parm%mnOfScinti)
  call awprmr(io,  'cf', parm%cf)
end subroutine epLightParamWrite

