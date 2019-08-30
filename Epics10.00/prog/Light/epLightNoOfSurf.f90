subroutine epLightNoOfSurf(cn, nsurf)
!         get no. of surfaces of the comp. # cn
  implicit none
  integer, intent(in):: cn ! component #
  integer, intent(out):: nsurf ! # of surfaces of the component
      
  character(len=16)::struc
  !///////////////////
   if( cn == 0 ) then
      write(0,*) ' cn =', cn, ' in NoOfSurf'
   endif
  !//////////////
  call epqstruc(cn, struc)
      
  if(struc(1:3) == "box") then
     nsurf = 6
  elseif(struc(1:7) == "octagon" ) then
     nsurf = 10
  elseif(struc(1:3) == "cyl") then
     nsurf = 3
  elseif(struc(1:4) == "ecyl") then
     nsurf = 3
  elseif(struc(1:4) == "pipe") then
     nsurf = 4
  elseif(struc(1:5) == "prism") then
     nsurf = 5
  else
     write(0,*) ' component #=',cn, ' struc=',struc
     write(0,*) ' not yet supported in epLightNoOfSurf'
     stop
  endif
end subroutine epLightNoOfSurf


