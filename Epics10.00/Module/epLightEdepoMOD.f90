module modepLightEdepo
  use modepLightMaxDef
!  energy deposit --> scintillation source 
 type scintiSrc
    real(4),allocatable::edepo(:,:,:)
 end type scintiSrc

 type(scintiSrc),allocatable::parts(:)
   !  for i-th light component
   !   if cylinder 
   !   allocate( parts(i)%edepo(-nx:nx, -ny:ny, 1:nz) )
   !      if box
   !   allocate( parts(i)%edepo(0:nx,0:ny,0:nz) )
end module modepLightEdepo
