!  correction of sub/ for init and ending of
!   one event, all event  

!  init for 1 event generation
!
!   
subroutine epLighti1ev
!  use modepLightMaxDef
  use modepLight
  use modepLightPty
  use modepLightCounter
  implicit none
  integer::i

  do i = 1, LightCompNo
     if( allocated( sensor(i)%pe ) ) then
        deallocate(sensor(i)%pe )
!        deallocate(sensor(i)%pet )
     endif
     if( allocated( exiting(i)%pc ) ) then
        deallocate( exiting(i)%pc )
        deallocate( exiting(i)%pct )
     endif
     if( allocated( entering(i)%pc ) ) then
        deallocate( entering(i)%pc  )
        deallocate( entering(i)%pct  )
     endif
!     if(allocated( lightcomp(i)%dE )) then
!        deallocate( lightcomp(i)%dE )
!     endif
  enddo
!          weight counter 
  sumni = 0.
  sumniwi = 0.

end subroutine epLighti1ev




!  end of 1evet
subroutine epLighte1ev
  use modepLight
  use modepLightPty
  use modepLightCounter
  implicit none
  integer::i, index
  real(8):: peakWL, Epeak, wl0
  character(12):: struc, matter  !>=MAX_MEDIANAMELENG=8 

  real(4):: lightcomp

!             alloc     not alloc
! sensor        pe:  0           
! sensor       pet:  0  
! lightcomp   edep:  

  do i = 1, LightCompNo
     if( Lcomp(i)%mn >= SensorNoMin ) then
!        if( allocated( lightcomp(i)%dE ) .or. allocated( sensor(i)%pe )) then
!           if(.not. allocated( lightcomp(i)%dE ) ) then
!              allocate(lightcomp(i)%dE )
!              lightcomp(i)%dE  = 0.
!           elseif( .not. allocated( sensor(i)%pe )) then
        if( .not. allocated( sensor(i)%pe )) then
           allocate( sensor(i)%pe(1:5))
!           allocate( sensor(i)%pet )
           sensor(i)%pe(:) = 0.
        endif
!           lightcomp(i)%dE = lightcomp(i)%dE * 1.e9   ! in eV
              ! convert direct hit effect into p.e.
              ! the sensor is attached to a scintillator of 
              ! which property can be  obtained via next index
!        //////////
!        call Lcompchk( ' new ',  Lcomp(i)%comInfoNo)
        !        ////////////
        call epLightMn2cominfoIdx( comInfo( Lcomp(i)%comInfoNo)%mnOfScinti,  index)
        if( index == 0 ) then
           write(0,*) ' warning: Sensor (comp #=',Lcomp(i)%compno,') seems to be attached'
           write(0,*) ' to a component which is not specified for light ray tracing.'
           write(0,*) ' Sensor mn=',  Lcomp(i)%mn
           write(0,*) ' This might be due to that you reset countDE for the scintillator '
           write(0,*) ' but forgot to reset sensor countDE for the sensor attached to it'
           call epqstruc(Lcomp(i)%compno, struc)
           call epqmat(Lcomp(i)%compno, matter)
           write(0,*) ' Sensor is ', trim(struc), ' media=',trim(matter)
        else
           peakWL = comInfo( index )%peakWL
           call epLightwl2E(peakWL, 1.d0, wl0, Epeak)
!              sensor(i)%pe(5) = lightcomp(i)%dE/Epeak * comInfo(index)%cf
           sensor(i)%pe(5) = lightcomp(i)*1.d9/Epeak * comInfo(index)%cf
        endif
        sensor(i)%pet = sum( sensor(i)%pe(:) )
     endif

!          photon # counter
     if(  Lcomp(i)%mn >= LightNoMin ) then
        if( allocated (entering(i)%pc ) ) then
           ! for each surface, sum of Scinti and Ceren
           entering(i)%pct(:) = sum ( entering(i)%pc(:,:),1) 
        endif
        if( allocated (exiting(i)%pc) ) then
                    ! for each surface, sum of Scinti and Ceren
           exiting(i)%pct(:) = sum ( exiting(i)%pc(:,:),1) 
        endif
     endif
  enddo

end subroutine epLighte1ev

!  return p.e (ingredient and total)
subroutine epqLightSensor(cn, edepo,  ing,  total, icon )
 ! use modepLightMaxDef
!  use modepLight
  use modepLightPty
  use modepLightCounter
  implicit none
  integer, intent(in):: cn  ! component number of the sensor ( PWO etc )
  real(4), intent(out):: edepo  ! total energy deposit  in the component (GeV)
                 ! has meaning if  (B,mn > 0). else 0. 
  real(4), intent(out):: ing(5) ! p.e by scinti, Ceren, syc, tran, hit
                            ! if this is not a sensor, all 0's are returned
  real(4), intent(out):: total  ! sum of the above;   if this is not a
                                ! sensor,  will be 0.
  integer, intent(out):: icon   ! 0--> obtained, 1--> no data for this cn.

  real(4)::lightcomp

  integer:: Lno 
!  cn to Lcom number
  call epLightCn2LcompNo(cn, Lno)
  if(Lno == 0 ) then
     ing(:) = 0.
     total = 0.
     edepo = 0.
     icon = 1
  elseif(  Lcomp(Lno)%mn < SensorNoMin) then  
     ing(:) = 0.
     total = 0.
     edepo = 0.
     icon = 1
  else
     if( allocated(sensor(Lno)%pe ) ) then
        ing(:) = sensor(Lno)%pe(:)
        total =  sensor(Lno)%pet
!        edepo =  lightcomp(Lno)%dE
        edepo =  lightcomp(Lno)*1.d9
        icon = 0
     else
        icon = 1
     endif
  endif
end  subroutine epqLightSensor


!  return photon counter
subroutine epqLightPC(cn, info,  pc,  pct, icon)
  use modepLightMaxDef
  use modepLightCounter
  implicit none
  integer, intent(in):: cn  ! component number of the sensor
  integer, intent(in):: info  ! 1--> get entering photons
                              ! 0--> get exiting photons 
  real(4), intent(out):: pc(2,*) ! p.c by scinti, Ceren at each  surface
                         ! :=  6 for box, 3 cyl and ecyl, 4 pipe
  real(4), intent(out):: pct(*)  ! sum of the above
  integer, intent(out):: icon    ! 0==> obtained, 1--> no data for this cn





  integer:: Lno, ns


  call epLightCn2LcompNo(cn, Lno)

  if(Lno == 0 ) then
!      in this case, pc and  pct are undefined.
     icon = 1
  else
     call epLightNoOfSurf(cn, ns)  ! # of surfaces
     if( info <= 0 ) then
        !  exiting
        if( allocated( exiting(Lno)%pc  ) ) then
           pc(1:2,1:ns) = exiting(Lno)%pc(1:2,1:ns)
           pct(1:ns) = exiting(Lno)%pct(1:ns)
           icon = 0
        else
           pc(1:2,1:ns) = 0.
           pct(1:ns) = 0.
           icon = 1
        endif
     else
        !  entering
        if( allocated( entering(Lno)%pc  ) ) then
           pc(1:2,1:ns) = entering(Lno)%pc(1:2,1:ns)
           pct(1:ns) = entering(Lno)%pct(1:ns)
           icon = 0
        else
           pc(1:2,1:ns) = 0.
           pct(1:ns) = 0.
           icon = 1
        endif
     endif
  endif

end  subroutine epqLightPC

subroutine epqLightAveW(avew, tracedP)
  use modepLight
  implicit none
!  returns average weight of traced photons
  real(4),intent(out)::avew  !  sum(NiWi)/sum(Ni)
  real(4),intent(out)::tracedP  ! Ni

  if( sumni > 0. ) then
     avew = sumniwi / sumni
     tracedP = sumni
  else
     avew = 0.
     tracedP = 0.
  endif
end subroutine epqLightAveW

