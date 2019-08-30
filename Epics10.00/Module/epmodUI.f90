module modUI
       ! In what follows, array index means the component number.
  real(4),save,allocatable::ElossT(:)  ! full (true) energy loss in GeV
                    ! available for all components
  integer,save:: minQComp=0    ! min and
  integer,save:: maxQComp=0    ! max comp. #  
                   ! of those components which have quenching property.
                   ! (Note, the components within this number range 
                   ! are not necessarily such ones).
       ! The size of the folllowing arrays will be defined as
       !      (minQComp:maxQComp)
  real(4),save,allocatable::ElossE(:)  ! effective energy loss 
                    ! for photon emission. (in GeV).  computed 
                    ! for components with quenching property.
  real(4),save,allocatable::RdEH(:) !sum of dE by heavy ion's restriced 
                            ! enegy loss
  real(4),save,allocatable::RdES(:) !same by single charge (i.e, by
                        ! delta rays)
  real(4),save,allocatable::dED(:) ! dE by K.E of dying ptcl.
                          ! this one is included in RdEH and/or RdES.
!  real(4),save,allocatable::qfactorCff(:) ! Quenching 
                  ! factor for dEdxf. In case of Tarle's formula
                  !  (1-c2)/(1+(1-c2)c1*dEdxf) + c2
!!  real(4),save,allocatable::qfactorQf(:) ! qf of 
                              ! dEdx of heavy is quenched
!!  real(4),save,allocatable::qfactorC0(:) ! Tarle's quenching
                         ! factor (1-c2)/(1+(1-c2)c1*dEdxf) 
!  real(4),save,allocatable::qfactorCf(:) !cf in  qf*R*cf
         !Above one is not used if cqbycf is not used
         !  (i.e, epUI.f-woTarle and epgetCq.f-woTarle are
         ! employed
!  real(4),save::qcf  !   qfactorCf of the current Cn
!!  real(4),save,allocatable::qfactorC2(:) !c2 
end module modUI

!  these are intentinally not contained in modUI (to keep user 
!  program unchaged)
  subroutine epAlloc(compno)
    use modUI
    implicit none
    integer,intent(in):: compno  ! total number of components
    allocate(ElossT(compno) )  ! true  loss for

    call epMinMaxQuenchComp(MinQComp,MaxQComp) ! get min/max comp# which has quench
                   !effect
    if( MinQComp > 0 ) then
       allocate(ElossE(MinQComp:MaxQComp) )
       allocate(RdEH(MinQComp:MaxQComp) )
       allocate(RdES(MinQComp:MaxQcomp) )
       allocate(dED(MinQComp:MaxQComp) )
!       allocate(qfactorCff(MinQComp:MaxQComp) )
!       allocate(qfactorC0(MinQComp:MaxQComp) )
!       allocate(qfactorQf(MinQComp:MaxQComp) )
!       allocate(qfactorC2(MinQComp:MaxQComp) )
!       allocate(qfactorCf(MinQComp:MaxQComp) )
    endif
  end subroutine epAlloc
 
!  subroutine epqEloss(i, dET, dEeff, dEcff)
  subroutine epqEloss(i, dET, dEeff)
    use modUI
    implicit none
    integer,intent(in)::i  ! comp. no
    real(4),intent(out)::dET  ! total energy deposit in that comp (GeV)
    real(4),intent(out)::dEeff  ! effective dE with quenching (GeV)
                      ! if comp. is not quenching one or
                      ! no  heavy ion passed it, same as dET
!    real(4),intent(out)::dEcff ! dET * qfactorCff. 
                     ! if comp. is not quenching one or
                     ! no  heavy ion passed it, same as dET
                     ! If two or more heavy ions passed, or
                     ! the heavy ion interacted inside, 
                     ! value is not reliable.
                     ! Otherwise, the most probalbe value
                     ! should be ~ dEeff but distribution would
                     ! be smaller than dEeff.
    dET = ElossT(i)
    if(i >= minQComp .and. i<= maxQComp) then
       if( ElossE(i) > 0. ) then
          dEeff = ElossE(i)
 !         dEcff = qfactorCff(i)*dET
       else
          dEeff = dET
 !        dEcff = dET
       endif
    else
       dEeff = ElossT(i)
!       dEcff = dET
    endif
  end subroutine epqEloss
  subroutine epqElossItem(i, dEbyH, dEbyS, dEbyD)
    use modUI
    implicit none
    integer,intent(in)::i  ! comp. no
    real(4),intent(out):: dEbyH 
    real(4),intent(out):: dEbyS
    real(4),intent(out):: dEbyD
    if( i>= MinQComp .and. i <= MaxQComp) then 
       dEbyH = RdEH(i)
       dEbyS = RdES(i)
       dEbyD = dED(i)
    else
       dEbyH = 0.
       dEbyS = 0.
       dEbyD = 0.
    endif
  end subroutine epqElossItem
! clear Eloss counter
  subroutine epcEloss
    use modUI
    implicit none
    
    ElossT(:) = 0.
    if( MinQComp > 0 ) then
       ElossE(MinQComp:MaxQComp) = 0.
       RdEH(MinQComp:MaxQComp) = 0
       RdES(MinQComp:MaxQComp) = 0
       dED(MinQComp:MaxQComp) = 0
!       qfactorCff(:) = 0
    endif
  end subroutine epcEloss

  subroutine epdEItemUpdt(i, charge, dE)
!    count energy deposit by itemizing 
!  it; by heavy ion or by  single charge
!  This should be called for comp. with quenching effect
    use modUI
    implicit none
    integer,intent(in):: i  ! component #
    integer(2),intent(in):: charge
    real(8),intent(in):: dE ! energy deposit by ionization
                     ! (from restricted energy  part)
    if( charge > 1 ) then ! heavy
       RdEH(i) = RdEH(i) + dE
    elseif( charge /= 0 ) then
       RdES(i) = RdES(i) + dE
    else
       write(0,*) " charge=", charge, " is invalid "
       write(0,*) " for epdEItemUpdate"
       stop
    endif
  end subroutine epdEItemUpdt

  subroutine epdEItemUpdt2(i, dE)
!    count energy deposit by dyeing  ptcl.
!  This should be called for comp. with quenching effect
    use modUI
    implicit none
    integer,intent(in):: i  ! component #
    real(8),intent(in):: dE ! energy deposit by dyeing 
       ! ptcl. (=K.E)
    dED(i) = dED(i) + dE
  end subroutine epdEItemUpdt2

!!  subroutine epqqfactor(i, cff)
!!    use modUI
!!    implicit none
!!    integer,intent(in):: i  ! component #
!!    real(8),intent(out):: cff ! Quenching factor
!!                            ! for the heavy primary of 
!!                            ! current event. if 1,
!!                            ! no heavy ion passed the comp. i
!!    if( i < MinQComp .or. i > MaxQComp ) then
!!       cff = 1.
!!    else
!! !      cff = qfactorCff(i)
!!    endif
!!  end subroutine epqqfactor
