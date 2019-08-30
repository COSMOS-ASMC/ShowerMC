module modMCSparam
  implicit none
!     multi socattering model
  character(100),save::MCSmodel='"Mol"'
  !      next few lines    COSMOSTOP /EPICSTOP=> LIBLOFT
  character(80),save::MCSdir = ' ' ! default: if from cosmos
  ! replaced by $COSMOSTOP/Data/MCS
  ! if from  Epics, replaced by $EPICSTOP/Data/MCS
  
  character(20),save::MCSparam = 'param_1.9-0.05'
  integer,parameter::MCSxyzCondmx=5
  ! a) suppose want to use MCSmodel only for  components in Z=(Z1,Z2)
  !  give 
  !       MCSzRange(1)=(Z1,Z2)
  !   If another range (Z3,Z4) is also target,
  !       MCSRange(2)= (Z3,Z4) may be given. (Z1< Z2 < Z3 < Z4)
  ! b) suppose want to use MCSmodel only for  component with comp #
  !             n1~n2
  !   give 
  !       MCSnumRange(1)=(n1, n2)
  !   For more compoents with n3~n4
  !   give
  !       MCSnumRange(2)=(n3,n4)  (n1< n2  < n3  < n4)
  !  If a component falls in the category, the actual model
  !  is determined by referring to energy and MCSmodel.
  !  If both a) and b) are specified, MCSandor is referrd.
  !  If it is 'and' (defualt), "and" is taken.    
  !  If it is 'or'  
  !  If (0,0), all compoents are target.
  !
  !  complex(8)::MCSxRange(MCSxyzCondmx)
  !  data  MCSxRange(:) /MCSxyzCondmx*(0.,0.)/
  !  complex(8)::MCSyRange(MCSxyzCondmx)
  !  data  MCSyRange(:) /MCSxyzCondmx*(0.,0.)/
  complex(8)::MCSzRange(MCSxyzCondmx)
  data  MCSzRange(:) /MCSxyzCondmx*(0.,0.)/

  integer,parameter:: MCSnumCondmx=5
  complex(8),save:: MCSnumRange(MCSnumCondmx)
  data  MCSnumRange(:)  /MCSnumCondmx*(0,0)/
  

  character(3),save::  MCSandor='and'
!  integer,save:: MCSxCond = 0   ! <= MCSxyzCondmx; actual cond.
!  integer,save:: MCSyCond = 0   ! <= MCSxyzCondmx; actual
  integer,save:: MCSzCond = 0   ! <= MCSxyzCondmx; actual # of conds
                                  !    x,y,z cond
  integer,save:: MCSnumCond = 0
  logical,save:: MCSrevert = .false.
  logical,save:: MCSdebug = .false.  ! if t, output world pos
  ! points with Moliere scat and El_hin scat on stderr
  ! format:  mol: x y z
  !          hin: x y z
  ! Be carefull not to output to many points. For 10GeV e incident
  ! 10 events is enough
  !
!  example of MCSmodel:
!    MCSmodel='"El_hin",10.0d-3,"El_con",100.0d-3, "Mol"'
end module modMCSparam

module modMCScontrol
  use modMCSparam, MCS=>MCSmodel
  implicit none
                !     max list.    Now 3 models
  integer,parameter::MaxHowMCS=6, NoOfUsableMCS=3
  real(8),parameter::MinErgMCS=100.0d-9  ! GeV ->100eV
  real(8),parameter::MaxErgMCS=1.0d0     ! GeV.  > this Moliere
                                   ! is referred.
  character(len=8),save::HowMCSList(MaxHowMCS)=' '
  real(8),save::MCSErg(MaxHowMCS)=1.d0
  integer,save:: NoOfHowMCS   ! # of actually specified List
  character(len=8),parameter::UsableMCS(NoOfUsableMCS)=  &
    (/"El_hin  ", "El_con  ", "Mol     "/)  !  embed " " gfortr

  real(8),save:: KEGeV, KEeV, mfpHardgr, lHardrl, lHardcm, lHardgr
  real(8),save:: tMCScm, pathScm
  logical,save:: doNewMCS   ! becomes F if always "Mol"
  character(len=4),save::MCSmode ! A1,A2,B1,B2,B3 used for El_hin 
  character(len=8),save::ActiveMCS  ! current MCS model; One of UsableMCS
  integer,save::MCSnumRangeMin(MCSnumCondmx) ! integer copy of
  integer,save::MCSnumRangeMax(MCSnumCondmx) ! MCSnumRangne
  data MCSnumRangeMin(:) /MCSnumCondmx*0/
  data MCSnumRangeMax(:) /MCSnumCondmx*0/

  !      length must be same as HowMCSList

end   module modMCScontrol

!    
subroutine cHowMCS
!
!      Establish how to treat Multiple Coulomb Scattering
!      below 1 GeV: any cobination of
!        El_hin, El_con, Mol  (hin means hinge, con means condense
!        Mol means to use old Moliere parameter).
!      is permitted. Format: E.g
!
!      MCS=' "El_hin" 1e-3 "El_con" 100e-3 "Mol"'
!  
!        100eV ~ 1MeV:  El_hinge 
!            ~ 100MeV:  El_condense
!            >       :  Mol (parameter Moliere is referred)
!      MCS='"Mol"'
!                 In the entire region, Moliere is referred.
!      
!      MCS=' "El_con" 10e-6 "El_hin" 100e-6 "El_con" '
!        100 eV ~ 10 keV: El_condense 
!               ~ 10MeV:  El_hinge
!               ~ 1GeV:   El_condense 
!               >         Mol
!  
!     This may used for test. If the result shows
!   almost the same as
!      MCS='"Mol"' or MCS='"El_hin"'
!   you may employ latter, since cpu time is much shorter.
!   (In stead of "El_con", "Mol" may be used in the last example).
!
!   The minimum energy used in Epics is normally 10keV to 100 keV
!   but, in the case of e+, it must be followed down to 0 energy
!   since it annihilates to produce 2Me absorbable energy.
!   In normal situations, it anhihilates before it's energy becomes
!   lower than 10 keV or so. However, sometimes we encounter
!   very low energy e+ of which energy is even below 100 eV.
!   Epics assumes thick  detectors and, treatment of
!   e-/e+ below *10 keV need not be so accurate because they
!   stop soon. So any Mol treatment is ok for them where no
!   data is available from  Elsepa.
! 
  use modMCScontrol
  implicit none

  integer i, j, l    

  l=len(MCS)
  MCS(l:l) = '/'
  read(MCS, *, Err=100) &
       (HowMCSList(i), MCSErg(i), i=1, MaxHowMCS)

  NoOfHowMCS = MaxHowMCS
  do i = 1, MaxHowMCS
     if(HowMCSList(i) .eq. ' ') then
        NoOfHowMCS = i-1
        exit
     endif
  enddo

  do i = 1, NoOfHowMCS-1
     if(MCSErg(i) >= MCSErg(i+1)) then !
!                 energy region invalid
        write(0,*)'MCSmodle is invalid; energy not ascending'
        write(0,*)' MCSmodel=', MCS
        stop
     endif
  enddo

  !    convert num ragen itno integer
  if( MCSnumCond > 0 ) then
     MCSnumRangeMin(1:MCSnumCond) = real(MCSnumRange(1:MCSnumCond))
     MCSnumRangeMax(1:MCSnumCond) = imag(MCSnumRange(1:MCSnumCond))
  endif
!         exam models
  do i = 1, NoOfHowMCS
     do j = 1, NoOfUsableMCS
        if( HowMCSList(i) ==  UsableMCS(j) ) goto 25
     enddo
     call cerrorMsg( HowMCSList(i),  1)
     call cerrorMsg('above MCSmodel is not yet registered', 0)
25   continue
  enddo
!!!!!!!!!!
  write(0,*) 'MCS model @KE(GeV)<'
  do i = 1, NoOfHowMCS
     write(0,*) HowMCSList(i), MCSErg(i)
  enddo
!!!!!!!!!!
  return
100 continue
  write(0,*) ' MCSmodel=',MCS, ' is invalid'
  stop
end subroutine cHowMCS

subroutine ciniMCS
  use modMCScontrol
  use modMCSparam
  implicit none
!    From Epics:
!        This must be called after config file has been read
!        and  cHowMCS has been called.  only once.
!    From Cosmos:  must be called after cHowMCS and cixsec 
!            has been read.
!
  
  integer::i
  doNeWMCS =.false.

  do i = 1, NoOfHowMCS
     if( index( HowMCSList(i), 'El_hin' )  > 0  .or.  &
         index( HowMCSList(i), 'El_con' )  > 0 ) then
        call ciniMixedMS
        doNewMCS =.true.
        exit
     endif
  enddo

  write(0,*) 'doNewMCS=',  doNewMCS
  if( doNewMCS ) then
     write(0,*) ' MCSparam is ', trim(MCSparam)
     if( MCSzCond > 0 ) then
        write(0,*) 'MCSmodel will be used for Z in:'
        do i = 1, MCSzCond
           write(0,*)  MCSzRange(i)
        enddo
     endif
     if( MCSnumCond > 0 ) then
        write(0,*) 'MCSmodel will be used for cn in:'
        do i = 1, MCSnumCond
           write(0,*) "(",MCSnumRangeMin(i),",", MCSnumRangeMax(i),")"
        enddo
     endif
     if( MCSnumCond > 0 .and. MCSzCond > 0 ) then
        write(0,*) &
         'both condistions are treated by: ', MCSandor
     endif

     if( MCSrevert ) then
        write(0,*) ' MCSrevert=', MCSrevert, ' so  we treat'
        write(0,*) ' El_hin part ==> Mol'
        write(0,*) ' Mol part    ==> El_hin'
     endif
  endif
end subroutine ciniMCS
subroutine csetActiveMCS
  use modMCScontrol
  implicit none
  ActiveMCS="Mol"
end subroutine csetActiveMCS
