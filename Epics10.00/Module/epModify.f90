!       Modify file format
!  before ---------------------- line, all are comment
!  a line starting form the 2nd col. is data.  Each data must
!  be followed by "/".   After "/", all are neglected.
!  others are  all commnet
! 123456789*123456789*
! ------------------------------
!  1  /    entry # 1.  this corresponds to the number in the modifier
!                      field of configdata.  then folloiwng data, if given,
!                      is for this component. The number need not be
!                      consequtive nor from small to large. But must be
!                      kept as small as possible to save the memory size
!                      maximum # must be < 2^15
!  Quench   a b c B /  Birks's quenching factor (only first 'a' is used but
!                      b c must be put.           
!                      Other format is
!                      Quench  a b  T /    Tarle's formula is assumed
!                      Quench  a b c L /   Log type formula is assumed
!  Emin    Eg Ee Er /  or Eg Ee Er Ek En / 
!                      minumum energy for the component.
!                      ***Energy is same as in epicsfile 
!                      all in GeV.  Ee is total energy.
!                      If Eg /, then only Eg is changed.
!                      Eg: for gamma
!                      Ee: for electron
!                      Er: recoil minimum energy above which random
!                          sampling of delta-ray is performed
!                      Ek: kinetic energy of other charged particles
!                      En: //      neutron. (typical one is ~20000 --> 20 MeV)
!   
!       *** If the value of some of Eg, Ee, ...is negative,  its absolute
!           value is used; however,
!           any  particle, whose energy beomces lower than the minimum 
!           specified here, is discarded instantly (That is, residual
!           range consideration by  AutoEmin=2 or AutoEmin=4  is
!           invalidated. ***

!
!               
module epModify
  implicit none
#if defined (Solaris) || (jaxa) || (jaxaflat)
#else
  private
#endif
!  next assignment is to avoid possilbe error when testCnf1.f is used
!  (which dose not read epicsfile and EminGsave etc is still not given
!  a value)
!
  real(8),save,public:: EminGsave = 100.d-6
  real(8),save,public:: EminEsave = 611.d-6
  real(8),save,public:: RecoilEsave = 100.d-6
  real(8),save,public:: KEminsave = 20.d-3      ! 20 MeV
  real(8),save,public:: Enminsave = 20.d-3      ! 20 MeV
  real(8),save,public:: ElecMass =  511.d-6    
  integer,public,parameter::bitEmin=1  ! 2nd bit
  integer,public,parameter::bitQuench=0  ! 1st bit

  type quench
     real(8):: a
     real(8):: b
     real(8):: c
     character*1:: id
  end type quench
  type minErg   
     real(8):: Egmin
     real(8):: Eemin
     real(8):: Recoilmin
     real(8):: KEmin
     real(8):: Enmin
     logical:: imperative =.false.  
             ! if some Emin < 0,  abso is used
             !  but  this will become .true. and
             !  any ptcl< Emin is discarded instantly.
             !  (No residual range considereation).    
  end type minErg

  type  modifier
     integer::kind   !   to specify what modification is included 
                     !   by using bit; yz-->  z: quench  y: emin 
     type(quench):: q
!     real(4)::rho    ! not yet used
     type(minErg):: Em
  end type modifier

  type(modifier),public,allocatable::modify(:) 
  public::epfixModifier   !  only this may be used by the user
  integer,save,public:: maxModifyNum=0
!               -------------------------
  logical,save::first=.true.  
  integer,private,parameter:: iowk=11  ! this must be the same as in ZepManager.h


  character*256,save:: line
  character*64,save::  temp

  integer,parameter::nfmax=8
  character*16,private,save:: Field(nfmax)
  integer nf
  integer,save:: n, maxn, pos, pos2
  integer,save:: ncomp ! how many components have modification
  integer,save:: error
  character*1  id


contains
  subroutine epfixModifier(modfile)
!    read Modify file and fix the modifiation list
!          this is called after epicsfile is read
!       and if ModifyMile is non blank. 
!     
    implicit none
    character(*),intent(in):: modfile !input ModifyFile 
    integer:: icon
!      If modfile is blank, this is not called but
!      check it here. 
    if( first .and. modfile /= " " ) then
       call epcountModifier(modfile)   ! count number of entries first
                                       ! and allocate needed memories. 
       rewind iowk  !      rewind file
       call epreadModifier(modfile)          !  actually read 
       first = .false.
    endif
  end subroutine epfixModifier

  subroutine epcountModifier(modfile)
    implicit none
    character(*),intent(in):: modfile
    integer::icon
    logical kalpha
!           read 
    call copenf(iowk, modfile, icon)
    if(icon /= 0) then
       write(0,*) ' ModifyFile=',trim(modfile)
       write(0,*) ' could not be opened'
       stop 87623
    endif
    call afsep(iowk)  ! skip until "---------------" line
!             see how many data
    maxn = 0  ! max index number 
    ncomp = 0 ! # of comp. counter
    do while (.true.)
       read(iowk, '(a)', end=100) temp
       if( temp(1:1) /= ' ' ) cycle
       if( temp(2:2) == ' ' ) cycle
       if( kalpha( temp(2:2) ) ) cycle
       call kgetField(temp(1:10), Field, 1,  nf)

       read(Field(1),*)  n
       if( n > maxn) maxn = n
       ncomp = ncomp + 1
    enddo
100 continue
    if(maxn > 2**15-1 )  then
       write(0,*) ' Too large entry # in the Modify file=', maxn
       write(0,*) ' must be <',  2**15-1
       write(0,*) ' File name is ', trim(modfile)
       stop
    else
       write(0,*) ' Modifyfile=',trim(modfile), ' has data for', &
            ncomp, ' components, max index=', maxn 
    endif

!          allocate size
    allocate( modify(maxn))
    maxModifyNum = maxn
  end subroutine epcountModifier
             

  subroutine epreadModifier(modfile)
    implicit none
    character(*),intent(in):: modfile  ! used only to show the file name
                       ! when error takes place
    integer::i

    logical:: kalpha
    logical:: first=.true.

    call afsep(iowk)
    error = 0

    modify(:)%kind = 0   ! clear  


    do while (.true.)
       read(iowk,'(a)', end=200)  line

       if( line(1:1) /= ' ') cycle
       if( line(2:2) == ' ') cycle
       call kgetField(line(1:10), Field, 1,  nf)
       if(nf <= 0) cycle
       if( trim( Field(1) ) == "#" ) cycle

       if(  .not. kalpha( line(2:2) ) ) then
          read(Field(1),*)  n          
          if( modify(n)%kind /= 0 ) then
             ! double definition
             error = error +1
             write(0,*) ' ModifyFile data error; double definition'
             write(0,*) ' entry # =', n, ' already defined'
             write(0,*) ' modify(n)%kind =', modify(n)%kind 
          endif
          first = .false.
       elseif( first ) then
            ! without entry #, data appeared; error
          write(0,*) 'In the modfication File ', trim(modfile)
          write(0,*) 'data appeared without entry # /'
          write(0,*) 'the line is ', trim(line)
          stop
       else
          !  data entry.
          if( Field(1) == "Quench" ) then
             call epquenchModify
          elseif( Field(1) == "Emin" ) then
             call epEminModify
          else
             write(0,*) ' data entry name =', Field(1), 'undefined'
             error = error + 1
          endif
       endif
    end do
200 continue
    if(error > 0 ) then
       write(0,*) 'ModifyFile', trim(modfile),  &
            '  has ', error, '  errors'
       stop 9999
    endif
    close(iowk)
  end subroutine epreadModifier

  subroutine epquenchModify
    implicit none
    real(8)::a, b, c
    ! get quenching field
    pos = index(line, "/")
    if( pos == 0 ) then
       write(0,*) ' Quench line for entry #=',n
       write(0,*) ' has no /'
       error = error + 1
       return
    endif

    temp = line(1:pos-1)  ! original would be below,but temp has no /
                       ! Quench  a b c B/ 
                       ! Quench  a b c B /
                       ! Quench  a b c  /
                       ! Quench  a b c/
                       ! Quench  a b c L /
                       ! Quench  a b T /
                       ! Quench  a b T/
!                       12345678 
    c =0.  ! dummy

    call kgetField(temp(8:), Field, 4,  nf)

    id = ' '
    if( nf == 4 ) then
       read(temp(8:), *)  a, b, c, id
       if(id == "B" .or. id =="L" ) then  ! ok
       else
          write(0,*) 'quenching filed ID =', id , ' undef, entry=', n 
          write(0,*) ' the line is ', trim(line)
          error  = error + 1
       endif
       modify(n)%q%c = c
       modify(n)%kind = IBSET( modify(n)%kind, bitQuench)  
    elseif( nf == 3 ) then
       read(temp(8:),*)  a, b, id
       if(id /= "T") then
          read(temp(8:),*) a, b, c   ! assumed to be Birks  type
          id = "B"
       endif
       modify(n)%kind = IBSET( modify(n)%kind, bitQuench)  
    else
!             no quenching modification
!       write(0,*) ' quenching filed=', temp, ' for index =', n, ' invalid'
!       error = error + 1
       return  !!!!!
    endif
    modify(n)%q%id = id             
    modify(n)%q%a = a
    modify(n)%q%b = b 

  end subroutine epquenchModify

  subroutine epEminModify 
    implicit none
    real(8):: Egmin, Eemin, Recoilmin, KEmin, Enmin
    ! get Emin field
    pos = index(line, "/")
    if( pos == 0 ) then
       write(0,*) ' Emin line for entry #=',n
       write(0,*) ' has no /'
       error = error + 1
       return
    endif
 ! energy is alwasy GeV.  Electron : total energy
    Egmin = EminGsave     ! GeV
    Eemin = EminEsave     ! total  Energy GeV
    Recoilmin = RecoilEsave
    KEmin = KEminSave
    Enmin = Enminsave

    read(line(7:pos), *)  Egmin, Eemin, Recoilmin, KEmin, Enmin
!          Emin
!          23456
    modify(n)%Em%Egmin = abs(Egmin)   ! in GeV
    modify(n)%Em%Eemin = abs(Eemin)  ! Total GeV
    modify(n)%Em%Recoilmin =abs(Recoilmin)
    modify(n)%Em%KEmin = abs(KEmin)
    modify(n)%Em%Enmin = abs(Enmin)
         ! if some Emin< 0, imperative Emin, i.e, no residual range 
         ! is considered.
    modify(n)%Em%imperative = Egmin< 0. .or. Eemin < 0. .or. Recoilmin < 0.   &
                         .or.   KEmin < 0. .or. Enmin < 0.
    
    modify(n)%kind = IBSET(modify(n)%kind,bitEmin)

  end subroutine epEminModify

end module epModify
