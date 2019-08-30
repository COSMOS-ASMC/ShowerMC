!
!       manages angle of gamma at brems  
!       manages angle of pair for pair creation
!   ref: Rev. Mod. Phys. Tsai, vol.46,1974
module modBremPairAng
  implicit none
!   at test time
!  private 
!  public cBremAng, cPairAng, f1,f2,f3,f4,FF1,FF2,Me
!    at normal run
  private
  public cBremAng, cPairAng, cBremPairAngInit
  type brempairAng
     real(8):: A, Z
     real(8):: G2
     real(8):: Z2CC
     real(8):: sd 
     integer:: used =0
  end type brempairAng
  integer,parameter:: maxmedia=20   !  MAX_MEDIA in ZepMaxdef.h
  type(brempairAng),save:: bpaCnst(maxmedia)  !save due to gfort
  
  real(8),save:: XmCC
  real(8),save:: f1,f2,f3,f4, X
  real(8),save:: FF1, FF2
!
!       tmin is insensitive to l=(Ee/m teta)**2 as far as
!       electron pair createion or X value (via sc).
!       So we regard it 1 if  force  is .true. ,
!        irrespectively of other conditions
!
  logical,save:: force=.true.
  real(8),parameter:: Me=0.511d-3  ! Me in GeV
  real(8),parameter:: pi=3.14159265358979d0, hpi=pi/2
  real(8),parameter:: smallA=0.35  ! sin(smallA) ~ smallA ; ~2% error
  integer,save:: case,  reject  ! only at test time

contains
  function cCoulombC(z) result(ans)
    implicit none
    real(8),intent(in):: z ! ( Z/137 )**2

    real(8):: ans
  ! looks diff. from Tsai's formula but bit more accurate
    ans =  z * (1./(1. + z) +0.20206 -0.0369*z + 0.0083*z**2 - 0.002*z**3)
  end function cCoulombC

  function cBremtprim(Ee, Eg, lin) result(ans)
    ! eq. 3.81 of Tsai.  We regard l to be 1, since
    !     (theta~Me/Ee, l = theta**2 (Ee/Me)**2)
    ! and t' is not important factor.
    real(8),intent(in):: Ee, Eg  ! in GeV
    real(8),intent(in):: lin  ! if < 0, l=1 is used.
    real(8):: ans

    real(8),save::l 

    if( lin < 0. .or. force ) then
       l = 1.
    else
       l = lin
    endif
!      Eq. 3.81 (1+l)**2 is N.G
    ans =( Eg*Me**2*(1+l)/(2*Ee*(Ee-Eg) ))**2
  end function cBremtprim

  subroutine cBremPairAngInit(mno, Ain, Zin)
    implicit none
    integer,intent(in):: mno ! media # unique to this A,Z
    real(8),intent(in):: Ain, Zin
    if( mno > maxmedia)   then
       write(0,*) ' too many  media for cBremPairAngInit'
       write(0,*) ' mno must be <= ',maxmedia
       write(0,*) ' Enlarge maxmedia there'
       stop
    endif
    if( bpaCnst(mno)%used == 0 ) then
       bpaCnst(mno)%A = Ain
       bpaCnst(mno)%Z = Zin
       bpaCnst(mno)%G2 = Zin*(Zin+1)
       bpaCnst(mno)%Z2CC=2* Zin*Zin*cCoulombC( (Zin/137)**2 )
       bpaCnst(mno)%sd = 0.164/Ain**0.66666d0
       bpaCnst(mno)%used = 1
    endif
  end subroutine cBremPairAngInit


  function cBremAngX(mno, Ee,  Eg, lin) result(ans)
!       X of Eq. 3.76 in Tsai
    implicit none
    integer,intent(in):: mno
    real(8),intent(in):: Ee, Eg, lin
!       if lin < 0, it is neglected and l=1 is
!  used
    real(8):: ans
    real(8):: sc, sb
    real(8),save::l=1.0  ! we can regard this const.
    !  for sc, we also regard l to be 1 as in cBremtprim

    if( bpaCnst(mno)%used ==0 ) then
       write(0,*)  ' media no=',mno, ' for cBremAng/cPairAng'
       write(0,*) ' invalid '
       stop
    endif

    if(lin < 0. .or. force) then
       l =1
    else
       l = lin
    endif
    sc =(Me*(1+l))**2 /bpaCnst(mno)%sd
    sb = cBremtprim(Ee, Eg, l)/bpaCnst(mno)%sd
    ans =bpaCnst(mno)%Z**2* ( (1+2*sb)*log( (1+1/sb)/(1+1/sc)) -  &
         (1+sb/sc)*(1+2*sc)/(1+sc) )
!!!!
!    write(0,*) ' l is ', l, ' ans', ans, ' sb, sc',sb,sc
!!!!

  end function cBremAngX

  function cBremsAFunc(mno, Ee, Eg, l) result(ans)
    implicit none
    integer,intent(in):: mno  !one of  media no used at ...Init( )
    real(8),intent(in):: Ee, Eg, l
    real(8):: ans
    real(8):: temp, y

    y = Eg/Ee
    temp = (1-y)*2
    f1 = -temp/(1+l)**2
    
    f2 = 6*temp*l/(1+l)**4

    f3 = (temp+y**2)/(1+l)**2

    f4 = -2*temp*l/(1+l)**4

    FF1 = (f1+f2)*bpaCnst(mno)%G2
    X = cBremAngX(mno, Ee, Eg, l)
    XmCC = X -bpaCnst(mno)%Z2CC
    FF2 = (f3+f4)*XmCC
    ans = FF1 +FF2
  end function cBremsAFunc

  function cPairAFunc(mno, Eg, Ee, l) result(ans)
    implicit none
    integer,intent(in):: mno  !one of  media no used at ...Init( )
    real(8),intent(in):: Ee, Eg, l
    real(8):: ans
    real(8):: temp, x

    x = Ee/Eg
    temp = x*(1-x)*2
    f1 = temp/(1+l)**2
    
    f2 = -6*temp*l/(1+l)**4

    f3 =(1-temp)/(1+l)**2

    f4 = -2*temp*l/(1+l)**4

    FF1 = (f1+f2)*bpaCnst(mno)%G2
    X = cBremAngX(mno, Ee, Eg, l)
    XmCC = X -bpaCnst(mno)%Z2CC
    FF2 = (f3+f4)*XmCC
    ans = FF1 +FF2
  end function cPairAFunc

  subroutine cBremAng(mno, Ee, Eg, teta)
      implicit none
      integer,intent(in):: mno ! one of mno given to the ..Init routine
      real(8),intent(in)::Ee ! Energy of electron/positron GeV  
                         ! may be in other unit.
      real(8),intent(in)::Eg ! Brems gamma energy.  same unit as Ee
      real(8),intent(out):: teta ! sampled angle of photons relative to
                   !   the inciedent electron

      real(8):: y  !  Eg/Ee
      real(8):: C1, C2,  temp, u, el
      y = Eg/Ee
      X = cBremAngX(mno, Ee, Eg, -1.d0 )
      XmCC = X - bpaCnst(mno)%Z2CC
!!!
!      write(0,*) ' X=',X, ' XmCC=',XmCC
!!!
      temp = (1-y)*2
      ! coeff. for 1/(1+l)**2
      !  include integration of 1/(1+l)**2 dl = 1
      C1 = -temp*bpaCnst(mno)%G2 + (temp +y*y)*XmCC  ! 3.80
      ! coeff. for l/(1+l)**4
      C2 =( 6*bpaCnst(mno)%G2 -2*XmCC )*temp   ! 3.80
      !  include integration of l/(1+l)**4 dl =  1/6
      C2 = C2/6.
!!!!!
!      write(*,*) 'C1 C2=', C1, C2
!!!!
      do 
         call rndc(u)
         if(C2 <= 0. .and. C1 > 0. ) then 
                 ! most frequent.  Ek=0.01-->98.6%
                 !          Ek=1e-4   100 %
             ! For Pb, Ek=1,case1:  9998097 among 1e7 events
             !                  2:  5
             !                  3:  1783
             !                  4:  115
            call ksampForBPA1(el)  
            case=1
!               write(0,*) ' C2*el/(1+el)**2/C1=', C2*el/(1+el)**2/C1
!         we may safely neglect next. (time is almost not affected)
!         C2*el... is order of  -1.e-3. for Ek=1e-4
            if( u > 1.+ C2*el/(1+el)**2/C1 ) continue
         elseif( C1 > 0. .and. C2 > 0. ) then
            if( C1/(C1+C2) < u ) then
               ! select 1/(1+l)**2 dl
               ! 0.07 %  Ek=0.01.  v various
               !  Ek=1  v> 0.9999
               call ksampForBPA1(el) 
               case=2
            else
                ! 2nd frequent Ek=0.01  1.25 %
                !  Ek =1.  v>0.99
               call ksampForBPA2(el)
               case =3         
            endif
         elseif(C2 > 0. ) then
            call ksampForBPA2(el)  ! not happen yet
            case = 4
            if( u > 1. + C1*(1+el)**2/el/C2) continue
         else
            write(0,*) ' strange C1,C2',C1,C2
            write(0,*) 'Ee, Eg=', Ee, Eg
            stop
         endif
         !  el =( E/m theta)**2
         teta = sqrt( el )*Me/Ee 
         if( teta < smallA) exit
         if(teta > hpi) continue
!         this   distribution is in dl= teta dteta  but true one must be
!            sin(teta) dteta ;  teta dteta = sin(teta) teta/sin(teta) dteta
!            so sampled one must be accepted only fraction of  sin(teta)/teta
         call rndc(u)         
         if( u < sin(teta)/teta ) exit
      enddo
    end subroutine cBremAng
      
    subroutine cPairAng(mno, Eg, Ee,  teta)
      implicit none
!          samples polar angle of an electron or positron 
!        from  pair creation.  It must be a smaller eneregy
!        one.
      integer,intent(in):: mno ! one of mno given to ..Init
      real(8),intent(in)::Eg ! parent photon energy in GeV
      real(8),intent(in)::Ee ! Energy of electron or  
           !        positron GeV (smaller one)

      real(8),intent(out):: teta
       ! sampled angle of e-/or e+ 0 relative to 
       ! the inciedent photon . radian (<pi/2)

      real(8):: sx  !  Ee/Eg
      real(8):: C1, C2,  temp, u, el


      sx = Ee/Eg
      X = cBremAngX(mno, Ee, Eg, -1.d0)
      XmCC = X - bpaCnst(mno)%Z2CC
      temp = sx*(1-sx)*2
      ! coeff. for 1/(1+l)**2
      !  include integration of 1/(1+l)**2 dl = 1
      C1 = temp*bpaCnst(mno)%G2  +  (1-temp)*XmCC  ! Eq.3.5
!      C1 = temp*bpaCnst(mno)%G2  -  (1+temp)*XmCC  ! ??
      ! coeff. for l/(1+l)**4
      C2 =(-6*bpaCnst(mno)%G2 -2*XmCC )*temp   ! Eq.3.5
!      C2 =(6*G2 +2*XmCC )*temp   !  ??
      !  include integration of l/(1+l)**4 dl =  1/6
      C2 = C2/6.
!!!!!
!      write(*,*) 'C1 C2=', C1, C2
!!!!
      if(C1 <= 0. ) then
!       both of C1 and C2 < 0. for Eg<1.75MeV or so.
!        work around 
         C1=0.1
      endif
!!!!
      reject = -1
      do 
         call rndc(u)
         if(C2 <= 0. .and. C1 > 0. ) then
            call ksampForBPA1(el)
            case =1
            reject = reject + 1  ! this is almost 0, 4 is rare
            if( u > 1.+ C2*el/(1+el)**2/C1 ) continue
         elseif( C1 > 0. .and. C2 > 0. ) then  ! never happen
            if( C1/(C1+C2) < u ) then
               ! select 1/(1+l)**2 dl
               case=2
               call ksampForBPA1(el)
            else
               case=3
               call ksampForBPA2(el)
            endif
         elseif(C2 > 0. ) then  ! never happen
            case = 4
            call ksampForBPA2(el)
            if( u > 1. + C1*(1+el)**2/el/C2) continue
         else  ! work around prevents to come here
            write(0,*) ' strange C1,C2',C1,C2
            write(0,*) 'Ee, Eg=', Ee, Eg
            stop
         endif
         !  el =( E/m theta)**2
         teta = sqrt( el )*Me/Ee 
         if( teta < smallA) exit
         if(teta > hpi) continue
!         this   distribution is in dl= teta dteta  but true one must be
!            sin(teta) dteta ;  teta dteta = sin(teta) teta/sin(teta) dteta
!            so sampled one must be accepted only fraction of  sin(teta) / teta
         call rndc(u)         
         if( u < sin(teta)/teta ) exit
      enddo
    end subroutine cPairAng

  end module modBremPairAng

! !   test  cPairAFunc. draw function
! program main
!   use modBremPairAng
!   implicit none
!   real(8):: Ee, Ek, Eg, v, vmin, vmax
!   real(8):: teta, l, fval, u
!   real(8):: Ain=14.55, Zin=7.266
! 
!   call cBremPairAngInit(1, Ain, Zin)
! 
! !  do 
!   write(0,*) 'Enter Gamma E in GeV (>2Me=1.023e-3)'
!   read(*,*) Eg
!   vmin = Me/Eg
!   vmax = (Eg-Me)/Eg
!   v = vmin
!   do while( v < vmax ) 
!      Ee = v*Eg
!      if( Ee > Me ) then
!         l = 1d-5
!         do while  (l < 40.) 
!            fval = cPairAFunc(1, Eg, Ee, l)
!            if( fval > 0 ) then
!               write(*,'(1p,9g14.4)') &
!                 l, fval, v, f1, f2, f3, f4,FF1, FF2
!            endif
!            l = l * 10.d0**0.01
!         enddo
!      endif
!      v = v +0.0025
!   enddo
!end program main

!!   test  cBremAng
!program main
!  use modBremPairAng
!  implicit none
!  real(8):: Ee, Ek,  Eg
!  real(8):: vmin=1.d-5, vmax, v, u
!  real(8):: teta
!  integer::i
!  real(8):: Ain=14.55, Zin=7.266
!!  real(8):: Ain=207, Zin=82
!!
!  call cBremPairAngInit(1, Ain, Zin)
!  call cBremPairAngInit(2, 207.d0, 82.d0)
!
!  Ek = 1.
!!!  do 
!  write(0,*) 'Enter Ek '
!  read(*,*) Ek
!  Ee = Ek + Me
!  vmax =Ek/Ee
!  do i = 1, 1000000 
!!!          simple test by dv/v gamma spectrun
!     call rndc(u)
!     v = vmin*(vmax/vmin)**u
!!!!!!!!!!!
!!     v = 0.01  !!!!!!!!!!
!!!!!!!!!!!!!
!     Eg = v*Ee
!     if(Ee > Me .and. Eg < Ee-Me) then
!!           write(0,*)'i Ee,Eg v=',i, Ee, Eg, v
!        call cBremAng(1,  Ee, Eg, teta) 
!        write(*,'(1p,4g15.6,i3)') Ee, Eg, teta, (teta*Ee/Me)**2,1, case
!        call cBremAng(2,  Ee, Eg, teta) 
!        write(*,'(1p,4g15.6,i3)') Ee, Eg, teta, (teta*Ee/Me)**2,2, case
!!        write(*,'(1p,4g15.6)') Ee, Eg, teta, (teta*Ee/Me)**2
!!           !        write(0,*) ' test=', teta
!     else
!        stop
!     endif
!  enddo
!!!  end do
!end program main

!   test    cPairAng  
!program main  
!  use modBremPairAng
!  implicit none
!  real(8):: Ee,  Eg
!  real(8)::  v, u, vmax, vmin
!  real(8):: Ain, Zin
!  real(8):: teta
!!        at test time when "case" and "reject" are used
!!        next maybe comment outed
!  real(8):: Me=0.511d-3
!  integer::i
!  Ain = 14.55
!  Zin = 7.266
!  Eg = 1.
!  do 
!     write(0,*) 'Enter Eg > 1.022e-3'
!     read(*,*) Eg
!     if(Eg <= 0 ) then
!        stop
!     endif
!     vmax = (Eg-Me)/Eg
!     vmin = Me/Eg
!     call cBremPairAngInit(1, Ain, Zin)
!     do i = 1, 1000000 
!!          sipmle test by vdv
!        call rndc(u)
!        v = (vmax - vmin)*u + vmin
!        Ee = v*Eg     ! inc mass
!        Ee = min(Ee, Eg-Ee)
!!           write(0,*)' Ee,Eg=',Ee, Eg
!        call cPairAng(1, Eg, Ee, teta)
!!        write(*,'(1p,4g15.6,i3,i5)')  &
!!             Eg, v, teta, (teta*Me/Ee)**2, case, reject
!        write(*,'(1p,4g15.6)')  &
!             Eg, v, teta, (teta*Me/Ee)**2
!     enddo
!  end do
! end program main
