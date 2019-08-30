module modBremAngTsai
  implicit none
  logical,save:: first=.true.
  integer,parameter:: nv=21, nE=6
  real(8),save:: v(nv)
  real(8),save:: Ek(nE), Log10Ek(nE)
  real(8),save:: p(nv, nE), q(nv, nE), r(nv, nE), s(nv, nE)
  real(8),save:: ps, qs, rs, ss  ! save fixed p,q,r,s
contains
  subroutine epReadTsaipqrs(io)
!    implicit none
    integer,intent(in)::io ! temporary io number for disk input
    integer:: icon
    character(100):: line
    integer:: nvc, nEc
    integer,external:: kgetenv2

    if( kgetenv2("EMDATA",line) == 0 ) then
       ! before EMDDATA is made
       call copenf(io, "$EPICSTOP/Data/BremPair/TsaiBrem.coef",icon)
       if(icon /= 0 ) then
          write(0,*) 'The file "TsaiBrem.coef" could not be opened '
          write(0,*) 'It should be in "$EPICSTOP/Data/BremPair/"'
          stop
       endif
    else
       call copenf(io, "$EMDATA/BremPair/TsaiBrem.coef", icon)
       if( icon /= 0  ) then
          write(0,*) 'The file "TsaiBrem.coef" could not be opened '
          write(0,*) 'It should be in "$EMDATA/BremPair/"'
          stop
       endif
    endif

    nvc = 1
    nEc = 1
    line = " "
    do
       read(io, '(a)') line

       if( line(1:1) == '#' ) cycle
       if( line(1:5) == '     ' ) cycle
          

       read(line,* ) Ek(nEc),  v(nvc), &
            p(nvc,nEc), q(nvc,nEc), r(nvc,nEc), s(nvc,nEc)

       if( nvc  == nv ) then
          nvc = 1
          if( nEc == nE ) exit
          nEc = nEc + 1
       else
          nvc = nvc + 1
       endif
    enddo
    Ek(:) = Ek(:) /1d3   ! to GeV
    v(1) = 0.            ! make exactly line
    v(nv) = 1.           !     thses
    close(io)
    
  end subroutine epReadTsaipqrs
end module modBremAngTsai

subroutine epSmpBremAngTsaiFE(Ee, Eg, teta)
  use modBremAngTsai
  implicit none
#include "Zmass90.h"
  
  real(8),intent(in):: Ee ! electron total energy in GeV.  >0.9d-3.
  real(8),intent(in):: Eg ! emitted brems photon energy in GeV, < Ee-masse

  real(8),intent(out):: teta ! sampled photon polar angle relative to
  ! electron direction before brems.  in rad  < pi.
  real(8):: Ekv, vv, g, gteta
  real(8):: pv, qv, rv, sv
  real(8),parameter:: pi=asin(1.d0)*2

  Ekv= Ee - masele
  vv = Eg/Ekv
  g = Ee/masele
  call epGetTsaipqrs(Ekv, vv, pv, qv, rv,sv)
  do 
     call epSmpBARF_m(pv, qv, rv, sv, gteta)
     teta = gteta/g
     if( teta < pi) exit
  enddo
end subroutine epSmpBremAngTsaiFE

subroutine epGetTsaipqrs(Ekv, vv, pv, qv, rv,sv)
  use modBremAngTsai
  implicit none
  real(8),intent(in):: Ekv ! electron kinetic energy  in GeV
  real(8),intent(in):: vv  ! Brems photon fractional  energg.
  real(8),intent(out):: pv, qv, rv, sv ! coefficient for (Ekv, vv) air
  
    !  photon energy is Ekv*vv

  real(8):: EkvLog  ! log10(Ekv)

  if( first ) then
     call  epReadTsaipqrs( 11 )
     first = .false.
  endif
       
    
  if( Ekv > 100.d-3) then
     EkvLog=-3
  elseif( Ekv > 0.3162d-3) then 
     EkvLog = log10(Ekv)
  else
     write(0,*) 'Ek =', Ekv, ' < 0.4 MeV: too small for Tsai'
     stop
  endif
!!!!!!!!
  
  call k4ptdi(p, nv, nE, nv, 0.d0, -3.5d0,  0.05d0, 0.5d0, vv, EkvLog, ps)
  call k4ptdi(q, nv, nE, nv, 0.d0, -3.5d0,  0.05d0, 0.5d0, vv, EkvLog, qs)
  call k4ptdi(r, nv, nE, nv, 0.d0, -3.5d0,  0.05d0, 0.5d0, vv, EkvLog, rs)
  call k4ptdi(s, nv, nE, nv, 0.d0, -3.5d0,  0.05d0, 0.5d0, vv, EkvLog, ss)
  pv = ps
  qv = qs
  rv = rs
  sv = ss
end subroutine epGetTsaipqrs

subroutine epqTsaipqrs(pv,qv,rv,sv)
  !  returns currently used p,q,r,s values
  use modBremAngTsai
  implicit none
  real(8),intent(out):: pv, qv, rv, sv
  pv = ps
  qv = qs
  rv = rs
  sv = ss
end subroutine epqTsaipqrs
  
subroutine epTsaiFuncValues(maxra)
  use modBremAngTsai
  ! using currently assigned p,q,r,s,
  ! output x, m(x), m1(x),m2(x) function values upto  x= maxra
  ! where x is reduced andgle (teta*Ee/me)

  real(8),intent(in) ::  maxra !max reduced angle . if 0, 6 is used.

  real(8),save:: dx = 0.01
  real(8)::x, peak, m1, m2, m, maxv

  real(8),external:: epBARF_m1, epBARF_m2
  real(8):: limitx

  if( maxra <= 0. ) then
     limitx = 6.0
  else
     imitx = maxra
  endif
  
  x=0.
  maxv = 0.
  do while( x <= limitx )
     m1 = epBARF_m1(ps,qs, x)
     m2 = epBARF_m2(rs,ss, x)
     m = m1  + m2
     if(m > maxv ) then
        maxv = m
        peak = x
     endif
     write(*,'(1p, 4g14.4)') x, m, m1, m2
     x = x + dx
  enddo
  write(*,'(a, 1p, 2g14.4)') '# (peak pos, peak value) ', peak, maxv
end subroutine epTsaiFuncValues




