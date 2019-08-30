!   epSmpBremAng2BN:  sample brems photon angle when the kinetic energy
!               of electron is Ek < 5.6 MeV.  (actual use should be < ...MeV)
!    Internally, the routines for this are classified further as
!        follows: (v= Eg/Ek)
!                        
module modBremAng2BN
  implicit none
  logical,save:: first=.true.
  !    v=0, 0.05, 0.1...   0.90 0.95  1.0
  !    Ek=0.01 MeV to 5.623 MeV step 10^0.25
  integer,parameter::  nE=12, nv=21
  real(8),save:: p(nv, nE), q(nv, nE), r(nv, nE), s(nv, nE)
  real(8),save:: ps, qs, rs, ss  ! save fixed p,q,r,s
  real(8),save:: EkvLog  ! log10(Ekv)
  character(100):: line
contains
  subroutine epRead2BNcoef(io)
!    implicit none
    integer,intent(in)::io ! temporary io number for disk input

    integer:: icon    

    real(8):: Ekv, vv 
    integer:: nvc, nEc
    integer,external:: kgetenv2

    if( kgetenv2("EMDATA",line) == 0 ) then
       ! before EMDDATA is made
       call copenf(io, "$EPICSTOP/Data/BremPair/2BN.coef",icon)
       if(icon /= 0 ) then
          write(0,*) 'The file "2BN.coef" could not be opened '
          write(0,*) 'It should be in "$EPICSTOP/Data/BremPair/"'
          stop
       endif
    else
       call copenf(io, "$EMDATA/BremPair/2BN.coef", icon)
       if( icon /= 0  ) then
          write(0,*) 'The file "2BN.coef" could not be opened '
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
          

       read(line,* ) Ekv, vv, &
            p(nvc,nEc), q(nvc,nEc), r(nvc,nEc), s(nvc,nEc)
       if( nvc  == nv ) then
          nvc = 1
          if( nEc == nE ) exit
          nEc = nEc + 1
       else
          nvc = nvc + 1
       endif
    enddo
    close(io)
  end subroutine epRead2BNcoef
  
end module modBremAng2BN

subroutine epSmpBremAng2BN(Ee, Eg, teta)
!  use modBremAng2BN
  implicit none
#include "Zmass90.h"
  
  real(8),intent(in):: Ee ! electron total energy in GeV. <= 5.6d-3.
  real(8),intent(in):: Eg ! emitted brems photon energy in GeV, < Ee-masse

  real(8),intent(out):: teta ! sampled photon polar angle relative to
  ! electron direction before brems.  in rad  < pi.

  real(8),parameter:: pi=asin(1.0d0)*2

  real(8):: Ekv
  real(8):: pv, qv, rv, sv, v, u

  Ekv=max( Ee - masele, 10.0d-6)  ! for small Ekv us 10keV data

  if( Ekv > 5.623d-3) then
     write(0,*) ' Ek=', Ekv, ' is too high for epSmpBremAng2BN '
     write(0,*) ' use epSmpBremAngTsai '
     stop
  else
     v = Eg/Ekv
     call epGet2BNpqrs(Ekv, v, pv, qv, rv,sv)
     do
        call epSmpBARF_m(pv, qv, rv, sv, teta)
        if( teta > pi) cycle
        if( Ekv > 316.270d-6) exit
           ! weight factor (1-(teta/pi)**8)
        call rndc(u)
        if( u < (1-(teta/pi)**8) ) exit
     enddo
  endif
end subroutine epSmpBremAng2BN


subroutine epGet2BNpqrs(Ekv, v, pv, qv, rv,sv)
  use modBremAng2BN

  implicit none
  real(8),intent(in):: Ekv ! electron kinetic energy  in GeV
  real(8),intent(in):: v   ! Brems photon fractional  energg.  Eg/Ekv
  real(8),intent(out):: pv, qv, rv, sv ! coefficient for (Ekv, vv) air
  
  !  photon energy is Ekv*v
  
  real(8),parameter:: Elogmin=log10(0.01d-3), Elogstep=0.25d0
  real(8),parameter:: vmin= 0.d0, vstep=0.05d0
  
  if( first ) then
     call  epRead2BNcoef( 11 )
     first = .false.
  endif
       
  
  if( Ekv < 0.01d-3) then
     write(0,*) ' should not come here: E=', Ekv,' < 0.01d-3 '
     stop
  else
     EkvLog = log10(Ekv)
     call k4ptdi(p, nv, nE, nv, vmin, Elogmin, vstep, Elogstep, &
          v, EkvLog, ps)
     call k4ptdi(q, nv, nE, nv, vmin, Elogmin, vstep, Elogstep, &
          v, EkvLog, qs)
     call k4ptdi(r, nv, nE, nv, vmin, Elogmin, vstep, Elogstep, &
          v,  EkvLog, rs)
     call k4ptdi(s, nv, nE, nv, vmin, Elogmin, vstep, Elogstep, &
          v, EkvLog, ss)
     pv = ps
     qv = qs
     rv = rs
     sv = ss
  endif
end subroutine epGet2BNpqrs

subroutine epq2BNpqrs(pv,qv,rv,sv)
  !  returns ccurrently used p,q,r,s values
  use modBremAng2BN
  implicit none
  real(8),intent(out):: pv, qv, rv, sv
  pv = ps
  qv = qs
  rv = rs
  sv = ss
end subroutine epq2BNpqrs
  
subroutine ep2BNFuncValues
  use modBremAng2BN
  ! must be used for 1d-3> Ek>0.056d-3
  ! using currently assigned p,q,r,s,
  ! output x=teta, m(x), m1(x), m2(x) function values; x in 0~pi

  real(8),save:: dx = 0.01
  real(8)::x, peak, m1, m2, m, maxv

  real(8),external:: epBARF_m1, epBARF_m2
  real(8),parameter:: pi=asin(1.0d0)*2
  real(8)::wf
  x=0.
  maxv = 0.
  do while( x < pi )
     wf = (1.-(x/pi)**8) 
     m1 = epBARF_m1(ps,qs, x) * wf 
     m2 = epBARF_m2(rs,ss, x) * wf
     m = m1  + m2 
     if(m > maxv ) then
        maxv = m
        peak = x
     endif
     write(*,'(1p, 4g14.4)') x, m, m1, m2
     x = x + dx
  enddo
  write(*,'(a, 1p, 2g14.4)') '# (peak pos, peak value) ', peak, maxv
end subroutine ep2BNFuncValues





