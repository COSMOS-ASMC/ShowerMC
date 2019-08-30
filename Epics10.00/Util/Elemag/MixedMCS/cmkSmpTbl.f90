! Fix  c1forHardScat below (D=0.1)
! If you need to print reduced pdf, make  
!     print_redpdf=.true.
! In usual run for mainkg sampling table, avoid print it.
!
!
!  When we make sampling tables, we work on grid points as
!  to energy and mu.  (some parts are
!  usable for non grid case, but they are comment out)
module modcDCSf
  use modTPXS
  use modDCS
  implicit none
  real(8),save:: c1forHardScat=0.1   ! 0.05 ~ 0.2 (for H.E,could be
                            !  made to be  >0.2 ?
  real(8),save:: maxHCSmfprl=0.05  ! hard cs m.f.p is kept smaller than this
  logical,save:: print_redpdf=.false. 
 
  type(DCSconst),pointer:: dcsname  ! assign some real DCSconst name like
  ! 
  !    dcsname => example
  ! where 
  !    type(DCSconst),target:: example 
  ! must have been declared.
  !   see an example  mkSampTbl/test.f90
  ! 
  real(8):: Ein  ! Energy of e-/e+ in eV
  integer:: Eindex ! index of KEele(:)
  real(8),parameter::fpi = 4*3.14159265358979323d0
  real(8),save:: adjfac=1.0   ! not used now
  real(8):: S0, S1, S0p
  real(8):: A,  tcs, logKE, dcsval, dcsW, intedsdmu
  real(8):: cm2topgrm
  character(len=2),save:: negaOrposi  ! "e-" or "e+" should be put
 ! to get better integration value of xs
 ! before integration sset adjfac= 1.0d-15 and
 ! f=xs/adjfac is integrated and later
 ! it is recovered to 1.0
contains
function credPDF(u, Ee) result(ans)
!  compute reduced pdf.  
!  get mu(= (1-cos)/2), which gives correct angle distribution
!  of Wenztel (=dcsW). Then, dcs/dcsW becomes almost flat
!  (eaaxct numerical dcs should be normalzied).
!  This gives  dcs/norm/dcsW.
!  (norm = intedcs)
  implicit none
  real(8),intent(in):: u   ! any value in  (0,1)
     !  but is related to a uniform random number.
  real(8),optional,intent(in):: Ee ! if given, e+/e- energy
                   ! in eV.   Otherwise, energy is on the grid
                   ! specified by index, Eindex 
  real(8):: ans

  real(8)::dcsv, muin  ! ds/dmu at current u

  muin = A*u/ (A+1-u)
  dcsW = A*(A+1)/(A+muin)**2

  if( PRESENT(Ee) )  then
     call cDCS(dcsname, muin, Ee, dcsv)
  else
     if(u == 0.) then
        dcsv = dcsname%dcs(1,Eindex)
     elseif(u == 1.0d0) then
        dcsv = dcsname%dcs(nmu,Eindex)
     elseif( u<1.0d0 .and. u>0.d0 ) then
        call cDCSEgrid(dcsname, muin, Eindex, dcsv)
     else
        write(0,*) 'u =',u,' invalid for credPDF'
        stop
     endif
  endif
  ans = dcsv/intedsdmu/dcsW
!!!!!!!!!!!!
  if( print_redpdf ) then
     write(*,'(a,a,  1p, 6E15.7)') &
          'redpdf ',negaOrposi, Ein,  u,  ans,  dcsv, intedsdmu, dcsW 
  endif
!!!!!!!!!!!!
end function credPDF
end module modcDCSf

module modSampRedPDF
  use modTPXS
  use modMCS0
  implicit none
  real(8):: sampledV(usize, nEneg)
  real(8):: redpdf(usize)
  real(8):: uNorm
  real(8):: uuNorm
  real(8):: lambdah(nEneg), lambdas1(nEneg), lambdas2(nEneg)    
  real(8):: muc(nEneg)
  real(8):: minNon0mucE
  integer:: minNon0mucEindex
end module modSampRedPDF

!  make sampling table of dcs.
!  the table can be used for all mu as well as
!  for hard scattering with mu >muc
subroutine cmkSmpTab(Dname, Tname, norp)
  use modcDCSf
  use  modSampRedPDF
  implicit none
  type(TPXSconst),intent(in):: Tname
  type(DCSconst),target,intent(in):: Dname
  character(len=2),intent(in):: norp  ! "e-" or  "e+"
  integer::n
  real(8),external:: cInteDCS
  integer:: iu, ie
  integer::icon

  n = Tname%n   ! same as nEneg but may be smaller for e+=nEpos  
            ! (not used yet).
  negaOrposi = norp  ! set in modcDCSf
!      dcsnmae is used as a general name

  dcsname => Dname
!     fix E independent consts: get uarray and cm2topgrm
!     (cm2 ==> g/cm2 conversion factor)
  call cinimkSmpTab

  minNon0mucE = 0.
  do ie= 1, nEneg
     Eindex = ie
     Ein = KEele(ie)
     logKE = log(Ein)
!      tcs: total xs.(int ds/dmu from mu=0 to 1).  To compare with
!      S0, 4pi must be multiplied.

     intedsdmu =   cInteDCS(1.0d0)
     tcs = intedsdmu*fpi
     S0 = Tname%S0(ie)  
!         if non grid energy, use next to get S0 by interpolation
!        (as Elsepa value).
!El     call kcsplIntp(logKEele, Tname%logS0, n, Tname%coefS0, n-1,logKE, S0p)
!El     S0p=exp(S0p)


!    tcs and S0 (S0p) should be the same but Elsepa given value (S0)
!    and tcs sometimes differ by < 1%. due to numerical integration
!    method.  At 1GeV or near 1GeV, the diff. is largest.

!    write(*,'(a,1p,3E14.4)') &
!          'S0 ',S0, S0p,  tcs

!        If difference is > 2 %, issue warning 
     if(abs( S0/tcs - 1.d0) > 0.02 ) then
        write(0,*)'*****************************warning'
        write(0,*) ' E=', Ein
        write(0,'(a,1p,3E14.4)') &
          'S0 ',S0, S0p,  tcs
        write(0,*) 'abs( S0/tcs - 1.d0)='
        write(0,*) abs( S0/tcs - 1.d0)
     endif


!        we employ A0 so that average of
!           g(mu)dmu=A0(A0+1)/(A0+mu)^2dmu (=<mu>) (dcsW; Wentzel
!           dcs). 
!        is numerical DCS's <> = S1/S0/2
!    This has been  already gvien by name%A0
!    if we use Ein not on grid, we may employ next to get A
!!El    call kcsplIntp(KEele, name%A0, n, name%coefA0, n-1, Ein, A)
!!E     call kcsplIntp(KEele, name%S1, n, name%coefS1, n-1, Ein, S1)

     A =  Tname%A0(ie)
!     A =  Tname%A(ie)  ! this makes some improvemnt ? buf
!        not so obvious. so use A0 as A.
!     write(*,'( a, 1p , 4g14.4)') &
!        'EA ', Ein, Tname%A0(ie), Tname%A(ie), Tname%S0(ie)



!    Inte 0 to mu of g(mu)dmu =
!          (A+1) mu/(A+mu)
!   To sample mu following g(mu), we may solve
!         (A+1)mu/(A+mu) = u
!  where u is  uniform random unumber in (0~1).  I.e
!       mu = A*u/(A+1-u)
!
! In analog to g(mu),      
! now we want to find mu for the f(mu) instead of g(mu) where
! f(mu) is numerial DCS.  The approximate solution would be mu
! for g(mu). 
!   
! We now convert variable from mu to u where 
!    u=  (A+1)mu/(A+mu) . As a result normalzed NPDF is
!   coverted to new pdf:   NPDF/tcs/dcsW
!  of which variation with energy is much smaller than NPDF.
!  Manipulate such reduced pdf.
         ! make sampling table in sampleV
     call cManipRedPDF
!       for hinge sampling, get lambda(h) --> stored in module modSampRedPDF
     call cgetLambdah(Tname, icon)
     if(icon == 0 .and. minNon0mucE==0.) then
        if( Ein > 1.0d3 ) then  ! we disregard above cond. at low E
           minNon0mucE = KEele(Eindex)
           minNon0mucEindex = Eindex
        else        
           icon = 1  ! and force soft 
        endif
     endif
!      and get muc; stored in modSampRedPDF
     call cgetMuc(icon)
!      using muc, get lambda1s, lambda2s -->  stored in module modSampRedPDF
     call cgetResTPXS
  enddo

  write(*,'(i3, 1p,E15.5,a,a)') &
    minNon0mucEindex, minNon0mucE, ' Min Non 0 muc idx& Energy ----- ',negaOrposi
  write(*,'(a,a)') 'E muc lh ls1 ls2 -------------- ',negaOrposi
  do ie= 1, nEneg
     write(*,'(1p, 5E16.7, 0p, i2)') &
          KEele(ie), muc(ie), lambdah(ie), lambdas1(ie), lambdas2(ie)
  enddo
  write(*,*)

  call coutSmpTabRedPDF
end subroutine cmkSmpTab


subroutine cinimkSmpTab
  use modSampRedPDF
  implicit none

  call ciniSmpTab

  call ciniGetMuc  ! obtain cm2topgram.

!  call coutSmpConst

end subroutine cinimkSmpTab

subroutine coutSmpConst
  use modSampRedPDF
  implicit none
  integer::i
  write(*,'(a)') 'uarray ---------- '
  write(*,'(10f11.7)')  (uarray(i), i = 1, usize)
  write(*,*) 
end subroutine coutSmpConst

subroutine cmkRedPDFTab
! npdf's variable mu is changed to u = (A+1)mu/(A+mu)
!  This leads to the reduced pdf: 
!      pdf= NPDF/Wpdf  (NPDF here is normalized
!  one:  ds/dmu/tcs)
!  
  use modcDCSf
  use modSampRedPDF
  implicit  none
  integer::iu

  real(8), external:: cinteRedPDF  ! integral of reduced pdf

  do iu = 1, usize
     redpdf(iu) = credPDF(uarray(iu))
  enddo

  uNorm = cinteRedPDF(1.d0)
end subroutine cmkRedPDFTab


function cinteRedPDF(u) result(ans)
!     integral of reduced pdf from 0 to u.
!     use table redpdf and trapezoidal integration
  use modcDCSf

  use modSampRedPDF
  implicit none
  real(8),intent(in):: u    ! integration from 0 to u
  real(8):: ans
  real(8):: error
  integer:: icon
  real(8),parameter:: ub=1.1d0
  real(8):: u1, ans1, ans2

  call ktrpzIntT2(redpdf, 1, usize, uarray, 1,  0.0d0, u, ans)

end function cinteRedPDF

subroutine  cManipRedPDF
  implicit none
!     make reduced  pdf table corresponding to uarray
  call cmkRedPDFTab

!     make sampledV(i,j) for i= uarray index,
!                              j= KEele index
!       v= sampledV(i,j)  is  (1-cos)/2
  call cmkSmpTabRedPDF
end subroutine cManipRedPDF

     

subroutine cwritendcs
  use modcDCSf
  implicit  none
!        compute N.PDF of DCS; cDCSf is ore an one argument function
!        while cDCS/cDCSEgrit  is subroutine with argments. 
!        dcsname must have been fixed.
  integer:: i
  real(8):: dcsv
  real(8),external:: cDCSf

  do i= 1, nmu
     dcsv = cDCSf(muval(i))
     write(*,'(a, a, 1p,3E14.4)') &
          'dcsv ', negaOrposi, Ein, muval(i), dcsv
  enddo
end subroutine cwritendcs

function cDCSf(mu) result(ans)
!        see comment in modcDCSf ; other 
!      preparation by calling cmkDCSconst must have been done.
  use modcDCSf
  implicit none
  real(8),intent(in):: mu(2)
  real(8):: ans
      
!  call cDCS(dcsname, mu(1),  Ein,  ans)
  call cDCSEgrid(dcsname, mu(1),  Eindex,  ans)
  ans = ans/adjfac  ! adjfac is != 1 when used from cInteDCS
end function  cDCSf

function cInteDCS(mumax) result(ans)
  use modcDCSf
  implicit none
  real(8),intent(in):: mumax    ! integration from 0 to mumax
  real(8):: ans
  real(8),external:: cDCSf
  real(8):: error
  integer:: icon
!  k16pGaussLeg is too bad
!                                    1.d-8 is still not enough
!  call kdexpIntF(cDCSf, 0.d0, mumax, 1.d-10, ans, error, icon)
!
!  with fine mesh mu's, next one seems accurate and  faster.
  call ktrpzIntT2(dcsname%dcs(1,Eindex), 1, nmu, muval, 1,  &
       0.0d0, mumax, ans)
  ans = ans*adjfac  ! adjfac is always 1 now.
  adjfac =1.0
end function cInteDCS


subroutine cmkSmpTabRedPDF
  use modcDCSf
  use modSampRedPDF
  implicit none
  real(8):: vl=0., vr=1.0
  integer:: icon
  real(8),save:: epsilon=1.d-13
  integer::i
  real(8)::u, v, vin
  real(8),external::  cSolveRedPDFf, cinteRedPDF

  vin = uarray(2)
  do i = 2, usize-1
     u = uarray(i)
     ! solve u*uNorm = cinteRedPDF(v)
     uuNorm =u*uNorm
     vl = 0.d0
     if(csolveRedPDFf(vin) > 0.) then
        vr = vin
     else
        vr = min(vin*2, 1.d0)
        do while(csolveRedPDFf(vr) < 0. ) 
           vr = min(vr + 1.d-3, 1.d0)
        enddo
     endif
     call kbinChop(csolveRedPDFf, vl, vr, vin, epsilon, v, icon)
     if( icon /= 0 ) then
        write(0,*) 'ERR csolveRedPDF; icon= ',icon 
        do while (vl < 1.0 )
           write(0,*) vl, csolveRedPDFf(vl)
           vl = vl + 1.e-3
        enddo
        stop  
     endif
     vin =  v
     sampledV(i, Eindex) = v
  enddo
  sampledV(1,Eindex) = 0.
  sampledV(usize,Eindex) = 1.
end subroutine cmkSmpTabRedPDF


function cSolveRedPDFf(v) result(ans)
  use modSampRedPDF
  implicit none
  real(8),external:: cInteRedPDF
  real(8),intent(in):: v
  real(8):: ans
  
  ans = cInteRedPDF(v)/uuNorm - 1.0d0

end function cSolveRedPDFf

subroutine coutSmpTabRedPDF
  use modTPXS
  use modcDCSf
  use modSampRedPDF
  implicit none
  integer::i,j

!  write(*,'(1p, 10E14.5)') (KEele(i), i = 1, nEneg)
!  write(*,*)

  write(*,'(a,a)') 'sampling tab ------------- ', negaOrposi
  do j= 1, nEneg
     write(*,'(1p,10E14.6)') &
         (sampledV(i, j), i = 1, usize)
  enddo
  write(*,*) 
end subroutine coutSmpTabRedPDF
