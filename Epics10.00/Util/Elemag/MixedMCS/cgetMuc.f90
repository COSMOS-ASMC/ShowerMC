  ! get muc for a gvien lambda(el,h) so that
  !  N*lambda(el,h)*sigma(el,h)=1; where sigma(el,h) is 
  !  integration of 4pids/dmu from muc to 1.
  !  ie  tcs - inte(0,muc); tcs = inte(0,1) of 4pids/dmu
  ! dcsname must be fixed before calling this
subroutine cgetMuc(icon)
  use modcDCSf
  use modSampRedPDF
  implicit none
  integer,intent(in)::icon !  if 0, muc>0 must be
              ! obtained here 
              ! else muc = 0
!  real(8),intent(out)::muc  ! obtaned mu (= (1-cos)/2)
!   above which we regard hard collision.  This is stored in
! muc array in modSampRedPDF
! To get muc, we solve (tcs-4pi*InteDCS(muc) )*N= 1/lambdah
! i.e,  
!     (tcs- 4pi*InteDCS(muc)*N*lambdah -1 =0.
!xs is in cm^2 and lambdah is g/cm^2. 
!  xs can be converted to g/cm2 by xs*cm2topgrm 
!  so the equation above becomes
!     (tcs- 4pi*InteDCS(muc)*cm2topgrm = lambdah
!  i.e
!     (tcs- 4pi*InteDCS(muc)*cm2topgrm * lambdah -1.0 =0. 
  real(8),external:: cLambdah2muc

  if(icon == 0) then
     call csolvehard(cLambdah2muc,  muc(Eindex))
  else
     muc(Eindex) = 0
  endif
end subroutine cgetMuc

function cLambdah2muc(mu) result(ans)
  use modcDCSf
  use modSampRedPDF
  implicit none
  real(8),intent(in):: mu
  real(8):: ans

  real(8),external:: cInteDCS
!!!!!!!!!
!  real(8):: temp
!  if(Eindex == 106 ) then
!     temp = cInteDCS(mu)
!     write(0,*) ' tcs=', tcs, ' mu=',mu,  ' cintedcs=', temp
!     write(0,*) 'cm2topgrm=', cm2topgrm
!     write(0,*) ' lambdah=', lambdah(Eindex)
!     write(0,*) ' tcs - fpi*cInteDCs=', tcs -fpi*temp
!     write(0,*) ' // x cm2topgrm==', (tcs -fpi*temp)*cm2topgrm
!  endif
!!!!

  ans = (tcs-fpi*cInteDCS(mu))*cm2topgrm * lambdah(Eindex) - 1.0
!!!!
!  if(Eindex == 106 ) then
!     write(0,*) ' ans =',ans 
!  endif
!!!
end function cLambdah2muc


subroutine cgetLambdah(Tname, icon)
  use modcDCSf  
  use modSampRedPDF
  implicit none
  type(TPXSconst),intent(in):: Tname
!     where lambdael is 1/(N*S0), lambda1 is 1/(N*S1)
!     c1=c1forHardScat which is ~ 0.1-0.05.
!     lambdah = max( lambda1*c1, lambdael)
  integer,intent(out):: icon ! 0 if lambda1*c1 is > lambdael
                         !   else 1 .
!     lambdah is put in moudlue  modSampRedPDF  (in g/cm2)

  real(8):: temp, temp1, temp2, g2rl, g2cm


  call getConvFac(g2rl, g2cm)
!   write(0,*) ' g2rl=',g2rl, ' g2cm=',g2cm
  temp1 = 1.d0/(Tname%S0(Eindex)*cm2topgrm)  !  g/cm2
  temp2 =  c1forHardScat/(Tname%S1(Eindex)*cm2topgrm) 
  temp2 = min( temp2, maxHCSmfprl/g2rl)    ! maximum is 0.05rl --> g/cm2

  if(temp2 >  temp1 ) then
     icon  = 0
     temp = temp2
  else
     icon = 1
     temp = temp1
  endif

  lambdah(Eindex) = temp

!  write(0,*) ' icon =',icon, ' cm2topgrm=',cm2topgrm
!  write(0,*) ' Eindex=',Eindex, ' lambdah=', lambdah 
!!!!!!!!!
!  if(Eindex == 106)  then
!     write(0,*) ' L(el)=',temp1/cm2topgrm
!  endif
!!!!
end subroutine cgetLambdah

subroutine cgetResTPXS
!  get TPXS for soft scattering. (restricted TPXS)
!   Ss1= 4pi int(0,muc) of 2mu DCS dmu
!   Ss2= 4pi int(0,muc) of (1-P2(cos)) DCS dmu
!     where P2(cos) = (3(cos)^2 -1 )/2 = (3(1-2mu)^2 -1)/2
!  lambdas1 = 1/(Ss1*N);  lambdas2=1/(Ss2*N)
!   N = cm2topgrm
  use modcDCSf    
  use modDCS
  use modSampRedPDF
  implicit none
  real(8):: mux
!  real(8),intent(out):: lambdas1, lambdas2  ! in g/cm2

  real(8):: cmu1sInte(nmu), cmu2sInte(nmu)
  real(8):: ans, ss1, ss2
!     make table of  (1-Pl) dcs(mu)  (l=1,2)
! for muval (E dependent)
  call cmkIntegrandRTPXS(cmu1sInte, cmu2sInte)
!   and integrate them
  mux = muc(Eindex)
  if(mux > 0.) then
     call ktrpzIntT2(cmu1sInte, 1, nmu, muval, 1,  0.0d0, mux, ans)
     ss1 = fpi*ans
     lambdas1(Eindex) = 1.d0/(ss1*cm2topgrm)
     call ktrpzIntT2(cmu2sInte, 1, nmu, muval, 1,  0.0d0, mux, ans)
     ss2 = fpi*ans
     lambdas2(Eindex) = 1.d0/(ss2*cm2topgrm)
  else
!       ss1,ss2 are 0. so lambda is infty.   We put 0 instead
!      and don't use them.
     lambdas1(Eindex) = 0.
     lambdas2(Eindex) = 0.
  endif
end subroutine cgetResTPXS

subroutine cmkIntegrandRTPXS( k1, k2)
!  make k1(:)=(1-P1(cos)*ds/dmu for mu's in uarray of modcDCSf
!       k2(:)=(1-P2(cos)* ds/dmu
!  dcs must be ready
  use modcDCSf  

  implicit none
  real(8),intent(out):: k1(nmu)  ! 
  real(8),intent(out):: k2(nmu)
  real(8):: p2, cs
  integer:: i

     k1(:) = 2*muval(:)*dcsname%dcs(:,Eindex)
     do i = 1, nmu
        cs = (1.0d0-2*muval(i))
        p2 = (3.d0*cs**2 - 1.d0)/2.d0
        k2(i) =(1.d0-p2)* dcsname%dcs(i,Eindex)
     enddo
end subroutine cmkIntegrandRTPXS


subroutine ciniGetMuc
  use modcDCSf
  use modXsecMedia
  implicit none
  cm2topgrm =  media(1)%mbtoPkgrm/ 1d-27 *10.0
end subroutine ciniGetMuc

subroutine csolvehard(f,  mu)
  use modcDCSf
  implicit none
  real(8),external::f
  real(8),intent(out):: mu
  real(8):: vl=0., vr=1.0
  integer:: icon
  real(8),save:: epsilon=1.d-10

  vl = 0.d0
  vr = 1.0

  if( f(vl) * f(vr) > 0. ) then
     mu = 0.
  else
     call kbinChop(f, vl, vr, 0.5d0, epsilon, mu, icon)
     if( icon /= 0 ) then
        write(0,*) 'eRRR; icon= ',icon 
        do while (vl < 1.0 )
           write(0,*) vl, f(vl)
           vl = vl + 1.e-3
        enddo
        stop  
     endif
  endif
end subroutine csolvehard

