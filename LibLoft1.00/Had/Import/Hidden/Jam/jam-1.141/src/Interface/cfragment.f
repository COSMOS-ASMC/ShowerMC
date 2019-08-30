!   test program is at end. use next
!  cpp -C -P -traditional  cfragment.f > temp.f
!  ifort temp.f -L$LIBLOFT/lib/$ARCH -lloft

! This is for Jam (original jam) which breaks up all spectator nucleons  into individul
!    nucleons; to avoid such, we employ sibyll method of making fragments.
!   This is to treat projectile spectator fragmentation.
!   Both need data in Cosmos/Data/Fragment/
#undef FragMom
!!!       #define FragMom
!       If above is undef, 
!       we don't use Pppcom which should contain 3 momentum
!       of evaporated nucleus (& nucleons) at the projectile
!       sysetm.  The sampling method in the program seems
!       wrong and Pz is not given. (Pz= p cos(theta); 
!       Px=p sin(theta) cos(fai); Py = p sin(theta)sin(fai)
!       where theta is uniform in dcos(theta)
!       But, in the program, it is uniform in d(theta)
!       (like  fai is uniform in dfai).
!       Since  Jam code gives momentum so we use them
!       and need not use Pppcom. (origian is ppp)

      module modfragment
!       This is based on sibyll program. Some modifications
!     have been done for use with fortran90 which is used
!     to avoid name collision etc.
!
      integer,parameter:: maxfrag=60
      integer,parameter:: nheavy=10  ! original is 10 upto Fe      
      integer, parameter:: neaa=10
!      common /fragmod/a(10,10,20),ae(10,10,20),eres(10,10),nflagg(10,10)
      real(8),save:: Acom(nheavy,neaa,20)
#if defined (FragMom)
      real(8),save:: Aecom(nheavy,neaa,20)
      real(8),save:: Erescom(nheavy,neaa)
#endif
      integer,save:: Nflaggcom(nheavy,neaa)
!      common /fragments/ ppp(3,60)
#if defined (FragMom)      
      real(8),save:: Pppcom(3,maxfrag)
#endif
      contains
!====================================================================
!...code of fragmentation  of spectator nucleons
!.  based on jon engel  abrasion-ablation algorithms
!====================================================================
      subroutine creadfragdata(io)
      implicit none
      integer,intent(in)::io  ! logical dev. # for input
      integer:: i, j, k

      do k =1, 20
         do j = 1, 10
            read(io,*) Acom(:,j,k)
         enddo
      enddo

      read(io,*)

      do k =1, 20
         do j = 1, 10
#if defined (FragMom)            
            read(io,*) Aecom(:,j,k)
#else
            read(io,*)
#endif
         enddo
      enddo

      read(io,*)

      do j = 1, 10
#if defined (FragMom)            
         read(io,*) Erescom(:,j)
#else
         read(io,*)
#endif
      enddo

      read(io,*) 

      do j = 1, 10
         read(io,*) Nflaggcom(:,j)
      enddo

      end       subroutine creadfragdata

      subroutine cwritefragdata
!      common /fragmod/a(10,10,20),ae(10,10,20),eres(10,10),nflagg(10,10)
      integer  i, j, k
      do k =1, 20
         do j = 1, 10
            write(0,'(1p,10g10.3)') Acom(:,j,k)
         enddo
      enddo

      write(0,*) ' '
#if defined (FramMom)
      do k =1, 20
         do j = 1, 10
            write(0,'(10f7.2)') Aecom(:,j,k)
         enddo
      enddo

      write(0,*) ' '

      do j = 1, 10
         write(0,'(1p,10g10.4)') Erescom(:,j)
      enddo

      write(0,*) ' '
#endif
      do j = 1, 10
         write(0,'(10i2)') Nflaggcom(:,j)
      enddo
      end subroutine cwritefragdata

      function s_rndm(x) result(u)
      implicit none
      integer:: x ! dummy
      real(8):: u
      call rndc(u)
      end      function s_rndm

      subroutine fragm(iat,iap, nw, b, nf, iaf, Z)
      implicit none
!...nuclear fragmentation, abrasion-ablation model, 
!...based on jon engel's routines abrabl 
!...this most recent version adds for all prefragment
!...masses > 10 the model calculation for the fragment
!...mass distribution and the energy carried by the fragment
!...of W. Friedmann
!...the average values are used to implement the model
!...in the montecarlo fashion / tss, dec '91
!.
      integer,intent(in)::  iap ! mass of incident nucleus; 60 is ok 
      integer,intent(in)::  iat ! mass of target   nucleus
      integer,intent(in)::  nw  ! number of wounded nucleons in the beam nucleus
                                ! (not spectator) (nfn in cTargetFrag.f)
                                ! spectator = iap - nw
      real(8),intent(in)::  b   ! impact parameter in the interaction (fm ?)
      
      integer,intent(out):: nf !  number of fragments  of the spectator nucleus
                               ! <= maxfrag
      integer,intent(out):: iaf(:) ! iaf(1:nf) = mass number of each fragment
      integer,intent(out):: Z(:)   ! charge of each fragment
!.          Pppcom(3,maxfrag) in modfragment contains
!.          the three momentum components (MeV/c) of each
!.          fragment in the projectile frame
!..............................................................
!      common /fragmod/a(10,10,20),ae(10,10,20),eres(10,10),nflagg(10,10)
!             heavy mass A

      real(8),save::aa(nheavy)=
     *    (/10.,  15.,  20., 25.,  30.,  35.,  40.,  45.,  50.,  56./)
!             fragment mass 

      real(8),save::eaa(neaa)=(/1.,2.,4.,6.,8.,10.,12.,16.,20.,30./)

      real(8):: ap, at, eb
      integer:: npf
      real(8)::fk, sig, ppfx, ppfy, eps, etot, pp, pxe, pye, arat
      real(8):: esum, fr, fr1 , ekin
      real(8):: pprob
      real(8):: s1, c1, s2, c2
      real(8):: ebp, erat, f1
      integer:: nnuc, nalp, nuc, nsum
      integer:: i, j, k, ja, je, jf, n1
      logical,save::first=.true.
      integer,parameter:: io = 11
      integer:: icon
      real(8):: u

      if( first ) then
         call copenf(io,
     *    '$LIBLOFT/Data/Fragment/sibyllFrag.data', icon)
         if(icon /= 0 ) then
            write(0,*)
     *         '$LIBLOFT/Data/Fragment/sibyllFrag.data could not',
     *         ' be opened'
            stop
         endif
         call creadfragdata(io)
         close(io)
         first = .false.
      endif

      npf = iap - nw    ! 
      if (npf .eq. 0) then
         nf = 0
         return   ! ****************
      endif
      ap=iap
      at=iat

      eb = estar(ap,at, b)
      ebp = estarp (npf, nw)
! contribution to E* from energy deposited by secondaries
      eb = eb + ebp
! total E* is the sum of the two components

!.....prefragment transverse momentum (mev/nucleon)...
      fk = fermk(ap)

! Fermi momentum of the projectile nucleus
      if (nw .lt. iap) then
         sig = fk*sqrt(nw*npf/(ap-1.))/3.162

! gaussian sigma in all three direction
      else
         sig = fk/3.162
! this is not correct, too large !!!!!!!!!!!!!!
      endif
      ppfx = sig*gasdev(0)/npf
      ppfy = sig*gasdev(0)/npf

! three momentum components per nucleon for the prefragment

!.............crude model for small prefragment mass .......
      if (npf <  10) then
         call evap(npf, eb, eps, nnuc, nalp)
         
!                  eps is the kinetic energy carried by the evaporated nucleons
         etot = 938. + eps
         pp = sqrt((etot*etot - 8.79844e5)/3.)

!                  average momentum of evaporated nucleons in each direction
         nuc = npf - nnuc - 4*nalp
         nf = 0
         if (nuc .gt. 0) then
            nf = nf + 1
            iaf(nf) = nuc
#if defined (FragMom)
            Pppcom(1,nf) = nuc*ppfx
            Pppcom(2,nf) = nuc*ppfy
!//////////
!            write(0,*) '1 Pppcom nf=', nf, Pppcom(1:2,nf)
!///////////
#endif
         endif
         if (nalp .ne. 0) then
            do i=1,nalp
               nf = nf + 1
               iaf(nf) = 4
#if defined (FramMom)
               call kcossn(s1,c1)
               call kcossn(s2,c2)
               pxe = 4.*pp*s1*s2
               pye = 4.*pp*s1*c2
               Pppcom(1,nf) = 4.*ppfx + pxe
               Pppcom(2,nf) = 4.*ppfy + pye
!//////////
!            write(0,*) '2 Pppcom nf=', nf, Pppcom(1:2,nf)
!///////////
               Pppcom(1,1) = Pppcom(1,1) - pxe
               Pppcom(2,1) = Pppcom(2,1) - pye
!//////////
!            write(0,*) '3 Pppcom nf=', 1, Pppcom(1:2,nf)
!///////////
#endif
            enddo
         endif
         if (nnuc .ne. 0) then
            do i=1,nnuc
               nf = nf + 1
               iaf(nf) = 1
#if defined (FramMom)
               call kcossn(s1,c1)
               call kcossn(s2,c2)
               pxe = pp*s1*s2
               pye = pp*s1*c2
               Pppcom(1,nf) = 4.*ppfx + pxe
               Pppcom(2,nf) = 4.*ppfy + pye
!//////////
!            write(0,*) '4 Pppcom nf=', nf, Pppcom(1:2,nf)
!///////////
               Pppcom(1,1) = Pppcom(1,1) - pxe
               Pppcom(2,1) = Pppcom(2,1) - pye
!//////////
!            write(0,*) '5 Pppcom nf=', 1, Pppcom(1:2,nf)
!///////////
#endif
            enddo
         endif
      else   !  npf >= 10
!.........more refined model calculation .............
         ja = npf/5 -1
         if (ja <  nheavy ) then
            if ((npf - aa(ja)) .gt. (aa(ja+1)-npf)) ja = ja + 1
         endif
         arat = float(npf)/aa(ja)

         do j=1,neaa
            if (eb .lt. eaa(j)) go to 29
         enddo
         je = neaa
         go to 39
 29      je = j
 39      if (je .gt. 1 .and. je .ne. neaa) then
            if ((eb - eaa(j-1)) .lt. (eaa(j)-eb))  then
               je = j - 1
            endif
         endif
         erat = eb/eaa(je)
         
         if (eb .lt. 1.) then
            erat = eb
         endif
! interpolate between eb=0. (nothing happens) and eb = 1. mev
      
         if (ja .eq. nheavy .and. je .gt. 6) then
            write(0,*)'wargning ja=',ja,',   je=',je, ' in fragm'
         endif
 43      esum = 0.
         nsum = 0
         jf = 0
         do j=20,1,-1
            fr =  Acom(ja, je, j)*arat*erat
            n1 = 1 + fr
            fr1 = fr/float(n1)
            do k=1, n1
               if (s_rndm(0) .lt. fr1) then
                  jf = jf + 1
                  iaf(jf) = j
                  nsum = nsum + j
#if defined (FramMom)
                  ekin = erat*Aecom(ja,je, j)
                  if (ekin .gt. 0.) then
                     esum = esum + ekin
                     etot = 938.*iaf(jf) + ekin
                     pp = sqrt(2.*(etot*etot - iaf(jf)**2*8.79844e5)/3.)
                     call kcossn(s1,c1)
                     call kcossn(s2,c2)
                     Pppcom(1,jf) = pp*s1*s2 + iaf(jf)*ppfx
                     Pppcom(2,jf) = pp*s1*c2 + iaf(jf)*ppfy
!//////////
!                     write(0,*) '6 Pppcom jf=', jf, Pppcom(1:2,jf)
!                     write(0,*) ' pp, ppfx, ppfy, s1, s2, iaf'
!                     write(0,*)  pp, ppfx, ppfy, s1, s2, iaf(jf)
!///////////
                  endif
#endif
                  if (nsum .gt. npf) then
!         write(*,*)' warning, nsum=', nsum,',  npf=',npf
!         write(*,*)'  arat =', arat
                     go to 43
                  else
                     if (nsum .eq. npf) then
                        go to 44
                     endif
                  endif
               endif
            enddo
         enddo
         if (Nflaggcom(ja,je) .eq. 0) then
!             'the residue' is a nuclear fragment
            jf = jf + 1
            iaf(jf) = npf - nsum
#if defined (FragMom)
            f1 = npf*eb - esum
            if (f1 .lt. 0.) f1 = 0.
! give the rest of eb to the fragment
            ekin = f1
            if (ekin .gt. 0.) then
               etot = 938.*iaf(jf) + ekin
               pp = sqrt(2.*(etot*etot - iaf(jf)**2*8.79844e5)/3.)
               call kcossn(s1,c1)
               call kcossn(s2,c2)
               Pppcom(1,jf) = pp*s1*s2 + iaf(jf)*ppfx
               Pppcom(2,jf) = pp*s1*c2 + iaf(jf)*ppfy
!//////////
!            write(0,*) '7 Pppcom jf=', jf, Pppcom(1:2,jf)
!///////////
            endif
#endif
         else
! 'the residue' consists of spectator nucleons
            n1 = npf - nsum
            do k=1,n1
               jf = jf + 1
               iaf(jf) = 1
#if defined (FragMom)
               ekin = erat*Erescom(ja,je)
               if (ekin .gt. 0.) then
                  etot = 938.*iaf(jf) + ekin
                  pp = sqrt(2.*(etot*etot - iaf(jf)**2*8.79844e5)/3.)
                  call kcossn(s1,c1)
                  call kcossn(s2,c2)
                  Pppcom(1,jf) = pp*s1*s2 + ppfx
                  Pppcom(2,jf) = pp*s1*c2 + ppfy
!//////////
!                  write(0,*) '8 Pppcom jf=', jf, Pppcom(1:2,jf)
!                  write(0,*) ' pp, ppfx, ppfy, s1, s2'
!                  write(0,*)  pp, ppfx, ppfy, s1, s2
!///////////
               endif
#endif
            enddo
         endif
 44      nf = jf
      endif

      do i = 1, nf
         call csetFragChg( iap, iaf(i), Z(i) ) 
      enddo
      end  subroutine fragm 

      function estarp (npf, nw) result(ans)
      implicit none
      integer,intent(in):: npf ! #  of spectator nucleons
      integer,intent(in):: nw  ! # of woonded projective nucleons
! contribution to E* from energy deposited by secondaries
! very naive version incorporating hueffner's ideas
      real(8):: ans

      real(8):: apf, f1,  f2
      integer:: i

      apf = npf
      f1 = 15.3/apf**0.666666666
! average kinetic energy/nucleon in prefragment (mev)
! per pathlength equal to the prefragment radius
      ans = 0.
      do i=1, nw
         if (s_rndm(0) > 0.5) then
            f2 = f1*rdis(0)
            ans = ans + f2
         endif
      enddo
      end function estarp
      
      function rdis(idum) result(ans)
      implicit none
      integer,intent(in):: idum ! dummy
      integer,parameter::nprob=20
      real(4),save:: probr(0:nprob)=(/0.,
     *      0.10000, 0.15748, 0.21778, 0.28605, 0.36060,
     *      0.43815, 0.51892, 0.60631, 0.70002, 0.79325,
     *      0.88863, 0.98686, 1.10129, 1.21202, 1.32932,
     *      1.44890, 1.57048, 1.70139, 1.83417, 2.00000/)
      real(8)::ans

      integer nr
      real(8)::f1, dr

      nr = nprob*s_rndm(0) + 1

!      if (nr .eq. 1) then
!      f1 = 0.
!      else
      f1 = probr(nr-1)
      dr = probr(nr) - f1
      ans = f1 + dr * s_rndm(0)
      end      function rdis


      function estar(ap,at,b) result(e1b)
      implicit none
      real(8),intent(in):: ap
      real(8),intent(in):: at
      real(8),intent(in):: b
      real(8):: e1b
      
      real(8):: sigma, rt, rp, alpha, beta, f, alf, alalf
      real(8):: gfac, gfac1, s1, s2, s3, g0, g1, g2, g3,pb
      real(8):: f11
      integer:: ii, n

!      sigma=4.5                 ! total n-n cross section in fm**2  45 mb
      sigma = 4.0            !    40 mb
      rt=.82*at**.3333 !target radius
      rp=.82*ap**.3333 !projectile radius
      alpha=(rt/rp)**2  
      beta=(b/rt)**2
      f=at*sigma/(3.14159*rt**2)
      alf = log(f)
      alalf = log(alpha)
      gfac=0
      gfac1=0
      s1=0.
      s2=0.
      s3=0.      
      ii=1
      do n = 0,  10 ! this limit may not need to be so high.
         if(n  >= 2) then
            gfac1=gfac
            gfac=gfac+log(float(n)) 
         endif

         g0=n*alf -n*beta*alpha/(n+alpha)+alalf
         g1=g0-log(alpha+n)-gfac
         g2=(n+2)*log(f)-(n+2)*beta*alpha/(n+2+alpha) 
     *      +log(n+2+alpha+beta*alpha**2)-3*log(n+2+alpha)-gfac
         g3=g0-2*log(n+alpha)-gfac1

         ii=-ii
         s1=s1+ii*exp(g1)
         s2=s2+ii*exp(g2)
         if( n >= 1) s3=s3+ii*exp(g3)
      enddo

      pb = s1
      e1b = 197.**2/(2*938.*rp**2*pb) *s2

      end      function estar

      subroutine evap(npf,eb,eps,nnuc,nalp) 
      implicit none
      integer,intent(in):: npf
      real(8),intent(in):: eb

      real(8),intent(out):: eps      
      integer,intent(out):: nnuc
      integer,intent(out):: nalp

      integer:: n 

      eps=7.5+sqrt(8*eb)
      n=min(npf*int(eb/eps),npf)
      nalp=n/5
      nnuc=n-4*nalp
      end      subroutine evap
!->
      function fermk(a) result(ans)
      implicit none
      real(8),intent(in):: a  ! A

      real(8)::ans
      real(8):: f11, f12, f13, f21, f22, f23
      real(8),save::aa(6)=(
     */4., 6., 12., 24., 40., 57./)
      real(8),save:: fk(6)=(
     * /130.,169.,221.,235.,251.,260./ )
      integer i
!!        A vs fk
!!        12  184
!!        20  230
!!        40  240
!!        93  255
!!       197  265     
      do i=2,4
         if (a .lt. aa(i)) go to 25
      enddo
      i = 5
 25   f11 = aa(i-1)

      f12 = aa(i)
      f13 = aa(i+1)
      f21 = fk(i-1)
      f22 = fk(i)
      f23 = fk(i+1)

      ans = quad_int(a,f11,f12,f13, f21,f22,f23)
      end      function fermk


      function quad_int (r,x0,x1,x2,v0,v1,v2) result(ans)
!...quadratic interpolation
      implicit none
      real(8),intent(in):: r   !  at which interpolated value needed
      real(8),intent(in):: x0,x1, x2  !  3 x values
      real(8),intent(in):: v0,v1, v2  ! function values
      real(8):: ans ! interpolated value
      real(8):: s0, s1, s2
      real(8):: r0, r1, r2

      r0=r-x0
      r1=r-x1
      r2=r-x2
      s0=x0-x1
      s1=x0-x2
      s2=x1-x2
      ans = v0*r1*r2/(s0*s1)-v1*r0*r2/(s0*s2)+v2*r0*r1/(s1*s2)
      end      function quad_int

      function gasdev(idum) result(ans)
!...gaussian deviation
      implicit none
      integer,intent(in):: idum ! dummy
      real(8):: ans
!! The method to utilize two variances is N.G for
!!       skeleton-flesh  method of Cosmos.
!!  So we use kgauss
!!      save gset
!!      data iset/0/
!!      if (iset.eq.0) then
!!1       v1=2.*s_rndm(0)-1.
!!        v2=2.*s_rndm(0)-1.
!!        r=v1**2+v2**2
!!        if(r.ge.1.)go to 1
!!        fac=sqrt(-2.*log(r)/r)
!!        gset=v1*fac
!!        gasdev=v2*fac
!!        iset=1
!!      else
!!        gasdev=gset
!!        iset=0
!!      endif
      call kgauss(0.d0, 1.d0, ans)
      end      function gasdev

      end module modfragment
!
!      program main
!      use modfragment
!      implicit none
!      integer:: iat, iap, nw 
!      real(8):: b
!      integer:: nf
!      integer:: iaf(maxfrag)
!      
!      integer:: i,j 
!      iat = 12
!      iap = 56  !  must be <=56
!      nw = 1
!      b = 5.0
!      do i = 1, 10000000 
!         call  fragm(iat, iap, nw,b, nf, iaf)
!         if(i< 100) then
!         write(*,*) ' nf =', nf
!         write(*,*) ' iaf =',iaf(1:nf)
!#if defined (FramMom)
!         do j = 1, nf
!            write(*,'(1p,3g12.3)') Pppcom(1:3,i)
!         enddo
!#endif
!         endif
!      enddo
!      end       program main

      subroutine csetFragChg(ia,  fm, fc)
!      fix   projectile fragment charge
      implicit none
      integer,intent(in):: ia  ! heavy proj. mass #
      integer,intent(in):: fm !  fragment mass #
      integer,intent(out):: fc !  fragment charge

      real(8):: pprob, u

      if( fm  == 1 ) then
         if( ia > 29 ) then
            pprob = 0.4         ! proton prob.
         else
            pprob = 0.5
         endif
         call rndc(u)
         if(u < pprob ) then
            fc = 1
         else
            fc = 0
         endif
      elseif( fm == 2 ) then
         fc = 1
      elseif( fm <= 4 ) then
         fc = 2
      elseif( fm < 29 ) then
         if( fm == 27) then
            fc = 13
         elseif( fm == 19 ) then
            fc = 9
         else
            fc = fm*0.5
         endif
      else
         if( fm == 56 ) then
            fc = 26
         elseif( fm == 40 ) then
            fc = 18
         elseif( fm == 48 )  then
            fc = 22
         elseif( fm < 70 ) then
            fc = fm*0.47
         else 
            fc =fm*0.4
         endif
      endif
      end   subroutine csetFragChg


      
      
      

