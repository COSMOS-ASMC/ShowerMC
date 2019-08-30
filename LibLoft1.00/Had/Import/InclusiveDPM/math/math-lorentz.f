!      program testlorentz
!      implicit real*8 (a-h,o-z)
!      dimension telement(5), p0(5), pl(5), p1(5), p2(5),pmu(5)
!      dimension spin0(5), spin1(5), spin2(5), zeta(3)
!
!      data spin0/0,0,1,0,0/
!      data  p0/0., 0., 0., .1056, .1056/
!c      data pmu/0.0, 0.0,0.030173714, 0.010982629, .1056/
!      data pmu/0.0, 0.0,0.0, 0.1056, .1056/
!      pmu(4) =sqrt(pmu(1)**2+pmu(2)**2+pmu(3)**2 + pmu(5)**2)
!
!c         write(*,*) pmu
!c         write(*,*) 'polarization of muon ?'
!c         read(*,*) pol
!
!c         call samplespin(pol, pmu, spin0)
!         write(*,*) 'spin0',spin0
!         call spin2zeta(spin0, pmu, zeta)
!         write(*,*) 'zeta0',zeta
!
!c         call mom2telement(pmu,telement)
!c         write(*,*) 'T0',telement
!c         do 10 j=1, 3
!c            telement(j) = -telement(j)
!c 10      continue
!c         call lboost(telement, pmu, p1)
!c         write(*,*) 'mu at motion', pmu
!c         write(*,*) 'mu at rest', p1
!
!c         call lboost(telement, spin0, spin1)
!c         write(*,*) 'spin at rest', spin1
!c         call spin2zeta(spin1,p1,zeta)
!c         write(*,*) zeta
!      do 1 i=1, 10
!         write(*,*) ' boost system (in puon momentum)'
!         read(*,*) (pl(k),k=1,3)
!
!         pl(5) = 1.0d0
!         pl(4) = sqrt(pl(1)**2 + pl(2)**2 + pl(3)**2 + pl(5)**2) 
!         write(*,*) 'P as proton',pl
!
!         call mom2telement(pl,telement)
!         write(*,*) 'T',i,telement
!c     &    ,    telement(1)**2+telement(2)**2+telement(3)**2
!
!         call lboost(telement, pmu, p1)
!         write(*,*) 'mu-boosted',p1
!         call lboost(telement, spin0, spin1)
!         write(*,*) 'spin-boosted',spin1
!
!         call spin2zeta(spin1,p1,zeta)
!         write(*,*) 'zeta boosted',zeta
!         do 11 k=1,5
!            spin0(k)=spin1(k)
!            pmu(k) = p1(k)
! 11      continue
!
!         call mom2telement(p1, telement)
!         write(*,*) 'T1',telement
!         do 12 j=1, 3
!            telement(j) = -telement(j)
! 12      continue
!
!         call lboost(telement, p1, p2)
!         write(*,*) 'muon at rest',p2
!
!         call lboost(telement, spin1, spin2)
!         write(*,*) 'spin at rest',spin2
!
!         call spin2zeta(spin2,p2,zeta)
!         write(*,*) 'zeta2 at rest',zeta
! 1    continue
!      end

      subroutine lboost(telement, p0, p1)
      implicit real*8 (a-h, o-z)
! Lorentz transformation by lorentz TRansformation ELEMENTs: 
! telement(1-3) : velosity direction
! telement(4) : gamma
! telement(5) : beta
      dimension telement(5), p0(5), p1(5)

      pdirg = p0(1)*telement(1)+p0(2)*telement(2)+p0(3)*telement(3)
      gamma = telement(4)
      beta  = telement(5)

      p12 = 0. 
      do 2 i=1, 3
         p1(i) = ((gamma-1.d0)*pdirg + gamma*beta*p0(4))*telement(i)
     &        + p0(i)
         p12 = p12 + p1(i)**2
 2    continue

!     Calculate E from E = sqrt(P^2 + M^2), since original
!     p1(4) = gamma*(p0(4) + beta*pdirg) sufferes from comput. errors.
!     For space-like 4-vectors i.e. spins, however, 
!     we calculate M^2 = -p(5)**2 to avoid the imaginary numbers, and 
!     determine the sign of p1(4) from original expression.

      if(p0(5).lt.0) then       
         p4 = p0(4) + beta*pdirg
         p1(4) = sign(sqrt(p12 - p0(5)**2), p4)
      else
         p1(4) =      sqrt(p12 + p0(5)**2)
      end if

      p1(5) = p0(5)
      end

      subroutine mom2telement(p,telement)
      implicit real*8 (a-h,o-z)
! Calulate lorentz TRansformation ELEMENTs: 
! telement(1-3) : velosity direction
! telement(4) : gamma
! telement(5) : beta
! 
      dimension p(5), telement(5)
      ap2 = p(1)**2 + p(2)**2 + p(3)**2
      ap = sqrt(ap2)

      gamma = p(4)/p(5)
      if(gamma.gt.1.000001d0) then ! Ek > 1eV 
         beta = sqrt(1.d0 - 1.d0/gamma**2)
         do 1 i=1, 3
            telement(i) = p(i)/ap
 1       continue
      else
         beta = 0.
         gamma = 1.d0
         do 2 i=1, 3
            telement(i) = 0.0
 2       continue
      end if
      telement(4) = gamma
      telement(5) = beta
      end

      subroutine unit3vec(a, dira, rr)
      implicit real*8 (a-h, o-z)
      dimension a(*),dira(*),u(3)
      rr=sqrt(a(1)**2 + a(2)**2 + a(3)**2)

      if(rr.gt.0.) then         ! normalize
         do 1 i=1, 3
            dira(i) = a(i)/rr
 1       continue
      else                      ! random unit 3-vector
         rr=2.0
         do 3 while (rr.gt.1.0d0 .or. rr.lt.1.0d-2)
            rr=0.0
            do 13 i=1,3
               call rndc(u1)
               u(i) = 2*u1 - 1.0d0
               rr = rr + u(i)**2
 13         continue
 3       continue
         rr = sqrt(rr)

         do 2 i=1,3
            dira(i) = u(i)/rr
 2       continue
      end if
      end

      subroutine samplespin(pol, p, spin)
      implicit real*8 (a-h, o-z)
!     construct a (muon) 4-vector spin for given p and polarization
      dimension p(5), dirp(3), zeta(3), u(3), spin(5)

      call sample_angvec(pol, p, zeta)
      call zeta2spin(zeta, p, spin)
      end

      subroutine sample_angvec(pol, p, zeta)
      implicit doubleprecision (a-h, o-z)
!     **Sample a unit vector zeta, satisfies pol = (zeta*p)/|p|
      dimension p(5), dirp(3), zeta(3), u(3), spin(5)

!     0:Make the direction of p
      call unit3vec(p, dirp, ap)

      if(pol.ge.0.99999d0) then
!     0-1: i.e. polarization = 1 

         do 14 i = 1, 3
            zeta(i) = dirp(i)
 14      continue
      else
!     1 polarization < 1, (Sample a vector perpendicular to dirp)

         ru =0.0
         do 1 while(ru.lt.1e-6)
!     1-1: sample a vector inside of a unit sphere

            rr=2.0
            do 3 while (rr.gt.1.0d0 .or. rr.lt.1.0d-2)
               rr=0.0
               do 13 i=1,3
                  call rndc(u1)
                  u(i) = 2*u1 - 1.0d0
                  rr = rr + u(i)**2
 13            continue
 3          continue

!     1-2: substruct the component parallel to dirp
            pinner = dirp(1)*u(1) + dirp(2)*u(2) + dirp(3)*u(3)

            ru = 0.0
            do 2 i=1,3
               u(i) = u(i) - dirp(i)*pinner
               ru   = ru + u(i)**2
 2          continue
!     1:Resulted vector shouldn't be zero (see condition in 'do 1 while')
 1       continue

!     2:sin(teta)* normalization
         tfact = sqrt((1.d0 - pol**2)/ru)
         do 4 i = 1, 3
            zeta(i) = dirp(i)*pol + u(i)*tfact
 4       continue
      end if
      end

      subroutine zeta2spin(zeta, p, spin)
      implicit real*8 (a-h,o-z)
! 3:Construct spin 4-vector (space-like) from spin direction zeta
      dimension zeta(3), p(5), spin(5)

      zp = zeta(1)*p(1) + zeta(2)*p(2) + zeta(3)*p(3)
      pcomp = zp/p(5)/(p(4)+p(5))
      ss = 0.
      do 5 i=1,3
         spin(i) = zeta(i) + pcomp*p(i)
         ss = ss + spin(i)**2
 5    continue
!     Because  s4^2 - (s1^2+s2^2+s3^2) = -1,
!     spin(4) = zp/p(5) is replaced by
      spin(4) = sqrt(ss - 1.d0)
      spin(5) = -1.d0
      end

      subroutine spin2zeta(spin, p, zeta)
      implicit real*8 (a-h,o-z)
      dimension spin(5), p(5), zeta(3)
! 4: Make spin direction zeta from boosted spin 4-vector.
      do 2 i=1, 3
         zeta(i) = spin(i) - spin(4)*p(i)/(p(4)+p(5))
 2    continue
      end
