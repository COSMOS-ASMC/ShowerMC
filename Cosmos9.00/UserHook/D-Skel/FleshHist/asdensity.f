!c       test asdensity
!      real r, s(10), den
!      real asdensity
!      integer code
!      integer ir, i, pw 
!      
!      do while (.true.)
!         write(0,*) 'enter code '
!         read(*,*) code
!         do i= 1, 10
!            s(i) = 0.
!         enddo
!         if(code  .eq. 3) then
!            write(0,*) 'enter depth/cog upto 10 with /'
!            pw=2
!         else
!            write(0,*) 'enter ages upto 10 with /'
!            pw=3
!         endif
!         read(*,*)  s      
!
!         do i = 1, 10
!            if(s(i) .gt. 0.) then
!               r= 0.01
!               do ir = 1, 42
!                  den=asdensity(1.e9, code, s(i), r)
!                  write(*,*)
!     *             r, den*2*3.1415*r**pw
!                  r= r*10.**0.1
!               enddo
!               write(*,*)
!            endif
!         enddo
!      enddo
!      end
      real function asdensity(E0, code, age, r)
      implicit none
!
!        air shower particle density /m.u^2 for vertical shower
!
      real E0 ! input. primary proton total E in GeV. not used

      integer code ! input particle code. 1---> gamma
                   !                      2---> electron
                   !                      3---> muon
                   !                      4---> hadrons

      real  age  !  input  age of the shower.  
!
      real r   !  input.  core distance in Moliere unit.
               !        measured at 2 r.l above the observation
               !        point along the shower axis.
!
      real  a, b, c, d, x
      real ag1, bg1, cg1, dg1
      real ag2, bg2, cg2, dg2
      real ae1, be1, ce1, de1
      real ae2, be2, ce2, de2
      real am, bm, cm, dm
      real s, f, rr
      real twopi
      real y, s2cogdep 
      parameter (twopi=2*3.141592)
!       for gamma   all r
!        a                  b               c
!
      ag1(s)= 8.610*s**(-1.585)
      bg1(s) =3.4716*s**(-0.640)
      cg1(s) = -0.41456*s*s + 0.9951*s -0.1680515
      dg1(s) = 0.47598*s*s -1.5426*s + 0.798
!       e r<2
      ae2(s) = 82.66*s**2.433
      be2(s) = 6.130*s**0.1376
      ce2(s) = 0.39671*s**(-0.29685)
      de2(s) = -(0.50456*s**2.05154)
!      e r>1
      ae1(s) =0.18857*s**2.2397
      be1(s) =0.1144*s**2.1224
      ce1(s) =0.94314*s**(-0.59357)
      de1(s) = 2.0624*s**(-0.69654)
!       mu  all r
      am(y) = 0.77492* y**(-1.804)
      bm(y) = 1.24209* y**(-0.7329)
      cm(y) = 0.56881* y**0.156110
      dm(y) = -0.4031*y** 0.781031
!
      x = age
      if(code .eq. 1) then
         a=ag1(x)
         b=bg1(x)
         c=cg1(x)
         d=dg1(x)
      elseif(code .eq. 2) then
         if(r .gt. 1.5) then
            a=ae1(x)
            b=be1(x)
            c=ce1(x)
            d=de1(x)
         else
            a=ae2(x)
            b=be2(x)
            c=ce2(x)
            d=de2(x)
         endif
      elseif(code .le. 6) then
         y = s2cogdep(E0, x)
         a=am(y)
         b=bm(y)
         c=cm(y)
         d=dm(y)
      else
         write(0,*) ' exor code for asdensity'
         stop 9753
      endif
      rr=max(r, 0.01)
      f = a*exp(-b*rr**c) / rr**d
      asdensity = f/twopi/rr
      end
