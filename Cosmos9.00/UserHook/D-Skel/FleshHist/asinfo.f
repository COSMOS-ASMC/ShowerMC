!      real depth, s, E0, sizex
!      integer code, id
!      do while (.true.)
!         write(0,*) ' enter code, E0'
!         E0=1.e8
!         code=2
!         read(*,*) code, E0
!         do id = 200, 1200, 100
!            depth = id
!            s= dep2s(E0, 1, depth)
!            write(*,*) depth, sizex(E0, 1, code, s)
!         enddo
!      enddo
!      end

      real function sizex(E0in, NN,  code, s)
      implicit none
      real E0in ! input primary total energy in GeV
      integer NN ! input. incident nucleon number 0--> photon, 1-->proton, 4-->He ...56-->Fe
      integer code ! input code
      real s ! input age

      real s2Ng, s2Ne, s2Nmu, s2Nh, E0, fac
      if(NN .ge. 1) then
         E0=E0in/NN
         fac=NN
      else
         E0=E0in
         fac =1
      endif
      if(E0 .gt. 1.e11 .or. E0 .lt. 1.e7) then
         write(0,*)
     *   "*******warning: E0/n=",E0, " is out of range ******"
      endif
         
      if(code .eq. 1) then
         sizex=s2Ng(E0, s)
         if(NN .eq. 0) fac = 1.2
      elseif(code .eq. 2) then
         sizex=s2Ne(E0, s)
         if(NN .eq. 0) fac = 1.2
      elseif(code .eq. 3) then
         sizex=s2Nmu(E0, s)
         if(NN .eq. 0) fac = 0.01
      elseif(code .le. 6) then
         sizex=s2Nh(E0, s)
         if(NN .eq. 0) fac = 0.01
      else
         write(0,*) ' code error=',code, ' to sizex'
         stop 6776
      endif
      sizex=sizex*fac
      end
!     ************************
      real function dep2cogdep(E0in, NN, dep)
!     ************************
      implicit none
      real E0in  ! input primary nucleus total energy in GeV.
      integer NN ! input incident nucleon number. 
      real dep ! input slant depth in g/cm^2

      real E0, eps
      eps = 3.e-3
      if(NN .eq. 0) then
         E0=E0in
      else
         E0=E0in/NN
      endif
!       function value is dep/cog. 
      if(abs(1.0- E0/1.e8) .lt. eps) then
         dep2cogdep=0.001442*dep
      elseif( abs(1.0-E0/1.e9) .lt. eps) then
         dep2cogdep=0.001247*dep
      elseif( abs(1.0-E0/1.e10) .lt. eps ) then
         dep2cogdep=0.001094*dep
      else
         dep2cogdep=(-7.556e-05*log(E0) + 0.002826999)* dep      
      endif
      end
!     ************************
      real function cogdep2Nmu(E0, cogdep)
!     ************************
      implicit none
      real E0 ! input GeV
      real cogdep ! input dep/cog
      real a, b, c, d
      a=3.2e8 
      if(E0 .eq. 1e8) then
        b=4.8316
        c=5.9692 
        d=0.709
      elseif(E0 .eq. 1e9) then
         b=4.6796
         c=3.7835
         d=1.091
      elseif(E0 .eq. 1e10) then
         b=4.1251
         c=1.8791
         d=1.818
      else
         b = -0.1534*log(E0) + 7.7246
         c = -0.8881*log(E0) + 22.282
         d = 0.24081*log(E0) -3.7845
      endif   
      cogdep2Nmu=a*cogdep**b * exp(-c*cogdep**d)
      end
!     ***********************
      real function cogdep2Nh(E0, cogdep)
!     ***********************
      implicit none
      real E0 ! input GeV
      real cogdep ! input dep/cog
!
      real cogdep2Nmu
!        use mu at present
      cogdep2Nh= cogdep2Nmu(E0,cogdep)
      end
!     ***********************
      real function s2Nh(E0, s)
!     ***********************
      implicit none
      real E0 ! input GeV
      real s  ! input age
!        use mu at present
      real s2Nmu
      s2Nh = s2Nmu(E0, s)
      end
!     ************************
      real function cogdep2Ne(E0, cogdep)
!     ************************
      implicit none
      real E0 ! input GeV
      real cogdep ! input dep/cog
      real a, b, c, d
      real eps
      eps = 3.e-3

      a=1.e11
      if(abs(1.- 1e8/E0) .lt. eps ) then
         b=8.1044
         c=7.4068
         d=1.12022
      elseif( abs(1.- 1e9/E0) .lt.eps)  then
         b= 7.278 
         c= 5.063 
         d= 1.555
      elseif( abs(1.-1.e10/E0) .lt. eps) then
         b=6.3661
         c=2.7725
         d=2.2605
      else
         b =  -0.3774*log(E0) +  15.07
         c = -1.0063*log(E0) +  25.93
         d = 0.24760*log(E0) -3.48601
      endif   
      cogdep2Ne=a*cogdep**b * exp(-c*cogdep**d)
      end
!     ************************
      real function cogdep2Ng(E0, cogdep)
!     ************************
      implicit none
      real E0 ! input GeV
      real cogdep ! input dep/cog
      real a, b, c, d

      a=1.e12
      if(E0 .eq. 1e8) then
         b=8.00  
         c=8.034 
         d=0.941 
      elseif(E0 .eq. 1e9) then
         b=7.71  
         c=5.652 
         d=1.383
      elseif(E0 .eq. 1e10) then
         b= 6.91 
         c= 3.271 
         d= 2.022
      else
         b=-0.23669*log(E0)+12.445
         c=-1.03427*log(E0)+ 27.085
         d=0.2347*log(E0) -3.41583
      endif   
      cogdep2Ng=a*cogdep**b * exp(-c*cogdep**d)
      end
!     ***************
      real function cogdep2s(E0,NN, cogdep)
!     ***************
!      for all energy      
      implicit none
      real E0   !  input not used now
      integer NN ! input incident nucleon number. not used now
      real cogdep ! input dep/cog
      cogdep2s = 10.28*cogdep**0.0407-9.228
      end
!     ***************
      real function dep2s(E0, NN, dep)
!     ***************
      implicit none 
      real E0  !input GeV
      integer NN ! input  incident nucleon number. not used now
      real dep ! input dep g/cm^2
      
      real dep2cogdep, cogdep2s
      dep2s=cogdep2s(E0,NN, dep2cogdep(E0,NN, dep))
      end
!     ***************
      real function s2Ng(E0, s)
!     ***************
      implicit none 
      real E0  !input total incident energy/n  GeV
      real s ! input  age
      
      real cogdep2Ng, s2cogdep
      s2Ng=cogdep2Ng(E0,s2cogdep(E0,s))
      end
!     ***************
      real function s2Ne(E0, s)
!     ***************
      implicit none 
      real E0  !input GeV
      real s    ! input age
      
      real s2cogdep, cogdep2Ne

      s2Ne=cogdep2Ne(E0, s2cogdep(E0, s))
      end
!     ***************
      real function s2Nmu(E0, s)
!     ***************
      implicit none 
      real E0  !input GeV
      real s    ! input age
      
      real s2cogdep, cogdep2Nmu

      s2Nmu=cogdep2Nmu(E0,s2cogdep(E0, s))
      end
!     *************
      real function s2cogdep(E0, s)
!     *************
      real E0    ! not used now
      real s     ! input. age
      s2cogdep= ( (s+9.228)/10.28 )**(1./0.0407)
      end
