! 
!       This program gets numerical values for the coefficients
!    in formulas expessing the air density as a function of the height.
!                   
!                 10/hn exp( (hb-z)/hn )          at z > hc
!      rho(z) = 
!                 10/hl/hm (  (ha -z)/hl ) **( 1/hm -1)   at z < hc
!
!     where   hn = 7000 (m) ; shoud be in 6700 to 7300.
!             hc = 11.1d3 (m);
!
!     Other coefficients are determined so that 
!    
!    1)  at z=0, rho(0) = rho0 = 1.225 or little bit less.  (kg/m3)
!    2)  at z=hc, rho(hc) =rho1 = 0.35932d0 (kg/m3)
!    3)  at z=hc, rho is continus and rho'(=d rho/dz) is also continus.
!
      real*8  z,  rhop1,temp, hmi, rho, hlhmi, z0
      real*8  rho0, rho1, ha, hb, hm, hn, hc, f1, f2
      integer ios
      real*8 solv, x
!         to fix ha.
      solv(x) =( (x-hc)/(x-z0) )**(-rhop1/rho1 *  (x-hc))/(rho1/rho0) -1.d0

      z0 = 0.d3
      hn =6.4d3
      hc = 11.1d3
      rho1 = 0.35932d0
!     rho0 = .73643
!      rho0 = 1.225
      rho0 = 1.2      ! make rho0 bit smaller to have better fit < 11.1km


      hb = log(rho1*hn/10.d0) *hn + hc
      rhop1 = -10.d0/hn/hn * exp( (hb-hc)/hn)
!         binary chop for getting ha
      ha1 = 13000.d0
      ha2 = 45000.d0
      f1 = solv(ha1)
      f2= solv(ha2)
      if(f1 *f 2 .gt. 0.) stop 9876
      do while(.true.)
         ha = (ha1 + ha2)/2
         temp = solv(ha)
         if(ha .eq. ha1 .or. ha .eq. ha2) goto 10
         if(temp * f1 .le. 0.) then
            f2 = temp
            ha2 = ha
         else
            f1= temp
            ha1 = ha
         endif
!         write(0, *) temp, ha
      enddo
 10   continue
            
      hmi = -  rhop1/rho1 * (ha - hc) + 1.d0
      hm = 1.d0/hmi
      hl = (10.d0/hm/rho0)**hm  * (ha-z0)**(1.d0 - hm)
      hlhmi = 1.d0/hl/hm
!      write(0, *) ha, hb, hl,  hm, hn,  hmi, hlhmi
      write(*, *) '       data ha/', ha,'/, hb/',hb,'/,',
     *        '     * hl/',hl,'/, hm/',hm,'/,',
     *        '     * hn/',hn,'/, hmi/',hmi,'/,',
     *        '     * hlhmi/',hlhmi,'/, hc/',hc,'/'
!
!        read z vs rho data for comparison 
      open(12,file='stdatmos3.d', status='old', form='formatted',
     *    access='sequential')
      do while(.true.)
         read(12, *, iostat = ios) z, rho
         if(ios .lt. 0) stop  9999
         
         if(z .le. hc) then
            rhox = 10.d0/hl/hm * (( ha-z )/hl)**(hmi -1.d0)
         else
            rhox = 10.d0/hn*exp( (hb-z)/hn)
         endif
         write(*, *) sngl(z), sngl(rho), sngl(rhox)
      enddo
      end

