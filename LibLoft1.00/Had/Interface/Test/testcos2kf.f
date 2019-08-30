      implicit none
      integer kf, code, subcode, charge
      
      do while(.true.)
         write(*,*) ' Enter code subcode charge  or 0 / to stop'
         read(*,*) code, subcode, charge
         if(code .eq. 0) stop
         call ccos2kf( code, subcode, charge, kf)
         write(*,*)   ' c,s,cg=', code, subcode, charge,
     *     ' kf=',kf
      enddo
      end
