      integer kf, code, subcode, charge
      
      do while(.true.)
         write(*,*) ' Enter KF code or 0 to stop'
         read(*,*)  kf
         if(kf .eq. 0) stop
         call ckf2cos(kf, code, subcode, charge)
         write(*,*)  'kf=',kf, ' c,s,cg=', code, subcode, charge
      enddo
      end
