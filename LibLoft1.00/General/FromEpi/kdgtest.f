!      integer dgpat, dgpos
!      logical kdgtest
!      integer kdgset, kdgclr
!      do while(.true.)
!         write(*,*) ' enter digit pattern and digit pos'
!         read(*,*)  dgpat, dgpos
!         write(*,*) ' dgpat=', dgpat, ' dgpos=',dgpos
!         write(*,*) ' digit is  on/off=',kdgtest(dgpat, dgpos)
!         write(*,*) ' set dgpos=', kdgset(dgpat, dgpos)
!         write(*,*) ' clear dgpos=',kdgclr(dgpat, dgpos)
!      enddo
!      end
      
!     ****************************************************************
!     *                                                              *
!     * kdgtest:examine if specified digit is on                     *
!     * kdgset:  set a specified digit on                            *
!     * kdgclr:  set a specified digit off                           *
!     *                                                              *
!
!
      logical function kdgtest(dgpat, dgpos)
      implicit none
!
      integer dgpat
      integer dgpos

      integer la

      la = dgpat/10**(dgpos-1)
      kdgtest = .not.( (la - (la/10)*10) .eq. 0 )

      end
!
!
      integer function kdgset(dgpat, dgpos)
      implicit none
      integer dgpat
      integer dgpos
      
      integer  k, amari, pow

      pow =  10**(dgpos-1)
      k = dgpat/pow
      amari = dgpat - k*pow
      kdgset = ((k/10)*10 + 1)*pow + amari
      end
!
      integer function kdgclr(dgpat, dgpos)
      implicit none
      integer dgpat
      integer dgpos
      
      integer  k, amari, pow

      pow =  10**(dgpos-1)
      k = dgpat/pow
      amari = dgpat - k*pow
      kdgclr =( (k/10)*10) *pow + amari
      end
