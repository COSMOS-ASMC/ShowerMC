!  Manipulte transport XSections.
! Usage:
! 0)  use modTPXS must be given before implicit none
!     In the  eprogame, define a type TPXSconst
!     as
!         type(TPXSconst),pointer:: some
!     or
!         type(TPXSconst),pointer:: some(n)  n is one of  1,2...
! 1)  open a data file in which E,s0,s1,s2 are contained with 
!     logical dev. # io
! 2)  
!      allocate(some)
!        or
!      allocate( some(i) )    (here i= one of n)
!
!     and
!     call cReadTPXS(io, pm,some) 
!       or
!     call cReadTPXS(io, pm, some(i))  
!where integer  pm >0 or < 0 for e+ or e-
!    At this point some%s0(j) etc are not in log.  (j=1~some%n)
!   Energy is KEele(j), and its log is logKEele(j).
!  To prepare for cubic spline interpolation in log-log scale 
!     call cPrepIntpTPXS(some) etc shoule be called. now
!   some%s0(j), some%s1(j), some%s2(j) are log of s0,s1,s2
!    (log ; natural log)
!   soome%coefS0(j), some%coefS1(j), some%coefS2(j) (j=1, some%n=1)
!   are ready now.
!!

    subroutine cCountTPXSdata(io, pm, n)     ! now normally not used 
!        count # of data ( E, S0, S1, S2 )  in file specified by io
!   if a line in the file start with "#" or blank, it is regarded as
!   comment
      use modTPXS
      implicit none
      integer,intent(in):: pm  ! <0 for elec, >0 for e+
      integer,intent(in):: io  ! logical dev. #  file  is not
                     ! closed but rewind is performed
      integer,intent(out):: n  ! # of data (E S0, S1, S2) 
      
      integer:: nf

      n = 0
      do
         read(io,'(a)', end=100) line
         call ksplit(line, 1, 1, field, nf)  ! get first character
         if( field(1)=="#" .or. nf==0 ) cycle
         n = n + 1
      enddo
100   continue
      if( pm < 0  ) then
         if( nEneg /= n) then
            write(0,*)  &
                 '# of data by cCountTPXSdata != previous one=',nEneg
            write(0,*) ' present one is ',n
            write(0,*) ' filename =',trim(filename)
            stop
         endif
      elseif( pm > 0 ) then
         if( nEpos /= n ) then
            write(0,*)  &
                 '# of data by cCountTPXSdata != previous one=',nEneg
            write(0,*) ' present one is ',n
            write(0,*) ' filename =',trim(filename)
            stop
         endif
      else
         write(0,*) ' pm=', pm, ' invalid to cCountTPXSdata'
         stop
      endif
            
      rewind io
    end subroutine cCountTPXSdata
    
    subroutine cTPXSalloc(pm,  name)   ! now normally not used 
      use modTPXS
      implicit none
      integer,intent(in):: pm ! <0 for e-, >0 for e+
      type(TPXSconst),intent(out):: name

      integer:: n  ! size of allocation 
      if( pm < 0) then
         n = nEneg
      elseif( pm > 0 ) then
         n = nEpos
      else
         write(0,*) ' error pm=',pm, ' to cTPXSalloc'
         stop
      endif
      name%n = n
!!!!!!!!!!!
!      write(0,*) ' nEneg, nEpos=', nEneg, nEpos
!!!!!!!!
!      if( .not. allocated( KEele ) ) then 
!         allocate( KEele(n))
!         KEele(:) = 0.
!         allocate( logKEele(n))
!      endif



!      if( .not. allocated( name%S0 ) ) then
!         allocate(name%S0(n))
!         allocate(name%S1(n))
!         allocate(name%S2(n))
!         allocate(name%A0(n))
!         allocate(name%A(n))
!         allocate(name%B(n))
!      !       cubic spline coef.
!         allocate(name%coefS0(n-1, 3))
!         allocate(name%coefS1(n-1, 3))
!         allocate(name%coefS2(n-1, 3))
!
!         allocate(name%coefA0(n-1, 3))
!         allocate(name%coefA(n-1, 3))
!         allocate(name%coefB(n-1, 3))

!      endif
    end subroutine cTPXSalloc

    subroutine cReadTPXS(io, pm, name)
      use modTPXS
      implicit none
      integer,intent(in):: io  ! logical dev # of the file. file will be closed on return
      integer,intent(in):: pm  ! >0 for e+ <0 for e-
      type(TPXSconst),intent(inout):: name

      integer:: n 
      integer:: i, nf
      real(8),save:: prevE
      i = 0
      n = nEneg
      if( pm > 0 ) n = nEpos

!!         next two not used now
!!      call cCountTPXSdata(io, pm, n)
!!      call cTPXSalloc(pm, name)
      do
         read(io,'(a)', end=100) line
         call ksplit(line, 1, 1, field, nf)  ! get first character
         if( field(1)=="#" .or. nf==0  ) cycle
         i = i + 1
         prevE= KEele(i)
         read(line,*) KEele(i), name%S0(i), name%S1(i), name%S2(i) 
         ! check energy
         if( prevE > 0. ) then
            if( abs( prevE/KEele(i) - 1.) > eps)  then
               write(0,*) ' Energy read=',KEele(i)
               write(0,*) ' previously stored =',prevE
               write(0,*) ' diff. large; cReadTPXS'
               stop
            endif
         else
!            logKEele(i) = log(KEele(i))
         endif
         if( i == n ) exit
      enddo
      name%n = n
!      name%S0(:) = name%S0(:)/1.0d-27  ! in mb
!      name%S1(:) = name%S1(:)/1.0d-27  ! in mb
!      name%S2(:) = name%S2(:)/1.0d-27  ! in mb
      close(io)
      return
100   continue
      write(0,*) &
      ' EOF is reached during data reading in cReadTPXS'
      write(0,*) ' last line read is ', line
      write(0,*) ' # of data already read is ', i
      write(0,*) ' # of data requested is ', n
      write(0,*) ' file may be ', trim(filename)
      stop
    end subroutine cReadTPXS
    
    subroutine cPrepIntpTPXS(name)
!         prepare for cubic spline interpolation
      use modTPXS
      implicit none
      type(TPXSconst),intent(inout):: name
      
      integer:: n, nc 
      if( logKEele(1) == 0. .and. logKEele(2) == 0.) then
         logKEele(:) = log( KEele(:))
      endif
!         take log of s0,s1,s2 
      name%logS0(:) = log( name%S0(:) )
      name%logS1(:) = log( name%S1(:) )
      name%logS2(:) = log( name%S2(:) )

      n = name%n
      nc = n-1
!         make cubic spline intp.
      call kcsplCoef(logKEele, name%logS0, n, name%coefS0, nc)
      call kcsplCoef(logKEele, name%logS1, n, name%coefS1, nc)
      call kcsplCoef(logKEele, name%logS2, n, name%coefS2, nc)
    end subroutine cPrepIntpTPXS

    subroutine cTPXS(name, cond, Ein, S0, S1, S2)
      use modTPXS
      implicit none
      type(TPXSconst),intent(in):: name ! 
      integer,intent(in):: cond  ! <0 log value of Si
                                 ! >0 non log value of Si
      real(8),intent(in):: Ein  ! e-/e+ energy in eV
      real(8),intent(out):: S0,S1,S2 !  transport XS (in  cm2) (or
                                    ! its  natural log

      integer::n
      real(8):: logKE

                             ! log of eV value if Ein is in GeV
      logKE = log(Ein) !      + 2.07232658369464111562d1

      n = name%n

      call kcsplIntp(logKEele, name%logS0, n, name%coefS0, n-1, logKE,&
           S0)
!      call kcsplIntp(KEele, name%S0, n, name%coefS0, n-1, Ein, S0)
     call kcsplIntp(logKEele, name%logS1, n, name%coefS1, n-1, logKE,&
           S1)
!!      call kcsplIntp(KEele, name%S1, n, name%coefS1, n-1, Ein, S1)
     call kcsplIntp(logKEele, name%logS2, n, name%coefS2, n-1, logKE,&
           S2)
!!      call kcsplIntp(KEele, name%S2, n, name%coefS2, n-1, Ein, S2)

      if(cond > 0 ) then
         S0 = exp(S0)
         S1 = exp(S1)
         S2 = exp(S2)
      elseif( cond == 0 ) then
         write(0,*) " error cond =0 to cTPXs"
         stop
      endif
    end subroutine cTPXS
    
