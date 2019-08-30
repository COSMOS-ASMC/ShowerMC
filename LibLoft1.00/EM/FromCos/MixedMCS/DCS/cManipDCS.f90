!  Manipulte transport XSections. ! not used now
subroutine cCountDCSdata(io, n)
!        count # of data ( mu, dcs )
!   if a line in the file start with "#" or blank, it is regarded as
!   comment
  use modTPXS
  use modDCS
  implicit none
  integer,intent(in):: io  ! logical dev. #  file  is not
  ! closed but rewind is performed
  integer,intent(out):: n  ! # of data lines (mu, dcs)
      
  integer:: nf

  n = 0
  do
     read(io,'(a)', end=100) line
     call ksplit(line, 1, 1, field, nf)  ! get first character
     if( field(1)=="#" .or. nf==0 ) cycle
     n = n + 1
  enddo
100 continue
!  if( nmu == 0 ) then
!     nmu = n
!  elseif( nmu /= n) then
!     write(0,*)  &
!          '# of data lines by cCountDCSdata != previous one=',nmu
!     write(0,*) ' present one is ',n
!     write(0,*) ' filename =',trim(filename)
!     stop
!  endif
!  rewind io
end subroutine cCountDCSdata

subroutine cDCSalloc(pm,  name)  ! now not used
  use modTPXS
  use modDCS
  implicit none
  integer,intent(in):: pm ! <0 for e-, >0 for e+
  type(DCSconst),intent(out):: name
  
  integer:: n  ! size of allocation 
  ! to fix 
  if( pm < 0) then
     n = nEneg
  elseif( pm > 0 ) then
     n = nEpos
  else
     write(0,*) ' error pm=',pm, ' to cTPXSalloc'
     stop
  endif
!  if(.not. allocated( muval )) then
!     allocate( muval(nmu))
!     muval(:) = -1.
!  endif
!  if( .not. allocated( name%dcs ) ) then
!     allocate(name%dcs(nmu, n))
!  endif
end subroutine cDCSalloc
! now not used
subroutine cPreReadDCS(io)
!  this is to  KEele
!  For KEele, we must allocate dummy TPXSconst

  use modTPXS
  implicit none
  integer,intent(in):: io

  integer:: icon, n
  type(TPXSconst),allocatable:: dummy(:)
!       use  Z= 1 and e-
  call cTPXS_Z2file(1, -1, filename)
  call copenf(io, filename, icon)
  if(icon /= 0 ) then
     write(0,*) trim(filename), ' caanot be opened in cPreReadDCS'
     stop
  endif
  allocate(dummy(1))
!  call cCountTPXSdata(io, -1, n)
!  allocate(dummy)
  call cReadTPXS(io,-1, dummy(1))
  close(io)
  deallocate(dummy)
!    e+
!  call cTPXS_Z2file(1, 1, filename)
!  call copenf(io, filename, icon)
!  if(icon /= 0 ) then
!     write(0,*) trim(filename), ' caanot be opened in cPreReadDCS'
!     stop
!  endif
!  call cCountTPXSdata(io, 1, n)
!  close(io)
end subroutine cPreReadDCS

subroutine cReadDCS(io, pm, ntherg, name)
  use modTPXS
  use modDCS
  implicit none
  integer,intent(in):: io  ! logical dev # of the file. file will be closed on return
  integer,intent(in):: pm  !  -1 for e-, 1 for e+
  integer,intent(in):: ntherg  !  dcs(:, ntherg) is target
  type(DCSconst),intent(inout):: name ! name%dcs is target

  real(8):: dummy
  integer:: n 
  integer:: i, nf
  real(8),save:: prevmu
  real(8):: val8
  
!  if(nmu == 0 ) then
!     call cCountDCSdata(io, n)  ! n should be nmu
!  endif

!  call cDCSalloc(pm, name)
  i = 0
  do
     read(io,'(a)', end=100) line
     call ksplit(line, 1, 1, field, nf)  ! get first character
     if( field(1)=="#" .or. nf==0  ) cycle
     i = i + 1
     prevmu = muval(i)
     read(line,*)  muval(i), val8
     name%dcs(i,ntherg) = val8
         ! check energy
         if( prevmu  >= 0. ) then
            if( abs( prevmu - muval(i) ) > eps)  then
               write(0,*) ' mu read=',muval(i)
               write(0,*) ' previously stored =',prevmu
               write(0,*) ' diff. large; cReadDCS'
               stop
            endif
         endif
         if( i == nmu ) exit
      enddo
      close(io)

      return
100   continue
      write(0,*) &
      ' EOF is reached during data reading in cReadDCS'
      write(0,*) ' last line read is ', line
      write(0,*) ' # of data already read is ', i
      write(0,*) ' # of data requested is ', nmu
      write(0,*) ' file may be ', trim(filename)
      stop
    end subroutine cReadDCS
    
    subroutine cPrepIntpDCS(name)
!         prepare for cubic spline interpolation
      use modTPXS
      use modDCS
      implicit none
      type(DCSconst):: name
      

      integer:: n, nc,  ie
      n = nmu
      nc = n-1
      name%logdcs(:,:) = log( name%dcs(:,:))
      logmuval(2:n)  = log(muval(2:n))
      logmuval(1) = log( muval(2)/5) ! not used
      logKEele(:) = log(KEele(:))
!      write(0,*) ' nmu =', nmu, ' nEneg=', nEneg
!      write(0,*) ' log muval'

!!cs      do ie =1, nEneg 
!!cs         call kcsplCoef(logmuval(2), name%logdcs(2,ie), nmu-1, &
!!cs           name%coefdcs(2,1,ie), nmu-2) 
!         call kcsplCoef(logmuval, name%logdcs(1,ie), nmu, &
!           coef, nmu-1)
!         write(0,*) ' moving coef '
!         name%coefdcs(:,:,ie) =coef(:,:)
!         coefAll(:,:, ie) =  coef(:,:)
!!cs      enddo
    end subroutine cPrepIntpDCS

    subroutine cDCS(name, mu,  Ein,  dcsval)
      use modTPXS
      use modDCS
      implicit none

      type(DCSconst),intent(in):: name ! 
      real(8),intent(in):: mu   ! (1-cos)/2
      real(8),intent(in):: Ein  ! e-/e+ energy in eV

      real(8),intent(out):: dcsval  ! ds/domega= 2pi ds/dmu/2 = 4pi ds/dmu   .   mb/sr

      real(8):: error, dcsval1, dcsval2
      real(8):: logKE, logmu
      integer:: ie

      logKE = log(Ein)
      if( mu < muval(2) ) then
         call kdwhereis(Ein, nEneg, KEele, 1, ie)
         if(ie >= nEneg) then
            dcsval = (name%dcs(2,ie)-name%dcs(1,ie))/(muval(2)-muval(1)) *mu &
            + name%dcs(1,ie) 
         else
!            write(0,*) ' Ein=', Ein , ' ie=', ie
!            dcsval = dcsval1
!         else
            dcsval1 = (name%dcs(2,ie)-name%dcs(1,ie))/(muval(2)-muval(1)) *mu &
            + name%dcs(1,ie) 
 
            dcsval2 = (name%dcs(2,ie+1)-name%dcs(1,ie+1))/(muval(2)-muval(1)) *mu &
            + name%dcs(1,ie+1) 
            dcsval = (dcsval2 - dcsval1)/(logKEele(ie+1)-logKEele(ie)) * &
              (logkE-logKEele(ie)) + dcsval1
         endif
      else
         logmu = log(mu)
!               next:     something wrong
!         call  kpolintp2(logmuval(2), 1, 0.d0, logKEele, 1, 0.d0, &
!              name%logdcs(2,1), nmu-1, nmu-1, nEneg, 3, 3, logmu, logKE, dcsval, error)
!         dcsval = exp(dcsval)
!               next:    something wrong
!!         call  kpolintp2(logmuval(2), 1, 0.d0, logKEele, 1, 0.d0, &
!!              name%dcs(2,1), nmu-1, nmu-1, nEneg, 3, 3, logmu, logKE, dcsval, error)
!               next: some zigzag but globally OK
!!         call  kpolintp2(muval, 1, 0.d0, KEele, 1, 0.d0, &
!!              name%dcs, nmu, nmu, nEneg, 2, 2, mu, Ein, dcsval, error)
!              next : much better than above.    
!         call  kpolintp2(logmuval, 1, 0.d0, logKEele, 1, 0.d0, &
!              name%logdcs,  nmu, nmu, nEneg, 2, 2, logmu, logKE, dcsval, error)
!         dcsval = exp(dcsval)
!              next : almost same as above. seems better than
!                   above?  max diff is half of the line thickness at
!                   u = arouund  0.1 or 0.2 .
         call  kpolintp2(logmuval, 1, 0.d0, logKEele, 1, 0.d0, &
              name%logdcs,  nmu, nmu, nEneg, 3, 3, logmu, logKE, dcsval, error)
         dcsval = exp(dcsval)
      endif

!        name%dcs, nmu, nmu, nEneg, 3, 3, mu, Ein, dcsval, error)  !NG
    end subroutine cDCS

    subroutine cDCSEgrid(name, mu,  ie,  dcsval)
!         same as cDCS but  energy is assumed to be on grid 
      use modTPXS
      use modDCS
      implicit none

      type(DCSconst),intent(in):: name ! 
      real(8),intent(in):: mu   ! (1-cos)/2  dcos = 2 dmu
      integer,intent(in):: ie   ! e-/e+ energy grid index

      real(8),intent(out):: dcsval  ! ds/domega= 2pi ds/dmu/2 =4pi ds/ dmu   .   mb/sr

      real(8):: error
      real(8):: logmu
!         without taking log, only linear intp. is ok. otherwise
!        negative dcs may appear for non grid point. (# of point is 2 )
      if( mu < muval(2) ) then
         dcsval = (name%dcs(2,ie)-name%dcs(1,ie))/(muval(2)-muval(1)) *mu &
            + name%dcs(1,ie) 
      else
         logmu = log(mu)
!!cs         call kcsplIntp(logmuval(2), name%logdcs(2,ie), nmu-1, &
!!cs              name%coefdcs(2,1,ie), nmu-2,logmu, dcsval)
!!cs           xa, xstep, ya, ystep, nt, m,  x, y, error)
         call kpolintpFE(logmuval(2), 1, name%logdcs(2,ie), 1, &
             nmu-1,  3, logmu, dcsval, error)
         dcsval = exp(dcsval)
     endif

!        name%dcs, nmu, nmu, nEneg, 3, 3, mu, Ein, dcsval, error)  !NG
    end subroutine cDCSEgrid

