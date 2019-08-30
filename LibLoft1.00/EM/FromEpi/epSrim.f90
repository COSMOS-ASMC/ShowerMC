! this manages SRIM data when it is availabe for dE/dx of heavy ion
! (proton may be included)
#include <ZepMaxdef.h>
module srimdata
  type stpw
     integer::size   ! size of Ekt, dEdx
     real(8),allocatable::Ekt(:)   ! total kinetic energy of heavy ion in GeV
     real(8),allocatable::dEdx(:)  ! dE/dx (GeV/(g/cm2) for that Ekt
  end type stpw

  type srimmedia
     integer::exist
     character(8)::name
     character(150)::path
     type(stpw)::c(MAXHEAVYCHG)
  end type srimmedia

  type(srimmedia):: srimno(MAX_SRIMMEDIA)

  contains
    subroutine epSrimChk(iowk, dir, name, index)
!         check if given media has SRIM data, if not  
!         index is given -1
!         else if this is n-th srim data, give n.
      implicit none
      integer,intent(in)::iowk !# to be used for  working disk 
      character(*),intent(in)::dir ! Media directory
      character(*),intent(in)::name ! Mddia name
      integer,intent(out)::index   !  media.srim. if SRIM
                    ! data exist, sequence no. is given else -1
      

      integer,save::seqno=0
      integer::ios, ios2
      character(150)::path
      character(2)::Z
      character(8)::nameX=" "   ! if this value is not give, some strange
                                ! thing will happen at character comparison
      integer::i


      index = -1
!         Z=2  is mandatory ion, so first check its existence
      path = dir//"/"//trim(name)//"_srim/"//"srim.02"
      
      call copenf(iowk, path,ios)
      if( ios /= 0 ) then
!            try lower case name
         call c2lowerCase(trim(name), nameX)
         if( trim(nameX) == trim(name) ) then
            call c2upperCase(name, nameX)
         endif
         
         path = dir//"/"//trim(nameX)//"_srim/"//"srim.02"
         
         call  copenf(iowk, path, ios)
      endif

      if(ios == 0 ) then
         close(iowk)
         if( seqno >= MAX_SRIMMEDIA) then
            write(0,*) 'SRIM data for media=',name
            write(0,*) 'will not be used. Too many srim data for'
            write(0,*) 'different media.  Increase max_srimmedia'
            write(0,*) 'in ZepMaxdef%h'
            stop 1111
         elseif(index == -1) then
            seqno = seqno + 1
            index = seqno
            srimno(index)%name = name
            srimno(index)%path = dir
         endif

         do i = 1, MAXHEAVYCHG
            path = dir//"/"//trim(name)//"_srim/"//"srim."    
            write(Z,'(i2.2)') i
            path = trim(path)//Z 
            call copenf(iowk, path, ios2)
            if(ios2 == 0 ) then
               close(iowk) 
               srimno(index)%c(i)%size = 0  ! later if dE/dx is requested
                          ! this is made to be actul size of aray
                          !  and  stpw data is  allocated and read.
            else
               srimno(index)%c(i)%size = -1  ! not data for Z=i
            endif
         enddo
      else
         !   no srim data is assumed
      endif
    end subroutine epSrimChk
    
    subroutine epSrimRead(iowk, smedia, chg)
      implicit none
      ! read srim data if data for charge Z=chg exists
      integer,intent(in):: chg   ! charge Z of heavy ion
      integer,intent(in):: iowk  ! io # to be used
      type(srimmedia),intent(inout)::smedia  ! srim media
      
      integer:: i, ios
      character(2)::Z

      integer::nline
      real(8)::x, y

      if( smedia%c(chg)%size == 0) then
         write(Z,'(i2.2)') chg
         call copenf(iowk,  &
          trim( smedia%path )//"/"//trim( smedia%name)//"_srim/"//"srim."//Z, &
          ios)

         if(ios /= 0 ) then
            write(0,*) ' media=',smedia%name, 'should have srim data'
            write(0,*) ' for charge =',chg, ' but not found'
            stop 2222
         endif
         

         nline = 0
         do while( .true. )
            read(iowk, *, end=100) x, y
            nline = nline + 1
         enddo
100      continue
         smedia%c(chg)%size = nline

         allocate(  smedia%c(chg)%Ekt(nline) )
         allocate(  smedia%c(chg)%dEdx(nline) )
         rewind(iowk)
 
         do i = 1, nline
            read(iowk, *, end=200) x, y
            smedia%c(chg)%Ekt(i) = x
            smedia%c(chg)%dEdx(i) = y
         enddo
200      continue
         close(iowk)
      endif
    end subroutine epSrimRead
    
    subroutine epSrimdEdx(smedia, Ekt, Z, dedt)
      implicit none
      type(srimmedia),intent(in)::smedia  ! srim media
      real(8),intent(in):: Ekt ! total KE of heavy GeV
      integer,intent(in):: Z   ! charge Z of heavy ion
      real(8),intent(out):: dedt  ! GeV/(g/cm2)

      real(8):: error
!///////////
!      if(Z>30) then
!         write(0,*) Z, Ekt
!      endif
!/////////////

      call kpolintpFE(smedia%c(Z)%Ekt, 1,  smedia%c(Z)%dEdx, 1,  &
              smedia%c(Z)%size, 3,  Ekt, dedt, error)
!///////////
!      write(*,'(a, i3, 1p, 2g12.3)') 'srim ', Z, Ekt,  dedt 
!///////////////
    end subroutine epSrimdEdx

  end module srimdata
