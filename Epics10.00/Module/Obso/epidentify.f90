!  test program is in prog/Light/Test/listAndPointer.f90
!  see top of the test prog.

!      call epiniIdentify
!      call epinsertItem(id, comp)
!      .....
!      call epid2cnList(id, pcomp, n)
!      call epcn2idList(cn, ids, n)
!  When config data is read, cn2idList and id2cnList are made
!  to be available. The use may access only these two pointer
!  arrays. 
!           cn2idList(i)
!           i          cn id(:)                                nid
!           1           5 TT33                                 1  
!           2          22 S1          trigger     layer1       3
!           3          23 S2                                   1     
!           4          34 S1          layer1                   2
!           5          66 TT33                                 1
!           6          90 S1                                   1
! -----------------------------------------
!          id2cnList(i)
!           i id                    cn(:)                         ncn
!           1 S1                    22          34          90     3
!           2 S2                    23                             1
!           3 TT33                   5          66                 2
!           4 layer1                22          34                 2
!           5 trigger               22                             1
!
!  
module  epidentify
  integer,parameter::idlen=12
  integer,parameter::snamelen=16
  type idlist
     !id may be such as  layer1 trigger S1 S2 pmt APD1 pd22 anti
     
         ! sideanti scifix1 scifiy2 trigscintiA trigscinitB
         !  
     character(len=idlen)   :: id
     integer             :: compNO
     character(len=snamelen) :: sname  ! struc name or subd name
     type(idlist), pointer :: next
  end type idlist

  type idblock
     type(idlist), pointer :: entrance, new, p
     type(idlist), pointer :: entrancecn,  newcn, pcn
  end type idblock


  type cnVsId
     integer::cn    ! comp. #
     integer::nid   ! # of id's belonging to this copm.
     character(len=idlen),pointer::id(:) !  : = 1, nid
  end type cnVsId

  type idVsCn
     character(len=idlen),pointer::id  ! id
     integer::ncn                      ! # of cn belonging to this id
     integer,pointer::cn(:)            ! : = 1, ncn
  end type idVsCn

  character(len=idlen),pointer::allid(:)
  integer,pointer::allcn(:)

  type(idVsCn),pointer::id2cnList(:)
  integer :: id2cnListSize      !  size of id2cnList

  type(cnVsId),pointer::cn2idList(:)
  integer :: cn2idListSize

  type(idblock),pointer::subdids(:)

contains

subroutine epinsertItem(block, id, comp, sname)
  implicit  none

#if defined (Solaris)
  type(idblock),pointer::block 
#else
  type(idblock),intent(inout),pointer::block 
#endif


  character(len=idlen),intent(in)::id
  integer,intent(in)::comp  ! component # associated to the id
  character(len=*),intent(in)::sname

  allocate(block%new)

  block%p => block%entrance
  do while( associated(block%p%next))            
     if( id  < block%p%next%id ) exit       
!     if( comp < block%p%next%compNo ) exit       
     block%p => block%p%next
  end do
  block%new = idlist(id, comp, sname,  block%p%next)         
  block%p%next => block%new                        
!   entrance  next0    next1    next2
!     ??       apd          pd      S1
!     ??        20          12        3  
!   next0     next1       next2    null
!-------------
!   if (id,comp)=(aa,15);  exit with p=entrance.  new=(aa, 15,next0)
!                                               entrance%next=new
!   
!   if (id,comp)=(b,17);   exit with p=next0.  new=(aa, 15,next1)
!   if (id,comp)=(S2,24);  exit with p=null.  new=(S2, 24, null)

  allocate(block%newcn)   ! for cn sorted

  block%pcn => block%entrancecn
  do while( associated(block%pcn%next))            
     if( comp < block%pcn%next%compNo ) exit       
     block%pcn => block%pcn%next
  end do
  block%newcn = idlist(id, comp, sname, block%pcn%next)         
  block%pcn%next => block%newcn                        

end subroutine epinsertItem

subroutine epiniIdentify(block)
  implicit none

#if defined (Solaris)
  type(idblock),pointer::block 
#else
  type(idblock),intent(inout),pointer::block 
#endif

  if( associated( block%entrance)) then
     nullify( block%entrance)
  endif
  allocate(block%entrance)
  nullify(block%entrance%next)   

  allocate(block%entrancecn)
  nullify(block%entrancecn%next)   
end subroutine epiniIdentify

subroutine epidChgNum(block, nnow, nnew)
  implicit none

#if defined (Solaris)
  type(idblock),pointer::block
#else
  type(idblock),intent(inout),pointer::block
#endif

  integer,intent(in):: nnow ! cn in the list
  integer,intent(in):: nnew ! # to replace nnow

  block%pcn => block%entrancecn
  do while( associated(block%pcn%next))            
     if( nnow < block%pcn%next%compNo ) exit       
     if( nnow == block%pcn%next%compNo ) then
        block%pcn%next%compNo = nnew
     endif
     block%pcn => block%pcn%next
  end do

  block%p => block%entrance
  do while( associated(block%p%next))            
     if( nnow == block%p%next%compNo )  then
        block%p%next%compNo = nnew
     endif
     block%p => block%p%next
  end do
end subroutine epidChgNum


!     count # of components with given id
!     and return the # and component list
subroutine epid2cnList(block, id, pcomp, n)
  implicit none

#if defined (Solaris)
  type(idblock),pointer::block
#else
  type(idblock),intent(in),pointer::block
#endif

  character(len=idlen),intent(in)::id

#if defined (Solaris)
  integer,pointer::pcomp(:)   ! comp # list
#else
  integer,intent(out), pointer::pcomp(:)   ! comp # list
#endif

  integer,intent(out):: n           !  # of comp.
!    how many list with give ID
  n = 0
  block%p => block%entrance%next
  do while( associated(block%p) )     
     if(block%p%id == id) then
        n = n + 1
     endif
     block%p => block%p%next
  end do

!      allocate pcomp
  if(n> 0) then
     allocate( pcomp(1:n) )  
     n = 0
     block%p => block%entrance%next
     do while( associated(block%p) )     
        if(block%p%id == id) then
           n = n + 1
           pcomp(n) = block%p%compNo
        endif
        block%p => block%p%next
     end do
  endif
end subroutine epid2cnList
!     count # of ID with given comp#
!     and return the # and ID list
subroutine epcn2idList(block, cn, ids, n)
  implicit none

#if defined (Solaris)
  type(idblock),pointer::block
#else
  type(idblock),intent(in),pointer::block
#endif

  integer,intent(in):: cn  ! comp. #

#if defined (Solaris)
  character(len=idlen),pointer::ids(:)   ! id list
#else
  character(len=idlen),intent(out), pointer::ids(:)   ! id list
#endif

  integer,intent(out):: n           !  # of id's

  n = 0
  block%pcn => block%entrancecn%next
  do while( associated(block%pcn) )     
     if(block%pcn%compNo == cn) then
        n = n + 1
     endif
     block%pcn => block%pcn%next
  end do

!      allocate ids
  if(n> 0) then
     allocate( ids(1:n) )
     n = 0
     block%pcn => block%entrancecn%next
     do while( associated(block%pcn) )     
        if(block%pcn%compNo == cn) then
           n = n + 1
           ids(n) = block%pcn%id
        endif
        block%pcn => block%pcn%next
     end do
  endif
end subroutine epcn2idList

subroutine epgetAllid(block, m)
  implicit none

#if defined (Solaris)
  type(idblock),pointer::block 
#else
  type(idblock),intent(in),pointer::block 
#endif

  integer,intent(out):: m  ! possible  max number
  m = 0
  if(associated(allid)) then
     nullify(allid)
  endif
  block%p => block%entrance%next
  do while( associated(block%p) )     
     m = m + 1
     block%p => block%p%next
  end do
  allocate( allid(m) )

  m = 0
  block%p => block%entrance%next
  do while( associated(block%p) )     
     m = m + 1
     allid(m) = block%p%id
     block%p => block%p%next
  end do
end subroutine epgetAllid

subroutine epgetAllcn(block, m)
  implicit none

#if defined (Solaris)
  type(idblock),pointer::block
#else
  type(idblock),intent(in),pointer::block
#endif

  integer,intent(out):: m  ! possible  max number
  m = 0
  if( associated(allcn)) then
     nullify(allcn)
  endif

  block%pcn => block%entrancecn%next
  do while( associated(block%pcn) )     
     m = m + 1
     block%pcn => block%pcn%next
  end do
  allocate( allcn(m) )

  m = 0
  block%pcn => block%entrancecn%next
  do while( associated(block%pcn) )     
     m = m + 1
     allcn(m) = block%pcn%compNo

     block%pcn => block%pcn%next
  end do

end subroutine epgetAllcn

!   get list of ID
subroutine epgetIdList(block, pid, n)
  implicit none

#if defined (Solaris)
  type(idblock),pointer::block
#else
  type(idblock),intent(in),pointer::block
#endif


#if defined (Solaris)
  character(len=idlen), pointer::pid(:)   
#else
  character(len=idlen),intent(out), pointer::pid(:)   
#endif

  integer,intent(out):: n           !  # of diff. id's 
!    how many diff. id

  character(len=idlen),pointer::idx
  integer:: m, i
  !        a a b b b c d e e f 
  
  call epgetAllid(block,m)
  
  if( m > 0 ) then
     n = 1
     i = n+1
     idx => allid(n)
     do while (i <= m )
        if( idx /= allid(i) ) then
           n = n + 1
           idx => allid(i)
           allid(n) = idx
        endif
        i = i + 1
     enddo
  else
     n= 0
  endif


  if(n> 0) then
     allocate( pid(1:n) )
     pid(1:n) = allid(1:n)
     deallocate( allid )
     allocate( allid(1:n) )
     allid(:) =  pid(:)
  endif
end subroutine epgetIdList

!   get list of cn
subroutine epgetCnList(block, pcn, n)
  implicit none

#if defined (Solaris)
  type(idblock),pointer::block
#else
  type(idblock),intent(in),pointer::block
#endif


#if defined (Solaris)
  integer, pointer::pcn(:)   
#else
  integer,intent(out), pointer::pcn(:)   
#endif

  integer,intent(out):: n           !  # of diff. cn's
!    how many diff. cn

  integer,pointer::idx
  integer:: m, i
  !        a a b b b c d e e f 
  
  call epgetAllcn(block, m)

  if( m > 0 ) then
     n = 1
     i = n+1
     idx => allcn(n)
     do while (i <= m )
        if( idx /= allcn(i) ) then
           n = n + 1
           idx => allcn(i)
           allcn(n) = idx
        endif
        i = i + 1
     enddo
  else
     n= 0
  endif


  if(n> 0) then
     allocate( pcn(1:n) )
     pcn(1:n) = allcn(1:n)
     deallocate( allcn )
     allocate( allcn(1:n) )
     allcn(:)=pcn(:)
  endif
end subroutine epgetCnList


subroutine epgetIdVsCn(block)
  implicit none
!   get all ID list and cn's for each of the ID
!   id2cnList(i)%id    id
!   id2cnList(i)%ncn   # of cn with this id
!   id2cnList(i)%cn(j). j = 1, ncn
!   i=1,n (n= id2cnListSize)


#if defined (Solaris)
  type(idblock),pointer::block
#else
  type(idblock),intent(in),pointer::block
#endif

  character(len=idlen),pointer::pid(:)
  integer:: m, i
  integer,pointer:: cns(:)         ! cn list

  integer::n 

  call epgetIdList(block, pid, n)  

  if( associated(id2cnList) ) then
     nullify(id2cnList)
  endif
  allocate(id2cnList(1:n))
  do i = 1, n
     call epid2cnList(block, pid(i), cns, m)
     allocate(id2cnList(i)%cn(1:m))
     id2cnList(i)%id => pid(i)
     id2cnList(i)%ncn = m
     id2cnList(i)%cn(:) = cns(:)
  enddo
  id2cnListSize = n
end subroutine epgetIdVsCn

subroutine epgetNoVsId(block)
  implicit none
!   get all comp. # list and id's for each of the comp.
!   cn2idList(i)%cn      comp. # 
!   cn2idList(i)%nid    # of id's
!   cn2idList(i)%id(j). j = 1, nid, i=1,n (n=cn2idListSize)


#if defined (Solaris)
  type(idblock),pointer::block
#else
  type(idblock),intent(in),pointer::block
#endif

  integer :: n           !  # of diff. cn's


  integer,pointer::pcn(:)
  integer:: m, i
  character(len=idlen), pointer::ids(:)   ! id list  

  call epgetCnList(block, pcn, n)  
  if(associated(cn2idList)) then
     nullify(cn2idList)
  endif
  allocate(cn2idList(1:n))
  do i = 1, n
     call epcn2idList(block, pcn(i), ids, m)
     allocate(cn2idList(i)%id(1:m))
     cn2idList(i)%cn=pcn(i)
     cn2idList(i)%nid = m
     cn2idList(i)%id(:) = ids(:)
  enddo
  cn2idListSize = n
end subroutine epgetNoVsId


end module epidentify
