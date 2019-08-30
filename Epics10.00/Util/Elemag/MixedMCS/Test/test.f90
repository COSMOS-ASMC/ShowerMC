module modA
  implicit none
  integer:: ijk
  type abc
     real(8):: xyz(20,30)
  end type abc
end module modA
module bridge
  use modA
  implicit none
  type(abc),pointer::pname
end module bridge
  
program test
  use bridge
  implicit none
  type(abc),target:: name1
  type(abc),target:: name2
  name1%xyz(2,3)= 6
  name2%xyz(2,3) = -6
  write(0,*) ' pass 1'
  pname => name1
  write(0,*) ' pass 2'
  write(0,*) ' pname(2,3)=', pname%xyz(2,3)
!  if(associated(pname)) then 
!     deallocate(pname) 
     pname => name2
!  endif
  write(0,*) ' pname(2,3)=', pname%xyz(2,3)
  call test2(pname)
  call test3(name1)
  write(0,*) ' aft sub 3'
  write(0,*) pname%xyz(2,3)
end program test

subroutine test2(name)
  use modA
  implicit none
  type(abc):: name
  write(0,*) ' in sub 2 '
  write(0,*) name%xyz(2,3)
end subroutine test2
subroutine test3(name)
  use bridge
  implicit none
  type(abc):: name
  pname = name
end subroutine test3

  

  
  
