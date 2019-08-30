      subroutine ddd
      character*30  string
      string = "abc  xyz"
      call xxx(string)
      write(*,*) 'after xxx; in Fortran ',string 
      end
