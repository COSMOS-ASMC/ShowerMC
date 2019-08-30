      subroutine iofortran
      integer kk
      write(*,*)  " this is output from fortran"
      write(0,*)  " this is  error output from fortran"
      read(*, *) kk   
      write(0,*)  " the value of k=", kk
      end

