    program main
      use NucMassFormula
      implicit none
      integer:: A, Z
      write(0,*) 'Enter A,Z'
      read(*,*) A, Z
      write(0,*) ' A=', A, ' Z=',Z
      write(0,*) ' Enter formula'
      read(*,*)  formula
      write(0,*) 'mass excecss in MeV',  cNucMassExcessMeV(A,Z)
      write(0,*) 'mass excess in A.M.U', cNucMassExcess(A,Z)
      write(0,*) 'binding Energy in MeV', cNucMassBE(A,Z)
      write(0,*) 'binding Energy in MeV/A', cNucMassBE(A,Z)/A
      write(0,*) 'mass defect in MeV', cNucMassDefect(A,Z)
      write(0,*) 'Gp ',  cProtonPhotoDissoci(A,Z)
    end program main
