!         magnetic field component expressions 
!    this should be placed at last part of the declaration
      integer max_mag_types
      parameter (max_mag_types = 3 )
      character*4 mag_types(max_mag_types)
      data mag_types(1) /'xyz'/,
     *   mag_types(2)/'hva'/, mag_types(3)/'ned'/
