      integer max_coord_types
      parameter (max_coord_types = 4 )
      character*4 coord_types(max_coord_types)
      data coord_types(1)/'xyz'/, 
     *   coord_types(2)/'llh'/, coord_types(3)/'sph'/,
     *   coord_types(4)/'det'/
