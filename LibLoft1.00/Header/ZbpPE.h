!            photo electric effect const.
       integer maxNoOfShells
       parameter(maxNoOfShells=20)
       type photoE 
       sequence
          real*8 b0        
          real*8 b1
          real*8 b2
          real*8 fa
          real*8 a
          real*8 b    
          real*8 p
          real*8 l
          real*8 ek
          real*8 cr     ! correction factor for low Z
                        ! cross-section at Eg=0.2 MeV
                        ! this is multiplied to old cross-section
	  real*4 shellE(maxNoOfShells)  ! from v8.70
	  integer noOfShells
       end type photoE
  
