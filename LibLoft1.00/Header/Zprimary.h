#include "Zmaxdef.h"

!                this is Primary class data type definition
!         Note that primary angle information is not here.
!
         integer NoOfSymbols,   ! max # of primary type symbols
     *           maxSegments,    ! max # of segments in each primary
     *           maxNoOfComps,  ! max # of components usable
!                                        at a time 
     *           maxErgUnit       ! max # of energy unit symbols
         parameter (NoOfSymbols = 67,
#ifdef MAX_SEGMENTS
     *              maxSegments = MAX_SEGMENTS,
#else
     *              maxSegments = 40, 
#endif
#ifdef MAX_NO_OF_COMPS
     *              maxNoOfComps = MAX_NO_OF_COMPS,
#else
     *              maxNoOfComps = 8,
#endif
     *              maxErgUnit=7 )
!
         type component     ! 1 component of 1ry 
           sequence
             integer label          ! composition label number
             character*16 symb        ! 'P', 'gamma' etc  12->16 (2019/Aug/29)
             character*3 eunit        ! 'GeV' etc
             character*4 etype        ! 'KE/n' etc
             character*1 diff_or_inte   ! 'd' or 'i'
             real*8 flatterer               ! dI/dE*E**flatterer
             real*8 cut                  !  lower cut off.
             real*8 cut2                 ! upper cut off.
             real*8 energy(maxSegments+1)  ! segment left energy
             real*8 flux(maxSegments+1)    ! input flux
!                                            dI/dE * E**flatterer
!                          above: from  input table directly
!                          below: made by subroutines
             integer code, subcode, charge  ! particle code
             real*8 togev

!!!!!        real*8 diff(maxSegments+1)      ! diff. flux
             real*8 norm_inte(maxSegments+1) ! normalized 
!                      integral flux > E at segment left value
             real*8 beta(maxSegments+1)    ! dI/dE=(true flux)=
!                                            const*E**(-beta)
             integer no_of_seg         ! no of segments given
             real*8 inte_value       !  integral flux from min. E
             real*8 emin, emax     ! min and max energy defined
!            integer histnbin
!            parameter (histnbin = 30)
!            type(histgl):: comphist
         end type component    
!         *******************************************
         type primaries   ! 1 set of primaries.
           sequence
              type(component):: each(maxNoOfComps)
              real*8 cummInteFlux(maxNoOfComps)
              integer no_of_comps    ! how many diff. compositions
              integer NoOfSamplings  ! total number of samplings including ones discarded by cutoff
              integer NoOfSampComp(maxNoOfComps, 2)  ! 1 is for number of sampling including
                                                      ! discarded ones due to cutoff. 2 is only for
                                                      ! employed ones.
!                        after a sampling of a 1ry, the following 
!                        is fixed.
              integer label            ! sampled primary label
              real*8 sampled_e     ! sampled energy(or momentum) as
!                                    defined in etype. If this is
!                                    in total energy in GeV, it is the
!                                    same as particle.fm.e
               type(ptcl):: particle
         end type primaries


