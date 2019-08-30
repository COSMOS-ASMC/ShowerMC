!         
!               variables generated during incident angle sampling
          complex*16 Obsvhour      !  observation hour for point source
          real*8 Cspsmx, Cspsmn    ! cos(zenith) of point source
! 
          type(coord):: AngleAtObsCopy, DcAtObsXyz
!        
!           AngleAtObsCopy  :  direction cos at the deepest Obssite
!                              in 'det' system
!           DcAtObsXyz :       its transfroamtion to "xyz" system
!
          type(track):: IncidentCopy
          common /Zincidentv/ IncidentCopy, AngleAtObsCopy,
     *    DcAtObsXyz,  Obsvhour, Cspsmx, Cspsmn
