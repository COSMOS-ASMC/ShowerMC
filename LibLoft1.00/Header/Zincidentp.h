!                parameters for primary angle sampling
!  (->	---------------------------------------------------

           complex*16  CosZenith   !1  Range of cos(zenith angle). Say, (0.5, 1.0). Used when Za2ry is 'is' 
                                   !   If ObsPlane=3 (spherical), real(CosZenith) must be $>$0, and means
                                   !   the zenith angle range at the incident point (not in Exyz system).
                                   !   In that case, azimuth is 0 to 2pi.
           complex*16  Azimuth     !1  Range of azimuthal angle in deg. Say, (0, 45). Default is (0,360).
                                   !   Can be such as (300., 390.).  Used when Za1ry is 'is'\newline
                                   !   If ObsPlane=3 (spherical), this is used to show the half opening angle
                                   !   range where the primary injection position is uniformly distributed 
                                   !   on a sphere.  The center of the opening angle is (Latit, Longit, HeightOfInj).
                                   !   In this case, for the upper opening angle,  min( Imag(Azimuth),180.) is used.
           character*4 Za1ry       !1  Specify the primary angle sampling method by one of 'is', 'ps' or 'aps'.\newline
                        !   "is" is isotropic. The range is by CosZenith.\newline
                        !   "ps" is for point source (See also SourceDec)\newline
                        !   "aps" is around point source (See also SourceDec and  Ddelta) \newline
                        !   If ObsPlane=3 (spherical), this must be "is".
           real*8  SourceDec    !1  Source declination of point source.(deg)
           real*8  Ddelta       !1  SourceDec $\pm$ Ddelta is the region for 'aps' (deg).
           real*8  HeightOfInj  !1  The vertical height of primary injection point (m).
                                !   If this is $<$ deepest obs. level and zeinth angle of primary is $< 0$, 
                                !   the primary is  assumed to be upgoing even if Reverse =0.
                                 !   NOTE: BorderHeightH must be given explicitly in this case.
           real*8  OffsetHeight !2  The vertical offset  height from  the deepest detector. 
                                ! The  primary is directed to this height above the detector.
                                ! If ObsPlane is  3 (spherical), not used.

!	<-)	----------------------------------------------------

           common /Zincident/ Azimuth, CosZenith, SourceDec, Ddelta,
     *      HeightOfInj,  OffsetHeight
           common /Zincidentc/ Za1ry
