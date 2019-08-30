!       Zobsvp.h---parameters to be given by input.
!       This must be preceded by Zobs.h

!	(->	---------------------------------------------------

         real*8  HeightList  !1  Height of observation levels in m. This is  made from DepthList internally. 
                            ! I.e., this one is usually not an input. However, if the DepthList values are 
                            ! negative, this is used as input and corresponding DepthList is computed internally.
        real*8  DepthList   !1	Depth List of Observation level in kg/m$^2$. If $< 0$, HeightList has priority. 
                            !  (See HeightList)
        real*8  ASHeightList	!1  This is HeightList for Air Shower observ.  Used only if Generate contains
                            !  "as". See  HeightList.
        real*8  ASDepthList     !1  This is DepthList for AS observation.  Used only if Generate contains 
                            ! "as". See DepthList.
        real*8  LatitOfSite     !1  Latitude of the deepest observation level in degree.  East is positive.
        real*8  LongitOfSite    !1  Longitude of the deepest observation level in degree.  North is positive.

     	real*8  DtGMT           !1  Difference of the local time of the observation place from GMT (hour).
         real*8  YearOfGeomag    !1  Like 1999.5. Year when Geomagnetic field is to be calculated.
         integer ObsPlane        !1    How to observe particles. \newline
                                !    0$ \Rightarrow $ no detector plane is used for observation. BorderHeightL
                                !    and BorderHeightH are used to detect particles. This is for, say, neutrino
                                !    observation. See BorderHeight{L,H}. However, the primary is directed to
                                !    the deepest depth.  \newline
                                !    1,-1$ \Rightarrow $ detector at the observation place is horizontal. Note 
                                !    that the horizontal means not tangential plane, but rather a spherical surface \newline
                                !    2,-2$ \Rightarrow $ detector is perpendicular to the primary. \newline
                                !    3$ \Rightarrow $ spherical observation. See text. \newline
                                !    For ObsPlane={1,2}, the user observation routine will receive coordinate values in
                                !    the corresponding detector system. However, if it is 0, 3 or negative, Exyz values
                                !    are obtained.
        integer NoOfSites2    !2   No of Sites for particle observation; not to be touched; for skeleton/flesh use.
         real*8 XaxisFromSouth   !2 Angle between the horizontal detector X-axis and the south(deg). + is counter
                                ! clockwise.  If $|$XaxisFromSouth$| > 360$, it is computed so that the direction is
                                ! to the magnetic east at the deepest observation point. Default is 361.
!	<-)	--------------------------------------------

   

        common /Zobsc/
     *	 HeightList(maxNoOfSites),
     *   DepthList(maxNoOfSites),
     *   ASHeightList(maxNoOfASSites),
     *   ASDepthList(maxNoOfASSites),
     *   LatitOfSite, LongitOfSite, DtGMT,
     *   XaxisFromSouth, YearOfGeomag,
     *   ObsPlane, NoOfSites2

