#ifndef  EP_MAX_DEF
#define  EP_MAX_DEF
#define  EPMAX_STACK 3000
#define MAX_MEDIA  20

       !    max number of different media (such as Pb, Si, BGO)

       !         #define MAX_MATRESKA  1000 obsolete

       !    max number of matreska, i.e, component that containes other
       !    components                        

#define MAX_COMPONENT  100000

       !     max number of components definable in Epics

#define MAX_MEDIA_DIR  3

       !     max number of directories where media file is  serached for.
       !     if this is increased, some changes in EPPARM is needed
#define MAX_ATTR   150000
       !     max total number of attributes of the all the components
       !    = max component nubmer x average attributes for 1 component.
       !     ( box has 3, sphere 1,  cyl 2, pipe 3...)
#define MAX_IATTR  50

      !  max number of attributes of an individual component (including
      !  derived attribues in the program).

#define MAX_VERTEX  40
       !     max number of vertexes of a component treatable.
       !     If a volume has a  curved surface (convex), vertex
       !     information is not used; NVTX=0 
#define MAX_TOTCHA 2048
       !    max number of total characters in a line.
       !    A line includes concatinated lines by backslash
#define MAX_SUBD  2000  
       !     max number of sub-detectors definable
#define MAX_COMPONENT_IN_SUBD  200000
       !     max number of components in a sub detector
#define MAX_SUBD_REF  2000
       !     max number of sub detectors that can be referred

#define MAX_DISPLAY  800

       !    max number of components displayable with a moderate cpu time

#define MAX_ERR_DIR_COS  3.d-3
       !     max error permissilbe in sum of direction cos^2

#define EPS_FOR_90  1.d-5

       !     scaler prod. of orthogonal unit vector must be <  this.
#define MAX_NEW_STRUC  20

       !     max number of new structures definable by the user

#define MAX_ELEMENTS  20

       !     max number of elements(atoms) in a medium.

#define MAX_EQUATE  800
       !    max number of #equate definable
#define MAX_CN_AREA  500000
       !  max total number of comp. that appear as contained by
       !    other comp.s
#define EPMAX_COMMENT  200
       !   max number of comment lines keptable and printable later
#define EPMAX_UHOOKC 10
#define EPMAX_UHOOKI 10
#define EPMAX_UHOOKR 10
      !     number of user definable character, integer, real*8
      !     variables in sepicsfile
#define MAXERGNODE_MUPAIR 10
      !  max # of energy nodes at which exact cross-sections 
      !  are given for pair creation by muon at low eneries
      !  (3 GeV to 1500 GeV)
#define MAXHEAVYCHG  30
      !  max heavy charge treatable by SRIM data when dE/dx is calculated
#define MAX_SRIMMEDIA  3
      !  max # of  medias for which SRIM data can be applied. 
#define MAX_STRUCCHR  16
      !  from 12 to 16  v9.164
#define MAX_MEDIANAMELENG 8
#endif
