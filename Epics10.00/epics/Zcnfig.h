!   #include "Zep3Vec.h" must preceed this.
#include "Zepcondc.h"
#include "ZepMaxdef.h"
#include "Zcnfig2.h"

!           configuration of detector

        integer  maxCnArea   ! max size to memorize comp. numbers.
        parameter( maxCnArea = MAX_CN_AREA )

        integer maxcom        ! max comment lines printable
        parameter (maxcom = EPMAX_COMMENT)


        logical  Incused      ! #inc is used or not
	integer ncmax         ! max components
        parameter ( ncmax = MAX_COMPONENT )
!((((((((((
!        integer MaxMatreska   ! max number of matreska
!        parameter ( MaxMatreska = MAX_MATRESKA )
!))))))))))
        integer ncmaxInSubD   ! max number of components in a sub det.
        parameter ( ncmaxInSubD = MAX_COMPONENT_IN_SUBD )

        integer maxSubD       ! max number of sub detectors
        parameter ( maxSubD = MAX_SUBD )

        integer maxSubDRef    ! max number of sub detectors referrable
        parameter ( maxSubDRef = MAX_SUBD_REF )

!        integer NPreDefName   ! number of predefined structures
!        parameter ( NPreDefName = 5 )
  
!        integer MaxNewStruc   ! max number of new structures constructable
!        parameter ( MaxNewStruc = MAX_NEW_STRUC )

!           location of attribute characterizing a volume
        integer boxa, boxb, boxc  ! location for a,b,c edges of a box
        parameter (boxa =1, boxb = 2, boxc = 3)

        integer pipeir, pipeor, pipeh ! location for inner, outer and height 
                                      ! of a pipe
        parameter (pipeir = 1, pipeor = 2, pipeh = 3)
     
        integer cylr, cylh   ! location for radius and height of a cylinder
        parameter( cylr = 1, cylh = 2 )

        integer sphr        ! location for radius of a sphere
        parameter( sphr = 1 )

        integer prisma, prismb, prismc, prismh  ! a, b, c and h of a prism
        parameter( prisma=1, prismb=2, prismc=3, prismh=4 )

        integer maxattr      ! max total number of attributes of all 
                              ! the components
        parameter ( maxattr = MAX_ATTR ) 
 
        integer maxiattr    ! max number of attributes of an individul
                         ! comp.
        parameter (maxiattr = MAX_IATTR )

        integer maxvertex    ! max number of vertexes of a given volume
                             !  that can be treated in Epics
        parameter (maxvertex =MAX_VERTEX)
        integer maxtotcha    ! max characters in a line
	  !  (concatinated if back slash exists)
        parameter (maxtotcha = MAX_TOTCHA)
        real*8 MaxErrDirCos ! permissible error in sum of direction cos^2
        parameter( MaxErrDirCos = MAX_ERR_DIR_COS )

        real*8 EpsFor90     ! scaler prod. of orthogonal unit vectors must
                            ! be smaller than this.
        parameter( EpsFor90 = EPS_FOR_90 )
        integer maxEq       !  max num. of #equate  efinable
        parameter ( maxEq = MAX_EQUATE )
!  -------------------------------------------------------------------------

       type Component 
       sequence
         real*8 orgx, orgy, orgz,
     *          direc(9)
         real   offsetx, offsety, offsetz, MaxPathL
	  character(len=MAX_STRUCCHR):: struc
          character(len=MAX_MEDIANAMELENG):: matter
	 real  rhoc  ! rho calibration factor 
         real  EminG, EminE, RecoilE ! 
         real  LengGmin, LengChmin 
!            LengGmin: If length to the boundary is > this (cm)
!                      gamma is regarded absorbed safely.
!                 fixed by Ex = min(10keV, EminG) and attenuation
!                 length in the media.
!   
!            LengChmin: If length to the boundary is > this(cm)
!                   charged ptcl is regarded absorbed safely.
!                   (but, anti-ptcl, decayable ptcl is not).
!                 fixed by Ex=min(10keV, EminE-mass) and
!                  dE/dx of media at Ex. 
!            
!         integer  Nattributes,
!     *            NMatreska,
!     *            NContainer, vol,
         integer  NMatreska, vol

!               vol; vol+i is the i-th attribute of the component
!                    in Volat
! ((((((((((
!     *            Contains(MaxMatreska),
!     *            ContainsR(MaxMatreska),
         integer  Contains,
     *            ContainsR,
! )))))))))))
! (((((((((((
!     *            PContained(MaxMatreska),
     *            PContained, cn, chno
!            cn: component number assigned to this cmponent
! )))))))))))


       integer*2 NPContainer, CountIO,  CountDE, Modifier
!  >>>>>>>>>>>>>>>>>>light        
         integer*2  LightCompNo
!              LightCompNo; seq # given to light related comp
!   <<<<<<<<<<<<<<<     

         integer*2 level

         integer*2  Nattributes, NContainer
         integer*2 strucNo         ! number expressing the structure
         logical*2 rotation
         integer*2 subdidx         ! ^^^^^^  
         integer*2 fsubdc      ! counter for no. of subdetectors
                               ! included in the last main boday.
       end type Component

       type Detector 
       sequence
       type(Component):: cmp(ncmax)
          integer Cn2media(ncmax)
          integer nct
!          integer  nbox, ncyl, npip, nprs, nsph,           
!!!!          integer*2 nworld, nnew(MaxNewStruc)
!!!xxx          integer nworld, nnew(MaxNewStruc)
          integer nworld
       end type Detector

       type(Component):: SubdArea(ncmaxInSubD)

       type SubDetector 
       sequence
          integer  loc
          integer nct
!          integer  nbox, ncyl, npip, nprs, nsph,           
!!!!          integer*2 nworld, nnew(MaxNewStruc)
!!!xxx
          integer nworld
       end type SubDetector
       integer  cumsubdloc

       real*8  PlusF, EqualF, MinusF
       integer NotGiven
       parameter( NotGiven = -9999999 )
       integer*2  EqualFshort
       type(ep3Vec)::  XYZthick(maxSubD)  ! --> use  real(8)
       real(8)::ORGsubd(3,maxSubD)  

!       ********************
       type(Detector):: Det
       type(SubDetector):: SubD(maxSubD)
!       ********************



	integer Nfield, NsubD, mode, SubDUsed, coment,
     *          PosInDet(maxSubDRef), SubDNumb(maxSubDRef),
     *          SumOffset(maxSubDRef), wrtcom, comcounter,
     *          comflag(maxcom), comloc 

!
!   #if defined IBMAIX
!        character*500  confdata
!   #else
        character(maxtotcha):: confdata
!   #endif
        character*80 comarea(maxcom)
        character*24  Field(maxiattr), FieldAsItis(maxiattr)
        character*4 form, OrgEqBy
!        character*8 SubDName(maxSubD),  PreDefName(NPreDefName)
        character(len=MAX_STRUCCHR) SubDName(maxSubD)  ! 16--> MAX_STRUCCHR v9.164
        integer SubDUsedList(maxSubD)			      
        character*16 EqName(maxEq)
        real*8  EqValue(maxEq)
        integer Nequates
        integer CnArea(maxCnArea)
        integer CnCounter   ! CnArea already used. 
        integer AttrCounter ! Volat already used.
        real*8  Volat(maxattr)
!                next 2nd is v.8.10. to keep offset given as =0.45
        real*4  VolatEq(maxattr)
        real*8  VTXx(maxvertex), VTXy(maxvertex), VTXz(maxvertex)
        integer NVTX       ! number of points in VTX*  
!          above are given by calling epqenvlper or epenvlpAll
!         if NVTX = 0, volume is assumed to be contained in a box
!         
        common /Zcnfig/ Det, SubdArea,
     1  SubD,  XYZthick, ORGsubd, PlusF, EqualF, MinusF, 
     2  EqValue, Volat,  VTXx, VTXy, VTXz, VolatEq, 
     3  PosInDet, SubDNumb, SumOffset,
     4  CnArea, comflag, CnCounter, AttrCounter, NVTX,
     5  coment, Nequates, cumsubdloc, SubDUsedList,
     6  Nfield, NsubD,  mode, SubDUsed, wrtcom, comcounter,
     7  comloc, Incused, 
     *  EqualFshort

        common /Zcnfic/ confdata, comarea, Field, FieldAsItis,
     *    SubDName, EqName, form, 
     *   OrgEqBy
