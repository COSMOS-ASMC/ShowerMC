/*
  c         
  c               variables generated during incident angle sampling

  complex*16 Obsvhour      !  observation hour for point source
  real*8 Cspsmx, Cspsmn    ! cos(zenith) of point source
  record /coord/ AngleAtObsCopy, DcAtObsXyz
   c        
   c           AngleAtObsCopy  :  direction cos at the deepest Obssite
   c                              in 'det' system
   c           DcAtObsXyz :       its transfroamtion to "xyz" system
   c
	  record /track/ IncidentCopy
          common /Zincidentv/ IncidentCopy, AngleAtObsCopy,
     *    DcAtObsXyz,  Obsvhour, Cspsm,x Cspsmn
*/


extern struct zincidentv {
  struct track incidentcopy;
  struct coord angleatobscopy; 
  struct coord dcatobsxyz;
  double obsvhour[2];  // size of  complex*16 
  double cspsmx;
  double cspsmn; 
} zincidentv_;

#define  IncidentCopy    zincidentv_.incidentcopy
#define  AngleAtObsCopy  zincidentv_.angleatobscopy
#define  DcAtObsXyz      zincidentv_.dcatobsxyz
#define  Obsvhour        zincidentv_.obsvhour
#define  Cspsmx          zincidentv_.cspsmx
#define  Cspsmn          zincidentv_.cspsmn

