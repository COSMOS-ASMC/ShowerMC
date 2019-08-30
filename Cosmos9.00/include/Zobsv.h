// *************this may better come after Zobsp.h (see NoOfSites2)
/*
c           need Zcoord.h  Zobs.h  Zpos.h Zmagfield.h
	integer NoOfSites           ! No of particle observation sites
        integer NoOfASSites
c          
	real*8 CosLatitude          ! oos of Latitude of deepest obs. site
        real*8 SinLatitude          ! sin
        real*8 CosLongitude         ! cos of Longitude
        real*8 SinLongitude         ! sin of ..

	record /coord/ DetZaxis     ! detector's Z axis in 'xyz' system
        record /coord/ DetXaxis     !  //        X    // 

        record /coord/ Xprimary     ! primary system x axis in 'xyz'
        record /coord/ Yprimary     ! primary system y axis in 'xyz'
        record /coord/ Zprimary     ! primary system z axis in 'xyz'
                                    ! these are computed in cprimxyz in
                                    ! ciniTracking in ceventLoop

        record /coord/ PolarInjPos  ! polar angle of the injection point in xyz.

	record /magfield/ MagfieldNED     ! mag in 'ned' at deepest obs. site
        record /magfield/ MagfieldHVA     ! mag in 'hva' at //. both in T.
        record /magfield/ MagfieldXYZ     ! mag in 'xyz' at //. both in T.

         structure /site/
               record /position/pos
               real*8  zpl           ! z value in 1ry system
	       real*8  mu
               real*8  minitime
         end structure 
         structure /assite/
               record /position/pos
               real*8  zpl 
               real*8  mu             ! Moliere Unit
               real*8  esize          ! electron size
               real*8  age            ! size weighted age
         end structure

	 record /site/ ObsSites(0:maxNoOfSites+1)
         record /assite/ ASObsSites(maxNoOfASSites)
c            to store Ne, age of a component shower for an electron
         real*8 CompASNe(maxNoOfASSites), CompASAge(maxNoOfASSites)
*/

struct site {
  struct position pos;
  double   zpl;
  double   mu;
  double   minitime;
};

struct assite {
  struct position pos;
  double zpl; 
  double mu;
  double esize;
  double age;
};


extern struct zobsvc {
  struct assite  asobssites[maxnoofassites];
  struct site    obssites[maxnoofsites+2];
  struct magfield magfieldned;
  struct magfield magfieldhva;
  struct magfield magfieldxyz;
  double compasne[maxnoofassites];
  double compasage[maxnoofassites];
  struct coord detzaxis;
  struct coord detxaxis;
  struct coord xprimary;
  struct coord yprimary;
  struct coord zprimary;
  struct coord polarinjpos;
  double coslatitude;
  double sinlatitude;
  double coslongitude;
  double sinlongitude;
  int    noofsites;
  int    noofassites;
} zobsvc_;



#define ASObsSites   zobsvc_.asobssites
#define ObsSites     zobsvc_.obssites
#define MagfieldNED  zobsvc_.magfieldned
#define MagfieldHVA  zobsvc_.magfieldhva
#define MagfieldXYZ  zobsvc_.magfieldxyz
#define CompASNe     zobsvc_.compasne
#define CompASAge    zobsvc_.compasage
#define DetZaxis     zobsvc_.detzaxis
#define DetXaxis     zobsvc_.detxaxis
#define Xprimary     zobsvc_.xprimary
#define Yprimary     zobsvc_.yprimary
#define Zprimary     zobsvc_.zprimary
#define PolarInjPos  zobsvc_.polarinjpos
#define CosLatitude  zobsvc_.coslatitude	 
#define SinLatitude  zobsvc_.sinlatitude
#define CosLongitude zobsvc_.coslongitude
#define SinLongitude zobsvc_.sinlongitude
#define NoOfSites    zobsvc_.noofsites
#define NoOfASSites  zobsvc_.noofassites


      
      

        
        
