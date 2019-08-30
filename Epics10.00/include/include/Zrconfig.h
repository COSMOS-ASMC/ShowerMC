#ifndef Zrconfig_
#define Zrconfig_
#include "ZepMaxdef.h"
#include "ZepDirec.h"
#include "ZepDirec.h"
#include "ZepPos.h"
#include "ZepComp.h"

struct rconfig {
  int ncomp;
  int NoOfSciFiLayers;
  int NoOfSciFi;

  int NoOfBGOLayers;
  int NoOfBGO;


  int NoOfSCIN;

  int NoOfSi;


  double scifidE[MAX_NoOfSciFiLayers][MAX_NoOfSciFiPerLayer];
  double scificoord[MAX_NoOfSciFiPerLayer];
  double scifiz[MAX_NoOfSciFiLayers];
  int NoOfSciFiPerLayer[MAX_NoOfSciFiLayers];
  char scifiXorY[MAX_NoOfSciFiLayers];

  double BGOdE[MAX_NoOfBGOLayers][MAX_NoOfBGOPerLayer];
  double BGOcoord[MAX_NoOfBGOPerLayer];
  double BGOz[MAX_NoOfBGOLayers];
  int NoOfBGOPerLayer[MAX_NoOfBGOLayers];
  char BGOXorY[MAX_NoOfBGOLayers];


  double SCINdE[MAX_NoOfSCIN];
  struct epPos SCINcent[MAX_NoOfSCIN];

  struct eachc comp[MAX_COMPONENT];
};


#endif /*  Zrconfig */
