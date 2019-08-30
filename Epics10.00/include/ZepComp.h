#ifndef ZepComp_
#define ZepComp_
#include "ZepMaxdef.h"
#include "ZepDirec.h"
#include "ZepPos.h"

struct eachc {
  int cn;
  char shape[8];
  char media[8];
  char dirflag;
  struct epPos org;
  struct epPos cent;
  struct epDirec dir[3];
  int nattrib;
  int nthl;  // n-th layer
  int nthc;  // n-th component in the n-th layer
  double abc[MAX_IATTR];
};

#endif
