#ifndef Zfuncdef_
#define Zfuncdef_
#include "ZepPos.h"
#include "ZepDirec.h"
#include "Zrconfig.h"
#include "ZsimHeader.h"
#include <stdio.h>

struct epPos getCenter(struct eachc *comp); 
struct epPos epR2r( struct epDirec X, struct epDirec Y, struct epDirec Z, struct epPos R) ;
struct epPos epr2R( struct epDirec X, struct epDirec Y, struct epDirec Z, struct epPos r) ;
int readReformedConfig(FILE *fd, struct eachc data[]);
int prepareAnalysis(FILE *fd,  struct rconfig *calet);
int inspectRconfig(struct rconfig *data);
int getHeader(FILE *fd, struct simHeader *hd);
int getCompVSdE(FILE *fd, struct rconfig *calet);

#endif
