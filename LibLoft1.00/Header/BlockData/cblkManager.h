      block data  cblkManager

#include  "Zmanagerp.h"
#ifdef HP9000
       data 
     *  Errorout /7/
#else
       data 
     *  Errorout /0/   ! sun4 ; others.  if 0 is ng, give by namelist
#endif
       data  
     * DestEventNo /1,0/ ,
     * Hidden /.false./ ,
     * TempDev /11/ ,
     * PrimaryFile /' '/ ,
     * Job /' '/ ,
     * ContFile /'ContInfo'/ , 
     * SkeletonFile /'SkeletonParam'/ ,
     * SeedFile /'Seed'/ ,
     * PrevEventNo /0/ ,
     * EventNo /0/, 
     * SeedFileDev /22/ ,
     * DeadLine /' '/ ,
     * Within /99999/ , 
     * BaseTime /10./ ,
     * BaseErg /1000./ ,
     * BasePower /1./
       data
     * CutOffFile /' '/ 
       data
     * GeomagFile /' '/
       data
     * UserHookc /MAX_USERHOOKC*' '/
       data
     * UserHooki /MAX_USERHOOKI*99999999/
       data
     * UserHookr /MAX_USERHOOKR*1.d65/
       data
     * DpmFile /' '/
       data
     * PercentEnv /' '/
       data
     * SharpEnv /' '/
       data
     * AtEnv /' '/
       data
     * NRL_period /4*0/
       data
     * AtmosModel /1/
       data
     * AtmosFile /' '/,  ObjFile /' '/
      end
