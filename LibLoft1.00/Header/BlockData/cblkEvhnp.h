      block data cblkEvhnp

#include  "Zevhnp.h"
!            currently usable models
      data RegMdls/'phits', 'jam', 'dpmjet3', 'fritiof7.02',
     *     'fritiof1.6', 'epos','sibyll',
     *     'gheisha', 'nucrin', 'ad-hoc', 'incdpm3', 'qgsjet2',
     *     'special'/ 

      
       data  
     * Cepic0 /0.00d0/ , 
     * Cekaon /0.d0/ ,
     * Ceneuc /0.35d0/ ,
     * DoNPadjust /0/ ,
     * dpmRareDecay /1/ ,
     * Mudirp /1.d0/ , 
     * K0sSemiLD /12/ ,
     * Kpilog /0.0062/ ,
     * Kpicns /0.077/ ,
#ifdef NEXT486
     * Elund /4.99/ , 
     * Elund2 /4.99/ ,
     * Elund3 /4.99/ ,
#else
     * Elund /500./ , 
     * Elund2 /500./ ,
     * Elund3 /500./ ,
#endif
     * SucPw /1.5/,
     * MulLow /1/,      ! old one is 0. bef. v.5.0
     * LundPara /1, 2, 0, 0, 1, 5*0/,
     * IntModel /'"dpmjet3" 1.d8'/ ,
     * TotXSopt /2/,
     * SxAbySxpOpt /1/,
     * XsecModel /' '/,
     * InclusiveFile /'$LIBLOFT/Had/Import/DPM/Inclusive/inputdata'/ ,
     * Efermi /10.e-3/
       data 
     * Eta2Pi0/0.2/,
#ifdef IBMAIX
     * SucInt /0/
#else
     * SucInt /1/     ! old one was 1 < v5.0; 0 <=v5.10   1 >= v5.11
#endif
      end
