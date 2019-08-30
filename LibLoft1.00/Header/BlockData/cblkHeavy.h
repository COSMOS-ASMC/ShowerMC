!     ********************
      block data cblkHeavy
!     ********************

#include  "Zcode.h"
#include  "Zheavyp.h"
      integer i, j

      data
     *  HeavyG2massN /1,   4,    8,    14,    25,    35,     56/ ,
     *  HeavyG2charge /1,   2,    4,     7,    12,    17,     26/ ,
     *  HeavyG2symbol /'p', 'Alfa',  'L', 'CNO', 'H', 'VH', 'Fe'/ ,
     *  Charge2heavyG /1,   2,   3,3,  4,4,4,  5*5,   5*6,   9*7/

      data 
     *  HeavyG2code /knuc, kalfa, klibe, kcno, khvy, kvhvy, kiron/    
      data
     * Code2heavyG(knuc)/1/ ,
     * Code2heavyG(kalfa)/2/ ,
     * Code2heavyG(klibe)/3/ ,
     * Code2heavyG(kcno)/4/ , 
     * Code2heavyG(khvy)/5/ ,
     * Code2heavyG(kvhvy)/6/ ,
     * Code2heavyG(kiron)/7/
      data
     * Code2massN(knuc)/1/ ,
     * Code2massN(kalfa)/4/ , 
     * Code2massN(klibe)/8/ ,
     * Code2massN(kcno) /14/,
     * Code2massN(khvy)/25/,
     * Code2massN(kvhvy)/35/,
     * Code2massN(kiron)/56/

      data  HowIntNuc/1/
!
      data ( ( FragmentTbl(i,j), j=1,maxHeavyG), i=1,maxHeavyG)/
     1             1.,    0.,    0.,    0.,    0.,    0.,   0.,
     2             4.,    0.,    0.,    0.,    0.,    0.,   0.,
     3             5.81,  .61,   .11,   0.,    0.,    0.,   0.,
     4             6.62,  .72,   .24,   .17,   0.,    0.,   0.,
     5             10.5,  .77,   .21,   .39,   .16,   0.,   0.,
     6             12.8,  1.17,  .17,   .20,   .42,   .06,  0.,
     7             18.1,  1.71,  .24,   .17,   .22,   .20,  .17/

         data PtAvNonInteNuc/50.e-3/   ! 50 MeV
         data PtAVFrag/50.e-3/     ! 50 MeV
      end
