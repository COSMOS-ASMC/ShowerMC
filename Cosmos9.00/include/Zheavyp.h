/*
c
c	(->	------------------------------------------
*/
extern struct zheavyc {
  double ptavnonintenuc;
  double ptavfrag;
  double fragmenttbl[maxheavyg][maxheavyg];
  int    charge2heavyg[maxheavycharge];
  int    heavyg2massn[maxheavyg];
  int    heavyg2charge[maxheavyg];
  int    heavyg2code[maxheavyg];
  int    code2massn[khvymax];
  int    code2heavyg[khvymax];
  int    howintnuc;
} zheavyc_;

#define  PtAvNonInteNuc  zheavyc_.ptavnonintenuc
#define  PtAvFrag        zheavyc_.ptavfrag
#define  FragmentTbl     zheavyc_.fragmenttbl
#define  Charge2heavyG   zheavyc_.charge2heavyg
#define  HeavyG2massN    zheavyc_.heavyg2massn
#define  HeavyG2charge   zheavyc_.heavyg2charge
#define  HeavyG2code     zheavyc_.heavyg2code
#define  Code2massN      zheavyc_.code2massn
#define  Code2heavyG     zheavyc_.code2heavyg
#define  HowIntNuc       zheavyc_.howintnuc

/*
  integer Charge2heavyG	!2  charge of heavy $\rightarrow$  heavy group index conversion array.
  integer HeavyG2massN    !2  heavy group index $\rightarrow$     mass number conversion array.
  integer HeavyG2charge	!2  heavy group index $\rightarrow$     charge of heavy conversion array.
  integer HeavyG2code     !2  heavy group index $\rightarrow$     particle code conversion array.
  integer Code2massN      !2  particle code $\rightarrow$     mass number conversion array.
  integer Code2heavyG	!2  particle code $\rightarrow$     heavy group index conversion array.
  real*8  FragmentTbl	!2  tbl(i,j)=$<$Number$>$  of frag. j when a heavy of heavy group index i
  !    breaks up at air.
  real*8  PtAvNonInteNuc  !2  $<$Pt$>$  of non interacting nucleons.
  real*8  PtAvFrag        !2  $<$Pt$>$  of heavy fragments.
  character*4 HeavyG2symbol !2   heavy group index $\rightarrow$  'Fe' etc conversion array.
  integer HowIntNuc       !2 If 0, the  number of interacting nucleons among a projectile heavy nucleus is 
  !  determined as the number of first collision of each interacting nucleon inside 
  ! the  nucleus.  If 1, the number is determined as the total number of collisions 
  !   including successive interactions. Default is 1. (There is uncertaninity in
  !  interpretation of the formula; value 1 gives larger number of interacting
  !  nucleons.)
  
  c	<-)	--------------------------------------
*/        

extern struct zheavycc {
  char heavyg2symbol[maxheavyg][5];
  DUMMYCHAR
} zheavycc_;
   //        common /Zheavycc/ HeavyG2symbol(maxHeavyG)

#define HeavyG2symbol  zheavycc_.heavyg2symbol


	
