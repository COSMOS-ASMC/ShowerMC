const int NpMax  = 50000 ;
struct ob {
  short where, code, subcode, charge;
  float atime, erg, mass, x, y, wx, wy, wz, zenith;
};



struct parent{
  double  posx, posy, posz, coszenith;
  /*
    ! these must be double
    ! otherwise small shift with
    ! children may appear
  */
  float depth, colHeight;
  double  height,  atime;
  short  where, code, asflag;
  float erg;
};




struct child {
  short code, subcode, charge;
  float fm[4], mass;
};





  struct ob o[NpMax];
  struct child c;
  struct parent p;
  double SumehMin, SumegMin,  Ethresh, Cuteg, Cutneg;
  int Np, NoOfLowE, NLowCounter, NhMin, NgMin, HowFlesh;
// next ist integer in Fortran;
  fstream  Wdev, Mdev, Rdev;

  int Accepted, Where;
  logical TopOfNode, RealBegin, RealEnd, Copy;


const int skfl = 120;
char    Mskel[skfl], Wskel[skfl];





