#ifndef Z90histo_
#define Z90histo_
/*
c       fortran 90 version
      type histogram
         character*4 label
         character*8 unit
         real  inc
         real bin
         real xmin
         real xm
         real sumw 
         integer*2 nhist
         integer*2 cent
         integer*2 imin
         integer*2 imax
         integer*2 step
         logical*1 tklg
         logical*1 ufl
         logical*1 ofl
      end type
*/
struct  histogram {
  char label[4];
  char unit[8];
  float inc;
  float bin;
  float xmin;
  float xm;
  float sumw;
  short int nhist;
  short int cent;
  short int imin;
  short int imax;
  short int step;
  short int  tklg;
  short int  ufl;
  short int  ofl;
};

#endif
