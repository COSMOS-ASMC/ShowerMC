#ifndef Z90hist1_
#define Z90hist1_
/*
c       fortran 90 version  1D histogram
      type histogram1
         type(histoc) c
         type(histogram) x
         real, allocatable ::  xw(:)
         real, allocatable ::  dnw(:)
         real, allocatable ::  dndx(:)
         real, allocatable ::  mean(:)
      end type
*/
struct histogram1 {
  struct histoc c;
  struct histogram x;
  float *xw;
  float *dnw;
  float *dndx;
  float *mean;
  float *xforI;
  float *integ;
};

#endif
