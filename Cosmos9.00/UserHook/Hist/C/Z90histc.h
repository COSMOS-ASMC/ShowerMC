#ifndef Z90histc_
#define Z90histc_
/*
      type histoc
         real norm
         real pw
         integer*2 eventno
         logical*1 logv
         character*8  init
         character*128 title
         character*8  categ
         character*96 id 
         character*32 dNunit
         character*128 dir
      end type
*/
struct histoc {
  float norm;
  float pw;
  short int eventno;
  //  unsigned char  logv;
  int logv;
  char  init[8];
  char  title[128];
  char  categ[8];
  char  id[96];
  char  dNunit[32];
  char  dir[128];
  char  dummy;
};

#endif
