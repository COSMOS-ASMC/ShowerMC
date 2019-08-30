// ***************:: this must be placed in last part of #include
//c             primary class variables
//c
//      record /primaries/ Prim   !to store 
//c                         primary information of various components
//      common /Zprimryv/ Prim
extern struct zprimaryv {
  struct primaries prim; 
} zprimaryv_;



#define Prim  zprimaryv_.prim

















