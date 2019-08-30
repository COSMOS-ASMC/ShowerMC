#ifndef FNODATDEF 
#define FNODATDEF 0
#endif

#if FNODATDEF > 0
      integer bufc, inbuf,  fnodat
      parameter (fnodat = FNODATDEF)
      type buffer
         integer*2 ldep
         integer*2 code
         integer*2 subcode
         integer*2 charge 
         integer*2 ridx
         integer*2 faiidx
         real  rinmu
         real  fai
         real  Ek
         real  t
         real wx
         real wy
         real wz
  	 real wgt
      end type buffer
#endif
