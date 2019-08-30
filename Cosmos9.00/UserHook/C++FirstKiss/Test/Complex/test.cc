extern "C" {
  extern void fort1_(int *);
  extern void fort2_(int *);
  extern void dummy_();
#include <math.h>
#include <complex.h>
#include "ccom.h"  
  void main(){
    com_.xxx =13;
    com_.instabc[0].i8=10;
    fort1_( &com_.xxx);
    fort1_( &com_.instabc[0].i8);
    dummy_();
  }
}
