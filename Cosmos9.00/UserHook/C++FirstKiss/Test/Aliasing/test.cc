extern "C" {
  extern void fort1_(int *);
  extern void fort2_(int *);
  extern void confirm_();
  void subprog();
#include "ccom.h"  
  int main(){
    com_.xxx =13;
    fort1_( &com_.xxx);
    confirm_();
    subprog();
    confirm_();
    com  * x;
    x = &com_;
    fort1_( &x->xxx ); //  this is the same as fort1_( &(x->xxx) );
    confirm_();
    return 0;
  }
  void  subprog(){
    com_.xxx =14;
    fort1_( &com_.xxx);
    confirm_();
  }
}
