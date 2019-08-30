#include <iostream>
extern "C" {
  extern void fort1_(int *);
  extern void fort2_(int *, int *);
  extern void fort3_(int *);
  extern void dummy_();
}
#include "ccom.h"  
using namespace std;
int main(){
    int temp[3] = {1,2,3};
    fort1_( &temp[1] );
    com_.xxx =13;
    com_.instabc[0].i8=10;
    com_.instabc[2].i8=5;
    fort1_( &com_.xxx);
    fort1_( &com_.instabc[0].i8);
    int ppp=8;
    fort2_( &com_.instabc[2].i8,  & ppp);
    dummy_();
    int ar[3] = {-2,-3,-4};
    fort3_(&ar[0]);
    cout << ar[0] << endl;
    cout << ar[1] << endl;
    cout << ar[2] << endl;
    return 0;
}




