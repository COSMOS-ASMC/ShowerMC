#include <iostream>
#define  Abc_  zexp_.abc
#define  Dd   zexp_.dd
#define  integer  int
extern "C" {
  struct zexp {
    double  dd; 
    integer abc;
  } zexp_;
  extern void xyz_();
}


using namespace std;
int  main(){
  xyz_();
  cout << zexp_.abc << endl;
  cout << zexp_.dd <<  endl;
  Abc_  = 100;
  Dd = -1.023;
  cout << zexp_.abc << endl;
  cout << Dd << endl;
  cout << zexp_.dd << endl;
}
