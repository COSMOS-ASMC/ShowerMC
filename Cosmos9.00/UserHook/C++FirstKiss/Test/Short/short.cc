#include <iostream>
#include <math.h>
using namespace std;
void main(){
  unsigned short int  x;
  x = (int) (pow(2,16)-1);
  cout << x << "\n";
  cout << sizeof(x)<<  "\n";
}

