#include <iostream>
using namespace std;
const double  i=1000;
const double  j=1000e3;
const int  k  = (int) (j/i);
double  yyy;
double *xxx;

void main(){
  char  zz[k];
  zz[0] = 'a';
  zz[1]= '\0';
  cout<< zz << "\n" ;
  xxx  = &yyy; 
  *xxx = 10.0;
  cout << *xxx << endl;
}

