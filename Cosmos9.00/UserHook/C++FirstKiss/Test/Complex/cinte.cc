#include <iostream.h>
#include <math.h>
#include "Complex.h"
complex f(complex z);
double  r(double teta);
complex fxy(double teta);
main(){
  double pi;
  pi = asin(1.0)*2;  //compute pi
  cout <<"pi=" << pi << "\n"; // show it
  double twopi = 2*pi;  

  complex sum =complex(0., 0.); // clear sum
  double teta=0.0; 
  complex z1 = fxy(teta); // start point z
  int imax =  10000; // divide 0-2pi by imax
  double  dteta = twopi/imax;  // upto twopi-dteta
  for(int i=0; i< imax; ++i){
    complex z2 =   fxy(teta+dteta);
    complex zeta = fxy(teta+dteta/2); // middle point
    sum +=  f(zeta) * (z2 - z1);  // integral
    z1 = z2;
    teta += dteta; 	
  }
  cout << sum << "\n";  
}
double r(double teta){
  return (2+ sin(teta*2 + cos(teta)));
}
complex fxy(double teta){
  double rt=r(teta);
  return  complex(rt*cos(teta), rt*sin(teta));
}
complex f(complex z){
  return z/(z-.5);
}

