#include <iostream>
#include <stdio.h>
#include <sstream>
extern "C" {
  int main(){
    char x[]="   123.e-3     5";
    //    float z;
    double z;
    int    p;
    std::sscanf(x, "%lf%d", &z, &p);
    cout << z << endl;
    cout << p << endl;
    istringstream iss(x);
    iss >> z;
    cout << z << endl;
    iss >> p;
    cout << p << endl;
  }
}


  


