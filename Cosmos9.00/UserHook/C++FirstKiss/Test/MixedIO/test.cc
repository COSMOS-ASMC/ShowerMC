#include <iostream>
#include <stdio.h>
using namespace std;
extern "C" {
  extern void iofortran_();
}
int main(){
    ios::sync_with_stdio();

    iofortran_();
    printf("this is test io\n");
    int x;
    cout << "eneter the value of x\n";
    cin >> x;
    cout << " the value of x is " <<  x << "\n" ;
}


  
