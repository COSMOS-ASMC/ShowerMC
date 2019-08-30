extern "C" {
  struct xxxx {
    float z;
    int ar[5];
    int j;
    int tf; 
    char string[2][10];
  } ;
  struct yyyy {
    struct xxxx a1, a2;
    float a3;
  };

  extern struct {
    struct yyyy ppp;
  } abc_;
  
  int cleng(char *);
  void callfort_();
#include <iostream>
  int main() {
    callfort_();
  }
  void xyz_( yyyy * a,  xxxx * zz, int * id) {
    cout << " a3="<< a->a3 << "  a1.ar[0]=" <<  a->a1.ar[0] << "\n";
    cout << " z="<< zz->z << "\n";
    cout << abc_.ppp.a3 << "  " << abc_.ppp.a1.ar[0] << "\n";
    cout <<  abc_.ppp.a2.z << "\n";
    cout << *id <<" " << *id << "\n";
  }

}

      
