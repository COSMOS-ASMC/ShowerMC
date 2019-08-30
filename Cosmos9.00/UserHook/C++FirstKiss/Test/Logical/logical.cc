
#include <iostream>
typedef   int logical;
extern "C" {
  extern  void ddd_();

  struct abc {
    double dv;
    double dv2;
    logical    torf;
    logical  torf2;
    int    ii;
    char  dummy;
  } ;
  extern struct yyy {
     struct abc test1[2];
  } yyy_;
}
int main(){
    ddd_();
    for(int i =1; i<=2; i++) {
      std::cout << yyy_.test1[i-1].dv << std::endl;
      std::cout << yyy_.test1[i-1].dv2 << std::endl;
      std::cout << yyy_.test1[i-1].torf << std::endl;
      std::cout << yyy_.test1[i-1].torf2 << std::endl;
      std::cout << yyy_.test1[i-1].ii << std::endl;
  }
}




