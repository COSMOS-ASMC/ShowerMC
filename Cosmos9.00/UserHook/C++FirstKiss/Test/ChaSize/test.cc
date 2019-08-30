
#include <iostream>
#include <string.h>
struct abc{
  double   x;
  char ch10[10];
  int  ii;
  char  ch8[8];
  char  ch4[4];
  char  dummy;
  //  int kk;
};

struct coord {
 double r[3];
 char sys[4]; 
 char ca[5][1024];  
  char dummy;
};



extern "C" {
  extern struct zz {
    struct coord mycoord;
    struct abc   mix;
  } zz_;
};

extern "C" void fort_();
main(){
    fort_();
    std::cout << "  x " <<  std::endl;
    std::cout << zz_.mix.x << std::endl;
    std::cout << "  ch10" <<  std::endl;
    std::cout << zz_.mix.ch10 <<  std::endl;
    std::cout << "  ii" <<  std::endl;
    std::cout << zz_.mix.ii << std::endl;
    std::cout << "  ch8" <<  std::endl;
    std::cout << zz_.mix.ch8 << std::endl;
    std::cout << "  ch4" <<  std::endl;
    std::cout << zz_.mix.ch4 <<  std::endl;

    std::cout << "  r1 " <<  std::endl;
    std::cout << zz_.mycoord.r[0]<< std::endl;
    std::cout << "  sys " <<  std::endl;
    std::cout << zz_.mycoord.sys<< std::endl; 
    std::cout << "  mycord " <<  std::endl;
    std::cout << zz_.mycoord.ca[1] << std::endl; 
}
