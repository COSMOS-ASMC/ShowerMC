#include <iostream>
void sub() {
  static  int x = 0;
  x ++;
std::cout << x << std::endl;
}

int main(){
  sub();
  sub();
  sub();
}

    
