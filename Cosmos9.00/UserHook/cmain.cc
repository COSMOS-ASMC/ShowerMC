extern "C" {
  extern void cmymain_(); // fortran subroutine special for c++ interrace
  int main(){
    cmymain_();
  }
}

