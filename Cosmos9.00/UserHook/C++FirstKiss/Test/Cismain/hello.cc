#include <iostream>
using namespace std;
extern "C" /* this line should be omitted in C */
{
  void hello_(long *integ, char *in_str, long *xxx, int in_str_length) {
    // must take care of the string since fortran does not provide a
    // \0 at the end of the string.

    // normal C/C++ programming here
    cout << "Hello\n\n";
    cout << "You passed the following to me:\n";
    cout << "long integer: " << *integ << '\n';
    cout << "string length: " << in_str_length << "\n";
    cout << "xxx:"<< *xxx << "\n";

    char* str=new char[in_str_length+1];        // C++ memory allocation
    for (unsigned int i=0; i<in_str_length; i++)
      str[i]=in_str[i];                 // strncpy is probably nicer here
    str[in_str_length]='\0';

    cout << "string: " << str << '\n';


    delete[] str;                               // C++ memory cleanup
  }
}

/* Need a small wrapper for the fortran main program, since fortran
      code should not be used as the main program, see information
      elsewhere */
extern "C" int mymain_();
int main()
{
  return mymain_();
}
