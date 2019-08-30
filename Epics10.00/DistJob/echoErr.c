#include <stdio.h>
/*    echoErr something; something will appear on errout  */
main(argc, argv)
     short argc;
     char *argv[];
{
  short nl = 1;
  if(strcmp(argv[1], "-n") == 0) {
    nl = 0;
    argc--;
    argv++;
  }
  while (--argc) {
    fputs(*++argv,stderr);
    if(argc > 1)
      putc(' ', stderr);
  }
  if (nl)
    putc('\n', stderr);
}
