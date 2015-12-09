#include <stdio.h>
#include <string.h>

extern char *strncat();

void chmod_(const char *cfx)
{
  int gl;
  char doit[] = "chmod 00755 ";

  gl=3;
  strncat(doit,cfx,gl);  
  system(doit); 
}
