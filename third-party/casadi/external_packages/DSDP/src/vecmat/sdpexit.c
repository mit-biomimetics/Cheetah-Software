#include "numchol.h"
#include <string.h>

void ShutDown(){

  /*   sdpdat* sdt = sdat; */
  printf("\n Shutdown --  ");

} /* ShutDown */


int ExitProc(int  ccode,
             char *str)
{
  xcode code=(xcode)ccode;

  printf("\n Exit -- %d : ",ccode);
  
  switch (code) {
    case OptFound:
      printf("optimal solution found");
      return code;
    case OutOfSpc:
      printf("out of memory space");
      break;
    default:
      break;
  }
  if (str) printf(", %s",str);

  ShutDown();

  printf("\n Exiting --  ");

  return 1;
} /* ExitProc */

