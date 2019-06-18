#include "Utilities/Utilities_print.h"

void printf_color(PrintColor color, const char *fmt, ...){
  auto color_id = (uint32_t)color;
  if(color_id)
    printf("\033[1;%dm", (uint32_t)color + 30);
  va_list args;
  va_start(args, fmt);
  vprintf(fmt, args);
  va_end(args);
  printf("\033[0m");
}