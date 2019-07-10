/*!
 * Segfault handler to print stack trace on crash.
 */

#include "include/Utilities/SegfaultHandler.h"
#include <cstdio>
#include <execinfo.h>
#include <csignal>
#include <cstdlib>
#include <unistd.h>
#include <cstring>
#include <include/Utilities/Utilities_print.h>


static void segfault_handler(int sig) {
  void* stack_frames[200];
  int size = backtrace(stack_frames, 200);
  fprintf_color(PrintColor::Red, stderr, "CRASH: Caught %d (%s)\n",
      sig, strsignal(sig));
  backtrace_symbols_fd(stack_frames, size, STDERR_FILENO);
  exit(1);
}


void install_segfault_handler() {
  signal(SIGSEGV, segfault_handler);
}