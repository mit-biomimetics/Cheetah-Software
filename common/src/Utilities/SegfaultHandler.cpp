/*!
 * @file SegfaultHandler.cpp
 * @brief Handler for segfaults.
 * This will catch a segfault, print a stack trace, and put an error code in the shared memory
 * (if it is connected), so that the simulator can provide a reasonable error message for why the
 * robot code disappears.
 */


#include <cstdio>
#include <execinfo.h>
#include <csignal>
#include <cstdlib>
#include <unistd.h>
#include <cstring>

#include "Utilities/SegfaultHandler.h"
#include "Utilities/Utilities_print.h"



static char* error_message_buffer;

/*!
 * Called on segfault.  Prints stack trace, flushes output, and sends error code to sim
 * @param sig : signal (should be SIGSEGV = 11)
 */
static void segfault_handler(int sig) {
  void* stack_frames[200];
  int size = backtrace(stack_frames, 200);
  fprintf_color(PrintColor::Red, stderr, "CRASH: Caught %d (%s)\n",
      sig, strsignal(sig));
  backtrace_symbols_fd(stack_frames, size, STDERR_FILENO);

  fflush(stderr);
  fflush(stdout);

  if(error_message_buffer)
    strcpy(error_message_buffer, "Segfault!\nCheck the robot controller output for more information.");
  exit(1);
}

static void sigint_handler(int sig) {
  if(sig == SIGINT && error_message_buffer) {
    strcpy(error_message_buffer, "Robot program has been killed by SIGINT");
  }
  exit(0);
}

/*!
 * Install the segfault handler function
 */
void install_segfault_handler(char* error_message) {
  signal(SIGSEGV, segfault_handler);
  signal(SIGINT, sigint_handler);
  error_message_buffer = error_message;
}