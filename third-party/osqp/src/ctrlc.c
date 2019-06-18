/*
 * Implements signal handling (ctrl-c) for OSQP.
 *
 * Under Windows, we use SetConsoleCtrlHandler.
 * Under Unix systems, we use sigaction.
 * For Mex files, we use utSetInterruptEnabled/utIsInterruptPending.
 *
 */

#include "ctrlc.h"

#if defined MATLAB

static int istate;

void osqp_start_interrupt_listener(void) {
  istate = utSetInterruptEnabled(1);
}

void osqp_end_interrupt_listener(void) {
  utSetInterruptEnabled(istate);
}

int osqp_is_interrupted(void) {
  return utIsInterruptPending();
}

#elif defined IS_WINDOWS

static int int_detected;
static BOOL WINAPI handle_ctrlc(DWORD dwCtrlType) {
  if (dwCtrlType != CTRL_C_EVENT) return FALSE;

  int_detected = 1;
  return TRUE;
}

void osqp_start_interrupt_listener(void) {
  int_detected = 0;
  SetConsoleCtrlHandler(handle_ctrlc, TRUE);
}

void osqp_end_interrupt_listener(void) {
  SetConsoleCtrlHandler(handle_ctrlc, FALSE);
}

int osqp_is_interrupted(void) {
  return int_detected;
}

#else /* Unix */

# include <signal.h>
static int int_detected;
struct sigaction oact;
static void handle_ctrlc(int dummy) {
  int_detected = dummy ? dummy : -1;
}

void osqp_start_interrupt_listener(void) {
  struct sigaction act;

  int_detected = 0;
  act.sa_flags = 0;
  sigemptyset(&act.sa_mask);
  act.sa_handler = handle_ctrlc;
  sigaction(SIGINT, &act, &oact);
}

void osqp_end_interrupt_listener(void) {
  struct sigaction act;

  sigaction(SIGINT, &oact, &act);
}

int osqp_is_interrupted(void) {
  return int_detected;
}

#endif /* END IF IS_MATLAB / WINDOWS */
