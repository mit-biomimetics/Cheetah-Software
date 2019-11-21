/*!
 * @file SegfaultHandler.h
 * @brief Handler for segfaults.
 * This will catch a segfault, print a stack trace, and put an error code in the shared memory
 * (if it is connected), so that the simulator can provide a reasonable error message for why the
 * robot code disappears.
 */

#ifndef PROJECT_SEGFAULTHANDLER_H
#define PROJECT_SEGFAULTHANDLER_H

#include <cstdint>

void install_segfault_handler(char* error_message);

#endif //PROJECT_SEGFAULTHANDLER_H
