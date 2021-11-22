
#include <main_helper.h>
#include "ControllerBridge.h"

int main(int argc, char** argv) {
  main_helper(argc, argv, new ControllerBridge());
  return 0;
}
