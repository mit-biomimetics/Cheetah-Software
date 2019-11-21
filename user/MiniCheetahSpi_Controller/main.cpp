
#include <main_helper.h>
#include "MiniCheetahSpi_Controller.h"

int main(int argc, char** argv) {
  main_helper(argc, argv, new MiniCheetahSpi_Controller());
  return 0;
}