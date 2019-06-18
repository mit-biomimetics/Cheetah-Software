#!/bin/bash
# This script gets the C or C++ system include paths
# Works with Clang and GCC
# Usage (C): ./get_system_includes.sh clang c > csystem_paths.txt
# Usage (C++): ./get_system_includes.sh gcc-5 c++ > cxxsystem_paths.txt
# Joel Andersson

$1 -E -x $2 - -v < /dev/null >/dev/null 2> /tmp/compiler_settings.txt
list_begin='#include <...> search starts here:'
list_end='End of search list.'
IN_SECTION=0
while read p; do
  if (($IN_SECTION)); then
    if [ "$p" == "$list_end" ]; then
      let IN_SECTION=0
    else
      echo $p
    fi
  else
    if [ "$p" == "$list_begin" ]; then
      let IN_SECTION=1
    fi
  fi
done </tmp/compiler_settings.txt
