cmake_minimum_required(VERSION 2.8.6)

set(TINYXML_SRC
  tinystr.cpp
  tinyxml.cpp
  tinyxmlerror.cpp
  tinyxmlparser.cpp
)

add_library(casadi_tinyxml STATIC ${TINYXML_SRC})
