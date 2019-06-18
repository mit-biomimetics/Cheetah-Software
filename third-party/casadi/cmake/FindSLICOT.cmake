find_library(SLICOT_LIB
  slicot slicot64_pic slicot_pic
  PATHS $ENV{SLICOT_LIBRARY_DIR})

if(SLICOT_LIB)
  set(SLICOT_LIBRARIES ${SLICOT_LIB})
  set(SLICOT_FOUND TRUE)
else()
  message(STATUS "SLICOT: Cound not find library slicot. Try stetting SLICOT_LIBRARY_DIR env var.")
endif()
