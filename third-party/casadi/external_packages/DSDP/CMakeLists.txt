cmake_minimum_required(VERSION 2.8.6)

find_package(LAPACK REQUIRED)

include_directories(include src/bounds src/lp src/sdp src/solver src/sys src/vecmat)

file(GLOB DSDP_SRC "src/bounds/*.c"  "src/bounds/*.h" "src/lp/*.c" "src/lp/*.h" "src/sdp/*.c" "src/sdp/*.h" "src/solver/*.c" "src/solver/*.h" "src/sys/*.c" "src/sys/*.h" "src/vecmat/*.c" "src/vecmat/*.h")

set(DSDP_LIBRARIES ${DSDP_LIBRARIES})

add_library(casadi_dsdp STATIC ${DSDP_SRC})

target_link_libraries(casadi_dsdp ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
