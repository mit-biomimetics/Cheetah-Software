# This file was partly taken from Sarcasm/irony-mode and modified

if (OLD_LLVM)
# Legacy version, to be removed
message(STATUS "Looking for clang 3.4.2")

set(LLVM_INSTALL_PREFIX $ENV{CLANG})
set(CLANG_LLVM_LIB_DIR ${LLVM_INSTALL_PREFIX}/lib)

set(LLVM_DEP LLVMJIT LLVMInterpreter LLVMX86CodeGen LLVMBitWriter LLVMIRReader LLVMipo LLVMLinker LLVMRuntimeDyld LLVMExecutionEngine LLVMAsmPrinter LLVMSelectionDAG LLVMX86Desc LLVMAsmParser LLVMBitReader LLVMVectorize LLVMMCParser LLVMCodeGen LLVMX86AsmPrinter LLVMX86Info LLVMObjCARCOpts LLVMScalarOpts LLVMX86Utils LLVMInstCombine LLVMTransformUtils LLVMipa LLVMAnalysis LLVMTarget LLVMCore LLVMMC LLVMObject LLVMSupport LLVMBitWriter LLVMIRReader LLVMipo LLVMLinker LLVMVectorize LLVMInstrumentation LLVMBitReader LLVMOption LLVMX86CodeGen LLVMAsmParser LLVMAArch64AsmParser LLVMAArch64Disassembler LLVMARMCodeGen LLVMARMAsmParser LLVMARMDisassembler LLVMCppBackendCodeGen LLVMHexagonCodeGen LLVMMipsCodeGen LLVMMipsAsmParser LLVMMipsDisassembler LLVMMSP430CodeGen LLVMNVPTXCodeGen LLVMPowerPCCodeGen LLVMPowerPCAsmParser LLVMR600CodeGen LLVMSparcCodeGen LLVMSystemZCodeGen LLVMSystemZAsmParser LLVMSystemZDisassembler LLVMX86AsmParser LLVMX86Disassembler LLVMX86Desc LLVMX86AsmPrinter LLVMX86Utils LLVMX86Info LLVMXCoreCodeGen LLVMXCoreDisassembler LLVMAArch64CodeGen LLVMAsmPrinter LLVMMCParser LLVMSelectionDAG LLVMCodeGen LLVMObjCARCOpts LLVMScalarOpts LLVMInstCombine LLVMTransformUtils LLVMipa LLVMAnalysis LLVMARMDesc LLVMCppBackendInfo LLVMHexagonAsmPrinter LLVMMipsDesc LLVMMSP430Desc LLVMNVPTXDesc LLVMPowerPCDesc LLVMR600Desc LLVMSparcDesc LLVMSystemZDesc LLVMXCoreDesc LLVMAArch64Desc LLVMARMAsmPrinter LLVMARMInfo LLVMHexagonDesc LLVMMipsAsmPrinter LLVMMipsInfo LLVMMSP430AsmPrinter LLVMMSP430Info LLVMNVPTXAsmPrinter LLVMNVPTXInfo LLVMPowerPCAsmPrinter LLVMPowerPCInfo LLVMR600AsmPrinter LLVMR600Info LLVMSparcInfo LLVMSystemZAsmPrinter LLVMSystemZInfo LLVMXCoreAsmPrinter LLVMXCoreInfo LLVMAArch64AsmPrinter LLVMAArch64Info LLVMTarget LLVMHexagonInfo LLVMAArch64Utils LLVMCore LLVMMC LLVMObject LLVMSupport)

set(CLANG_DEFINITIONS "-DCLANG_ENABLE_OBJC_REWRITER" "-DGTEST_HAS_RTTI=0" "-D_GNU_SOURCE" "-D__STDC_CONSTANT_MACROS" "-D__STDC_FORMAT_MACROS" "-D__STDC_LIMIT_MACROS")
set(CLANG_CXX_FLAGS "-fPIC -fvisibility-inlines-hidden -ffunction-sections -fdata-sections -fno-common -fno-strict-aliasing")

# Clang shares include directory with llvm
set(CLANG_INCLUDE_DIR ${LLVM_INSTALL_PREFIX}/include)

# All clang libraries
set(CLANG_DEP clangFrontend clangDriver clangCodeGen clangRewriteFrontend clangSerialization clangParse clangSema clangAnalysis clangEdit clangAST clangLex clangBasic ${LLVM_DEP})

# Get libraries
set(CLANG_LIBRARIES)
foreach(D ${CLANG_DEP})
  find_library(CLANG_DEP_${D} ${D} ${LLVM_INSTALL_PREFIX}/lib)
  set(CLANG_LIBRARIES ${CLANG_LIBRARIES} ${CLANG_DEP_${D}})
endforeach()

else (OLD_LLVM)

# Return variables
set(CLANG_LIBRARIES)
set(CLANG_DEFINITIONS)
set(CLANG_CXX_FLAGS)

# most recent versions come first: http://llvm.org/apt/
set(CLANG_KNOWN_LLVM_VERSIONS 4.0.1 4.0.0
  3.7.0 3.7
  3.6.2 3.6.1 3.6.0 3.6
  3.5.2 3.5.1 3.5.0 3.5
  3.4.2 3.4.1 3.4
  3.3
  3.2
  3.1)

# List of likely locations of llvm
set(llvm_hints)
foreach (version ${CLANG_KNOWN_LLVM_VERSIONS})
  string(REPLACE "." "" undotted_version "${version}")
  if(APPLE)
    # LLVM Homebrew
    set(llvm_hints ${llvm_hints} "/usr/local/Cellar/llvm/${version}/bin")
    # LLVM MacPorts
    set(llvm_hints ${llvm_hints} "/opt/local/libexec/llvm-${version}/bin")
  elseif(UNIX)
    # FreeBSD ports versions
    set(llvm_hints ${llvm_hints} "/usr/local/llvm${undotted_version}/bin")
  endif()
endforeach()

# Locate the LLVM config script
find_program(CLANG_LLVM_CONFIG llvm-config HINTS $ENV{CLANG}/bin ${llvm_hints})

# LLVM version
execute_process(COMMAND ${CLANG_LLVM_CONFIG} --version
                OUTPUT_VARIABLE CLANG_LLVM_VERSION
                OUTPUT_STRIP_TRAILING_WHITESPACE)
message(STATUS "Found LLVM ${CLANG_LLVM_VERSION}")

# LLVM installation prefix
execute_process(COMMAND ${CLANG_LLVM_CONFIG} --prefix
                OUTPUT_VARIABLE CLANG_LLVM_PREFIX
                OUTPUT_STRIP_TRAILING_WHITESPACE)

# LLVM include directory
execute_process(COMMAND ${CLANG_LLVM_CONFIG} --includedir
                OUTPUT_VARIABLE CLANG_LLVM_INCLUDE_DIR
                OUTPUT_STRIP_TRAILING_WHITESPACE)

# LLVM library directory
execute_process(COMMAND ${CLANG_LLVM_CONFIG} --libdir
                OUTPUT_VARIABLE CLANG_LLVM_LIB_DIR
                OUTPUT_STRIP_TRAILING_WHITESPACE)

# Clang shares include directory with llvm
set(CLANG_INCLUDE_DIR ${CLANG_LLVM_INCLUDE_DIR})

# All clang libraries
set(CLANG_DEP clangFrontend clangDriver clangCodeGen clangRewriteFrontend clangSerialization
              clangParse clangSema clangAnalysis clangEdit clangAST clangLex clangBasic)

# Get libraries
foreach(D ${CLANG_DEP})
  find_library(CLANG_DEP_${D} ${D} HINTS ${CLANG_LLVM_LIB_DIR})
  set(CLANG_LIBRARIES ${CLANG_LIBRARIES} ${CLANG_DEP_${D}})
endforeach()

# Get the LLVM libraries
execute_process(COMMAND ${CLANG_LLVM_CONFIG} --libfiles
                OUTPUT_VARIABLE CLANG_LLVM_LIBRARIES
                OUTPUT_STRIP_TRAILING_WHITESPACE)
separate_arguments(CLANG_LLVM_LIBRARIES)
set(CLANG_LIBRARIES ${CLANG_LIBRARIES} ${CLANG_LLVM_LIBRARIES})

# Get system libs
execute_process(COMMAND ${CLANG_LLVM_CONFIG} --system-libs
                OUTPUT_VARIABLE CLANG_LLVM_SYSTEM_LIBS
                OUTPUT_STRIP_TRAILING_WHITESPACE)
separate_arguments(CLANG_LLVM_SYSTEM_LIBS)
foreach(D ${CLANG_LLVM_SYSTEM_LIBS})
  if(${D} STREQUAL "-lm")
    # Standard math
    set(CLANG_LIBRARIES ${CLANG_LIBRARIES} m)
  elseif(${D} STREQUAL "-lz")
    # LLVM needs to be linked with zlib
    find_package(ZLIB QUIET REQUIRED)
    set(CLANG_LIBRARIES ${CLANG_LIBRARIES} ${ZLIB_LIBRARIES})
  elseif(${D} STREQUAL "-lcurses")
    if(APPLE)
      # ??
      find_library(CLANG_CURSES ncurses)
    else()
      # ??
      find_library(CLANG_CURSES curses)
    endif()
    set(CLANG_LIBRARIES ${CLANG_LIBRARIES} ${CLANG_CURSES})
  elseif(${D} STREQUAL "-lpthread")
    find_package(Threads QUIET REQUIRED)
    set(CLANG_LIBRARIES ${CLANG_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
  else()
    message(WARNING "Ignoring LLVM dependency '${D}'")
  endif()
endforeach()

# C++ compiler flags
message(STATUS "CLANG_LLVM_CONFIG: ${CLANG_LLVM_CONFIG}")
execute_process(COMMAND ${CLANG_LLVM_CONFIG} --cxxflags
                OUTPUT_VARIABLE CLANG_LLVM_CXXFLAGS
                OUTPUT_STRIP_TRAILING_WHITESPACE)
separate_arguments(CLANG_LLVM_CXXFLAGS)
foreach(D ${CLANG_LLVM_CXXFLAGS})
  if(${D} MATCHES "-D")
    set(CLANG_DEFINITIONS ${CLANG_DEFINITIONS} ${D})
  else()
    message(STATUS "Ignoring LLVM C++ flag '${D}'")
  endif()
endforeach()

endif (OLD_LLVM)

if(MINGW)
  message("/usr/${PREFIX}/lib")
  find_library(IMAGEHLP imagehlp /usr/${PREFIX}/lib)
  if (IMAGEHLP)
    set(CLANG_LIBRARIES ${CLANG_LIBRARIES} ${IMAGEHLP})
  endif()
  find_library(SHLWAPI libshlwapi /usr/${PREFIX}/lib)
  if (SHLWAPI)
    set(CLANG_LIBRARIES ${CLANG_LIBRARIES} ${SHLWAPI})
  endif()
  message("${CLANG_LIBRARIES}")
endif()


set(CLANG_FOUND TRUE)
