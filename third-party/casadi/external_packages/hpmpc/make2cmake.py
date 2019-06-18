import subprocess
import re
import os

targets = set()
for L in file("src/Makefile","r"):
  m = re.match("ifeq \(\$\(TARGET\), ([A-Z0-9_]+)\)", L)
  if m:
    targets.add(m.group(1))

targets = sorted(list(targets))

combinations = [
("X64_AVX","GENERIC"),
("X64_AVX2","X64_INTEL_HASWELL"),
("X64_AVX","X64_INTEL_SANDY_BRIDGE"),
("X64_SSE3","X64_INTEL_CORE"),
("X64_SSE3","X64_AMD_BULLDOZER"),
("C99_4X4","GENERIC")]
#("C99_4X4_PREFETCH","GENERIC")]

#("CORTEX_A7","GENERIC"),
#("CORTEX_A9","GENERIC"),
#("CORTEX_A15","GENERIC"),
#("CORTEX_A57","GENERIC"),
#("C99_4X4","GENERIC"),
#("C99_4X4_PREFETCH","GENERIC")
#]

target_data = {}

this_pwd = os.getcwd()
for t in targets:
  p = subprocess.Popen(["make","CC=echo compiler_action `pwd`","TARGET=%s" % t],stdout=subprocess.PIPE,stderr=subprocess.PIPE,cwd="src")

  stdout,stderr = p.communicate()

  data = []

  for L in stdout.split("\n"):
    if L.startswith("compiler_action"):
      line = L.split(" ")
      pwd = line[1][len(this_pwd)+1:]
      cfile = line[-1]
      flags = [f for f in line[2:-4] if not f.startswith("-I") and not f.startswith("-DOS")]
      data.append((cfile,pwd,flags))

  for cfile,pwd,f in data:
    assert f==data[0][2] or len(f)==0
  target_data[t] = data

with file('CMakeLists.txt','w') as f:
  f.write("cmake_minimum_required(VERSION 2.8.6)\n")
  f.write("""
    if(WIN32)
    set(HPMPC_FLAGS -DOS_WINDOWS)
    elseif(APPLE)
     set(HPMPC_FLAGS -DOS_MAC)
    else()
     set(HPMPC_FLAGS -DOS_LINUX)
    endif()
    if()
    set(CMAKE_SHARED_LINKER_FLAGS ${CMAKE_SHARED_LINKER_FLAGS} -Wl,--unresolved-symbols=report-all)
    endif()
    """)
  f.write("include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../blasfeo/src/include)\n")
  f.write("include_directories(${CMAKE_CURRENT_SOURCE_DIR}/dummy_include)\n")
  f.write("include_directories(${CMAKE_CURRENT_SOURCE_DIR}/dummy_include/foo/bar)\n")
  #f.write("include_directories(/opt/blasfeo/include)\n")

  for t,blasfeo_target in combinations:
    data = target_data[t]

    libname = "casadi_hpmpc_%s_%s" % (t,blasfeo_target)
    blasfeo_libname = "casadi_blasfeo_%s" % blasfeo_target
    bitness64 = "-m64" in data[0][2]
    if bitness64:
      f.write("if( CMAKE_SIZEOF_VOID_P EQUAL 8 )\n")
    f.write("add_library(%s SHARED %s)\n" % (libname, "\n".join([ os.path.join(pwd,cfile) for cfile, pwd, _ in data])))
    f.write("SET_TARGET_PROPERTIES(%s PROPERTIES COMPILE_FLAGS \"-DLA_BLASFEO -DTARGET_%s -DTARGET_%s %s ${HPMPC_FLAGS}\")\n\n" % (libname, blasfeo_target,t , " ".join(data[0][2])))
    f.write("target_link_libraries(%s %s m)\n" % (libname, blasfeo_libname))
    if bitness64:
      f.write("endif()\n")
