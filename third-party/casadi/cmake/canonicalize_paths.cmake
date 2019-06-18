# The following code canonicalizes paths:
# e.g.  converts "foo/bar/..;foo/../foo"   to "foo;foo"
# When a library is present in "foo" directory, "foo/bar/..;foo/../foo" would make CMake
#  complain about this library being present in two different folders,
#  which would be a recipe for trouble where it not for
#  the fact that these directories are really identical.
#  It could be considered a workaround for a bug in CMake.
function(canonicalize_paths MYPATH)
set(MYPATH_NEW "")

foreach(arg ${${MYPATH}})
  # get_filename_component(CANONICAL_PATH ${arg} REALPATH)
  #  REALPATH is not defined in cmake 2.6.0
  #
  # The following is a workaround:
  set(CANONICAL_PATH "${arg}")
  foreach(i RANGE 10)
    string(REGEX REPLACE "/([^/]*)/\\.\\./" "/" CANONICAL_PATH "${CANONICAL_PATH}")
  endforeach()
  string(REGEX REPLACE "/([^/]*)/\\.\\.$" "" CANONICAL_PATH "${CANONICAL_PATH}")
  foreach(i RANGE 10)
    string(REGEX REPLACE "\\\\([^\\\\]*)\\\\\\.\\.\\\\" "\\\\" CANONICAL_PATH "${CANONICAL_PATH}")
  endforeach()
  string(REGEX REPLACE "\\\\([^\\\\]*)\\\\\\.\\.$" ""  CANONICAL_PATH "${CANONICAL_PATH}")
  set(MYPATH_NEW ${MYPATH_NEW} "${CANONICAL_PATH}")
endforeach()
set(${MYPATH} ${MYPATH_NEW} PARENT_SCOPE)
endfunction()
