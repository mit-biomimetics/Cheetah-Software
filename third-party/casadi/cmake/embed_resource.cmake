file(STRINGS "${INPUT}" FILE_CONTENTS)
file(READ "${LICENSE_FILE}" COPYRIGHT_HEADER)

file(WRITE ${OUTPUT}.hpp ${COPYRIGHT_HEADER})
file(APPEND ${OUTPUT}.hpp "#ifndef RESOURCE_${SYMBOL}_HPP\n#define RESOURCE_${SYMBOL}_HPP\nnamespace casadi {\nextern const char * resource_${SYMBOL};" )

file(APPEND ${OUTPUT}.hpp "\n} // namespace casadi\n#endif // RESOURCE_${SYMBOL}_HPP\n" )

file(WRITE ${OUTPUT}.cpp ${COPYRIGHT_HEADER})
file(APPEND ${OUTPUT}.cpp "#include \"resource_${SYMBOL}.hpp\"\n namespace casadi {\nconst char * resource_${SYMBOL} =" )

foreach(line ${FILE_CONTENTS})
  string(REGEX REPLACE "\\\\" "\\\\\\\\" line "${line}")
  string(REGEX REPLACE "\"" "\\\\\"" line "${line}")
  file(APPEND ${OUTPUT}.cpp "\n  \"${line}\\n\"" )
endforeach()

file(APPEND ${OUTPUT}.cpp ";\n} // namespace casadi\n" )
