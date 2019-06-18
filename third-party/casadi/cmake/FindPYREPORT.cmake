find_path(PYREPORT_BIN
pyreport
/usr/local/bin/
/usr/bin/
${PYREPORT_BIN_DIR}
)

if(PYREPORT_BIN)
set(FOUND_PYREPORT TRUE)
endif()