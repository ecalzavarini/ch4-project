set(XSLDIR ${CMAKE_CURRENT_SOURCE_DIR}/xsl)
set(XMLDIR ${CMAKE_CURRENT_BINARY_DIR}/xml)
set(SRCDIR ${CMAKE_CURRENT_SOURCE_DIR})

set(XMLHEADERS
	${XMLDIR}/parameter-xml.h
	${XMLDIR}/gvars-xml.h
	${XMLDIR}/common-xml.h
)

add_custom_command(
	OUTPUT ${XMLDIR}/parameter-xml.h
	COMMAND /usr/bin/xsltproc --param what 0 ${XSLDIR}/parameter.xsl ${SRCDIR}/parameter.xml | sed -n '/<?xml version/!p' > ${XMLDIR}/parameter-xml.h
	DEPENDS ${SRCDIR}/parameter.xml
)

add_custom_command(
	OUTPUT ${XMLDIR}/gvars-xml.h
	COMMAND /usr/bin/xsltproc --param what 1 ${XSLDIR}/parameter.xsl ${SRCDIR}/parameter.xml | sed -n '/<?xml version/!p' > ${XMLDIR}/gvars-xml.h
	DEPENDS ${SRCDIR}/parameter.xml
)

add_custom_command(
	OUTPUT ${XMLDIR}/common-xml.h
	COMMAND /usr/bin/xsltproc --param what 2 ${XSLDIR}/parameter.xsl ${SRCDIR}/parameter.xml | sed -n '/<?xml version/!p' > ${XMLDIR}/common-xml.h
	DEPENDS ${SRCDIR}/parameter.xml
)

add_custom_command(
	OUTPUT ${XMLDIR}/CMakeLists-xml.txt
	COMMAND /usr/bin/xsltproc --param what 3 ${XSLDIR}/parameter.xsl ${SRCDIR}/parameter.xml | sed -n '/<?xml version/!p' > ${XMLDIR}/CMakeLists-xml.txt
#	COMMAND cd ${CMAKE_BINARY_DIR} && cmake .
	COMMAND cd ${CMAKE_BINARY_DIR} && make rebuild_cache
	DEPENDS ${SRCDIR}/parameter.xml
)

add_custom_target( process_xml
	DEPENDS ${XMLHEADERS} ${XMLDIR}/CMakeLists-xml.txt
	SOURCES ${SRCDIR}/parameter.xml
)

set_source_files_properties(${XMLHEADERS} PROPERTIES GENERATED TRUE HEADER_FILE_ONLY TRUE)

set(XMLDIR ${CMAKE_CURRENT_BINARY_DIR}/xml)
if(NOT EXISTS ${XMLDIR}/CMakeLists-xml.txt)
	message(STATUS "Missing ${XMLDIR}/CMakeLists-xml.txt\n\tcmake will rerun on build")
else()
	message(STATUS "Found ${XMLDIR}/CMakeLists-xml.txt -- I am happy.")
	include(${XMLDIR}/CMakeLists-xml.txt)
endif()

include_directories("${XMLDIR}")

