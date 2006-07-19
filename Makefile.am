# the library we want to create and its source files
lib_LTLIBRARIES     = libtoontag.la
libtoontag_la_SOURCES   = src/absorient.cpp src/fourpointpose.cpp

# the header files to be installed
pkginclude_HEADERS  = \
	tag/absorient.h \
	tag/array.h \
	tag/fourpointpose.h \
	tag/kalmanfilter.h \
	tag/constantposition.h \
	tag/constantvelocity.h \
	tag/measurements.h \
	tag/ransac.h \
	tag/ransac_estimators.h  \
	tag/stdpp.h \
	tag/printf.h \
	tag/tuple.h \
	tag/fn.h \
	tag/intersection.h

docs:
	doxygen Doxyfile
