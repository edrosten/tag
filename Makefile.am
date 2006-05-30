# the library we want to create and its source files
lib_LTLIBRARIES     = libtag.la
libtag_la_SOURCES   = src/absorient.cpp src/fourpointpose.cpp

# the header files to be installed
pkginclude_HEADERS  = tag/absorient.h tag/fourpointpose.h tag/kalmanfilter.h tag/ransac.h tag/ransac_estimators.h  tag/stdpp.h tag/printf.h tag/tuple.h

docs:
	doxygen Doxyfile
