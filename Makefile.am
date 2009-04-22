# the library we want to create and its source files

lib_LTLIBRARIES     = libtoontag.la
libtoontag_la_SOURCES   = src/quartic.cpp src/threepointpose.cpp src/absorient.cpp \
                          src/handeye.cpp src/five_point.cpp src/five_point_buildmatrix.cpp src/plrt.c

# the header files to be installed
pkginclude_HEADERS  = \
	tag/quartic.h \
	tag/threepointpose.h \
	tag/absorient.h \
	tag/array.h \
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
	tag/intersection.h \
	tag/helpers.h \
	tag/unscented.h \
	tag/handeye.h \
	tag/five_point.h

docs:
	doxygen Doxyfile
