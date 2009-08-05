// Tag main documentation file

/**
@mainpage tag - TooN Algorithm library

@section sIntro Introduction

tag is a companion library to TooN, a C++ header file collection for numerics computation.
It provides a set of algorithms typically found in Computer Vision applications:

    - @ref kalmanfiltergroup 
    - @ref posegroup
    - @ref absorient
    - @ref handeye
    - @ref essentialgroup
    - @ref ransac
    - @ref unscented

Moreover it also provides general C++ standard library related helper classes and functions, dealing with:
    
    - @ref stdpp
    - @ref printf
    - @ref tuple

Why use this library ?
    - Based on TooN for efficient numerics calculations based on LAPACK and BLAS
    - Implementations tested in research prototypes and usually reasonable correct and stable

@section sDownload Getting the code and installing

To get the code from cvs use:

@code
cvs -z3 -d:pserver:anoncvs@cvs.savannah.nongnu.org:/cvsroot/toon co tag
@endcode

@subsection unix Unix

On a unix system, <code>./configure && make install </code> will compile the %tag
library and install it to the correct place. It creates dynamic and static libraries
with the library name @b toontag.

@subsection win32 Windows

For Win32 systems, the @c build directory contains project files for different versions
of Visual Studio. Currently the vc2008 solutin is supported and should work out of the box. There
are two projects, one for compiling tag and one for installing it into a common directory tree.
Both projects assume the existence of three environment variables describing the location of header, 
library and binary files (for DLLs). 
	
	- @c INCLUDEDIR contains the header files. tag headers will be copied into @c \%INCLUDEDIR%\\\tag
	- @c LIBDIR contains library files. tag static libraries (debug and release verions) will be copied into @c \%LIBDIR\%
	- @c BINDIR is not used for tag, but would be the default directory for DLLs

@section sLinks Links
    - TooN - http://mi.eng.cam.ac.uk/~twd20/TooN/html/index.html
    - CVS web interface - http://cvs.savannah.gnu.org/viewvc/tag/?root=toon

*/
