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

On a unix system, <code>./configure && make install </code> will compile the %tag
library and install it to the correct place. It creates dynamic and static libraries
with the library name @b toontag.

For Win32 systems, the @c build directory contains project files for different versions
of Visual Studio.

@section sLinks Links
    - TooN - http://mi.eng.cam.ac.uk/~twd20/TooN/html/index.html
    - CVS web interface - http://cvs.savannah.gnu.org/viewvc/tag/?root=toon

*/
