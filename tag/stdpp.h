#ifndef __PRINT_H_
#define __PRINT_H_

namespace tag {

/// @defgroup stdpp C++ std lib enhancements
/// This group contains enhancements to make the std library more useable

/// @defgroup print stream simplifications
/// The comma operator for streams is defined to allow a simple statement similar to print in Python
/// or other languages. The stream object and any values to be printed are written as a comma separated
/// list. Between the values printed, the stream's fill character is output to delimit values. Some examples:
/// @code
/// cout, 12, "hello", 13, endl;
/// cout, setfill('|');
/// cout, 12, 13, 14, endl;
/// @endcode
/// will result in the following output.
/// @code
/// 12 hello 13
/// 12|13|14
/// @endcode
/// endl is treated specially and no separator is produced before or after it. Other iomanips are not
/// treated specially and produce a fill character before their position (if they are not the first in the list).
/// @ingroup stdpp
/// @{

template <class T> struct NotFirst {
    inline NotFirst(T & d) : data(d) {};
    T & data;
};

template <class T, class O>
inline NotFirst<O> operator,(O  & stream, const T & data ){
    stream << data;
    return NotFirst<O>(stream);
}

template <class T, class O>
inline NotFirst<O> operator,(NotFirst<O> nf, const T & data ){
    nf.data << nf.data.fill() << data;
    return nf;
}

template <class O>
inline O & operator,(O  & stream, O & (*modifier)(O &)){
    stream << modifier;
    return stream;
}

template <class O>
inline O & operator,(NotFirst<O> nf, O & (*modifier)(O &)){
    nf.data << modifier;
    return nf.data;
}

/// @}

} // namespace tag

#endif // __PRINT_H_
