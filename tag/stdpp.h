#ifndef __PRINT_H_
#define __PRINT_H_

#include <iostream>

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


#ifndef DOXYGEN_IGNORE_INTERNAL
namespace Internal
{

	struct add_fill_s{};
	struct like_print_s{};
	struct no_space_s{};

	template<class S> struct add_fill_bound
	{
		add_fill_bound(S& os)
		:o(os),first(1) {}


		template<class C> add_fill_bound& operator<<(const C& c)
		{
				if(first == true)
					first=false;
				else
					o << o.fill();

				o << c;
				return *this;
		}
		
		add_fill_bound& operator<<(const no_space_s&)
		{
			first=true;
			return *this;
		}
		
		add_fill_bound& operator<<(S& (*fptr)(S&) )
		{
				o << fptr;

				if(fptr == static_cast<S&(*)(S&)>(std::endl))
					first=true;

				return *this;
		}

		private:
			S& o;
			bool first;
	};


	template<class S> struct like_print_bound:public add_fill_bound<S>
	{
		like_print_bound(S&os)
		:add_fill_bound<S>(os)
		{
		}
		
		~like_print_bound()
		{
			add_fill_bound<S>::operator<<(std::endl);
		}
	};

}

#endif 

/// Ostream modifier which puts spaces between elements using the fill character. endl is treated
/// specially, so a space is not put after an endl
/// @code
///	cout << add_fill << 1 << 2 << 3 << endl << 4  << 5 << endl;
/// @endcode
/// This will print:
/// @code
/// 1 2 3
/// 5 4
/// @endcode
/// @ingroup print
static struct Internal::add_fill_s add_fill;

/// Ostream modifier similar to @ref add_fill, except that an endl is added at the end automatically
/// @code
/// cout << print << 1 << 2 << 3;
/// cout << "hello" << endl;
/// @endcode
/// Will print:
/// @code
/// 1 2 3
/// hello
/// @endcode
/// @ingroup print
static struct Internal::like_print_s print;


/// An ostream modifier to use with @ref print and @ref add_fill which prevents a space being added
/// before the next element.
/// @code
/// cout << @ref add_fill  << 0 << 1 << no_space << 2 << 3 << endl;
/// @endcode
/// will print:
/// @code
/// 0 12 3
/// @endcode
/// @ingroup print
static struct Internal::no_space_s no_space;


#ifndef DOXYGEN_IGNORE_INTERNAL

template<class Char, class Traits> Internal::add_fill_bound<std::basic_ostream<Char,Traits> > operator<<(std::basic_ostream<Char,Traits>& o, const Internal::add_fill_s&)
{
		return Internal::add_fill_bound<std::basic_ostream<Char,Traits> >(o);
}
 
template<class Char, class Traits> Internal::like_print_bound<std::basic_ostream<Char,Traits> > operator<<(std::basic_ostream<Char,Traits>& o, const Internal::like_print_s&)
{
		return Internal::like_print_bound<std::basic_ostream<Char,Traits> >(o);
}

#endif

} // namespace tag

#endif // __PRINT_H_
