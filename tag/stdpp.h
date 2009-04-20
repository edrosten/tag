#ifndef TAG_STDPP_H
#define TAG_STDPP_H

#include <iostream>
#include <utility>

namespace tag {

/// @defgroup stdpp C++ std lib enhancements
/// This group contains enhancements to make the std library more useable

/**
@defgroup printgroup stream simplifications
This group contains two sets of enhancements for using the output streams in the standard library.

The first is a set of additional modifiers for use with the standard << operator. See @ref add_fill ,
@ref no_space and @ref print for more details.

The second set introduces a new syntax using the comma operator for streams.
It creates simple statements similar to print in Python or other languages.
The stream object and any values to be printed are written as a comma separated
list. Between the values printed, the stream's fill character is output to delimit values. Some examples:
@code
cout, 12, \"hello\", 13, endl;
cout, setfill('|');
cout, 12, 13, 14, endl;
@endcode
will result in the following output
@code
12 hello 13
12|13|14
@endcode
endl is treated specially and no separator is produced before or after it. Other iomanips are not
treated specially and produce a fill character before their position (if they are not the first in the list).
@ingroup stdpp
*/

#ifndef DOXYGEN_IGNORE_INTERNAL

static struct noendl_s {} noendl;

template <class T> struct NotFirst {
    inline NotFirst(T & d) : data(d), last(true) {}
    inline ~NotFirst() { if(last) data << std::endl; }
    T & data;
    bool last;
    
    template <class S>
    inline NotFirst<T> & operator,( const S & arg ){
        data << data.fill() << arg;
        return *this;
    }

    inline T & operator,(T & (*modifier)(T &)){
        data << modifier;
        last = false;
        return data;
    }

    inline NotFirst<T> & operator,( noendl_s & arg ){
        last = false;
        return *this;
    }
};

template <class T, class Char, class Traits> 
inline NotFirst<std::basic_ostream<Char,Traits> > operator,(std::basic_ostream<Char,Traits>  & stream, const T & arg ){
    stream << arg;
    return NotFirst<std::basic_ostream<Char,Traits> >(stream);
}

template <class Char, class Traits>
inline std::basic_ostream<Char,Traits> & operator,(std::basic_ostream<Char,Traits>  & stream, std::basic_ostream<Char,Traits> & (*modifier)(std::basic_ostream<Char,Traits> &)){
    stream << modifier;
    return stream;
}

template <class Char, class Traits>
inline std::basic_ostream<Char,Traits> & operator,(std::basic_ostream<Char,Traits>  & stream, noendl_s & arg ){
    return stream;
}

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

/**
Ostream modifier which puts spaces between elements using the fill character. endl is treated
specially, so a space is not put after an endl
@code
cout << add_fill << 1 << 2 << 3 << endl << 4<< 5 << endl;
@endcode
This will print:
@code
1 2 3
5 4
@endcode
@ingroup printgroup
*/
static struct Internal::add_fill_s add_fill;

/**
Ostream modifier similar to @ref add_fill, except that an endl is added at the end automatically
@code
cout << print << 1 << 2 << 3;
cout << "hello" << endl;
@endcode
Will print:
@code
1 2 3
hello
@endcode
@ingroup printgroup
*/
static struct Internal::like_print_s print;

/**
An ostream modifier to use with @ref print and @ref add_fill which prevents a space being added
before the next element.
@code
cout << @ref add_fill  << 0 << 1 << no_space << 2 << 3 << endl;
@endcode
will print:
@code
0 12 3
@endcode
@ingroup printgroup
*/
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

#ifndef DOXYGEN_IGNORE_INTERNAL
namespace Internal
{
	template<class A, class B> struct refpair
	{
		A& a;
		B& b;
		refpair(A& aa, B& bb)
		:a(aa),b(bb)
		{}

		void operator=(const std::pair<A,B>& p)
		{
			a=p.first;
			b=p.second;
		}
	};
}

#endif

/**
Similar to <code>std::make_pair</code>, but for references. This can be used
for multiple return values:
@code
float f;
int i;
rpair(f,i) = make_pair(2.2f, 1);
@endcode
@param aa first value
@param bb second value
**/
template<class A, class B> Internal::refpair<A,B> rpair(A&aa, B&bb)
{
	return Internal::refpair<A,B>(aa, bb);
}



} // namespace tag

#endif // __PRINT_H_
