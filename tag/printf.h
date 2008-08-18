#ifndef TAG_PRINTF_H
#define TAG_PRINTF_H

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cctype>

#include <tag/tuple.h>

namespace tag
{
	#ifndef DOXYGEN_IGNORE_INTERNAL
	namespace Internal
	{

		// Code for parsing and interpreting printf style format specifiers.
		//Code for groking format strings
		struct format
		{
			enum
			{
				ALT = 1, ZP=2, LEFT=4, SPACE=8, SIGN=16, PERCENT=32, BAD=64,
				NO_PRECISION=-1, NO_WIDTH=-1
			};

			int flags;
			int  width, precision;
			char conversion;

			inline int parse(const std::string& fmt, int pos)
			{
				flags=0;
				width=NO_WIDTH;
				precision=NO_PRECISION;
				conversion=0;

				int s = fmt.size();

				if(pos == s)
				{
					flags |= BAD;
					return 0;
				}

				//Check for literal
				if(fmt[pos] == '%')
				{
					flags = PERCENT;
					return pos+1;
				}

				//Get flags (if any)
				for(; pos < s; pos++)
				{
					char c = fmt[pos];

					if(c == '#')
						flags |= ALT;
					else if(c == '0')
						flags |= ZP;
					else if(c == ' ')
						flags |= SPACE;
					else if(c == '-')
						flags |= LEFT;
					else if(c == '+')
						flags |= SIGN;
					else
						break;
				}

				//Get width
				for(; pos < s; pos++)
				{
					char c = fmt[pos];
					if(isdigit(c))
						if(width == NO_WIDTH)
							width = c - '0';
						else
							width = width * 10 + (c-'0');
					else
						break;
				}

				//Check for a precison
				if(pos < s && fmt[pos] == '.')
				{
					precision = 0;
					pos++;

					//Get width
					for(; pos < s; pos++)
					{
						char c = fmt[pos];
						if(isdigit(c))
							precision = precision * 10 + (c-'0');
						else
							break;
					}
				}

				//Now, this should be the conversion
				if(pos < s && isalpha(fmt[pos]))
					conversion = fmt[pos++];
				else
					flags = BAD;

				return pos;
			}
		};

		//To make it looks like ostream << format works
		template<class Char, class Traits> struct bound_format
		{
			bound_format(std::basic_ostream<Char, Traits>& os, const format& ff)
			:o(os), f(ff)
			{}

			std::basic_ostream<Char, Traits>& o;
			const format& f;
		};


		template<class Char, class Traits> bound_format<Char, Traits> operator<<(std::basic_ostream<Char, Traits>& o, const format& f)
		{
			return bound_format<Char, Traits>(o, f);
		}

		//Evaluate the result of osteram << format << X
		template<class Char, class Traits, class C> std::basic_ostream<Char, Traits>& operator<<(bound_format<Char,Traits> f, const C& c)
		{
			using std::ios;
			using namespace std;

			bool precision_is_max_width=0;


			//Save old stream state.
			int old_p = f.o.precision();
			Char old_fill = f.o.fill();
			ios::fmtflags old_flags = f.o.flags();

			//Conversion specific tricks
			//Defaults
			f.o.unsetf(ios::floatfield | ios::boolalpha);
			f.o.setf(ios::dec);
			f.o.fill(' ');

			//Process conversion characters. These can affect the formatting
			//parameters below.
			switch(f.f.conversion)
			{
				case 'f':
				case 'F':
					f.o.setf(ios::fixed);
					break;

				case 'e':
				case 'E':
					f.o.setf(ios::scientific);
					break;

				case 'g':
				case 'G':
					f.o.unsetf(ios::floatfield);
					break;

				case 'x':
				case 'X':
					f.o << hex;
					break;

				case 'o':
				case 'O':
					f.o << oct;
					break;

				case 'b':
				case 'B':
					f.o << boolalpha;
					break;

				case 's':
				case 'S':
					precision_is_max_width=1;
					break;

				case 'k':
				case 'K':
					return f.o;
					break;
			}



			if(f.f.width != format::NO_WIDTH)
				f.o.width(f.f.width);

			if(f.f.flags & format::ZP && !(f.f.flags & format::LEFT))
				f.o.fill('0');

			if(f.f.flags & format::SIGN)
				f.o.setf(ios::showpos);
			else
				f.o.unsetf(ios::showpos);

			if(f.f.flags & format::LEFT)
				f.o.setf(ios::left);
			else
				f.o.setf(ios::internal);

			if(f.f.flags * format::ALT)
				f.o.setf(ios::showbase | ios::showpoint);
			else
				f.o.unsetf(ios::showbase | ios::showpoint);


			if(isupper(f.f.conversion))
				f.o.setf(ios::uppercase);
			else
				f.o.unsetf(ios::uppercase);


			if(f.f.precision != format::NO_PRECISION && ! precision_is_max_width)
				f.o.precision(f.f.precision);

			if(precision_is_max_width && f.f.precision != format::NO_PRECISION)
			{
				ostringstream tmp;
				tmp.copyfmt(f.o);

				//Since we're doing the truncation by hand, then there should be
				//no width specification
				tmp << setw(0) << c;

				f.o << tmp.str().substr(0, f.f.precision);
			}
			else
				f.o << c;


			//Reset values
			f.o.precision(old_p);
			f.o.fill(old_fill);
			f.o.setf(old_flags);

			return f.o;
		}

		//This parses a format tring and uses it to print a typelist. The typelist has to be
		//accessed in reverse order, since the last element added to the typelist is the
		//one at the start of the list. Since lists are produced with (Fmt, a, b, c), that means
		//that c would be accessed first. To make it look like printf, a needs to be accessed first.
		template<class C, class D, int i, int max> struct print_typelist
		{
			static void print(std::ostream& o, const std::string& fmt, int fpos, const T_list<C,D>& l)
			{
				size_t ppos;

				while(1)
				{
					ppos = fmt.find('%', fpos);

					if(ppos == fmt.npos)
					{
						//No more format specifiers, so output the rest of the string
						o <<  &fmt[fpos];
						return;
					}

					//else output the strung up to the specifier
					o  << fmt.substr(fpos, ppos - fpos);

					//Parse the format string
					format f;
					int pos = f.parse(fmt, ppos+1);

					if(f.flags & format::PERCENT)
					{
						o << '%';
						fpos = pos;
						continue;
					}
					else if(f.flags & format::BAD)
					{
						o  << "<Malformed format>" << &fmt[ppos];
						return ;
					}
					else
					{
						o << f << l.template index<i>();
						print_typelist<C,D, i+1, max>::print(o, fmt, pos, l);
						return;
					}
				}
			}

		};

		template<class C, class D, int max> struct print_typelist<C, D, max, max>
		{
			static void print(std::ostream& o, const std::string& fmt, int fpos, const T_list<C,D>&)
			{
				size_t ppos;

				while(1)
				{
					ppos = fmt.find('%', fpos);

					if(ppos == fmt.npos)
					{
						//No more format specifiers, so output the rest of the string
						o <<  &fmt[fpos];
						return;
					}

					//else output the strung up to the specifier
					o  << fmt.substr(fpos, ppos - fpos);

					//Parse the format string
					format f;
					int pos = f.parse(fmt, ppos+1);

					if(f.flags & format::PERCENT)
					{
						o << '%';
						fpos = pos;
						continue;
					}
					else if(f.flags & format::BAD)
					{
						o  << "<Malformed format>" << &fmt[ppos];
						return ;
					}
					else
					{
						o << "<Missing value>";
						fpos = pos;
						continue;
					}
				}
			}
		};
	}
	#endif


	///@defgroup printf Typesafe, variadic printf
	///This group provides a reimplementation of the C standard library
	///printf family of functions in a typesafe manner, using the @ref tag::T_list tuple type.
	///@ingroup stdpp

	/**
	 This function provides the mose generic implementation, where a typelist and
	 an ostream can be provided.


	 @param o ostream to print to
	 @param fmt format string
	 @param l typelist of arguments

	 The format specifiers work much like the standard printf ones, but with some
	 notable differences. The format is:

	 @code
	 %[flags][width][.[precision]]<conversion>
	 @endcode

	 <h5>Flags</h5>

	 Flags contains zero or more of:
	 <table>
	 <tr><td>#</td>     <td>   Alternate form: Include trailing zeros for float, use 0x or 0 prefix for
	                     hex and octal ints </td></tr>
	 <tr><td>0</td>     <td>   Zero pad everything, including (oddly) strings. </td></tr>
	 <tr><td>' '</td>   <td> Space. Does nothing, I haven't figured out how to do this one.  </td></tr>
	 <tr><td>-</td>     <td>   Left justify. Overrides 0 </td></tr>
	 <tr><td>+</td>     <td>   Always show sign for numbers </td></tr>
	 </table>

	 <h5>Width</h5>
	 The width is optional and must start with a non-zero digit. It seems to work
	 even with overloaded operator<<

	 <h5>Precision</h5>
	  @code
	  .[precision]
	  @endcode

	  This can mean one of two things:
	  - The precision for floating point numbers (affects all numbers in the following operator<<)
	  - The maximum width for "string" types. Truncation happens from the right.

	  <h5>Conversion</h5>
	  This is a single alphabetic character.  This is the largets departure
	  from the C style printf. This flag affects the "type" of the datum being
	  formatted. Type safety still applies. Uppercase and lowercase conversions
	  behave in exactly the same way, except uppercase ones cause uppercased
	  output (where applicable). If a float looks like 1E6 intead of 1e6
	  <table>
	  <tr><td>x</td>	<td>Make ints come out in hexadecimal</td></tr>
	  <tr><td>o</td>	<td>Make ints come out in octal</td></tr>
	  <tr><td>f</td>    <td>Make floats come out as fixed</td></tr>
	  <tr><td>g</td>    <td>Make floats come out in the general format</td></tr>
	  <tr><td>e</td>    <td>Make floats come out in exponential (scientific) format</td></tr>
	  <tr><td>b</td>     <td>Make bools coma out as "true" or "false"</td></tr>
	  <tr><td>s</td>     <td>Treat next output as string: precision means max width.</td></tr>
	  <tr><td>k</td>     <td>Kill (ignore) the data</td></tr>
	  </table>


	  <h5>Differences campared to C-stype printf</h5>

	  Unlike C, the C++ sytle printf can output values of any type for which
	  operator<< is defined. Further, printf itself will not cause any undefined
	  behaviour. Access beyond the end of the vargs list produces <tt>\<Missing value\></tt>
	  in the output, instead of a segfault. A malformed format string causes
	  <tt>\<Malformed format\></tt> to appear in the output followed by the remaining string
	  including the malformed part. No further conversions are attempted after
	  a malformed format.

	  Note: length modifiers are not supported.

	  See also other members of the Printf family:
	  - @ref vPrintf
	  - @ref vsPrintf
	  - @ref Printf
	  - @ref fPrintf
	  - @ref sPrintf
	  And list making facilities
	  - @ref Fmt
	  - @ref T_list

	  <h5>Usage</h5>
	  @code
	  vfPrintf(cout, "%s, %s\n", (Fmt, "Hello", "world"));
	  @endcode
	  @ingroup printf
	**/
	template<class A, class B, class C, class D> void vfPrintf(std::basic_ostream<A, B>& o, const std::string fmt, const T_list<C,D>& l)
	{
		Internal::print_typelist<C, D, 0, T_list<C,D>::elements >::print(o, fmt, 0, l);
	}

	///This prints to cout. See @ref vfPrintf for details.
	///@ingroup printf
	///@param fmt format string
	///@param l typelist of arguments
	template<class C, class D> void vPrintf(const std::string fmt, const T_list<C,D>& l)
	{
		vfPrintf(std::cout, fmt, l);
	}

	///This prints to a string. See @ref vfPrintf for details.
	///@ingroup printf
	///@param fmt format string
	///@param l typelist of arguments
	///\return The resulting application of fmt to l, returned as a string
	template<class C, class D> std::string vsPrintf(const std::string& fmt, const T_list<C,D>& l)
	{
		std::ostringstream o;
		vfPrintf(o, fmt, l);
		return o.str();
	}

	///List head for argument lists. This is a synonym for @ref TupleHead
	///
	/// Argument lists can be made by doing:
	/// @code
	/// (Fmt, arg1, arg2, ...)
	/// @endcode
	///
	///@ingroup printf
	static const T_ListEnd Fmt=TupleHead;

	//Old cvs version has by-hand variable argument list up to 10 args
}

///This is the equivalent to the C-style printf: it provides a variadic
///interface and prints to cout. See @ref vfPrintf for details.
///@param A the format string
///@param ... the arguments
///@ingroup printf
#define Printf(A, ...) vPrintf(A, (tag::Fmt,## __VA_ARGS__))

///This is the equivalent to the C-style fprintf: it provides a variadic
///interface and prints to give ostream. See @ref vfPrintf for details.
///@ingroup printf
///@param A the ostream to write to
///@param B the format string
///@param ... the arguments
#define fPrintf(A,B, ...) vfPrintf(A,B, (tag::Fmt,## __VA_ARGS__))

///This is the equivalent to the C-style sprintf: it provides a variadic
///interface returns the string. See @ref vfPrintf for details.
///@param A the format string
///@param ... the arguments
///@return The retulting application of fmt to the arguments
///@ingroup printf
#define sPrintf(A, ...) vsPrintf(A, (tag::Fmt,## __VA_ARGS__))
#endif
