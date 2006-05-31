#ifndef TAG_TUPLE_H
#define TAG_TUPLE_H

namespace tag
{
	////////////////////////////////////////////////////////////////////////////////	
	//
	// Typesafe, varadic arguments are provided by a typelist
	//

	template<class C, class D> class T_list;

	#ifndef  DOXYGEN_IGNORE_INTERNAL
	namespace Internal
	{
		struct null{};
		
		//Index the list, retrieving the value and the type
		//with a good optimizer, this should be constant time at run-time.
		template<class C, class D, int i> struct index_l
		{
			typedef typename index_l<typename D::val_type, typename D::next_type, i-1>::val_type val_type;	

			static const val_type& value(const T_list<C,D>& l)
			{	
				return index_l<typename D::val_type, typename D::next_type, i-1>::value(l.next);
			}
		};

		template<class C, class D> struct index_l<C,D,0>
		{
			typedef C val_type;
			static const val_type& value(const T_list<C,D>& l)
			{
				return l.val;
			}
		};

		//Compute the length of a list
		template<class C, class D> struct length
		{
			enum
			{
				len = length<typename D::val_type, typename D::next_type>::len + 1
			};
		};

		template<> struct length<null, null>
		{
			enum
			{
				len = 0
			};
		};
	

		template<class C, class D, int i> struct T_index_forward
		{
			enum
			{
				len = length<C,D>::len
			};

			typedef typename index_l<C,D,len-i-1>::val_type val_type;

			static const val_type& value(const T_list<C,D>& l)
			{	
				return index_l<C, D, len-i-1>::value(l);
			}
		};
	}
	#endif


	/// @defgroup tuple Tuple types of arbitrary length
	/// This group containg classes for dealing with tuples (or typelists)
	/// of arbitrary length. These are used for making variadic functions.

	/** Tuple/ typelist class.
	@ingroup tuple
	
	 A typelist containing an int, float and char* represented by the class:
	@code
	T_list<int, T_list<float, T_list< char*, @ref T_ListEnd> > >
	@endcode

	The following code gives a simple example of the usage. The code simply
	prints out a list of values.
	@code
	#include <tag/tuple.h>
	#include <iostream>

	using namespace EdUtil;
	using namespace std;

	template<class C, int i, int max> struct dump_list
	{
		static void dump(const C& l)
		{
			cout << l.template index<i>()<< endl;
			dump_list<C, i+1,  max>::dump(l);
		}
	};

	template<class C, int max> struct dump_list<C, max, max>
	{
		static void dump(const C&)
		{}
	};

	template<class C, class D> void print_typelist(const T_list<C, D>& l)
	{
		dump_list<T_list<C,D>,0,T_list<C,D>::elements>::dump(l);
	}

	#define print(...)  print_typelist((TupleHead,##__VA_ARGS__))

	int main()
	{
		print(1, "hello", 3.3, (void*)"aa", "world", cin.rdbuf(), 99.9);
	}
	@endcode

	This will print
	@code
	1
	hello
	3.3
	0xdeadbeefwhatever
	world
	<whatever the user types here>
	99.0
	@endcode
	
	**/
	
	template<class C, class D> struct T_list
	{
		///The value of the current element (the end of the list)
		const C& val;
		/// The type of the current element
		typedef C val_type;

		///The rest of the list
		const D& next;
		/// The type of the rest of the list
		typedef D next_type;
		
		
		///Construct a typelist.
		///@param c The value of the end of the list
		///@param d The rest of the list
		T_list(const C& c, const D& d)
		:val(c),next(d){}

		///This operator can be used to build a typelist from values in code
		///
		///The list:
		///@code
		///T_list<int, T_list<float, T_list< char*, T_ListEnd> > >
		///@endcode
		///containing 1, 2.0 and 3.3 can be built up by:
		///@code
		///TupleHead, 1, 2.0, 3.3
		///@endcode
		///Though it is important to note that the list is built in reverse, so that 
		///3.3 is the most accessible element.
		///@param c List to which an element should be appended.
		template<class X> T_list<X, T_list<C, D> > operator,(const X& c) const
		{
			 return T_list<X, T_list<C, D> >(c, *this);
		}

		///Index the list in logical order and return the value. Index 0 is the first element added to 
		///the list (the one closest to @ref T_ListEnd)
		template<int i> const typename Internal::T_index_forward<C,D,i>::val_type& index() const
		{
			return Internal::T_index_forward<C,D,i>::value(*this);
		}

		///Index the list in logical order and return the type. Index 0 is the first element added to 
		///the list (the one closest to @ref T_ListEnd)
		template<int i> class index_type
		{	
			typedef typename Internal::T_index_forward<C,D,i>::val_type type;
		};

		enum
		{
			///The number of elements in the list.
			elements = (Internal::length<C,D>::len)
		};
		
	};

	///Guard type to mark the end of the list	
	///@ingroup tuple
	typedef T_list<Internal::null,Internal::null> T_ListEnd;
	///The value of the end of a list
	///@ingroup tuple
	static const T_ListEnd TupleHead=T_ListEnd(Internal::null(), Internal::null());
}

///This is a variadic ``function'' which creates and returns a tuple.
///@ingroup tuple
#define make_tuple(...) ((TupleHead, ## __VA_ARGS__))

#endif


