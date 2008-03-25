#ifndef TAG_TUPLE_H
#define TAG_TUPLE_H

namespace tag
{
	////////////////////////////////////////////////////////////////////////////////	
	//
	// Typesafe, varadic arguments are provided by a typelist
	//

	#ifndef  DOXYGEN_IGNORE_INTERNAL
	namespace Internal
	{
		struct null{};
		
		static null nv;

		//Index the list, retrieving the value and the type
		//with a good optimizer, this should be constant time at run-time.
		template<template<class,class> class List, class C, class D, int i> struct index_l
		{
			typedef typename index_l<List, typename D::val_type, typename D::next_type, i-1>::val_type val_type;	

			static const val_type& value(const List<C,D>& l)
			{	
				return index_l<List, typename D::val_type, typename D::next_type, i-1>::value(l.next);
			}
		};

		template<template<class,class> class List, class C, class D> struct index_l<List, C,D,0>
		{
			typedef C val_type;
			static const val_type& value(const List<C,D>& l)
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
	

		template<template<class,class> class List, class C, class D, int i> struct T_index_forward
		{
			enum
			{
				len = length<C,D>::len
			};

			typedef typename index_l<List, C,D,len-i-1>::val_type val_type;

			static const val_type& value(const List<C,D>& l)
			{	
				return index_l<List, C, D, len-i-1>::value(l);
			}
		};
	}
	#endif


	/** @defgroup tuple Tuple types of arbitrary length
	    This group containg classes for dealing with tuples (or typelists)
		of arbitrary length. 
		
		There are two types:
		- tag::T_list: Typelist with <code>const&</code> mebers.
			- Fully functional type list which supports static indexing of both values and types.
			- Used for typesafe variadic functions such as @ref Printf
			- See tag::T_list for a usage example.

		- tag::V_list: Typelist with value members.
			- Used for returning tuples from functions, since references can not be used.
			- tag::V_tuple is a convinience class for generating the types.
			- Supports basic functionality, but not indexing.

		V_list is designed primarily to allow multiple
		return values from functions (see also tag::rpair) and can be uses as follows:
		@code
			V_tuple<int, float, char>::type func()
			{
				int i;
				float f;
				char c;

				//Function body goes here
				return make_vtuple(i, f, c);
			}

		    void foo()
			{
				int a;
				float b;
				char c;

				//Call func() and collect all return values
				make_rtuple(a,b,c) = func();
				
				//a,b,c are now set.
			}
		@endcode


	*/

	/** Tuple/ typelist class. Passed by const reference.
	@ingroup tuple
	
	 A typelist containing an int, float and char* represented by the class:
	@code
	T_list<int, T_list<float, T_list< char*, T_ListEnd> > >
	@endcode

	Typelists can be conviniently be built through use of the <code>,</code>
	operator. The following code will have the type given above:
	@code
	(TupleHead, "hello", 2.2f, 1)
	@endcode
	or more convieniently with a macro:
	@code
	make_tuple("hello", 2.2f, 1)
	@endcode
	Note that the order is reversed between how the type is written and
	how a list is written using <code>operator,()</code>. The lists can be 
	statically indexed, and the indexing matches the order written using 
	<code>,</code>. That is index 0 will refer to <code>(char*) "hello"</code>.
	Indexing can be performed by value:
	@code
	char* ptr = some_tuple.index<0>();
	@endcode
	and by type:
	@code
	some_tuples_type::index_type<0>::type a_value = some_tuple.index<0>();
	@endcode
	The length of the list is a static member, so the last element can be eccessed
	by:
	@code
		some_tuple.index<some_tuples_type::elements-1>()
	@endcode





	The following code gives a simple example of the usage. The code simply
	prints out a list of values.
	@code
	#include <tag/tuple.h>
	#include <iostream>

	using namespace tag;
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
		template<int i> const typename Internal::T_index_forward<tag::T_list, C,D,i>::val_type& index() const
		{
			return Internal::T_index_forward<tag::T_list, C,D,i>::value(*this);
		}

		///Index the list in logical order and return the type. Index 0 is the first element added to 
		///the list (the one closest to @ref T_ListEnd)
		template<int i> class index_type
		{	
			typedef typename Internal::T_index_forward<tag::T_list, C,D,i>::val_type type;
		};

		enum
		{
			///The number of elements in the list.
			elements = (Internal::length<C,D>::len)
		};
		
	};
	

	///This is a Tuple/typelist class where members are stored by value.
	///It is designed primarily allowing multiple return values from 
	///functions. See ``@ref tuple'', @ref make_vtuple and @ref make_rtuple. You probably
	///do not waht to use this class directly, since it is rather cumbersome. See
	///tag::V_list.
	///@ingroup tuple
	template<class C, class D> struct V_list
	{
		#ifndef DOXYGEN_IGNORE_INTERNAL
		C val;
		typedef C val_type;

		D next;
		typedef D next_type;
		
		
		V_list(const C& c, const D& d)
		:val(c),next(d){}

		template<class X> V_list<X, V_list<C, D> > operator,(const X& c) 
		{
			 return V_list<X, V_list<C, D> >(c, *this);
		}
		
		template<int i> const typename Internal::T_index_forward<tag::V_list, C,D,i>::val_type& index() const
		{
			return Internal::T_index_forward<tag::V_list, C,D,i>::value(*this);
		}

		template<int i> class index_type
		{	
			typedef typename Internal::T_index_forward<tag::V_list, C,D,i>::val_type type;
		};

		enum
		{
			///The number of elements in the list.
			elements = (Internal::length<C,D>::len)
		};
		
		#endif
	};


	#ifndef DOXYGEN_IGNORE_INTERNAL
	//This class converts a typelist in to an equivalent V_list
	//This allows R_list to figure out the T_list it should allow
	//assignment from.
	namespace Internal
	{
		template<class X, class Y> struct r2v
		{
			typedef V_list<X, typename r2v<typename Y::val_type, typename Y::next_type>::v> v;
		};

		template<> struct r2v<Internal::null, Internal::null>
		{
			typedef V_list<null, null> v;
		};
	}
	#endif

	#ifndef DOXYGEN_IGNORE_INTERNAL
	template<class C, class D> struct R_list
	{
		private:
			typedef typename Internal::r2v<C,D>::v equivalent_V_list;
			R_list(const R_list& );
			void operator=(const R_list&);


		public:
			C& val;
			typedef C val_type;

			D& next;
			typedef D next_type;
			
			R_list(C& c, D& d)
			:val(c),next(d){}

			template<class X> R_list<X, R_list<C, D> > operator,(X& c) 
			{
				 return R_list<X, R_list<C, D> >(c, *this);
			}

			template<int i> const typename Internal::T_index_forward<tag::R_list, C,D,i>::val_type& index() const
			{
				return Internal::T_index_forward<tag::R_list, C,D,i>::value(*this);
			}

			template<int i> class index_type
			{	
				typedef typename Internal::T_index_forward<tag::R_list, C,D,i>::val_type type;
			};

			enum
			{
				///The number of elements in the list.
				elements = (Internal::length<C,D>::len)
			};
		
	
			template<class X, class Y> friend class R_list;
			void operator=(const equivalent_V_list& vlist)
			{
				val = vlist.val;
				next = vlist.next;
			}
	};
	#endif
	
	#ifndef DOXYGEN_IGNORE_INTERNAL
	typedef R_list<Internal::null,Internal::null> R_ListEnd;
	typedef V_list<Internal::null,Internal::null> V_ListEnd;


	static  R_ListEnd R_TupleHead(Internal::nv, Internal::nv);
	static  V_ListEnd V_TupleHead(Internal::nv, Internal::nv);
	template<class A, class B, class C=Internal::null, class D=Internal::null, class E=Internal::null, class F=Internal::null, class G=Internal::null, class H=Internal::null, class I=Internal::null, class J=Internal::null> struct V_tuple
	{
		typedef V_list<J,V_list<I, V_list<H, V_list<G, V_list<F, V_list<E, V_list<D, V_list<C, V_list<B, V_list<A, V_ListEnd> > > > > > > > > > type;
	};

	template<class A, class B, class C, class D, class E, class F, class G, class H, class I> struct V_tuple<A,B,C,D,E,F,G,H,I,Internal::null>
	{
		typedef V_list<I, V_list<H, V_list<G, V_list<F, V_list<E, V_list<D, V_list<C, V_list<B, V_list<A, V_ListEnd> > > > > > > > > type;
	};

	template<class A, class B, class C, class D, class E, class F, class G, class H> struct V_tuple<A,B,C,D,E,F,G,H,Internal::null,Internal::null>
	{
		typedef V_list<H, V_list<G, V_list<F, V_list<E, V_list<D, V_list<C, V_list<B, V_list<A, V_ListEnd> > > > > > > > type;
	};

	template<class A, class B, class C, class D, class E, class F, class G> struct V_tuple<A,B,C,D,E,F,G,Internal::null,Internal::null,Internal::null>
	{
		typedef V_list<G, V_list<F, V_list<E, V_list<D, V_list<C, V_list<B, V_list<A, V_ListEnd> > > > > > > type;
	};

	template<class A, class B, class C, class D, class E, class F> struct V_tuple<A,B,C,D,E,F,Internal::null,Internal::null,Internal::null,Internal::null>
	{
		typedef V_list<F, V_list<E, V_list<D, V_list<C, V_list<B, V_list<A, V_ListEnd> > > > > > type;
	};

	template<class A, class B, class C, class D, class E> struct V_tuple<A,B,C,D,E,Internal::null,Internal::null,Internal::null,Internal::null,Internal::null>
	{
		typedef V_list<E, V_list<D, V_list<C, V_list<B, V_list<A, V_ListEnd> > > > > type;
	};

	template<class A, class B, class C, class D> struct V_tuple<A,B,C,D,Internal::null,Internal::null,Internal::null,Internal::null,Internal::null,Internal::null>
	{
		typedef V_list<D, V_list<C, V_list<B, V_list<A, V_ListEnd> > > > type;
	};

	template<class A, class B, class C> struct V_tuple<A,B,C,Internal::null,Internal::null,Internal::null,Internal::null,Internal::null,Internal::null,Internal::null>
	{
		typedef V_list<C, V_list<B, V_list<A, V_ListEnd> > > type;
	};

	template<class A, class B> struct V_tuple<A,B,Internal::null,Internal::null,Internal::null,Internal::null,Internal::null,Internal::null,Internal::null>
	{
		typedef V_list<B, V_list<A, V_ListEnd> > type;
	};
	#else
	///A convinience type for generating tag::V_list  types. Note that 
	///the arguments are written in the order which matches use of <code>operator,()</code>.
	///This can be used to generate V_list types with up to 10 elements.
	///@ingroup tuple
	template<class A, class B, class C, ...> struct V_tuple
	{
		///The V_list type with elements A, B, C, ...
		typedef V_list<..., V_list<C, V_list<B, V_list<A, V_ListEnd> > > > type;
	};
	#endif





	///Guard type to mark the end of the list	
	///@ingroup tuple
	typedef T_list<Internal::null,Internal::null> T_ListEnd;
	///The value of the end of a list
	///@ingroup tuple
	static const T_ListEnd TupleHead=T_ListEnd(Internal::null(), Internal::null());
}

///This is a variadic ``function'' which creates and returns a tuple.
///@ingroup tuple
#define make_tuple(...) ((tag::TupleHead, ## __VA_ARGS__))


///This macro is used to return multiple values from a function.
///@ingroup tuple
#define make_vtuple(...) ((tag::V_TupleHead, ## __VA_ARGS__))

///This macro is used to collect multiple return values from a function.
///@ingroup tuple
#define make_rtuple(...) ((tag::R_TupleHead, ## __VA_ARGS__))


#endif


