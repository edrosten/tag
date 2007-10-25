#ifndef TAG_INC_ARRAY_H
#define TAG_INC_ARRAY_H
#include <tag/tuple.h>

namespace tag
{

template<class C=int, int I=-1> class array;

template<> struct array<int, -1>
{
	struct Underfill{};
};


template<class C, int I> class array
{
	public:

		typedef C* iterator;
		typedef const C* const_iterator;

		iterator begin(){return data;}
		const_iterator begin()const {return data;}

		iterator end(){return data+I;}
		const_iterator end()const {return data+I;}
			
		operator C*() { return data; }
		operator const C*() const { return data; }

		C& operator*(){return *data;}
		const C& operator*() const {return *data;}

		C& operator[](int i){return data[i];}
		const C& operator[](int i)const {return data[i];}

		C* operator+(int i){return data+i;}
		const C* operator+(int i)const{return data+i;}


		int size()const
		{
			return I;
		}

		array(){}

		template<class D, class E> array(const T_list<D, E>& l)
		{
			//array<> is just used as a placeholder class to prevent template
			//instantiation happening too early. A proper class (instead of a
			//builtin) is required to generate the proper error messages.
			typedef typename  sizetoolong<array<>, (T_list<D,E>::elements > I) >::dummy  foo;
			typedef typename sizetooshort<array<>, (T_list<D,E>::elements < I) >::dummy  bar;
			array_filler<D, E, T_list<D,E>::elements -1>::fill(data, l);
		}

		template<class D, class E> array(const T_list<typename array<int,-1>::Underfill, T_list<D,E> >& l)
		{
			typedef typename  sizetoolong<array<>, (T_list<D,E>::elements > I) >::dummy  foo;
			array_filler<D, E, T_list<D,E>::elements -1>::fill(data, l.next);
		}

	private:
		C data[I];



		template<class D, class E, int i> struct array_filler
		{
			static void fill(C* data, const T_list<D, E>& l)
			{
				data[i] = l.val;
				array_filler<typename E::val_type, typename E::next_type, i-1>::fill(data, l.next);
			}
		};

		template<class D, class E> struct array_filler<D,E,-1>
		{
			static void fill(C*, const T_list<D, E>&){}
		};

		template<class D, int i> struct sizetoolong
		{
			typedef typename D::Error_initializer_list_is_too_long dummy;
		};

		template<class D> struct sizetoolong<D, 0>
		{
			typedef int dummy;
		};

		template<class D, int i> struct sizetooshort
		{
			typedef typename D::Error_initializer_list_is_too_short dummy;
		};

		template<class D> struct sizetooshort<D, 0>
		{
			typedef int dummy;
		};
};
}
#endif
