#ifndef TAG_FN_H
#define TAG_FN_H

#include <functional>
#include <iterator>

namespace tag {

///@defgroup functional additional functors for <functional>
///This group provides additional functors to complement the <functional> header of STL.
///@ingroup stdpp
//@{

template <typename A, typename m>
struct mem_data_ref_t : std::unary_function<A &, m &> {
    m A::*data;
    inline mem_data_ref_t( m A::*d ) : data(d) {};
    inline m & operator()(A & instance) const {
        return instance.*data;
    }
};

template <typename A, typename m>
struct const_mem_data_ref_t : std::unary_function<const A &, const m &> {
    const m A::*data;
    inline const_mem_data_ref_t( const m A::*d ) : data(d) {};
    inline const m & operator()(const A & instance) const {
        return instance.*data;
    }
};

template <typename A, typename m>
inline struct mem_data_ref_t<A,m> mem_data_ref( m A::*data ){
    return mem_data_ref_t<A,m>(data);
}

template <typename A, typename m>
inline struct const_mem_data_ref_t<A,m> const_mem_data_ref( const m A::*data ){
    return const_mem_data_ref_t<A,m>(data);
}

template <typename A, typename m>
struct mem_data_t : std::unary_function<A *, m &> {
    m A::*data;
    inline mem_data_t( m A::*d ) : data(d) {};
    inline m & operator()(A * instance) const {
        return instance->*data;
    }
};

template <typename A, typename m>
struct const_mem_data_t : std::unary_function<const A *, const m &> {
    const m A::*data;
    inline const_mem_data_t(const m A::*d ) : data(d) {};
    inline const m & operator()(const A * instance) const {
        return instance->*data;
    }
};

template <typename A, typename m>
inline struct mem_data_t<A,m> mem_data( m A::*data ){
    return mem_data_t<A,m>(data);
}

template <typename A, typename m>
inline struct const_mem_data_t<A,m> const_mem_data( const m A::*data ){
    return const_mem_data_t<A,m>(data);
}

template <typename G, typename F>
struct bind_t : std::unary_function<typename F::argument_type, typename G::result_type> {
    const F & f;
    const G & g;
    inline bind_t( const F & f_, const G & g_) :f(f_), g(g_) {};
    inline typename G::result_type operator()( typename F::argument_type a) const {
        return g(f(a));
    }
};

template <typename G, typename F>
inline struct bind_t<G,F> bind( const G & g, const F & f ){
    return bind_t<G,F>(f,g);
}

//@}

///@defgroup iterator additional iterators
///This group provides additional iterators to complement the <iterator> header of STL.
///@ingroup stdpp
//@{


/**
An iterator wrapper that returns a member of a struct the wrapped iterator would point to.
@code
struct simple { int a; float b; };
vector<simple> test;
member_iterator_t<vector<simple>::iterator, int> ita( &simple::a );
ita = test.begin();
cout << *ita; // prints the value of a
@endcode
*/
template <typename It, typename m>
struct member_iterator_t : public It {
    typedef typename std::iterator_traits<It>::value_type Value;
    m Value::*data;
    inline member_iterator_t( m Value::*d ) : data(d) {};
    inline member_iterator_t( const It & it,  m Value::*d ) : data(d) { *this = it; };
    template <typename Other> inline member_iterator_t & operator=(const Other & other) {
        It::operator=(other);
        return *this;
    }
    inline member_iterator_t & operator=( const member_iterator_t & other){
        data = other.data;
        It::operator=(other);
        return *this;
    }
    inline m & operator*(void){
        return It::operator*().*data;
    }
    inline m & operator->(void){
        return It::operator->()->*data;
    }
    inline const m & operator*(void) const {
        return It::operator*().*data;
    }
    inline const m & operator->(void) const {
        return It::operator->()->*data;
    }
};

/**
helper function to simplify the use of @ref member_iterator_t wrapper. This is useful for passing
member iterators as arguments.
@arg it the iterator to wrap, the new member_iterator_t returned will point to the same position
@arg d the member to wrap
*/
template <typename It, typename m>
inline struct member_iterator_t<It, m> member_iterator( const It & it, m std::iterator_traits<It>::value_type::*d ){
    return member_iterator_t<It, m>(it, d);
}

//@}

}

#endif // TAG_FN_H
