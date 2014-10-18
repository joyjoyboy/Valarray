/*
 * EPL_traits.h
 *
 *  Created on: Apr 30, 2013
 *      Author: chase
 */

#ifndef EPL_TRAITS_H_
#define EPL_TRAITS_H_


#include <complex>
#include <vector>
#include <string>
//#include <Valarray.h>


template <bool P, typename T, typename F>
struct IfThen {
	using Type = T;
};

template <typename T, typename F>
struct IfThen<false, T, F> {
	using Type = F;
};

/* the default is empty
 * in that way, if EPL_traits is used on an unsupported type, then we'll have an compile-time error
 */
template <typename T> struct EPL_traits { };

template <> struct EPL_traits<void> { // some useful constants, exported to the real metafunctions below
	static const int INT = 1;
	static const int FLOAT = 2;
	static const int DOUBLE = 3;

	static std::string baseTypeName(const int s) {
		switch (s) { // s will be the SRank
		case INT: return "int";
		case FLOAT: return "float";
		case DOUBLE: return "double";
		default: return "unrecognized scalar base type";
		}
	}
};

template <> struct EPL_traits<int> {
	static const int SRank = EPL_traits<void>::INT;
	static const bool CRank = false;
};

template <> struct EPL_traits<float> {
	static const int SRank = EPL_traits<void>::FLOAT;
	static const bool CRank = false;
};

template <> struct EPL_traits<double> {
	static const int SRank = EPL_traits<void>::DOUBLE;
	static const bool CRank = false;
};

template <typename T> struct EPL_traits<std::complex<T>> {
	static const int SRank = EPL_traits<T>::SRank;
	static const bool CRank = true;
};



template <int rank> struct SType;
template <> struct SType<EPL_traits<void>::INT> {  using Type = int; };
template <> struct SType<EPL_traits<void>::FLOAT> {  using Type = float; };
template <> struct SType<EPL_traits<void>::DOUBLE> {  using Type = double; };

constexpr uint64_t max(uint64_t x, uint64_t y) {
	return (x < y) ? y : x;
}

constexpr uint64_t min(uint64_t x, uint64_t y) {
	return (y < x) ? y : x;
}

template <typename L, typename R>
struct ChooseType {
	static const int lrank = EPL_traits<L>::SRank;
	static const int rrank = EPL_traits<R>::SRank;
	static const int rank = max(lrank, rrank);
	using scalar_value_type = typename SType<rank>::Type;

	static const bool lcomplex = EPL_traits<L>::CRank;
	static const bool rcomplex = EPL_traits<R>::CRank;
	static const bool comp = lcomplex || rcomplex;

	using value_type = typename IfThen<comp, std::complex<scalar_value_type>, scalar_value_type>::Type;
};


#endif /* EPL_TRAITS_H_ */
