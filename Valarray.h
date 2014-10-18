#include <iostream>
#include <vector>
#include <complex>
#include "EPL_traits.h"
#include <climits>
#include <cstddef>
#include <cmath>

template <typename T, typename R> class Valarray;

template <typename T, typename R> struct EPL_traits<Valarray<T, R> > {
	static const int SRank = EPL_traits<T>::SRank;
	static const bool CRank = EPL_traits<T>::CRank;
};



template <typename L, typename R, typename F>
class BinOp{
	public:
	L	left;
	R	right;
	F	func;

	BinOp(){}
	BinOp(const L& l, const R& r) : left(l), right(r){

	}
	BinOp(const BinOp& that) : left(that.left), right(that.right){
		
	}
	~BinOp(){}
	using value_type = typename ChooseType<typename L::value_type, typename R::value_type>::value_type;

	value_type operator[](int k) const{
		return func(left[k], right[k]);
	}

	unsigned int size(void) const{
		if(left.size() >= 0 && left.size() <= right.size()){
			return left.size();
		}
		else if(right.size() >=0 && right.size() <= left.size()){
			return right.size();
		}
		return 0;
	}

    class const_iterator{
        int seq = 0;
        const BinOp* b;
        using Same = const_iterator;

    public:
        typedef ptrdiff_t difference_type;
        typedef size_t size_type;

        typedef const value_type* pointer;
        typedef const value_type& reference;
        typedef std::random_access_iterator_tag iterator_category;
        const_iterator(int k, const BinOp* that) : seq(k) , b(that){
        }
	const_iterator(const const_iterator& that) : seq(that.seq) , b(that.b){}

        Same& operator=(const Same& that){
            seq=that.seq;
            b=that.b;
            return (*this);
        }

        Same& operator++(){
            seq++;
            return (*this);
        }

        Same operator++(int){
            Same temp=*this;
            seq--;
            return temp;
        }

        Same& operator--(){
            seq--;
            return (*this);
        }

        Same operator--(int){
            Same temp=*this;
            seq--;
            return temp;
        }

        bool operator==(const Same& that){
            return (seq==that.seq) && (b==that.b);
        }

        bool operator!=(const Same& that){
            return (seq!=that.seq) || (b!=that.b);
        }

        const value_type operator*(void) const{
            return (*b)[seq];
        }

        friend class BinOp<L, R, F>;

    };

    const_iterator begin(void) const{
        return const_iterator(0,this);
    }

    const_iterator end(void) const{
        return const_iterator(this->size(),this);
    }

};

template <typename T>
class Scalar{
	public:
	T val;
	Scalar(){}
	Scalar(const Scalar& that) : val(that.val){}
	Scalar(const T& arg) : val(arg){}
	~Scalar(){}

	typedef T value_type;
	value_type operator[](int k) const{
		return val;
	}

	unsigned int size (void) const{
		return UINT_MAX;
	}
};


template <typename T, typename F>
class UniOp{
	public:
	T val;
	F fun;
	UniOp(){}
	UniOp(const UniOp& that) : val(that.val){}
	UniOp(const T& arg) : val(arg){}
	~UniOp(){}

	typedef typename T::value_type value_type;
	value_type operator[](int k) const{
		return fun(val[k]);
	}

	unsigned int size (void) const{
		return val.size();
	}
};

// **********

	template <typename InputType>
	class squareRoot{
		public:
		
		static const bool iscomplex = EPL_traits<InputType>::CRank;
		typedef typename IfThen<iscomplex, std::complex<double>, double>::Type result_type;
		result_type operator() (const InputType& input) {return std::sqrt(input);}
	};

// **********

template <typename T, typename R=std::vector<T> >
class Valarray : public R{
	public:
	explicit Valarray(int k) : R(k){}	// Use k to initialize vector<T>
	explicit Valarray(const R& r) : R(r){}
	
	Valarray() : R(){}
	~Valarray(){}

	template<typename TT, typename RR>
	Valarray<T, R>& operator=(const Valarray<TT, RR>& that){
		std::vector<T>::resize(that.size());
		for(unsigned int i = 0; i < that.size(); i++){
			(*this)[i] = that[i];
		}
		return *this;
	}

	Valarray<T, R>& operator=(const T& val){
		for(unsigned int i = 0; i < this->size(); i++){
			(*this)[i] = val;
		}
		return *this;
	}

	T sum(void){
		typename R::const_iterator it(this->begin());
		T res = 0;
		for(it = this->begin(); it != this->end(); it++){
			res = res + (*it);
		}
		return res;
	}

	template <typename F>
	T accumulate(F calc, const T& ini){
		typename R::const_iterator it(this->begin());
		T accu = ini;
		for(it = this->begin(); it != this->end(); it++){
			accu = calc(accu, (*it));
		}
		return accu;
	}

	template <typename F>
	Valarray<typename F::result_type> apply(F calc){
		Valarray<typename F::result_type> temp(this->size());
		for(unsigned int i = 0; i < this->size(); i++){
			temp[i] = calc((*this)[i]);
		}
		return temp;
	}

	template <typename F = typename squareRoot<T>::result_type>
	Valarray<F> sqrt(void){
		return apply(squareRoot<T>());
	}
};




// Add
template <typename L, typename R>
Valarray<typename ChooseType<L, R>::value_type, BinOp<Valarray<L>, Valarray<R>, std::plus<typename ChooseType<L, R>::value_type> > >
	operator+(const Valarray<L>& left, const Valarray<R>& right) {
		BinOp<Valarray<L>, Valarray<R>, std::plus<typename ChooseType<L, R>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<L, R>::value_type, BinOp<Valarray<L>, Valarray<R>, std::plus<typename ChooseType<L, R>::value_type> > > (temp);
	}

template <typename T1, typename T, typename R>
Valarray<typename ChooseType<T1, T>::value_type, BinOp<Valarray<T1>, Valarray<T, R>, std::plus<typename ChooseType<T1, T>::value_type> > >
	operator+(const Valarray<T1>& left, const Valarray<T, R>& right) {
		BinOp<Valarray<T1>, Valarray<T, R>, std::plus<typename ChooseType<T1, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<T1, T>::value_type, BinOp<Valarray<T1>, Valarray<T, R>, std::plus<typename ChooseType<T1, T>::value_type> > > (temp);
	}

template <typename T1, typename T, typename R>
Valarray<typename ChooseType<T1, T>::value_type, BinOp<Valarray<T, R>, Valarray<T1>, std::plus<typename ChooseType<T1, T>::value_type> > >
	operator+(const Valarray<T, R>& left, const Valarray<T1>& right) {
		BinOp<Valarray<T, R>, Valarray<T1>, std::plus<typename ChooseType<T1, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<T1, T>::value_type, BinOp<Valarray<T, R>, Valarray<T1>, std::plus<typename ChooseType<T1, T>::value_type> > > (temp);
	}

template <typename T1, typename R1, typename T, typename R>
Valarray<typename ChooseType<T1, T>::value_type, BinOp<Valarray<T, R>, Valarray<T1, R1>, std::plus<typename ChooseType<T1, T>::value_type> > >
	operator+(const Valarray<T, R>& left, const Valarray<T1, R1>& right) {
		BinOp<Valarray<T, R>, Valarray<T1, R1>, std::plus<typename ChooseType<T1, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<T1, T>::value_type, BinOp<Valarray<T, R>, Valarray<T1, R1>, std::plus<typename ChooseType<T1, T>::value_type> > > (temp);
	}


// ********** Scalar **********

template <typename L>
Valarray<typename ChooseType<L, int>::value_type, BinOp<Valarray<L>, Scalar<int>, std::plus<typename ChooseType<L, int>::value_type> > >
	operator+(const Valarray<L>& left, const int& right) {
		BinOp<Valarray<L>, Scalar<int>, std::plus<typename ChooseType<L, int>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<L, int>::value_type, BinOp<Valarray<L>, Scalar<int>, std::plus<typename ChooseType<L, int>::value_type> > > (temp);
	}

template <typename R>
Valarray<typename ChooseType<int, R>::value_type, BinOp<Scalar<int>, Valarray<R>, std::plus<typename ChooseType<int, R>::value_type> > >
	operator+(const int& left, const Valarray<R>& right) {
		BinOp<Scalar<int>, Valarray<R>, std::plus<typename ChooseType<int, R>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<int, R>::value_type, BinOp<Scalar<int>, Valarray<R>, std::plus<typename ChooseType<int, R>::value_type> > > (temp);
	}

template <typename T, typename R>
Valarray<typename ChooseType<int, T>::value_type, BinOp<Scalar<int>, Valarray<T, R>, std::plus<typename ChooseType<int, T>::value_type> > >
	operator+(const int& left, const Valarray<T, R>& right) {
		BinOp<Scalar<int>, Valarray<T, R>, std::plus<typename ChooseType<int, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<int, T>::value_type, BinOp<Scalar<int>, Valarray<T, R>, std::plus<typename ChooseType<int, T>::value_type> > > (temp);
	}

template <typename T, typename R>
Valarray<typename ChooseType<int, T>::value_type, BinOp<Valarray<T, R>, Scalar<int>, std::plus<typename ChooseType<int, T>::value_type> > >
	operator+(const Valarray<T, R>& left, const int& right) {
		BinOp<Valarray<T, R>, Scalar<int>, std::plus<typename ChooseType<int, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<int, T>::value_type, BinOp<Valarray<T, R>, Scalar<int>, std::plus<typename ChooseType<int, T>::value_type> > > (temp);
	}

// **********

template <typename L>
Valarray<typename ChooseType<L, float>::value_type, BinOp<Valarray<L>, Scalar<float>, std::plus<typename ChooseType<L, float>::value_type> > >
	operator+(const Valarray<L>& left, const float& right) {
		BinOp<Valarray<L>, Scalar<float>, std::plus<typename ChooseType<L, float>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<L, float>::value_type, BinOp<Valarray<L>, Scalar<float>, std::plus<typename ChooseType<L, float>::value_type> > > (temp);
	}

template <typename R>
Valarray<typename ChooseType<float, R>::value_type, BinOp<Scalar<float>, Valarray<R>, std::plus<typename ChooseType<float, R>::value_type> > >
	operator+(const float& left, const Valarray<R>& right) {
		BinOp<Scalar<float>, Valarray<R>, std::plus<typename ChooseType<float, R>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<float, R>::value_type, BinOp<Scalar<float>, Valarray<R>, std::plus<typename ChooseType<float, R>::value_type> > > (temp);
	}

template <typename T, typename R>
Valarray<typename ChooseType<float, T>::value_type, BinOp<Scalar<float>, Valarray<T, R>, std::plus<typename ChooseType<float, T>::value_type> > >
	operator+(const float& left, const Valarray<T, R>& right) {
		BinOp<Scalar<float>, Valarray<T, R>, std::plus<typename ChooseType<float, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<float, T>::value_type, BinOp<Scalar<float>, Valarray<T, R>, std::plus<typename ChooseType<float, T>::value_type> > > (temp);
	}

template <typename T, typename R>
Valarray<typename ChooseType<float, T>::value_type, BinOp<Valarray<T, R>, Scalar<float>, std::plus<typename ChooseType<float, T>::value_type> > >
	operator+(const Valarray<T, R>& left, const float& right) {
		BinOp<Valarray<T, R>, Scalar<float>, std::plus<typename ChooseType<float, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<float, T>::value_type, BinOp<Valarray<T, R>, Scalar<float>, std::plus<typename ChooseType<float, T>::value_type> > > (temp);
	}

// **********

template <typename L>
Valarray<typename ChooseType<L, double>::value_type, BinOp<Valarray<L>, Scalar<double>, std::plus<typename ChooseType<L, double>::value_type> > >
	operator+(const Valarray<L>& left, const double& right) {
		BinOp<Valarray<L>, Scalar<double>, std::plus<typename ChooseType<L, double>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<L, double>::value_type, BinOp<Valarray<L>, Scalar<double>, std::plus<typename ChooseType<L, double>::value_type> > > (temp);
	}

template <typename R>
Valarray<typename ChooseType<double, R>::value_type, BinOp<Scalar<double>, Valarray<R>, std::plus<typename ChooseType<double, R>::value_type> > >
	operator+(const double& left, const Valarray<R>& right) {
		BinOp<Scalar<double>, Valarray<R>, std::plus<typename ChooseType<double, R>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<double, R>::value_type, BinOp<Scalar<double>, Valarray<R>, std::plus<typename ChooseType<double, R>::value_type> > > (temp);
	}

template <typename T, typename R>
Valarray<typename ChooseType<double, T>::value_type, BinOp<Scalar<double>, Valarray<T, R>, std::plus<typename ChooseType<double, T>::value_type> > >
	operator+(const double& left, const Valarray<T, R>& right) {
		BinOp<Scalar<double>, Valarray<T, R>, std::plus<typename ChooseType<double, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<double, T>::value_type, BinOp<Scalar<double>, Valarray<T, R>, std::plus<typename ChooseType<double, T>::value_type> > > (temp);
	}

template <typename T, typename R>
Valarray<typename ChooseType<double, T>::value_type, BinOp<Valarray<T, R>, Scalar<double>, std::plus<typename ChooseType<double, T>::value_type> > >
	operator+(const Valarray<T, R>& left, const double& right) {
		BinOp<Valarray<T, R>, Scalar<double>, std::plus<typename ChooseType<double, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<double, T>::value_type, BinOp<Valarray<T, R>, Scalar<double>, std::plus<typename ChooseType<double, T>::value_type> > > (temp);
	}


// **********

template <typename L, typename T>
Valarray<typename ChooseType<L, std::complex<T> >::value_type, BinOp<Valarray<L>, Scalar<std::complex<T> >, std::plus<typename ChooseType<L, std::complex<T> >::value_type> > >
	operator+(const Valarray<L>& left, const std::complex<T>& right) {
		BinOp<Valarray<L>, Scalar<std::complex<T> >, std::plus<typename ChooseType<L, std::complex<T> >::value_type> > temp(left, right);
                return Valarray<typename ChooseType<L, std::complex<T> >::value_type, BinOp<Valarray<L>, Scalar<std::complex<T> >, std::plus<typename ChooseType<L, std::complex<T> >::value_type> > > (temp);
	}

template <typename R, typename T>
Valarray<typename ChooseType<std::complex<T>, R>::value_type, BinOp<Scalar<std::complex<T> >, Valarray<R>, std::plus<typename ChooseType<std::complex<T>, R>::value_type> > >
	operator+(const std::complex<T>& left, const Valarray<R>& right) {
		BinOp<Scalar<std::complex<T> >, Valarray<R>, std::plus<typename ChooseType<std::complex<T>, R>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<std::complex<T>, R>::value_type, BinOp<Scalar<std::complex<T> >, Valarray<R>, std::plus<typename ChooseType<std::complex<T>, R>::value_type> > > (temp);
	}

template <typename T, typename R, typename C>
Valarray<typename ChooseType<std::complex<C>, T>::value_type, BinOp<Scalar<std::complex<C> >, Valarray<T, R>, std::plus<typename ChooseType<std::complex<C>, T>::value_type> > >
	operator+(const std::complex<C>& left, const Valarray<T, R>& right) {
		BinOp<Scalar<std::complex<C> >, Valarray<T, R>, std::plus<typename ChooseType<std::complex<C>, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<std::complex<C>, T>::value_type, BinOp<Scalar<std::complex<C> >, Valarray<T, R>, std::plus<typename ChooseType<std::complex<C>, T>::value_type> > > (temp);
	}

template <typename T, typename R, typename C>
Valarray<typename ChooseType<std::complex<C>, T>::value_type, BinOp<Valarray<T, R>, Scalar<std::complex<C> >, std::plus<typename ChooseType<std::complex<C>, T>::value_type> > >
	operator+(const Valarray<T, R>& left, const std::complex<C>& right) {
		BinOp<Valarray<T, R>, Scalar<std::complex<C> >, std::plus<typename ChooseType<std::complex<C>, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<std::complex<C>, T>::value_type, BinOp<Valarray<T, R>, Scalar<std::complex<C> >, std::plus<typename ChooseType<std::complex<C>, T>::value_type> > > (temp);
	}


// **********        **********

// Minus
template <typename L, typename R>
Valarray<typename ChooseType<L, R>::value_type, BinOp<Valarray<L>, Valarray<R>, std::minus<typename ChooseType<L, R>::value_type> > >
	operator-(const Valarray<L>& left, const Valarray<R>& right) {
		BinOp<Valarray<L>, Valarray<R>, std::minus<typename ChooseType<L, R>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<L, R>::value_type, BinOp<Valarray<L>, Valarray<R>, std::minus<typename ChooseType<L, R>::value_type> > > (temp);
	}

template <typename T1, typename T, typename R>
Valarray<typename ChooseType<T1, T>::value_type, BinOp<Valarray<T1>, Valarray<T, R>, std::minus<typename ChooseType<T1, T>::value_type> > >
	operator-(const Valarray<T1>& left, const Valarray<T, R>& right) {
		BinOp<Valarray<T1>, Valarray<T, R>, std::minus<typename ChooseType<T1, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<T1, T>::value_type, BinOp<Valarray<T1>, Valarray<T, R>, std::minus<typename ChooseType<T1, T>::value_type> > > (temp);
	}

template <typename T1, typename T, typename R>
Valarray<typename ChooseType<T1, T>::value_type, BinOp<Valarray<T, R>, Valarray<T1>, std::minus<typename ChooseType<T1, T>::value_type> > >
	operator-(const Valarray<T, R>& left, const Valarray<T1>& right) {
		BinOp<Valarray<T, R>, Valarray<T1>, std::minus<typename ChooseType<T1, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<T1, T>::value_type, BinOp<Valarray<T, R>, Valarray<T1>, std::minus<typename ChooseType<T1, T>::value_type> > > (temp);
	}

template <typename T1, typename R1, typename T, typename R>
Valarray<typename ChooseType<T1, T>::value_type, BinOp<Valarray<T, R>, Valarray<T1, R1>, std::minus<typename ChooseType<T1, T>::value_type> > >
	operator-(const Valarray<T, R>& left, const Valarray<T1, R1>& right) {
		BinOp<Valarray<T, R>, Valarray<T1, R1>, std::minus<typename ChooseType<T1, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<T1, T>::value_type, BinOp<Valarray<T, R>, Valarray<T1, R1>, std::minus<typename ChooseType<T1, T>::value_type> > > (temp);
	}

// ********** Scalar **********

template <typename L>
Valarray<typename ChooseType<L, int>::value_type, BinOp<Valarray<L>, Scalar<int>, std::minus<typename ChooseType<L, int>::value_type> > >
	operator-(const Valarray<L>& left, const int& right) {
		BinOp<Valarray<L>, Scalar<int>, std::minus<typename ChooseType<L, int>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<L, int>::value_type, BinOp<Valarray<L>, Scalar<int>, std::minus<typename ChooseType<L, int>::value_type> > > (temp);
	}

template <typename R>
Valarray<typename ChooseType<int, R>::value_type, BinOp<Scalar<int>, Valarray<R>, std::minus<typename ChooseType<int, R>::value_type> > >
	operator-(const int& left, const Valarray<R>& right) {
		BinOp<Scalar<int>, Valarray<R>, std::minus<typename ChooseType<int, R>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<int, R>::value_type, BinOp<Scalar<int>, Valarray<R>, std::minus<typename ChooseType<int, R>::value_type> > > (temp);
	}

template <typename T, typename R>
Valarray<typename ChooseType<int, T>::value_type, BinOp<Scalar<int>, Valarray<T, R>, std::minus<typename ChooseType<int, T>::value_type> > >
	operator-(const int& left, const Valarray<T, R>& right) {
		BinOp<Scalar<int>, Valarray<T, R>, std::minus<typename ChooseType<int, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<int, T>::value_type, BinOp<Scalar<int>, Valarray<T, R>, std::minus<typename ChooseType<int, T>::value_type> > > (temp);
	}

template <typename T, typename R>
Valarray<typename ChooseType<int, T>::value_type, BinOp<Valarray<T, R>, Scalar<int>, std::minus<typename ChooseType<int, T>::value_type> > >
	operator-(const Valarray<T, R>& left, const int& right) {
		BinOp<Valarray<T, R>, Scalar<int>, std::minus<typename ChooseType<int, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<int, T>::value_type, BinOp<Valarray<T, R>, Scalar<int>, std::minus<typename ChooseType<int, T>::value_type> > > (temp);
	}

// **********

template <typename L>
Valarray<typename ChooseType<L, float>::value_type, BinOp<Valarray<L>, Scalar<float>, std::minus<typename ChooseType<L, float>::value_type> > >
	operator-(const Valarray<L>& left, const float& right) {
		BinOp<Valarray<L>, Scalar<float>, std::minus<typename ChooseType<L, float>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<L, float>::value_type, BinOp<Valarray<L>, Scalar<float>, std::minus<typename ChooseType<L, float>::value_type> > > (temp);
	}

template <typename R>
Valarray<typename ChooseType<float, R>::value_type, BinOp<Scalar<float>, Valarray<R>, std::minus<typename ChooseType<float, R>::value_type> > >
	operator-(const float& left, const Valarray<R>& right) {
		BinOp<Scalar<float>, Valarray<R>, std::minus<typename ChooseType<float, R>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<float, R>::value_type, BinOp<Scalar<float>, Valarray<R>, std::minus<typename ChooseType<float, R>::value_type> > > (temp);
	}

template <typename T, typename R>
Valarray<typename ChooseType<float, T>::value_type, BinOp<Scalar<float>, Valarray<T, R>, std::minus<typename ChooseType<float, T>::value_type> > >
	operator-(const float& left, const Valarray<T, R>& right) {
		BinOp<Scalar<float>, Valarray<T, R>, std::minus<typename ChooseType<float, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<float, T>::value_type, BinOp<Scalar<float>, Valarray<T, R>, std::minus<typename ChooseType<float, T>::value_type> > > (temp);
	}

template <typename T, typename R>
Valarray<typename ChooseType<float, T>::value_type, BinOp<Valarray<T, R>, Scalar<float>, std::minus<typename ChooseType<float, T>::value_type> > >
	operator-(const Valarray<T, R>& left, const float& right) {
		BinOp<Valarray<T, R>, Scalar<float>, std::minus<typename ChooseType<float, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<float, T>::value_type, BinOp<Valarray<T, R>, Scalar<float>, std::minus<typename ChooseType<float, T>::value_type> > > (temp);
	}

// **********

template <typename L>
Valarray<typename ChooseType<L, double>::value_type, BinOp<Valarray<L>, Scalar<double>, std::minus<typename ChooseType<L, double>::value_type> > >
	operator-(const Valarray<L>& left, const double& right) {
		BinOp<Valarray<L>, Scalar<double>, std::minus<typename ChooseType<L, double>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<L, double>::value_type, BinOp<Valarray<L>, Scalar<double>, std::minus<typename ChooseType<L, double>::value_type> > > (temp);
	}

template <typename R>
Valarray<typename ChooseType<double, R>::value_type, BinOp<Scalar<double>, Valarray<R>, std::minus<typename ChooseType<double, R>::value_type> > >
	operator-(const double& left, const Valarray<R>& right) {
		BinOp<Scalar<double>, Valarray<R>, std::minus<typename ChooseType<double, R>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<double, R>::value_type, BinOp<Scalar<double>, Valarray<R>, std::minus<typename ChooseType<double, R>::value_type> > > (temp);
	}

template <typename T, typename R>
Valarray<typename ChooseType<double, T>::value_type, BinOp<Scalar<double>, Valarray<T, R>, std::minus<typename ChooseType<double, T>::value_type> > >
	operator-(const double& left, const Valarray<T, R>& right) {
		BinOp<Scalar<double>, Valarray<T, R>, std::minus<typename ChooseType<double, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<double, T>::value_type, BinOp<Scalar<double>, Valarray<T, R>, std::minus<typename ChooseType<double, T>::value_type> > > (temp);
	}

template <typename T, typename R>
Valarray<typename ChooseType<double, T>::value_type, BinOp<Valarray<T, R>, Scalar<double>, std::minus<typename ChooseType<double, T>::value_type> > >
	operator-(const Valarray<T, R>& left, const double& right) {
		BinOp<Valarray<T, R>, Scalar<double>, std::minus<typename ChooseType<double, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<double, T>::value_type, BinOp<Valarray<T, R>, Scalar<double>, std::minus<typename ChooseType<double, T>::value_type> > > (temp);
	}


// **********

template <typename L, typename T>
Valarray<typename ChooseType<L, std::complex<T> >::value_type, BinOp<Valarray<L>, Scalar<std::complex<T> >, std::minus<typename ChooseType<L, std::complex<T> >::value_type> > >
	operator-(const Valarray<L>& left, const std::complex<T>& right) {
		BinOp<Valarray<L>, Scalar<std::complex<T> >, std::minus<typename ChooseType<L, std::complex<T> >::value_type> > temp(left, right);
                return Valarray<typename ChooseType<L, std::complex<T> >::value_type, BinOp<Valarray<L>, Scalar<std::complex<T> >, std::minus<typename ChooseType<L, std::complex<T> >::value_type> > > (temp);
	}

template <typename R, typename T>
Valarray<typename ChooseType<std::complex<T>, R>::value_type, BinOp<Scalar<std::complex<T> >, Valarray<R>, std::minus<typename ChooseType<std::complex<T>, R>::value_type> > >
	operator-(const std::complex<T>& left, const Valarray<R>& right) {
		BinOp<Scalar<std::complex<T> >, Valarray<R>, std::minus<typename ChooseType<std::complex<T>, R>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<std::complex<T>, R>::value_type, BinOp<Scalar<std::complex<T> >, Valarray<R>, std::minus<typename ChooseType<std::complex<T>, R>::value_type> > > (temp);
	}

template <typename T, typename R, typename C>
Valarray<typename ChooseType<std::complex<C>, T>::value_type, BinOp<Scalar<std::complex<C> >, Valarray<T, R>, std::minus<typename ChooseType<std::complex<C>, T>::value_type> > >
	operator-(const std::complex<C>& left, const Valarray<T, R>& right) {
		BinOp<Scalar<std::complex<C> >, Valarray<T, R>, std::minus<typename ChooseType<std::complex<C>, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<std::complex<C>, T>::value_type, BinOp<Scalar<std::complex<C> >, Valarray<T, R>, std::minus<typename ChooseType<std::complex<C>, T>::value_type> > > (temp);
	}

template <typename T, typename R, typename C>
Valarray<typename ChooseType<std::complex<C>, T>::value_type, BinOp<Valarray<T, R>, Scalar<std::complex<C> >, std::minus<typename ChooseType<std::complex<C>, T>::value_type> > >
	operator-(const Valarray<T, R>& left, const std::complex<C>& right) {
		BinOp<Valarray<T, R>, Scalar<std::complex<C> >, std::minus<typename ChooseType<std::complex<C>, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<std::complex<C>, T>::value_type, BinOp<Valarray<T, R>, Scalar<std::complex<C> >, std::minus<typename ChooseType<std::complex<C>, T>::value_type> > > (temp);
	}


// **********        **********


// Multiply
template <typename L, typename R>
Valarray<typename ChooseType<L, R>::value_type, BinOp<Valarray<L>, Valarray<R>, std::multiplies<typename ChooseType<L, R>::value_type> > >
	operator*(const Valarray<L>& left, const Valarray<R>& right) {
		BinOp<Valarray<L>, Valarray<R>, std::multiplies<typename ChooseType<L, R>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<L, R>::value_type, BinOp<Valarray<L>, Valarray<R>, std::multiplies<typename ChooseType<L, R>::value_type> > > (temp);
	}

template <typename T1, typename T, typename R>
Valarray<typename ChooseType<T1, T>::value_type, BinOp<Valarray<T1>, Valarray<T, R>, std::multiplies<typename ChooseType<T1, T>::value_type> > >
	operator*(const Valarray<T1>& left, const Valarray<T, R>& right) {
		BinOp<Valarray<T1>, Valarray<T, R>, std::multiplies<typename ChooseType<T1, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<T1, T>::value_type, BinOp<Valarray<T1>, Valarray<T, R>, std::multiplies<typename ChooseType<T1, T>::value_type> > > (temp);
	}

template <typename T1, typename T, typename R>
Valarray<typename ChooseType<T1, T>::value_type, BinOp<Valarray<T, R>, Valarray<T1>, std::multiplies<typename ChooseType<T1, T>::value_type> > >
	operator*(const Valarray<T, R>& left, const Valarray<T1>& right) {
		BinOp<Valarray<T, R>, Valarray<T1>, std::multiplies<typename ChooseType<T1, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<T1, T>::value_type, BinOp<Valarray<T, R>, Valarray<T1>, std::multiplies<typename ChooseType<T1, T>::value_type> > > (temp);
	}

template <typename T1, typename R1, typename T, typename R>
Valarray<typename ChooseType<T1, T>::value_type, BinOp<Valarray<T, R>, Valarray<T1, R1>, std::multiplies<typename ChooseType<T1, T>::value_type> > >
	operator*(const Valarray<T, R>& left, const Valarray<T1, R1>& right) {
		BinOp<Valarray<T, R>, Valarray<T1, R1>, std::multiplies<typename ChooseType<T1, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<T1, T>::value_type, BinOp<Valarray<T, R>, Valarray<T1, R1>, std::multiplies<typename ChooseType<T1, T>::value_type> > > (temp);
	}

// ********** Scalar **********

template <typename L>
Valarray<typename ChooseType<L, int>::value_type, BinOp<Valarray<L>, Scalar<int>, std::multiplies<typename ChooseType<L, int>::value_type> > >
	operator*(const Valarray<L>& left, const int& right) {
		BinOp<Valarray<L>, Scalar<int>, std::multiplies<typename ChooseType<L, int>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<L, int>::value_type, BinOp<Valarray<L>, Scalar<int>, std::multiplies<typename ChooseType<L, int>::value_type> > > (temp);
	}

template <typename R>
Valarray<typename ChooseType<int, R>::value_type, BinOp<Scalar<int>, Valarray<R>, std::multiplies<typename ChooseType<int, R>::value_type> > >
	operator*(const int& left, const Valarray<R>& right) {
		BinOp<Scalar<int>, Valarray<R>, std::multiplies<typename ChooseType<int, R>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<int, R>::value_type, BinOp<Scalar<int>, Valarray<R>, std::multiplies<typename ChooseType<int, R>::value_type> > > (temp);
	}

template <typename T, typename R>
Valarray<typename ChooseType<int, T>::value_type, BinOp<Scalar<int>, Valarray<T, R>, std::multiplies<typename ChooseType<int, T>::value_type> > >
	operator*(const int& left, const Valarray<T, R>& right) {
		BinOp<Scalar<int>, Valarray<T, R>, std::multiplies<typename ChooseType<int, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<int, T>::value_type, BinOp<Scalar<int>, Valarray<T, R>, std::multiplies<typename ChooseType<int, T>::value_type> > > (temp);
	}

template <typename T, typename R>
Valarray<typename ChooseType<int, T>::value_type, BinOp<Valarray<T, R>, Scalar<int>, std::multiplies<typename ChooseType<int, T>::value_type> > >
	operator*(const Valarray<T, R>& left, const int& right) {
		BinOp<Valarray<T, R>, Scalar<int>, std::multiplies<typename ChooseType<int, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<int, T>::value_type, BinOp<Valarray<T, R>, Scalar<int>, std::multiplies<typename ChooseType<int, T>::value_type> > > (temp);
	}

// **********

template <typename L>
Valarray<typename ChooseType<L, float>::value_type, BinOp<Valarray<L>, Scalar<float>, std::multiplies<typename ChooseType<L, float>::value_type> > >
	operator*(const Valarray<L>& left, const float& right) {
		BinOp<Valarray<L>, Scalar<float>, std::multiplies<typename ChooseType<L, float>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<L, float>::value_type, BinOp<Valarray<L>, Scalar<float>, std::multiplies<typename ChooseType<L, float>::value_type> > > (temp);
	}

template <typename R>
Valarray<typename ChooseType<float, R>::value_type, BinOp<Scalar<float>, Valarray<R>, std::multiplies<typename ChooseType<float, R>::value_type> > >
	operator*(const float& left, const Valarray<R>& right) {
		BinOp<Scalar<float>, Valarray<R>, std::multiplies<typename ChooseType<float, R>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<float, R>::value_type, BinOp<Scalar<float>, Valarray<R>, std::multiplies<typename ChooseType<float, R>::value_type> > > (temp);
	}

template <typename T, typename R>
Valarray<typename ChooseType<float, T>::value_type, BinOp<Scalar<float>, Valarray<T, R>, std::multiplies<typename ChooseType<float, T>::value_type> > >
	operator*(const float& left, const Valarray<T, R>& right) {
		BinOp<Scalar<float>, Valarray<T, R>, std::multiplies<typename ChooseType<float, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<float, T>::value_type, BinOp<Scalar<float>, Valarray<T, R>, std::multiplies<typename ChooseType<float, T>::value_type> > > (temp);
	}

template <typename T, typename R>
Valarray<typename ChooseType<float, T>::value_type, BinOp<Valarray<T, R>, Scalar<float>, std::multiplies<typename ChooseType<float, T>::value_type> > >
	operator*(const Valarray<T, R>& left, const float& right) {
		BinOp<Valarray<T, R>, Scalar<float>, std::multiplies<typename ChooseType<float, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<float, T>::value_type, BinOp<Valarray<T, R>, Scalar<float>, std::multiplies<typename ChooseType<float, T>::value_type> > > (temp);
	}

// **********

template <typename L>
Valarray<typename ChooseType<L, double>::value_type, BinOp<Valarray<L>, Scalar<double>, std::multiplies<typename ChooseType<L, double>::value_type> > >
	operator*(const Valarray<L>& left, const double& right) {
		BinOp<Valarray<L>, Scalar<double>, std::multiplies<typename ChooseType<L, double>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<L, double>::value_type, BinOp<Valarray<L>, Scalar<double>, std::multiplies<typename ChooseType<L, double>::value_type> > > (temp);
	}

template <typename R>
Valarray<typename ChooseType<double, R>::value_type, BinOp<Scalar<double>, Valarray<R>, std::multiplies<typename ChooseType<double, R>::value_type> > >
	operator*(const double& left, const Valarray<R>& right) {
		BinOp<Scalar<double>, Valarray<R>, std::multiplies<typename ChooseType<double, R>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<double, R>::value_type, BinOp<Scalar<double>, Valarray<R>, std::multiplies<typename ChooseType<double, R>::value_type> > > (temp);
	}

template <typename T, typename R>
Valarray<typename ChooseType<double, T>::value_type, BinOp<Scalar<double>, Valarray<T, R>, std::multiplies<typename ChooseType<double, T>::value_type> > >
	operator*(const double& left, const Valarray<T, R>& right) {
		BinOp<Scalar<double>, Valarray<T, R>, std::multiplies<typename ChooseType<double, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<double, T>::value_type, BinOp<Scalar<double>, Valarray<T, R>, std::multiplies<typename ChooseType<double, T>::value_type> > > (temp);
	}

template <typename T, typename R>
Valarray<typename ChooseType<double, T>::value_type, BinOp<Valarray<T, R>, Scalar<double>, std::multiplies<typename ChooseType<double, T>::value_type> > >
	operator*(const Valarray<T, R>& left, const double& right) {
		BinOp<Valarray<T, R>, Scalar<double>, std::multiplies<typename ChooseType<double, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<double, T>::value_type, BinOp<Valarray<T, R>, Scalar<double>, std::multiplies<typename ChooseType<double, T>::value_type> > > (temp);
	}


// **********

template <typename L, typename T>
Valarray<typename ChooseType<L, std::complex<T> >::value_type, BinOp<Valarray<L>, Scalar<std::complex<T> >, std::multiplies<typename ChooseType<L, std::complex<T> >::value_type> > >
	operator*(const Valarray<L>& left, const std::complex<T>& right) {
		BinOp<Valarray<L>, Scalar<std::complex<T> >, std::multiplies<typename ChooseType<L, std::complex<T> >::value_type> > temp(left, right);
                return Valarray<typename ChooseType<L, std::complex<T> >::value_type, BinOp<Valarray<L>, Scalar<std::complex<T> >, std::multiplies<typename ChooseType<L, std::complex<T> >::value_type> > > (temp);
	}

template <typename R, typename T>
Valarray<typename ChooseType<std::complex<T>, R>::value_type, BinOp<Scalar<std::complex<T> >, Valarray<R>, std::multiplies<typename ChooseType<std::complex<T>, R>::value_type> > >
	operator*(const std::complex<T>& left, const Valarray<R>& right) {
		BinOp<Scalar<std::complex<T> >, Valarray<R>, std::multiplies<typename ChooseType<std::complex<T>, R>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<std::complex<T>, R>::value_type, BinOp<Scalar<std::complex<T> >, Valarray<R>, std::multiplies<typename ChooseType<std::complex<T>, R>::value_type> > > (temp);
	}

template <typename T, typename R, typename C>
Valarray<typename ChooseType<std::complex<C>, T>::value_type, BinOp<Scalar<std::complex<C> >, Valarray<T, R>, std::multiplies<typename ChooseType<std::complex<C>, T>::value_type> > >
	operator*(const std::complex<C>& left, const Valarray<T, R>& right) {
		BinOp<Scalar<std::complex<C> >, Valarray<T, R>, std::multiplies<typename ChooseType<std::complex<C>, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<std::complex<C>, T>::value_type, BinOp<Scalar<std::complex<C> >, Valarray<T, R>, std::multiplies<typename ChooseType<std::complex<C>, T>::value_type> > > (temp);
	}

template <typename T, typename R, typename C>
Valarray<typename ChooseType<std::complex<C>, T>::value_type, BinOp<Valarray<T, R>, Scalar<std::complex<C> >, std::multiplies<typename ChooseType<std::complex<C>, T>::value_type> > >
	operator*(const Valarray<T, R>& left, const std::complex<C>& right) {
		BinOp<Valarray<T, R>, Scalar<std::complex<C> >, std::multiplies<typename ChooseType<std::complex<C>, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<std::complex<C>, T>::value_type, BinOp<Valarray<T, R>, Scalar<std::complex<C> >, std::multiplies<typename ChooseType<std::complex<C>, T>::value_type> > > (temp);
	}


// **********        **********


// Divide
template <typename L, typename R>
Valarray<typename ChooseType<L, R>::value_type, BinOp<Valarray<L>, Valarray<R>, std::divides<typename ChooseType<L, R>::value_type> > >
	operator/(const Valarray<L>& left, const Valarray<R>& right) {
		BinOp<Valarray<L>, Valarray<R>, std::divides<typename ChooseType<L, R>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<L, R>::value_type, BinOp<Valarray<L>, Valarray<R>, std::divides<typename ChooseType<L, R>::value_type> > > (temp);
	}

template <typename T1, typename T, typename R>
Valarray<typename ChooseType<T1, T>::value_type, BinOp<Valarray<T1>, Valarray<T, R>, std::divides<typename ChooseType<T1, T>::value_type> > >
	operator/(const Valarray<T1>& left, const Valarray<T, R>& right) {
		BinOp<Valarray<T1>, Valarray<T, R>, std::divides<typename ChooseType<T1, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<T1, T>::value_type, BinOp<Valarray<T1>, Valarray<T, R>, std::divides<typename ChooseType<T1, T>::value_type> > > (temp);
	}

template <typename T1, typename T, typename R>
Valarray<typename ChooseType<T1, T>::value_type, BinOp<Valarray<T, R>, Valarray<T1>, std::divides<typename ChooseType<T1, T>::value_type> > >
	operator/(const Valarray<T, R>& left, const Valarray<T1>& right) {
		BinOp<Valarray<T, R>, Valarray<T1>, std::divides<typename ChooseType<T1, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<T1, T>::value_type, BinOp<Valarray<T, R>, Valarray<T1>, std::divides<typename ChooseType<T1, T>::value_type> > > (temp);
	}

template <typename T1, typename R1, typename T, typename R>
Valarray<typename ChooseType<T1, T>::value_type, BinOp<Valarray<T, R>, Valarray<T1, R1>, std::divides<typename ChooseType<T1, T>::value_type> > >
	operator/(const Valarray<T, R>& left, const Valarray<T1, R1>& right) {
		BinOp<Valarray<T, R>, Valarray<T1, R1>, std::divides<typename ChooseType<T1, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<T1, T>::value_type, BinOp<Valarray<T, R>, Valarray<T1, R1>, std::divides<typename ChooseType<T1, T>::value_type> > > (temp);
	}

// ********** Scalar **********

template <typename L>
Valarray<typename ChooseType<L, int>::value_type, BinOp<Valarray<L>, Scalar<int>, std::divides<typename ChooseType<L, int>::value_type> > >
	operator/(const Valarray<L>& left, const int& right) {
		BinOp<Valarray<L>, Scalar<int>, std::divides<typename ChooseType<L, int>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<L, int>::value_type, BinOp<Valarray<L>, Scalar<int>, std::divides<typename ChooseType<L, int>::value_type> > > (temp);
	}

template <typename R>
Valarray<typename ChooseType<int, R>::value_type, BinOp<Scalar<int>, Valarray<R>, std::divides<typename ChooseType<int, R>::value_type> > >
	operator/(const int& left, const Valarray<R>& right) {
		BinOp<Scalar<int>, Valarray<R>, std::divides<typename ChooseType<int, R>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<int, R>::value_type, BinOp<Scalar<int>, Valarray<R>, std::divides<typename ChooseType<int, R>::value_type> > > (temp);
	}

template <typename T, typename R>
Valarray<typename ChooseType<int, T>::value_type, BinOp<Scalar<int>, Valarray<T, R>, std::divides<typename ChooseType<int, T>::value_type> > >
	operator/(const int& left, const Valarray<T, R>& right) {
		BinOp<Scalar<int>, Valarray<T, R>, std::divides<typename ChooseType<int, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<int, T>::value_type, BinOp<Scalar<int>, Valarray<T, R>, std::divides<typename ChooseType<int, T>::value_type> > > (temp);
	}

template <typename T, typename R>
Valarray<typename ChooseType<int, T>::value_type, BinOp<Valarray<T, R>, Scalar<int>, std::divides<typename ChooseType<int, T>::value_type> > >
	operator/(const Valarray<T, R>& left, const int& right) {
		BinOp<Valarray<T, R>, Scalar<int>, std::divides<typename ChooseType<int, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<int, T>::value_type, BinOp<Valarray<T, R>, Scalar<int>, std::divides<typename ChooseType<int, T>::value_type> > > (temp);
	}

// **********

template <typename L>
Valarray<typename ChooseType<L, float>::value_type, BinOp<Valarray<L>, Scalar<float>, std::divides<typename ChooseType<L, float>::value_type> > >
	operator/(const Valarray<L>& left, const float& right) {
		BinOp<Valarray<L>, Scalar<float>, std::divides<typename ChooseType<L, float>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<L, float>::value_type, BinOp<Valarray<L>, Scalar<float>, std::divides<typename ChooseType<L, float>::value_type> > > (temp);
	}

template <typename R>
Valarray<typename ChooseType<float, R>::value_type, BinOp<Scalar<float>, Valarray<R>, std::divides<typename ChooseType<float, R>::value_type> > >
	operator/(const float& left, const Valarray<R>& right) {
		BinOp<Scalar<float>, Valarray<R>, std::divides<typename ChooseType<float, R>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<float, R>::value_type, BinOp<Scalar<float>, Valarray<R>, std::divides<typename ChooseType<float, R>::value_type> > > (temp);
	}

template <typename T, typename R>
Valarray<typename ChooseType<float, T>::value_type, BinOp<Scalar<float>, Valarray<T, R>, std::divides<typename ChooseType<float, T>::value_type> > >
	operator/(const float& left, const Valarray<T, R>& right) {
		BinOp<Scalar<float>, Valarray<T, R>, std::divides<typename ChooseType<float, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<float, T>::value_type, BinOp<Scalar<float>, Valarray<T, R>, std::divides<typename ChooseType<float, T>::value_type> > > (temp);
	}

template <typename T, typename R>
Valarray<typename ChooseType<float, T>::value_type, BinOp<Valarray<T, R>, Scalar<float>, std::divides<typename ChooseType<float, T>::value_type> > >
	operator/(const Valarray<T, R>& left, const float& right) {
		BinOp<Valarray<T, R>, Scalar<float>, std::divides<typename ChooseType<float, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<float, T>::value_type, BinOp<Valarray<T, R>, Scalar<float>, std::divides<typename ChooseType<float, T>::value_type> > > (temp);
	}

// **********

template <typename L>
Valarray<typename ChooseType<L, double>::value_type, BinOp<Valarray<L>, Scalar<double>, std::divides<typename ChooseType<L, double>::value_type> > >
	operator/(const Valarray<L>& left, const double& right) {
		BinOp<Valarray<L>, Scalar<double>, std::divides<typename ChooseType<L, double>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<L, double>::value_type, BinOp<Valarray<L>, Scalar<double>, std::divides<typename ChooseType<L, double>::value_type> > > (temp);
	}

template <typename R>
Valarray<typename ChooseType<double, R>::value_type, BinOp<Scalar<double>, Valarray<R>, std::divides<typename ChooseType<double, R>::value_type> > >
	operator/(const double& left, const Valarray<R>& right) {
		BinOp<Scalar<double>, Valarray<R>, std::divides<typename ChooseType<double, R>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<double, R>::value_type, BinOp<Scalar<double>, Valarray<R>, std::divides<typename ChooseType<double, R>::value_type> > > (temp);
	}

template <typename T, typename R>
Valarray<typename ChooseType<double, T>::value_type, BinOp<Scalar<double>, Valarray<T, R>, std::divides<typename ChooseType<double, T>::value_type> > >
	operator/(const double& left, const Valarray<T, R>& right) {
		BinOp<Scalar<double>, Valarray<T, R>, std::divides<typename ChooseType<double, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<double, T>::value_type, BinOp<Scalar<double>, Valarray<T, R>, std::divides<typename ChooseType<double, T>::value_type> > > (temp);
	}

template <typename T, typename R>
Valarray<typename ChooseType<double, T>::value_type, BinOp<Valarray<T, R>, Scalar<double>, std::divides<typename ChooseType<double, T>::value_type> > >
	operator/(const Valarray<T, R>& left, const double& right) {
		BinOp<Valarray<T, R>, Scalar<double>, std::divides<typename ChooseType<double, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<double, T>::value_type, BinOp<Valarray<T, R>, Scalar<double>, std::divides<typename ChooseType<double, T>::value_type> > > (temp);
	}


// **********

template <typename L, typename T>
Valarray<typename ChooseType<L, std::complex<T> >::value_type, BinOp<Valarray<L>, Scalar<std::complex<T> >, std::divides<typename ChooseType<L, std::complex<T> >::value_type> > >
	operator/(const Valarray<L>& left, const std::complex<T>& right) {
		BinOp<Valarray<L>, Scalar<std::complex<T> >, std::divides<typename ChooseType<L, std::complex<T> >::value_type> > temp(left, right);
                return Valarray<typename ChooseType<L, std::complex<T> >::value_type, BinOp<Valarray<L>, Scalar<std::complex<T> >, std::divides<typename ChooseType<L, std::complex<T> >::value_type> > > (temp);
	}

template <typename R, typename T>
Valarray<typename ChooseType<std::complex<T>, R>::value_type, BinOp<Scalar<std::complex<T> >, Valarray<R>, std::divides<typename ChooseType<std::complex<T>, R>::value_type> > >
	operator/(const std::complex<T>& left, const Valarray<R>& right) {
		BinOp<Scalar<std::complex<T> >, Valarray<R>, std::divides<typename ChooseType<std::complex<T>, R>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<std::complex<T>, R>::value_type, BinOp<Scalar<std::complex<T> >, Valarray<R>, std::divides<typename ChooseType<std::complex<T>, R>::value_type> > > (temp);
	}

template <typename T, typename R, typename C>
Valarray<typename ChooseType<std::complex<C>, T>::value_type, BinOp<Scalar<std::complex<C> >, Valarray<T, R>, std::divides<typename ChooseType<std::complex<C>, T>::value_type> > >
	operator/(const std::complex<C>& left, const Valarray<T, R>& right) {
		BinOp<Scalar<std::complex<C> >, Valarray<T, R>, std::divides<typename ChooseType<std::complex<C>, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<std::complex<C>, T>::value_type, BinOp<Scalar<std::complex<C> >, Valarray<T, R>, std::divides<typename ChooseType<std::complex<C>, T>::value_type> > > (temp);
	}

template <typename T, typename R, typename C>
Valarray<typename ChooseType<std::complex<C>, T>::value_type, BinOp<Valarray<T, R>, Scalar<std::complex<C> >, std::divides<typename ChooseType<std::complex<C>, T>::value_type> > >
	operator/(const Valarray<T, R>& left, const std::complex<C>& right) {
		BinOp<Valarray<T, R>, Scalar<std::complex<C> >, std::divides<typename ChooseType<std::complex<C>, T>::value_type> > temp(left, right);
                return Valarray<typename ChooseType<std::complex<C>, T>::value_type, BinOp<Valarray<T, R>, Scalar<std::complex<C> >, std::divides<typename ChooseType<std::complex<C>, T>::value_type> > > (temp);
	}

// **********        **********

template <typename T>
Valarray<T, UniOp<Valarray<T>, std::negate<T> > >
operator-(const Valarray<T>& that){
	UniOp<Valarray<T>, std::negate<T> > temp(that);
	return Valarray<T, UniOp<Valarray<T>, std::negate<T> > > (temp);
}

template <typename T, typename R>
Valarray<T, UniOp<Valarray<T, R>, std::negate<T> > >
operator-(const Valarray<T, R>& that){
	UniOp<Valarray<T, R>, std::negate<T> > temp(that);
	return Valarray<T, UniOp<Valarray<T, R>, std::negate<T> > > (temp);
}


// **********        **********

template <typename T, typename R>
std::ostream& operator<<(std::ostream& out, const Valarray<T, R>& x) {
	const char* prefix = "";
	for (T v : x) {
		out << prefix << v;
		prefix = ", ";
	}
	return out;
}


