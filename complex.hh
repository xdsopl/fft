/*
complex - fast complex math
Written in 2015 by <Ahmet Inan> <xdsopl@googlemail.com>
To the extent possible under law, the author(s) have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide. This software is distributed without any warranty.
You should have received a copy of the CC0 Public Domain Dedication along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/

#ifndef COMPLEX_HH
#define COMPLEX_HH

template <typename T>
class Complex
{
	T re, im;
public:
	typedef T value_type;
	Complex() : re(0), im(0) {}
	Complex(T r, T i) : re(r), im(i) {}
	inline T &real() { return re; }
	inline T &imag() { return im; }
	inline Complex<T> operator = (T a)
	{
		real() = a;
		imag() = 0;
		return *this;
	}
	inline Complex<T> operator += (Complex<T> a)
	{
		return *this = a + *this;
	}
	inline Complex<T> operator *= (Complex<T> a)
	{
		return *this = a * *this;
	}
	inline Complex<T> operator *= (T a)
	{
		return *this = a * *this;
	}
	inline Complex<T> operator /= (T a)
	{
		return *this = *this / a;
	}
};

template <typename T>
static inline Complex<T> operator + (Complex<T> a, Complex<T> b)
{
	return Complex<T>(a.real() + b.real(), a.imag() + b.imag());
}

template <typename T>
static inline Complex<T> operator + (Complex<T> a)
{
	return a;
}

template <typename T>
static inline Complex<T> operator - (Complex<T> a, Complex<T> b)
{
	return Complex<T>(a.real() - b.real(), a.imag() - b.imag());
}

template <typename T>
static inline Complex<T> operator - (Complex<T> a)
{
	return Complex<T>(-a.real(), -a.imag());
}

template <typename T>
static inline Complex<T> operator * (T a, Complex<T> b)
{
	return Complex<T>(a * b.real(), a * b.imag());
}

template <typename T>
static inline Complex<T> operator / (Complex<T> a, T b)
{
	return Complex<T>(a.real() / b, a.imag() / b);
}

template <typename T>
static inline Complex<T> operator * (Complex<T> a, Complex<T> b)
{
	return Complex<T>(a.real() * b.real() - a.imag() * b.imag(), a.real() * b.imag() + a.imag() * b.real());
}

template <typename T>
static inline Complex<T> operator / (Complex<T> a, Complex<T> b)
{
	return Complex<T>((a.real() * b.real() + a.imag() * b.imag()) / (b.real() * b.real() + b.imag() * b.imag()),
			(a.imag() * b.real() - a.real() * b.imag()) / (b.real() * b.real() + b.imag() * b.imag()));
}

template <typename T>
static inline Complex<T> exp(Complex<T> a)
{
	return Complex<T>(exp(a.real()) * cos(a.imag()), exp(a.real()) * sin(a.imag()));
}

template <typename T>
static inline T abs(Complex<T> a)
{
	return hypot(a.real(), a.imag());
}

template <typename T>
static inline T arg(Complex<T> a)
{
	return atan2(a.imag(), a.real());
}

#endif
