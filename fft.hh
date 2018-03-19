/*
fft - mixed radix fft
Written in 2015 by <Ahmet Inan> <xdsopl@googlemail.com>
To the extent possible under law, the author(s) have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide. This software is distributed without any warranty.
You should have received a copy of the CC0 Public Domain Dedication along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/

#ifndef FFT_HH
#define FFT_HH

namespace FFT {

template <typename TYPE>
static inline TYPE half(TYPE a) { return typename TYPE::value_type(0.5) * a; }
template <typename TYPE>
static inline TYPE sqrt3(TYPE a) { return std::sqrt(typename TYPE::value_type(3)) * a; }
template <typename TYPE>
static inline TYPE rsqrt2(TYPE a) { return (std::sqrt(typename TYPE::value_type(2)) / typename TYPE::value_type(2)) * a; }
template <typename TYPE>
static inline TYPE cos2pi5(TYPE a) { return std::cos(typename TYPE::value_type(2 * M_PI / 5)) * a; }
template <typename TYPE>
static inline TYPE sin2pi5(TYPE a) { return std::sin(typename TYPE::value_type(2 * M_PI / 5)) * a; }
template <typename TYPE>
static inline TYPE cos4pi5(TYPE a) { return std::cos(typename TYPE::value_type(4 * M_PI / 5)) * a; }
template <typename TYPE>
static inline TYPE sin4pi5(TYPE a) { return std::sin(typename TYPE::value_type(4 * M_PI / 5)) * a; }
template <typename TYPE>
static inline TYPE cos2pi7(TYPE a) { return std::cos(typename TYPE::value_type(2 * M_PI / 7)) * a; }
template <typename TYPE>
static inline TYPE sin2pi7(TYPE a) { return std::sin(typename TYPE::value_type(2 * M_PI / 7)) * a; }
template <typename TYPE>
static inline TYPE cos4pi7(TYPE a) { return std::cos(typename TYPE::value_type(4 * M_PI / 7)) * a; }
template <typename TYPE>
static inline TYPE sin4pi7(TYPE a) { return std::sin(typename TYPE::value_type(4 * M_PI / 7)) * a; }
template <typename TYPE>
static inline TYPE cos6pi7(TYPE a) { return std::cos(typename TYPE::value_type(6 * M_PI / 7)) * a; }
template <typename TYPE>
static inline TYPE sin6pi7(TYPE a) { return std::sin(typename TYPE::value_type(6 * M_PI / 7)) * a; }

template <typename TYPE>
static inline TYPE fiddle(TYPE a, TYPE b)
{
	return TYPE(a.real() + a.imag() - b.real() + b.imag(), a.imag() - a.real() - b.real() - b.imag());
}

template <typename TYPE>
static inline TYPE twiddle(TYPE a, TYPE b)
{
	return TYPE(a.imag() - b.imag(), b.real() - a.real());
}

static constexpr int pow2(int N)
{
	return !(N & (N - 1));
}

static constexpr int pow4(int N)
{
	return pow2(N) && (N & 0x55555555);
}

static constexpr int pow8(int N)
{
	return pow2(N) && (N & 0xb6db6db6);
}

static constexpr int split(int N)
{
	return (!(N%7)) ? 7 : (!(N%5)) ? 5 : (!(N%3)) ? 3 : (!(N%8)&&pow8(N)) ? 8 : (!(N%4)&&pow4(N)) ? 4 : (!(N%2)) ? 2 : 1;
}

template <int RADIX, int BINS, int STRIDE, typename TYPE, int SIGN>
struct Dit {};

template <int STRIDE, typename TYPE, int SIGN>
struct Dit<1, 1, STRIDE, TYPE, SIGN>
{
	static inline void dit(TYPE *out, const TYPE *in, const TYPE *)
	{
		*out = *in;
	}
};

template <int STRIDE, typename TYPE, int SIGN>
struct Dit<2, 2, STRIDE, TYPE, SIGN>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE in0, TYPE in1)
	{
		*out0 = in0 + in1;
		*out1 = in0 - in1;
	}
	static inline void dit(TYPE *out, const TYPE *in, const TYPE *)
	{
		dft(out, out + 1, in[0], in[STRIDE]);
	}
};

template <int STRIDE, typename TYPE>
struct Dit<3, 3, STRIDE, TYPE, -1>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE *out2, TYPE in0, TYPE in1, TYPE in2)
	{
		TYPE a(in1 + in2), b(sqrt3(twiddle(in2, in1)));
		*out0 = in0 + a;
		*out1 = in0 - half(a + b);
		*out2 = in0 - half(a - b);
	}
	static inline void dit(TYPE *out, const TYPE *in, const TYPE *)
	{
		dft(out, out + 1, out + 2, in[0], in[STRIDE], in[2 * STRIDE]);
	}
};

template <int STRIDE, typename TYPE>
struct Dit<3, 3, STRIDE, TYPE, 1>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE *out2, TYPE in0, TYPE in1, TYPE in2)
	{
		TYPE a(in1 + in2), b(sqrt3(twiddle(in2, in1)));
		*out0 = in0 + a;
		*out1 = in0 - half(a - b);
		*out2 = in0 - half(a + b);
	}
	static inline void dit(TYPE *out, const TYPE *in, const TYPE *)
	{
		dft(out, out + 1, out + 2, in[0], in[STRIDE], in[2 * STRIDE]);
	}
};

template <int STRIDE, typename TYPE>
struct Dit<4, 4, STRIDE, TYPE, -1>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE *out2, TYPE *out3,
			TYPE in0, TYPE in1, TYPE in2, TYPE in3)
	{
		TYPE a(in0 + in2), b(in0 - in2), c(in1 + in3), d(twiddle(in1, in3));
		*out0 = a + c;
		*out1 = b + d;
		*out2 = a - c;
		*out3 = b - d;
	}
	static inline void dit(TYPE *out, const TYPE *in, const TYPE *)
	{
		dft(out, out + 1, out + 2, out + 3, in[0], in[STRIDE], in[2 * STRIDE], in[3 * STRIDE]);
	}
};

template <int STRIDE, typename TYPE>
struct Dit<4, 4, STRIDE, TYPE, 1>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE *out2, TYPE *out3,
			TYPE in0, TYPE in1, TYPE in2, TYPE in3)
	{
		TYPE a(in0 + in2), b(in0 - in2), c(in1 + in3), d(twiddle(in1, in3));
		*out0 = a + c;
		*out1 = b - d;
		*out2 = a - c;
		*out3 = b + d;
	}
	static inline void dit(TYPE *out, const TYPE *in, const TYPE *)
	{
		dft(out, out + 1, out + 2, out + 3, in[0], in[STRIDE], in[2 * STRIDE], in[3 * STRIDE]);
	}
};

template <int STRIDE, typename TYPE>
struct Dit<5, 5, STRIDE, TYPE, -1>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE *out2, TYPE *out3, TYPE *out4,
			TYPE in0, TYPE in1, TYPE in2, TYPE in3, TYPE in4)
	{
		TYPE a(in1 + in4), b(in2 + in3), c(twiddle(in1, in4)), d(twiddle(in2, in3));
		TYPE c1a(cos2pi5(a)), c1b(cos2pi5(b)), s1c(sin2pi5(c)), s1d(sin2pi5(d));
		TYPE c2a(cos4pi5(a)), c2b(cos4pi5(b)), s2c(sin4pi5(c)), s2d(sin4pi5(d));
		*out0 = in0 + a + b;
		*out1 = in0 + c1a + c2b + s1c + s2d;
		*out2 = in0 + c2a + c1b + s2c - s1d;
		*out3 = in0 + c2a + c1b - s2c + s1d;
		*out4 = in0 + c1a + c2b - s1c - s2d;
	}
	static inline void dit(TYPE *out, const TYPE *in, const TYPE *)
	{
		dft(out, out + 1, out + 2, out + 3, out + 4, in[0], in[STRIDE], in[2 * STRIDE], in[3 * STRIDE], in[4 * STRIDE]);
	}
};

template <int STRIDE, typename TYPE>
struct Dit<5, 5, STRIDE, TYPE, 1>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE *out2, TYPE *out3, TYPE *out4,
			TYPE in0, TYPE in1, TYPE in2, TYPE in3, TYPE in4)
	{
		TYPE a(in1 + in4), b(in2 + in3), c(twiddle(in1, in4)), d(twiddle(in2, in3));
		TYPE c1a(cos2pi5(a)), c1b(cos2pi5(b)), s1c(sin2pi5(c)), s1d(sin2pi5(d));
		TYPE c2a(cos4pi5(a)), c2b(cos4pi5(b)), s2c(sin4pi5(c)), s2d(sin4pi5(d));
		*out0 = in0 + a + b;
		*out1 = in0 + c1a + c2b - s1c - s2d;
		*out2 = in0 + c2a + c1b - s2c + s1d;
		*out3 = in0 + c2a + c1b + s2c - s1d;
		*out4 = in0 + c1a + c2b + s1c + s2d;
	}
	static inline void dit(TYPE *out, const TYPE *in, const TYPE *)
	{
		dft(out, out + 1, out + 2, out + 3, out + 4, in[0], in[STRIDE], in[2 * STRIDE], in[3 * STRIDE], in[4 * STRIDE]);
	}
};

template <int STRIDE, typename TYPE>
struct Dit<7, 7, STRIDE, TYPE, -1>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE *out2, TYPE *out3, TYPE *out4, TYPE *out5, TYPE *out6,
			TYPE in0, TYPE in1, TYPE in2, TYPE in3, TYPE in4, TYPE in5, TYPE in6)
	{
		TYPE a(in1 + in6), b(in2 + in5), c(in3 + in4), d(twiddle(in1, in6)), e(twiddle(in2, in5)), f(twiddle(in3, in4));
		TYPE c1a(cos2pi7(a)), c1b(cos2pi7(b)), c1c(cos2pi7(c)), s1d(sin2pi7(d)), s1e(sin2pi7(e)), s1f(sin2pi7(f));
		TYPE c2a(cos4pi7(a)), c2b(cos4pi7(b)), c2c(cos4pi7(c)), s2d(sin4pi7(d)), s2e(sin4pi7(e)), s2f(sin4pi7(f));
		TYPE c3a(cos6pi7(a)), c3b(cos6pi7(b)), c3c(cos6pi7(c)), s3d(sin6pi7(d)), s3e(sin6pi7(e)), s3f(sin6pi7(f));

		*out0 = in0 + a + b + c;
		*out1 = in0 + c1a + c2b + c3c + s1d + s2e + s3f;
		*out2 = in0 + c2a + c3b + c1c + s2d - s3e - s1f;
		*out3 = in0 + c3a + c1b + c2c + s3d - s1e + s2f;
		*out4 = in0 + c3a + c1b + c2c - s3d + s1e - s2f;
		*out5 = in0 + c2a + c3b + c1c - s2d + s3e + s1f;
		*out6 = in0 + c1a + c2b + c3c - s1d - s2e - s3f;
	}
	static inline void dit(TYPE *out, const TYPE *in, const TYPE *)
	{
		dft(out, out + 1, out + 2, out + 3, out + 4, out + 5, out + 6, in[0], in[STRIDE], in[2 * STRIDE], in[3 * STRIDE], in[4 * STRIDE], in[5 * STRIDE], in[6 * STRIDE]);
	}
};

template <int STRIDE, typename TYPE>
struct Dit<7, 7, STRIDE, TYPE, 1>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE *out2, TYPE *out3, TYPE *out4, TYPE *out5, TYPE *out6,
			TYPE in0, TYPE in1, TYPE in2, TYPE in3, TYPE in4, TYPE in5, TYPE in6)
	{
		TYPE a(in1 + in6), b(in2 + in5), c(in3 + in4), d(twiddle(in1, in6)), e(twiddle(in2, in5)), f(twiddle(in3, in4));
		TYPE c1a(cos2pi7(a)), c1b(cos2pi7(b)), c1c(cos2pi7(c)), s1d(sin2pi7(d)), s1e(sin2pi7(e)), s1f(sin2pi7(f));
		TYPE c2a(cos4pi7(a)), c2b(cos4pi7(b)), c2c(cos4pi7(c)), s2d(sin4pi7(d)), s2e(sin4pi7(e)), s2f(sin4pi7(f));
		TYPE c3a(cos6pi7(a)), c3b(cos6pi7(b)), c3c(cos6pi7(c)), s3d(sin6pi7(d)), s3e(sin6pi7(e)), s3f(sin6pi7(f));

		*out0 = in0 + a + b + c;
		*out1 = in0 + c1a + c2b + c3c - s1d - s2e - s3f;
		*out2 = in0 + c2a + c3b + c1c - s2d + s3e + s1f;
		*out3 = in0 + c3a + c1b + c2c - s3d + s1e - s2f;
		*out4 = in0 + c3a + c1b + c2c + s3d - s1e + s2f;
		*out5 = in0 + c2a + c3b + c1c + s2d - s3e - s1f;
		*out6 = in0 + c1a + c2b + c3c + s1d + s2e + s3f;
	}
	static inline void dit(TYPE *out, const TYPE *in, const TYPE *)
	{
		dft(out, out + 1, out + 2, out + 3, out + 4, out + 5, out + 6, in[0], in[STRIDE], in[2 * STRIDE], in[3 * STRIDE], in[4 * STRIDE], in[5 * STRIDE], in[6 * STRIDE]);
	}
};

template <int STRIDE, typename TYPE>
struct Dit<8, 8, STRIDE, TYPE, -1>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE *out2, TYPE *out3, TYPE *out4, TYPE *out5, TYPE *out6, TYPE *out7,
			TYPE in0, TYPE in1, TYPE in2, TYPE in3, TYPE in4, TYPE in5, TYPE in6, TYPE in7)
	{
		TYPE a(in0 + in4), b(in0 - in4), c(in1 + in5), d(in1 - in5), e(in2 + in6), f(twiddle(in2, in6)), g(in3 + in7), h(in3 - in7);
		TYPE cpg(c + g), tcg(twiddle(c, g)), fdh(rsqrt2(fiddle(d, h))), fhd(rsqrt2(fiddle(h, d)));
		*out0 = a + e + cpg;
		*out1 = b + f + fdh;
		*out2 = a - e + tcg;
		*out3 = b - f + fhd;
		*out4 = a + e - cpg;
		*out5 = b + f - fdh;
		*out6 = a - e - tcg;
		*out7 = b - f - fhd;
	}
	static inline void dit(TYPE *out, const TYPE *in, const TYPE *)
	{
		dft(out, out + 1, out + 2, out + 3, out + 4, out + 5, out + 6, out + 7, in[0], in[STRIDE], in[2 * STRIDE], in[3 * STRIDE], in[4 * STRIDE], in[5 * STRIDE], in[6 * STRIDE], in[7 * STRIDE]);
	}
};

template <int STRIDE, typename TYPE>
struct Dit<8, 8, STRIDE, TYPE, 1>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE *out2, TYPE *out3, TYPE *out4, TYPE *out5, TYPE *out6, TYPE *out7,
			TYPE in0, TYPE in1, TYPE in2, TYPE in3, TYPE in4, TYPE in5, TYPE in6, TYPE in7)
	{
		TYPE a(in0 + in4), b(in0 - in4), c(in1 + in5), d(in1 - in5), e(in2 + in6), f(twiddle(in2, in6)), g(in3 + in7), h(in3 - in7);
		TYPE cpg(c + g), tcg(twiddle(c, g)), fdh(rsqrt2(fiddle(d, h))), fhd(rsqrt2(fiddle(h, d)));
		*out0 = a + e + cpg;
		*out1 = b - f - fhd;
		*out2 = a - e - tcg;
		*out3 = b + f - fdh;
		*out4 = a + e - cpg;
		*out5 = b - f + fhd;
		*out6 = a - e + tcg;
		*out7 = b + f + fdh;
	}
	static inline void dit(TYPE *out, const TYPE *in, const TYPE *)
	{
		dft(out, out + 1, out + 2, out + 3, out + 4, out + 5, out + 6, out + 7, in[0], in[STRIDE], in[2 * STRIDE], in[3 * STRIDE], in[4 * STRIDE], in[5 * STRIDE], in[6 * STRIDE], in[7 * STRIDE]);
	}
};

template <int BINS, int STRIDE, typename TYPE, int SIGN>
struct Dit<2, BINS, STRIDE, TYPE, SIGN>
{
	static const int RADIX = 2;
	static const int QUOTIENT = BINS / RADIX;
	static void dit(TYPE *out, const TYPE *in, const TYPE *z)
	{
		for (int o = 0, i = 0; o < BINS; o += QUOTIENT, i += STRIDE)
			Dit<split(QUOTIENT), QUOTIENT, RADIX * STRIDE, TYPE, SIGN>::dit(out + o, in + i, z);
		for (int k0 = 0, k1 = QUOTIENT, l1 = 0; k0 < QUOTIENT; ++k0, ++k1, l1 += STRIDE)
			Dit<RADIX, RADIX, STRIDE, TYPE, SIGN>::dft(out + k0, out + k1, out[k0], z[l1] * out[k1]);
	}
};

template <int BINS, int STRIDE, typename TYPE, int SIGN>
struct Dit<3, BINS, STRIDE, TYPE, SIGN>
{
	static const int RADIX = 3;
	static const int QUOTIENT = BINS / RADIX;
	static void dit(TYPE *out, const TYPE *in, const TYPE *z)
	{
		for (int o = 0, i = 0; o < BINS; o += QUOTIENT, i += STRIDE)
			Dit<split(QUOTIENT), QUOTIENT, RADIX * STRIDE, TYPE, SIGN>::dit(out + o, in + i, z);
		for (int k0 = 0, k1 = QUOTIENT, k2 = 2 * QUOTIENT,
				l1 = 0, l2 = 0;
				k0 < QUOTIENT;
				++k0, ++k1, ++k2,
				l1 += STRIDE, l2 += 2 * STRIDE)
			Dit<RADIX, RADIX, STRIDE, TYPE, SIGN>::dft(out + k0, out + k1, out + k2,
				out[k0], z[l1] * out[k1], z[l2] * out[k2]);
	}
};

template <int BINS, int STRIDE, typename TYPE, int SIGN>
struct Dit<4, BINS, STRIDE, TYPE, SIGN>
{
	static const int RADIX = 4;
	static const int QUOTIENT = BINS / RADIX;
	static void dit(TYPE *out, const TYPE *in, const TYPE *z)
	{
		for (int o = 0, i = 0; o < BINS; o += QUOTIENT, i += STRIDE)
			Dit<split(QUOTIENT), QUOTIENT, RADIX * STRIDE, TYPE, SIGN>::dit(out + o, in + i, z);
		for (int k0 = 0, k1 = QUOTIENT, k2 = 2 * QUOTIENT, k3 = 3 * QUOTIENT,
				l1 = 0, l2 = 0, l3 = 0;
				k0 < QUOTIENT;
				++k0, ++k1, ++k2, ++k3,
				l1 += STRIDE, l2 += 2 * STRIDE, l3 += 3 * STRIDE)
			Dit<RADIX, RADIX, STRIDE, TYPE, SIGN>::dft(out + k0, out + k1, out + k2, out + k3,
				out[k0], z[l1] * out[k1], z[l2] * out[k2], z[l3] * out[k3]);
	}
};

template <int BINS, int STRIDE, typename TYPE, int SIGN>
struct Dit<5, BINS, STRIDE, TYPE, SIGN>
{
	static const int RADIX = 5;
	static const int QUOTIENT = BINS / RADIX;
	static void dit(TYPE *out, const TYPE *in, const TYPE *z)
	{
		for (int o = 0, i = 0; o < BINS; o += QUOTIENT, i += STRIDE)
			Dit<split(QUOTIENT), QUOTIENT, RADIX * STRIDE, TYPE, SIGN>::dit(out + o, in + i, z);
		for (int k0 = 0, k1 = QUOTIENT, k2 = 2 * QUOTIENT, k3 = 3 * QUOTIENT, k4 = 4 * QUOTIENT,
				l1 = 0, l2 = 0, l3 = 0, l4 = 0;
				k0 < QUOTIENT;
				++k0, ++k1, ++k2, ++k3, ++k4,
				l1 += STRIDE, l2 += 2 * STRIDE, l3 += 3 * STRIDE, l4 += 4 * STRIDE)
			Dit<RADIX, RADIX, STRIDE, TYPE, SIGN>::dft(out + k0, out + k1, out + k2, out + k3, out + k4,
				out[k0], z[l1] * out[k1], z[l2] * out[k2], z[l3] * out[k3], z[l4] * out[k4]);
	}
};

template <int BINS, int STRIDE, typename TYPE, int SIGN>
struct Dit<7, BINS, STRIDE, TYPE, SIGN>
{
	static const int RADIX = 7;
	static const int QUOTIENT = BINS / RADIX;
	static void dit(TYPE *out, const TYPE *in, const TYPE *z)
	{
		for (int o = 0, i = 0; o < BINS; o += QUOTIENT, i += STRIDE)
			Dit<split(QUOTIENT), QUOTIENT, RADIX * STRIDE, TYPE, SIGN>::dit(out + o, in + i, z);
		for (int k0 = 0, k1 = QUOTIENT, k2 = 2 * QUOTIENT, k3 = 3 * QUOTIENT, k4 = 4 * QUOTIENT, k5 = 5 * QUOTIENT, k6 = 6 * QUOTIENT,
				l1 = 0, l2 = 0, l3 = 0, l4 = 0, l5 = 0, l6 = 0;
				k0 < QUOTIENT;
				++k0, ++k1, ++k2, ++k3, ++k4, ++k5, ++k6,
				l1 += STRIDE, l2 += 2 * STRIDE, l3 += 3 * STRIDE, l4 += 4 * STRIDE, l5 += 5 * STRIDE, l6 += 6 * STRIDE)
			Dit<RADIX, RADIX, STRIDE, TYPE, SIGN>::dft(out + k0, out + k1, out + k2, out + k3, out + k4, out + k5, out + k6,
				out[k0], z[l1] * out[k1], z[l2] * out[k2], z[l3] * out[k3], z[l4] * out[k4], z[l5] * out[k5], z[l6] * out[k6]);
	}
};

template <int BINS, int STRIDE, typename TYPE, int SIGN>
struct Dit<8, BINS, STRIDE, TYPE, SIGN>
{
	static const int RADIX = 8;
	static const int QUOTIENT = BINS / RADIX;
	static void dit(TYPE *out, const TYPE *in, const TYPE *z)
	{
		for (int o = 0, i = 0; o < BINS; o += QUOTIENT, i += STRIDE)
			Dit<split(QUOTIENT), QUOTIENT, RADIX * STRIDE, TYPE, SIGN>::dit(out + o, in + i, z);
		for (int k0 = 0, k1 = QUOTIENT, k2 = 2 * QUOTIENT, k3 = 3 * QUOTIENT, k4 = 4 * QUOTIENT, k5 = 5 * QUOTIENT, k6 = 6 * QUOTIENT, k7 = 7 * QUOTIENT,
				l1 = 0, l2 = 0, l3 = 0, l4 = 0, l5 = 0, l6 = 0, l7 = 0;
				k0 < QUOTIENT;
				++k0, ++k1, ++k2, ++k3, ++k4, ++k5, ++k6, ++k7,
				l1 += STRIDE, l2 += 2 * STRIDE, l3 += 3 * STRIDE, l4 += 4 * STRIDE, l5 += 5 * STRIDE, l6 += 6 * STRIDE, l7 += 7 * STRIDE)
			Dit<RADIX, RADIX, STRIDE, TYPE, SIGN>::dft(out + k0, out + k1, out + k2, out + k3, out + k4, out + k5, out + k6, out + k7,
				out[k0], z[l1] * out[k1], z[l2] * out[k2], z[l3] * out[k3], z[l4] * out[k4], z[l5] * out[k5], z[l6] * out[k6], z[l7] * out[k7]);
	}
};

template <int BINS, typename TYPE, int SIGN>
class Factors
{
	TYPE z[BINS];
public:
	typedef typename TYPE::value_type value_type;
	Factors()
	{
		for (int n = 0; n < BINS; ++n)
			z[n] = exp(TYPE(0, value_type(SIGN * 2 * M_PI) * value_type(n) / value_type(BINS)));
	}
	inline operator const TYPE * () const
	{
		return z;
	}
};

template <int BINS, typename TYPE>
class Forward
{
	Factors<BINS, TYPE, -1> factors;
public:
	typedef typename TYPE::value_type value_type;
	inline void operator ()(TYPE *out, const TYPE *in)
	{
		Dit<split(BINS), BINS, 1, TYPE, -1>::dit(out, in, factors);
	}
};

template <int BINS, typename TYPE>
class Backward
{
	Factors<BINS, TYPE, 1> factors;
public:
	typedef typename TYPE::value_type value_type;
	inline void operator ()(TYPE *out, const TYPE *in)
	{
		Dit<split(BINS), BINS, 1, TYPE, 1>::dit(out, in, factors);
	}
};

template <int BINS, typename TYPE>
class Normalize
{
	typedef typename TYPE::value_type value_type;
	const value_type factor = sqrt(value_type(1) / value_type(BINS));
public:
	inline void operator ()(TYPE *io)
	{
		for (int n = 0; n < BINS; ++n)
			io[n] *= factor;
	}
};

}

#endif
