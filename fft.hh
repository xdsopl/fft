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
static inline TYPE half(TYPE a)
{
	return typename TYPE::value_type(0.5) * a;
}

template <typename TYPE>
static inline TYPE sqrt3(TYPE a)
{
	return std::sqrt(typename TYPE::value_type(3)) * a;
}

template <typename TYPE>
static inline TYPE rsqrt2(TYPE a)
{
	return (std::sqrt(typename TYPE::value_type(2)) / typename TYPE::value_type(2)) * a;
}

template <int n, int N, typename TYPE>
static inline TYPE cx(TYPE a)
{
	return std::cos(typename TYPE::value_type(n * 2 * M_PI / N)) * a;
}

template <int n, int N, typename TYPE>
static inline TYPE sx(TYPE a)
{
	return std::sin(typename TYPE::value_type(n * 2 * M_PI / N)) * a;
}

template <int n, int N, typename TYPE>
static inline TYPE ex(TYPE a)
{
	return exp(TYPE(0, n * 2 * M_PI / N)) * a;
}

template <typename TYPE>
static inline TYPE fiddle(TYPE a, TYPE b)
{
	TYPE c(a + b), d(a - b);
	return TYPE(d.real() + c.imag(), d.imag() - c.real());
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
	return
		(!(N % 29)) ? 29 :
		(!(N % 23)) ? 23 :
		(!(N % 19)) ? 19 :
		(!(N % 17)) ? 17 :
		(!(N % 13)) ? 13 :
		(!(N % 11)) ? 11 :
		(!(N % 7)) ? 7 :
		(!(N % 5)) ? 5 :
		(!(N % 3)) ? 3 :
		(!(N % 8) && pow8(N)) ? 8 :
		(!(N % 4) && pow4(N)) ? 4 :
		(!(N % 2)) ? 2 :
		1;
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

template <int STRIDE, typename TYPE>
struct Dit<3, 3, STRIDE, TYPE, -1>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE *out2, TYPE in0, TYPE in1, TYPE in2)
	{
		Dit<3, 3, STRIDE, TYPE, 1>::dft(out0, out2, out1, in0, in1, in2);
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
		TYPE a(in1 + in2), b(sqrt3(twiddle(in1, in2)));
		*out0 = in0 + a;
		*out1 = in0 - half(a + b);
		*out2 = in0 - half(a - b);
	}
	static inline void dit(TYPE *out, const TYPE *in, const TYPE *)
	{
		dft(out, out + 1, out + 2, in[0], in[STRIDE], in[2 * STRIDE]);
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

template <int STRIDE, typename TYPE>
struct Dit<4, 4, STRIDE, TYPE, -1>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE *out2, TYPE *out3,
			TYPE in0, TYPE in1, TYPE in2, TYPE in3)
	{
		Dit<4, 4, STRIDE, TYPE, 1>::dft(out0, out3, out2, out1,
			in0, in1, in2, in3);
	}
	static inline void dit(TYPE *out, const TYPE *in, const TYPE *)
	{
		dft(out, out + 1, out + 2, out + 3,
			in[0], in[STRIDE], in[2 * STRIDE], in[3 * STRIDE]);
	}
};

template <int STRIDE, typename TYPE>
struct Dit<4, 4, STRIDE, TYPE, 1>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE *out2, TYPE *out3,
			TYPE in0, TYPE in1, TYPE in2, TYPE in3)
	{
		TYPE a(in0 + in2), b(in0 - in2);
		TYPE c(in1 + in3), d(twiddle(in1, in3));
		*out0 = a + c;
		*out1 = b - d;
		*out2 = a - c;
		*out3 = b + d;
	}
	static inline void dit(TYPE *out, const TYPE *in, const TYPE *)
	{
		dft(out, out + 1, out + 2, out + 3,
			in[0], in[STRIDE], in[2 * STRIDE], in[3 * STRIDE]);
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

template <int STRIDE, typename TYPE>
struct Dit<5, 5, STRIDE, TYPE, -1>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE *out2, TYPE *out3, TYPE *out4,
			TYPE in0, TYPE in1, TYPE in2, TYPE in3, TYPE in4)
	{
		Dit<5, 5, STRIDE, TYPE, 1>::dft(out0, out4, out3, out2, out1,
			in0, in1, in2, in3, in4);
	}
	static inline void dit(TYPE *out, const TYPE *in, const TYPE *)
	{
		dft(out, out + 1, out + 2, out + 3, out + 4,
			in[0], in[STRIDE], in[2 * STRIDE], in[3 * STRIDE], in[4 * STRIDE]);
	}
};

template <int STRIDE, typename TYPE>
struct Dit<5, 5, STRIDE, TYPE, 1>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE *out2, TYPE *out3, TYPE *out4,
			TYPE in0, TYPE in1, TYPE in2, TYPE in3, TYPE in4)
	{
		TYPE a1(in1 + in4), t1(twiddle(in1, in4));
		TYPE a2(in2 + in3), t2(twiddle(in2, in3));
		TYPE c1(cx<1,5>(a1) + cx<2,5>(a2));
		TYPE s1(sx<1,5>(t1) + sx<2,5>(t2));
		TYPE c2(cx<2,5>(a1) + cx<1,5>(a2));
		TYPE s2(sx<2,5>(t1) - sx<1,5>(t2));
		*out0 = in0 + a1 + a2;
		*out1 = in0 + c1 - s1;
		*out2 = in0 + c2 - s2;
		*out3 = in0 + c2 + s2;
		*out4 = in0 + c1 + s1;
	}
	static inline void dit(TYPE *out, const TYPE *in, const TYPE *)
	{
		dft(out, out + 1, out + 2, out + 3, out + 4,
			in[0], in[STRIDE], in[2 * STRIDE], in[3 * STRIDE], in[4 * STRIDE]);
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

template <int STRIDE, typename TYPE>
struct Dit<7, 7, STRIDE, TYPE, -1>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE *out2, TYPE *out3, TYPE *out4, TYPE *out5, TYPE *out6,
			TYPE in0, TYPE in1, TYPE in2, TYPE in3, TYPE in4, TYPE in5, TYPE in6)
	{
		Dit<7, 7, STRIDE, TYPE, 1>::dft(out0, out6, out5, out4, out3, out2, out1,
			in0, in1, in2, in3, in4, in5, in6);
	}
	static inline void dit(TYPE *out, const TYPE *in, const TYPE *)
	{
		dft(out, out + 1, out + 2, out + 3, out + 4, out + 5, out + 6,
			in[0], in[STRIDE], in[2 * STRIDE], in[3 * STRIDE], in[4 * STRIDE], in[5 * STRIDE], in[6 * STRIDE]);
	}
};

template <int STRIDE, typename TYPE>
struct Dit<7, 7, STRIDE, TYPE, 1>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE *out2, TYPE *out3, TYPE *out4, TYPE *out5, TYPE *out6,
			TYPE in0, TYPE in1, TYPE in2, TYPE in3, TYPE in4, TYPE in5, TYPE in6)
	{
		TYPE a1(in1 + in6), t1(twiddle(in1, in6));
		TYPE a2(in2 + in5), t2(twiddle(in2, in5));
		TYPE a3(in3 + in4), t3(twiddle(in3, in4));
		TYPE c1(cx<1,7>(a1) + cx<2,7>(a2) + cx<3,7>(a3));
		TYPE s1(sx<1,7>(t1) + sx<2,7>(t2) + sx<3,7>(t3));
		TYPE c2(cx<2,7>(a1) + cx<3,7>(a2) + cx<1,7>(a3));
		TYPE s2(sx<2,7>(t1) - sx<3,7>(t2) - sx<1,7>(t3));
		TYPE c3(cx<3,7>(a1) + cx<1,7>(a2) + cx<2,7>(a3));
		TYPE s3(sx<3,7>(t1) - sx<1,7>(t2) + sx<2,7>(t3));
		*out0 = in0 + a1 + a2 + a3;
		*out1 = in0 + c1 - s1;
		*out2 = in0 + c2 - s2;
		*out3 = in0 + c3 - s3;
		*out4 = in0 + c3 + s3;
		*out5 = in0 + c2 + s2;
		*out6 = in0 + c1 + s1;
	}
	static inline void dit(TYPE *out, const TYPE *in, const TYPE *)
	{
		dft(out, out + 1, out + 2, out + 3, out + 4, out + 5, out + 6,
			in[0], in[STRIDE], in[2 * STRIDE], in[3 * STRIDE], in[4 * STRIDE], in[5 * STRIDE], in[6 * STRIDE]);
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

template <int STRIDE, typename TYPE>
struct Dit<8, 8, STRIDE, TYPE, -1>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE *out2, TYPE *out3, TYPE *out4, TYPE *out5, TYPE *out6, TYPE *out7,
			TYPE in0, TYPE in1, TYPE in2, TYPE in3, TYPE in4, TYPE in5, TYPE in6, TYPE in7)
	{
		Dit<8, 8, STRIDE, TYPE, 1>::dft(out0, out7, out6, out5, out4, out3, out2, out1,
			in0, in1, in2, in3, in4, in5, in6, in7);
	}
	static inline void dit(TYPE *out, const TYPE *in, const TYPE *)
	{
		dft(out, out + 1, out + 2, out + 3, out + 4, out + 5, out + 6, out + 7,
			in[0], in[STRIDE], in[2 * STRIDE], in[3 * STRIDE], in[4 * STRIDE], in[5 * STRIDE], in[6 * STRIDE], in[7 * STRIDE]);
	}
};

template <int STRIDE, typename TYPE>
struct Dit<8, 8, STRIDE, TYPE, 1>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE *out2, TYPE *out3, TYPE *out4, TYPE *out5, TYPE *out6, TYPE *out7,
			TYPE in0, TYPE in1, TYPE in2, TYPE in3, TYPE in4, TYPE in5, TYPE in6, TYPE in7)
	{
		TYPE a(in0 + in4), b(in0 - in4);
		TYPE c(in1 + in5), d(in1 - in5);
		TYPE e(in2 + in6), f(twiddle(in2, in6));
		TYPE g(in3 + in7), h(in3 - in7);
		TYPE cpg(c + g), tcg(twiddle(c, g));
		TYPE fdh(rsqrt2(fiddle(d, h))), fhd(rsqrt2(fiddle(h, d)));
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
		dft(out, out + 1, out + 2, out + 3, out + 4, out + 5, out + 6, out + 7,
			in[0], in[STRIDE], in[2 * STRIDE], in[3 * STRIDE], in[4 * STRIDE], in[5 * STRIDE], in[6 * STRIDE], in[7 * STRIDE]);
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

template <int STRIDE, typename TYPE>
struct Dit<11, 11, STRIDE, TYPE, -1>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE *out2, TYPE *out3, TYPE *out4, TYPE *out5, TYPE *out6, TYPE *out7, TYPE *out8, TYPE *out9, TYPE *out10,
			TYPE in0, TYPE in1, TYPE in2, TYPE in3, TYPE in4, TYPE in5, TYPE in6, TYPE in7, TYPE in8, TYPE in9, TYPE in10)
	{
		Dit<11, 11, STRIDE, TYPE, 1>::dft(out0, out10, out9, out8, out7, out6, out5, out4, out3, out2, out1,
			in0, in1, in2, in3, in4, in5, in6, in7, in8, in9, in10);
	}
	static inline void dit(TYPE *out, const TYPE *in, const TYPE *)
	{
		dft(out, out + 1, out + 2, out + 3, out + 4, out + 5, out + 6, out + 7, out + 8, out + 9, out + 10,
			in[0], in[STRIDE], in[2 * STRIDE], in[3 * STRIDE], in[4 * STRIDE], in[5 * STRIDE], in[6 * STRIDE], in[7 * STRIDE], in[8 * STRIDE], in[9 * STRIDE], in[10 * STRIDE]);
	}
};

template <int STRIDE, typename TYPE>
struct Dit<11, 11, STRIDE, TYPE, 1>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE *out2, TYPE *out3, TYPE *out4, TYPE *out5, TYPE *out6, TYPE *out7, TYPE *out8, TYPE *out9, TYPE *out10,
			TYPE in0, TYPE in1, TYPE in2, TYPE in3, TYPE in4, TYPE in5, TYPE in6, TYPE in7, TYPE in8, TYPE in9, TYPE in10)
	{
		TYPE a1(in1 + in10), t1(twiddle(in1, in10));
		TYPE a2(in2 + in9), t2(twiddle(in2, in9));
		TYPE a3(in3 + in8), t3(twiddle(in3, in8));
		TYPE a4(in4 + in7), t4(twiddle(in4, in7));
		TYPE a5(in5 + in6), t5(twiddle(in5, in6));
		TYPE c1(cx<1,11>(a1) + cx<2,11>(a2) + cx<3,11>(a3) + cx<4,11>(a4) + cx<5,11>(a5));
		TYPE s1(sx<1,11>(t1) + sx<2,11>(t2) + sx<3,11>(t3) + sx<4,11>(t4) + sx<5,11>(t5));
		TYPE c2(cx<2,11>(a1) + cx<4,11>(a2) + cx<5,11>(a3) + cx<3,11>(a4) + cx<1,11>(a5));
		TYPE s2(sx<2,11>(t1) + sx<4,11>(t2) - sx<5,11>(t3) - sx<3,11>(t4) - sx<1,11>(t5));
		TYPE c3(cx<3,11>(a1) + cx<5,11>(a2) + cx<2,11>(a3) + cx<1,11>(a4) + cx<4,11>(a5));
		TYPE s3(sx<3,11>(t1) - sx<5,11>(t2) - sx<2,11>(t3) + sx<1,11>(t4) + sx<4,11>(t5));
		TYPE c4(cx<4,11>(a1) + cx<3,11>(a2) + cx<1,11>(a3) + cx<5,11>(a4) + cx<2,11>(a5));
		TYPE s4(sx<4,11>(t1) - sx<3,11>(t2) + sx<1,11>(t3) + sx<5,11>(t4) - sx<2,11>(t5));
		TYPE c5(cx<5,11>(a1) + cx<1,11>(a2) + cx<4,11>(a3) + cx<2,11>(a4) + cx<3,11>(a5));
		TYPE s5(sx<5,11>(t1) - sx<1,11>(t2) + sx<4,11>(t3) - sx<2,11>(t4) + sx<3,11>(t5));
		*out0 = in0 + a1 + a2 + a3 + a4 + a5;
		*out1 = in0 + c1 - s1;
		*out2 = in0 + c2 - s2;
		*out3 = in0 + c3 - s3;
		*out4 = in0 + c4 - s4;
		*out5 = in0 + c5 - s5;
		*out6 = in0 + c5 + s5;
		*out7 = in0 + c4 + s4;
		*out8 = in0 + c3 + s3;
		*out9 = in0 + c2 + s2;
		*out10 = in0 + c1 + s1;
	}
	static inline void dit(TYPE *out, const TYPE *in, const TYPE *)
	{
		dft(out, out + 1, out + 2, out + 3, out + 4, out + 5, out + 6, out + 7, out + 8, out + 9, out + 10,
			in[0], in[STRIDE], in[2 * STRIDE], in[3 * STRIDE], in[4 * STRIDE], in[5 * STRIDE], in[6 * STRIDE], in[7 * STRIDE], in[8 * STRIDE], in[9 * STRIDE], in[10 * STRIDE]);
	}
};

template <int BINS, int STRIDE, typename TYPE, int SIGN>
struct Dit<11, BINS, STRIDE, TYPE, SIGN>
{
	static const int RADIX = 11;
	static const int QUOTIENT = BINS / RADIX;
	static void dit(TYPE *out, const TYPE *in, const TYPE *z)
	{
		for (int o = 0, i = 0; o < BINS; o += QUOTIENT, i += STRIDE)
			Dit<split(QUOTIENT), QUOTIENT, RADIX * STRIDE, TYPE, SIGN>::dit(out + o, in + i, z);
		for (int k0 = 0, k1 = QUOTIENT, k2 = 2 * QUOTIENT, k3 = 3 * QUOTIENT, k4 = 4 * QUOTIENT, k5 = 5 * QUOTIENT, k6 = 6 * QUOTIENT, k7 = 7 * QUOTIENT, k8 = 8 * QUOTIENT, k9 = 9 * QUOTIENT, k10 = 10 * QUOTIENT,
				l1 = 0, l2 = 0, l3 = 0, l4 = 0, l5 = 0, l6 = 0, l7 = 0, l8 = 0, l9 = 0, l10 = 0;
				k0 < QUOTIENT;
				++k0, ++k1, ++k2, ++k3, ++k4, ++k5, ++k6, ++k7, ++k8, ++k9, ++k10,
				l1 += STRIDE, l2 += 2 * STRIDE, l3 += 3 * STRIDE, l4 += 4 * STRIDE, l5 += 5 * STRIDE, l6 += 6 * STRIDE, l7 += 7 * STRIDE, l8 += 8 * STRIDE, l9 += 9 * STRIDE, l10 += 10 * STRIDE)
			Dit<RADIX, RADIX, STRIDE, TYPE, SIGN>::dft(out + k0, out + k1, out + k2, out + k3, out + k4, out + k5, out + k6, out + k7, out + k8, out + k9, out + k10,
				out[k0], z[l1] * out[k1], z[l2] * out[k2], z[l3] * out[k3], z[l4] * out[k4], z[l5] * out[k5], z[l6] * out[k6], z[l7] * out[k7], z[l8] * out[k8], z[l9] * out[k9], z[l10] * out[k10]);
	}
};

template <int STRIDE, typename TYPE>
struct Dit<13, 13, STRIDE, TYPE, -1>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE *out2, TYPE *out3, TYPE *out4, TYPE *out5, TYPE *out6, TYPE *out7, TYPE *out8, TYPE *out9, TYPE *out10, TYPE *out11, TYPE *out12,
			TYPE in0, TYPE in1, TYPE in2, TYPE in3, TYPE in4, TYPE in5, TYPE in6, TYPE in7, TYPE in8, TYPE in9, TYPE in10, TYPE in11, TYPE in12)
	{
		Dit<13, 13, STRIDE, TYPE, 1>::dft(out0, out12, out11, out10, out9, out8, out7, out6, out5, out4, out3, out2, out1,
			in0, in1, in2, in3, in4, in5, in6, in7, in8, in9, in10, in11, in12);
	}
	static inline void dit(TYPE *out, const TYPE *in, const TYPE *)
	{
		dft(out, out + 1, out + 2, out + 3, out + 4, out + 5, out + 6, out + 7, out + 8, out + 9, out + 10, out + 11, out + 12,
			in[0], in[STRIDE], in[2 * STRIDE], in[3 * STRIDE], in[4 * STRIDE], in[5 * STRIDE], in[6 * STRIDE], in[7 * STRIDE], in[8 * STRIDE], in[9 * STRIDE], in[10 * STRIDE], in[11 * STRIDE], in[12 * STRIDE]);
	}
};

template <int STRIDE, typename TYPE>
struct Dit<13, 13, STRIDE, TYPE, 1>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE *out2, TYPE *out3, TYPE *out4, TYPE *out5, TYPE *out6, TYPE *out7, TYPE *out8, TYPE *out9, TYPE *out10, TYPE *out11, TYPE *out12,
			TYPE in0, TYPE in1, TYPE in2, TYPE in3, TYPE in4, TYPE in5, TYPE in6, TYPE in7, TYPE in8, TYPE in9, TYPE in10, TYPE in11, TYPE in12)
	{
		TYPE a1(in1 + in12), t1(twiddle(in1, in12));
		TYPE a2(in2 + in11), t2(twiddle(in2, in11));
		TYPE a3(in3 + in10), t3(twiddle(in3, in10));
		TYPE a4(in4 + in9), t4(twiddle(in4, in9));
		TYPE a5(in5 + in8), t5(twiddle(in5, in8));
		TYPE a6(in6 + in7), t6(twiddle(in6, in7));
		TYPE c1(cx<1,13>(a1) + cx<2,13>(a2) + cx<3,13>(a3) + cx<4,13>(a4) + cx<5,13>(a5) + cx<6,13>(a6));
		TYPE s1(sx<1,13>(t1) + sx<2,13>(t2) + sx<3,13>(t3) + sx<4,13>(t4) + sx<5,13>(t5) + sx<6,13>(t6));
		TYPE c2(cx<2,13>(a1) + cx<4,13>(a2) + cx<6,13>(a3) + cx<5,13>(a4) + cx<3,13>(a5) + cx<1,13>(a6));
		TYPE s2(sx<2,13>(t1) + sx<4,13>(t2) + sx<6,13>(t3) - sx<5,13>(t4) - sx<3,13>(t5) - sx<1,13>(t6));
		TYPE c3(cx<3,13>(a1) + cx<6,13>(a2) + cx<4,13>(a3) + cx<1,13>(a4) + cx<2,13>(a5) + cx<5,13>(a6));
		TYPE s3(sx<3,13>(t1) + sx<6,13>(t2) - sx<4,13>(t3) - sx<1,13>(t4) + sx<2,13>(t5) + sx<5,13>(t6));
		TYPE c4(cx<4,13>(a1) + cx<5,13>(a2) + cx<1,13>(a3) + cx<3,13>(a4) + cx<6,13>(a5) + cx<2,13>(a6));
		TYPE s4(sx<4,13>(t1) - sx<5,13>(t2) - sx<1,13>(t3) + sx<3,13>(t4) - sx<6,13>(t5) - sx<2,13>(t6));
		TYPE c5(cx<5,13>(a1) + cx<3,13>(a2) + cx<2,13>(a3) + cx<6,13>(a4) + cx<1,13>(a5) + cx<4,13>(a6));
		TYPE s5(sx<5,13>(t1) - sx<3,13>(t2) + sx<2,13>(t3) - sx<6,13>(t4) - sx<1,13>(t5) + sx<4,13>(t6));
		TYPE c6(cx<6,13>(a1) + cx<1,13>(a2) + cx<5,13>(a3) + cx<2,13>(a4) + cx<4,13>(a5) + cx<3,13>(a6));
		TYPE s6(sx<6,13>(t1) - sx<1,13>(t2) + sx<5,13>(t3) - sx<2,13>(t4) + sx<4,13>(t5) - sx<3,13>(t6));
		*out0 = in0 + a1 + a2 + a3 + a4 + a5 + a6;
		*out1 = in0 + c1 - s1;
		*out2 = in0 + c2 - s2;
		*out3 = in0 + c3 - s3;
		*out4 = in0 + c4 - s4;
		*out5 = in0 + c5 - s5;
		*out6 = in0 + c6 - s6;
		*out7 = in0 + c6 + s6;
		*out8 = in0 + c5 + s5;
		*out9 = in0 + c4 + s4;
		*out10 = in0 + c3 + s3;
		*out11 = in0 + c2 + s2;
		*out12 = in0 + c1 + s1;
	}
	static inline void dit(TYPE *out, const TYPE *in, const TYPE *)
	{
		dft(out, out + 1, out + 2, out + 3, out + 4, out + 5, out + 6, out + 7, out + 8, out + 9, out + 10, out + 11, out + 12,
			in[0], in[STRIDE], in[2 * STRIDE], in[3 * STRIDE], in[4 * STRIDE], in[5 * STRIDE], in[6 * STRIDE], in[7 * STRIDE], in[8 * STRIDE], in[9 * STRIDE], in[10 * STRIDE], in[11 * STRIDE], in[12 * STRIDE]);
	}
};

template <int BINS, int STRIDE, typename TYPE, int SIGN>
struct Dit<13, BINS, STRIDE, TYPE, SIGN>
{
	static const int RADIX = 13;
	static const int QUOTIENT = BINS / RADIX;
	static void dit(TYPE *out, const TYPE *in, const TYPE *z)
	{
		for (int o = 0, i = 0; o < BINS; o += QUOTIENT, i += STRIDE)
			Dit<split(QUOTIENT), QUOTIENT, RADIX * STRIDE, TYPE, SIGN>::dit(out + o, in + i, z);
		for (int k0 = 0, k1 = QUOTIENT, k2 = 2 * QUOTIENT, k3 = 3 * QUOTIENT, k4 = 4 * QUOTIENT, k5 = 5 * QUOTIENT, k6 = 6 * QUOTIENT, k7 = 7 * QUOTIENT, k8 = 8 * QUOTIENT, k9 = 9 * QUOTIENT, k10 = 10 * QUOTIENT, k11 = 11 * QUOTIENT, k12 = 12 * QUOTIENT,
				l1 = 0, l2 = 0, l3 = 0, l4 = 0, l5 = 0, l6 = 0, l7 = 0, l8 = 0, l9 = 0, l10 = 0, l11 = 0, l12 = 0;
				k0 < QUOTIENT;
				++k0, ++k1, ++k2, ++k3, ++k4, ++k5, ++k6, ++k7, ++k8, ++k9, ++k10, ++k11, ++k12,
				l1 += STRIDE, l2 += 2 * STRIDE, l3 += 3 * STRIDE, l4 += 4 * STRIDE, l5 += 5 * STRIDE, l6 += 6 * STRIDE, l7 += 7 * STRIDE, l8 += 8 * STRIDE, l9 += 9 * STRIDE, l10 += 10 * STRIDE, l11 += 11 * STRIDE, l12 += 12 * STRIDE)
			Dit<RADIX, RADIX, STRIDE, TYPE, SIGN>::dft(out + k0, out + k1, out + k2, out + k3, out + k4, out + k5, out + k6, out + k7, out + k8, out + k9, out + k10, out + k11, out + k12,
				out[k0], z[l1] * out[k1], z[l2] * out[k2], z[l3] * out[k3], z[l4] * out[k4], z[l5] * out[k5], z[l6] * out[k6], z[l7] * out[k7], z[l8] * out[k8], z[l9] * out[k9], z[l10] * out[k10], z[l11] * out[k11], z[l12] * out[k12]);
	}
};

template <int STRIDE, typename TYPE>
struct Dit<17, 17, STRIDE, TYPE, -1>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE *out2, TYPE *out3, TYPE *out4, TYPE *out5, TYPE *out6, TYPE *out7, TYPE *out8, TYPE *out9, TYPE *out10, TYPE *out11, TYPE *out12, TYPE *out13, TYPE *out14, TYPE *out15, TYPE *out16,
			TYPE in0, TYPE in1, TYPE in2, TYPE in3, TYPE in4, TYPE in5, TYPE in6, TYPE in7, TYPE in8, TYPE in9, TYPE in10, TYPE in11, TYPE in12, TYPE in13, TYPE in14, TYPE in15, TYPE in16)
	{
		Dit<17, 17, STRIDE, TYPE, 1>::dft(out0, out16, out15, out14, out13, out12, out11, out10, out9, out8, out7, out6, out5, out4, out3, out2, out1,
			in0, in1, in2, in3, in4, in5, in6, in7, in8, in9, in10, in11, in12, in13, in14, in15, in16);
	}
	static inline void dit(TYPE *out, const TYPE *in, const TYPE *)
	{
		dft(out, out + 1, out + 2, out + 3, out + 4, out + 5, out + 6, out + 7, out + 8, out + 9, out + 10, out + 11, out + 12, out + 13, out + 14, out + 15, out + 16,
			in[0], in[STRIDE], in[2 * STRIDE], in[3 * STRIDE], in[4 * STRIDE], in[5 * STRIDE], in[6 * STRIDE], in[7 * STRIDE], in[8 * STRIDE], in[9 * STRIDE], in[10 * STRIDE], in[11 * STRIDE], in[12 * STRIDE], in[13 * STRIDE], in[14 * STRIDE], in[15 * STRIDE], in[16 * STRIDE]);
	}
};

template <int STRIDE, typename TYPE>
struct Dit<17, 17, STRIDE, TYPE, 1>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE *out2, TYPE *out3, TYPE *out4, TYPE *out5, TYPE *out6, TYPE *out7, TYPE *out8, TYPE *out9, TYPE *out10, TYPE *out11, TYPE *out12, TYPE *out13, TYPE *out14, TYPE *out15, TYPE *out16,
			TYPE in0, TYPE in1, TYPE in2, TYPE in3, TYPE in4, TYPE in5, TYPE in6, TYPE in7, TYPE in8, TYPE in9, TYPE in10, TYPE in11, TYPE in12, TYPE in13, TYPE in14, TYPE in15, TYPE in16)
	{
		TYPE a1(in1 + in16), t1(twiddle(in1, in16));
		TYPE a2(in2 + in15), t2(twiddle(in2, in15));
		TYPE a3(in3 + in14), t3(twiddle(in3, in14));
		TYPE a4(in4 + in13), t4(twiddle(in4, in13));
		TYPE a5(in5 + in12), t5(twiddle(in5, in12));
		TYPE a6(in6 + in11), t6(twiddle(in6, in11));
		TYPE a7(in7 + in10), t7(twiddle(in7, in10));
		TYPE a8(in8 + in9), t8(twiddle(in8, in9));
		TYPE c1(cx<1,17>(a1) + cx<2,17>(a2) + cx<3,17>(a3) + cx<4,17>(a4) + cx<5,17>(a5) + cx<6,17>(a6) + cx<7,17>(a7) + cx<8,17>(a8));
		TYPE s1(sx<1,17>(t1) + sx<2,17>(t2) + sx<3,17>(t3) + sx<4,17>(t4) + sx<5,17>(t5) + sx<6,17>(t6) + sx<7,17>(t7) + sx<8,17>(t8));
		TYPE c2(cx<2,17>(a1) + cx<4,17>(a2) + cx<6,17>(a3) + cx<8,17>(a4) + cx<7,17>(a5) + cx<5,17>(a6) + cx<3,17>(a7) + cx<1,17>(a8));
		TYPE s2(sx<2,17>(t1) + sx<4,17>(t2) + sx<6,17>(t3) + sx<8,17>(t4) - sx<7,17>(t5) - sx<5,17>(t6) - sx<3,17>(t7) - sx<1,17>(t8));
		TYPE c3(cx<3,17>(a1) + cx<6,17>(a2) + cx<8,17>(a3) + cx<5,17>(a4) + cx<2,17>(a5) + cx<1,17>(a6) + cx<4,17>(a7) + cx<7,17>(a8));
		TYPE s3(sx<3,17>(t1) + sx<6,17>(t2) - sx<8,17>(t3) - sx<5,17>(t4) - sx<2,17>(t5) + sx<1,17>(t6) + sx<4,17>(t7) + sx<7,17>(t8));
		TYPE c4(cx<4,17>(a1) + cx<8,17>(a2) + cx<5,17>(a3) + cx<1,17>(a4) + cx<3,17>(a5) + cx<7,17>(a6) + cx<6,17>(a7) + cx<2,17>(a8));
		TYPE s4(sx<4,17>(t1) + sx<8,17>(t2) - sx<5,17>(t3) - sx<1,17>(t4) + sx<3,17>(t5) + sx<7,17>(t6) - sx<6,17>(t7) - sx<2,17>(t8));
		TYPE c5(cx<5,17>(a1) + cx<7,17>(a2) + cx<2,17>(a3) + cx<3,17>(a4) + cx<8,17>(a5) + cx<4,17>(a6) + cx<1,17>(a7) + cx<6,17>(a8));
		TYPE s5(sx<5,17>(t1) - sx<7,17>(t2) - sx<2,17>(t3) + sx<3,17>(t4) + sx<8,17>(t5) - sx<4,17>(t6) + sx<1,17>(t7) + sx<6,17>(t8));
		TYPE c6(cx<6,17>(a1) + cx<5,17>(a2) + cx<1,17>(a3) + cx<7,17>(a4) + cx<4,17>(a5) + cx<2,17>(a6) + cx<8,17>(a7) + cx<3,17>(a8));
		TYPE s6(sx<6,17>(t1) - sx<5,17>(t2) + sx<1,17>(t3) + sx<7,17>(t4) - sx<4,17>(t5) + sx<2,17>(t6) + sx<8,17>(t7) - sx<3,17>(t8));
		TYPE c7(cx<7,17>(a1) + cx<3,17>(a2) + cx<4,17>(a3) + cx<6,17>(a4) + cx<1,17>(a5) + cx<8,17>(a6) + cx<2,17>(a7) + cx<5,17>(a8));
		TYPE s7(sx<7,17>(t1) - sx<3,17>(t2) + sx<4,17>(t3) - sx<6,17>(t4) + sx<1,17>(t5) + sx<8,17>(t6) - sx<2,17>(t7) + sx<5,17>(t8));
		TYPE c8(cx<8,17>(a1) + cx<1,17>(a2) + cx<7,17>(a3) + cx<2,17>(a4) + cx<6,17>(a5) + cx<3,17>(a6) + cx<5,17>(a7) + cx<4,17>(a8));
		TYPE s8(sx<8,17>(t1) - sx<1,17>(t2) + sx<7,17>(t3) - sx<2,17>(t4) + sx<6,17>(t5) - sx<3,17>(t6) + sx<5,17>(t7) - sx<4,17>(t8));
		*out0 = in0 + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8;
		*out1 = in0 + c1 - s1;
		*out2 = in0 + c2 - s2;
		*out3 = in0 + c3 - s3;
		*out4 = in0 + c4 - s4;
		*out5 = in0 + c5 - s5;
		*out6 = in0 + c6 - s6;
		*out7 = in0 + c7 - s7;
		*out8 = in0 + c8 - s8;
		*out9 = in0 + c8 + s8;
		*out10 = in0 + c7 + s7;
		*out11 = in0 + c6 + s6;
		*out12 = in0 + c5 + s5;
		*out13 = in0 + c4 + s4;
		*out14 = in0 + c3 + s3;
		*out15 = in0 + c2 + s2;
		*out16 = in0 + c1 + s1;
	}
	static inline void dit(TYPE *out, const TYPE *in, const TYPE *)
	{
		dft(out, out + 1, out + 2, out + 3, out + 4, out + 5, out + 6, out + 7, out + 8, out + 9, out + 10, out + 11, out + 12, out + 13, out + 14, out + 15, out + 16,
			in[0], in[STRIDE], in[2 * STRIDE], in[3 * STRIDE], in[4 * STRIDE], in[5 * STRIDE], in[6 * STRIDE], in[7 * STRIDE], in[8 * STRIDE], in[9 * STRIDE], in[10 * STRIDE], in[11 * STRIDE], in[12 * STRIDE], in[13 * STRIDE], in[14 * STRIDE], in[15 * STRIDE], in[16 * STRIDE]);
	}
};

template <int BINS, int STRIDE, typename TYPE, int SIGN>
struct Dit<17, BINS, STRIDE, TYPE, SIGN>
{
	static const int RADIX = 17;
	static const int QUOTIENT = BINS / RADIX;
	static void dit(TYPE *out, const TYPE *in, const TYPE *z)
	{
		for (int o = 0, i = 0; o < BINS; o += QUOTIENT, i += STRIDE)
			Dit<split(QUOTIENT), QUOTIENT, RADIX * STRIDE, TYPE, SIGN>::dit(out + o, in + i, z);
		for (int k0 = 0, k1 = QUOTIENT, k2 = 2 * QUOTIENT, k3 = 3 * QUOTIENT, k4 = 4 * QUOTIENT, k5 = 5 * QUOTIENT, k6 = 6 * QUOTIENT, k7 = 7 * QUOTIENT, k8 = 8 * QUOTIENT, k9 = 9 * QUOTIENT, k10 = 10 * QUOTIENT, k11 = 11 * QUOTIENT, k12 = 12 * QUOTIENT, k13 = 13 * QUOTIENT, k14 = 14 * QUOTIENT, k15 = 15 * QUOTIENT, k16 = 16 * QUOTIENT,
				l1 = 0, l2 = 0, l3 = 0, l4 = 0, l5 = 0, l6 = 0, l7 = 0, l8 = 0, l9 = 0, l10 = 0, l11 = 0, l12 = 0, l13 = 0, l14 = 0, l15 = 0, l16 = 0;
				k0 < QUOTIENT;
				++k0, ++k1, ++k2, ++k3, ++k4, ++k5, ++k6, ++k7, ++k8, ++k9, ++k10, ++k11, ++k12, ++k13, ++k14, ++k15, ++k16,
				l1 += STRIDE, l2 += 2 * STRIDE, l3 += 3 * STRIDE, l4 += 4 * STRIDE, l5 += 5 * STRIDE, l6 += 6 * STRIDE, l7 += 7 * STRIDE, l8 += 8 * STRIDE, l9 += 9 * STRIDE, l10 += 10 * STRIDE, l11 += 11 * STRIDE, l12 += 12 * STRIDE, l13 += 13 * STRIDE, l14 += 14 * STRIDE, l15 += 15 * STRIDE, l16 += 16 * STRIDE)
			Dit<RADIX, RADIX, STRIDE, TYPE, SIGN>::dft(out + k0, out + k1, out + k2, out + k3, out + k4, out + k5, out + k6, out + k7, out + k8, out + k9, out + k10, out + k11, out + k12, out + k13, out + k14, out + k15, out + k16,
				out[k0], z[l1] * out[k1], z[l2] * out[k2], z[l3] * out[k3], z[l4] * out[k4], z[l5] * out[k5], z[l6] * out[k6], z[l7] * out[k7], z[l8] * out[k8], z[l9] * out[k9], z[l10] * out[k10], z[l11] * out[k11], z[l12] * out[k12], z[l13] * out[k13], z[l14] * out[k14], z[l15] * out[k15], z[l16] * out[k16]);
	}
};

template <int STRIDE, typename TYPE>
struct Dit<19, 19, STRIDE, TYPE, -1>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE *out2, TYPE *out3, TYPE *out4, TYPE *out5, TYPE *out6, TYPE *out7, TYPE *out8, TYPE *out9, TYPE *out10, TYPE *out11, TYPE *out12, TYPE *out13, TYPE *out14, TYPE *out15, TYPE *out16, TYPE *out17, TYPE *out18,
			TYPE in0, TYPE in1, TYPE in2, TYPE in3, TYPE in4, TYPE in5, TYPE in6, TYPE in7, TYPE in8, TYPE in9, TYPE in10, TYPE in11, TYPE in12, TYPE in13, TYPE in14, TYPE in15, TYPE in16, TYPE in17, TYPE in18)
	{
		Dit<19, 19, STRIDE, TYPE, 1>::dft(out0, out18, out17, out16, out15, out14, out13, out12, out11, out10, out9, out8, out7, out6, out5, out4, out3, out2, out1,
			in0, in1, in2, in3, in4, in5, in6, in7, in8, in9, in10, in11, in12, in13, in14, in15, in16, in17, in18);
	}
	static inline void dit(TYPE *out, const TYPE *in, const TYPE *)
	{
		dft(out, out + 1, out + 2, out + 3, out + 4, out + 5, out + 6, out + 7, out + 8, out + 9, out + 10, out + 11, out + 12, out + 13, out + 14, out + 15, out + 16, out + 17, out + 18,
			in[0], in[STRIDE], in[2 * STRIDE], in[3 * STRIDE], in[4 * STRIDE], in[5 * STRIDE], in[6 * STRIDE], in[7 * STRIDE], in[8 * STRIDE], in[9 * STRIDE], in[10 * STRIDE], in[11 * STRIDE], in[12 * STRIDE], in[13 * STRIDE], in[14 * STRIDE], in[15 * STRIDE], in[16 * STRIDE], in[17 * STRIDE], in[18 * STRIDE]);
	}
};

template <int STRIDE, typename TYPE>
struct Dit<19, 19, STRIDE, TYPE, 1>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE *out2, TYPE *out3, TYPE *out4, TYPE *out5, TYPE *out6, TYPE *out7, TYPE *out8, TYPE *out9, TYPE *out10, TYPE *out11, TYPE *out12, TYPE *out13, TYPE *out14, TYPE *out15, TYPE *out16, TYPE *out17, TYPE *out18,
			TYPE in0, TYPE in1, TYPE in2, TYPE in3, TYPE in4, TYPE in5, TYPE in6, TYPE in7, TYPE in8, TYPE in9, TYPE in10, TYPE in11, TYPE in12, TYPE in13, TYPE in14, TYPE in15, TYPE in16, TYPE in17, TYPE in18)
	{
		TYPE a1(in1 + in18), t1(twiddle(in1, in18));
		TYPE a2(in2 + in17), t2(twiddle(in2, in17));
		TYPE a3(in3 + in16), t3(twiddle(in3, in16));
		TYPE a4(in4 + in15), t4(twiddle(in4, in15));
		TYPE a5(in5 + in14), t5(twiddle(in5, in14));
		TYPE a6(in6 + in13), t6(twiddle(in6, in13));
		TYPE a7(in7 + in12), t7(twiddle(in7, in12));
		TYPE a8(in8 + in11), t8(twiddle(in8, in11));
		TYPE a9(in9 + in10), t9(twiddle(in9, in10));
		TYPE c1(cx<1,19>(a1) + cx<2,19>(a2) + cx<3,19>(a3) + cx<4,19>(a4) + cx<5,19>(a5) + cx<6,19>(a6) + cx<7,19>(a7) + cx<8,19>(a8) + cx<9,19>(a9));
		TYPE s1(sx<1,19>(t1) + sx<2,19>(t2) + sx<3,19>(t3) + sx<4,19>(t4) + sx<5,19>(t5) + sx<6,19>(t6) + sx<7,19>(t7) + sx<8,19>(t8) + sx<9,19>(t9));
		TYPE c2(cx<2,19>(a1) + cx<4,19>(a2) + cx<6,19>(a3) + cx<8,19>(a4) + cx<9,19>(a5) + cx<7,19>(a6) + cx<5,19>(a7) + cx<3,19>(a8) + cx<1,19>(a9));
		TYPE s2(sx<2,19>(t1) + sx<4,19>(t2) + sx<6,19>(t3) + sx<8,19>(t4) - sx<9,19>(t5) - sx<7,19>(t6) - sx<5,19>(t7) - sx<3,19>(t8) - sx<1,19>(t9));
		TYPE c3(cx<3,19>(a1) + cx<6,19>(a2) + cx<9,19>(a3) + cx<7,19>(a4) + cx<4,19>(a5) + cx<1,19>(a6) + cx<2,19>(a7) + cx<5,19>(a8) + cx<8,19>(a9));
		TYPE s3(sx<3,19>(t1) + sx<6,19>(t2) + sx<9,19>(t3) - sx<7,19>(t4) - sx<4,19>(t5) - sx<1,19>(t6) + sx<2,19>(t7) + sx<5,19>(t8) + sx<8,19>(t9));
		TYPE c4(cx<4,19>(a1) + cx<8,19>(a2) + cx<7,19>(a3) + cx<3,19>(a4) + cx<1,19>(a5) + cx<5,19>(a6) + cx<9,19>(a7) + cx<6,19>(a8) + cx<2,19>(a9));
		TYPE s4(sx<4,19>(t1) + sx<8,19>(t2) - sx<7,19>(t3) - sx<3,19>(t4) + sx<1,19>(t5) + sx<5,19>(t6) + sx<9,19>(t7) - sx<6,19>(t8) - sx<2,19>(t9));
		TYPE c5(cx<5,19>(a1) + cx<9,19>(a2) + cx<4,19>(a3) + cx<1,19>(a4) + cx<6,19>(a5) + cx<8,19>(a6) + cx<3,19>(a7) + cx<2,19>(a8) + cx<7,19>(a9));
		TYPE s5(sx<5,19>(t1) - sx<9,19>(t2) - sx<4,19>(t3) + sx<1,19>(t4) + sx<6,19>(t5) - sx<8,19>(t6) - sx<3,19>(t7) + sx<2,19>(t8) + sx<7,19>(t9));
		TYPE c6(cx<6,19>(a1) + cx<7,19>(a2) + cx<1,19>(a3) + cx<5,19>(a4) + cx<8,19>(a5) + cx<2,19>(a6) + cx<4,19>(a7) + cx<9,19>(a8) + cx<3,19>(a9));
		TYPE s6(sx<6,19>(t1) - sx<7,19>(t2) - sx<1,19>(t3) + sx<5,19>(t4) - sx<8,19>(t5) - sx<2,19>(t6) + sx<4,19>(t7) - sx<9,19>(t8) - sx<3,19>(t9));
		TYPE c7(cx<7,19>(a1) + cx<5,19>(a2) + cx<2,19>(a3) + cx<9,19>(a4) + cx<3,19>(a5) + cx<4,19>(a6) + cx<8,19>(a7) + cx<1,19>(a8) + cx<6,19>(a9));
		TYPE s7(sx<7,19>(t1) - sx<5,19>(t2) + sx<2,19>(t3) + sx<9,19>(t4) - sx<3,19>(t5) + sx<4,19>(t6) - sx<8,19>(t7) - sx<1,19>(t8) + sx<6,19>(t9));
		TYPE c8(cx<8,19>(a1) + cx<3,19>(a2) + cx<5,19>(a3) + cx<6,19>(a4) + cx<2,19>(a5) + cx<9,19>(a6) + cx<1,19>(a7) + cx<7,19>(a8) + cx<4,19>(a9));
		TYPE s8(sx<8,19>(t1) - sx<3,19>(t2) + sx<5,19>(t3) - sx<6,19>(t4) + sx<2,19>(t5) - sx<9,19>(t6) - sx<1,19>(t7) + sx<7,19>(t8) - sx<4,19>(t9));
		TYPE c9(cx<9,19>(a1) + cx<1,19>(a2) + cx<8,19>(a3) + cx<2,19>(a4) + cx<7,19>(a5) + cx<3,19>(a6) + cx<6,19>(a7) + cx<4,19>(a8) + cx<5,19>(a9));
		TYPE s9(sx<9,19>(t1) - sx<1,19>(t2) + sx<8,19>(t3) - sx<2,19>(t4) + sx<7,19>(t5) - sx<3,19>(t6) + sx<6,19>(t7) - sx<4,19>(t8) + sx<5,19>(t9));
		*out0 = in0 + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9;
		*out1 = in0 + c1 - s1;
		*out2 = in0 + c2 - s2;
		*out3 = in0 + c3 - s3;
		*out4 = in0 + c4 - s4;
		*out5 = in0 + c5 - s5;
		*out6 = in0 + c6 - s6;
		*out7 = in0 + c7 - s7;
		*out8 = in0 + c8 - s8;
		*out9 = in0 + c9 - s9;
		*out10 = in0 + c9 + s9;
		*out11 = in0 + c8 + s8;
		*out12 = in0 + c7 + s7;
		*out13 = in0 + c6 + s6;
		*out14 = in0 + c5 + s5;
		*out15 = in0 + c4 + s4;
		*out16 = in0 + c3 + s3;
		*out17 = in0 + c2 + s2;
		*out18 = in0 + c1 + s1;
	}
	static inline void dit(TYPE *out, const TYPE *in, const TYPE *)
	{
		dft(out, out + 1, out + 2, out + 3, out + 4, out + 5, out + 6, out + 7, out + 8, out + 9, out + 10, out + 11, out + 12, out + 13, out + 14, out + 15, out + 16, out + 17, out + 18,
			in[0], in[STRIDE], in[2 * STRIDE], in[3 * STRIDE], in[4 * STRIDE], in[5 * STRIDE], in[6 * STRIDE], in[7 * STRIDE], in[8 * STRIDE], in[9 * STRIDE], in[10 * STRIDE], in[11 * STRIDE], in[12 * STRIDE], in[13 * STRIDE], in[14 * STRIDE], in[15 * STRIDE], in[16 * STRIDE], in[17 * STRIDE], in[18 * STRIDE]);
	}
};

template <int BINS, int STRIDE, typename TYPE, int SIGN>
struct Dit<19, BINS, STRIDE, TYPE, SIGN>
{
	static const int RADIX = 19;
	static const int QUOTIENT = BINS / RADIX;
	static void dit(TYPE *out, const TYPE *in, const TYPE *z)
	{
		for (int o = 0, i = 0; o < BINS; o += QUOTIENT, i += STRIDE)
			Dit<split(QUOTIENT), QUOTIENT, RADIX * STRIDE, TYPE, SIGN>::dit(out + o, in + i, z);
		for (int k0 = 0, k1 = QUOTIENT, k2 = 2 * QUOTIENT, k3 = 3 * QUOTIENT, k4 = 4 * QUOTIENT, k5 = 5 * QUOTIENT, k6 = 6 * QUOTIENT, k7 = 7 * QUOTIENT, k8 = 8 * QUOTIENT, k9 = 9 * QUOTIENT, k10 = 10 * QUOTIENT, k11 = 11 * QUOTIENT, k12 = 12 * QUOTIENT, k13 = 13 * QUOTIENT, k14 = 14 * QUOTIENT, k15 = 15 * QUOTIENT, k16 = 16 * QUOTIENT, k17 = 17 * QUOTIENT, k18 = 18 * QUOTIENT,
				l1 = 0, l2 = 0, l3 = 0, l4 = 0, l5 = 0, l6 = 0, l7 = 0, l8 = 0, l9 = 0, l10 = 0, l11 = 0, l12 = 0, l13 = 0, l14 = 0, l15 = 0, l16 = 0, l17 = 0, l18 = 0;
				k0 < QUOTIENT;
				++k0, ++k1, ++k2, ++k3, ++k4, ++k5, ++k6, ++k7, ++k8, ++k9, ++k10, ++k11, ++k12, ++k13, ++k14, ++k15, ++k16, ++k17, ++k18,
				l1 += STRIDE, l2 += 2 * STRIDE, l3 += 3 * STRIDE, l4 += 4 * STRIDE, l5 += 5 * STRIDE, l6 += 6 * STRIDE, l7 += 7 * STRIDE, l8 += 8 * STRIDE, l9 += 9 * STRIDE, l10 += 10 * STRIDE, l11 += 11 * STRIDE, l12 += 12 * STRIDE, l13 += 13 * STRIDE, l14 += 14 * STRIDE, l15 += 15 * STRIDE, l16 += 16 * STRIDE, l17 += 17 * STRIDE, l18 += 18 * STRIDE)
			Dit<RADIX, RADIX, STRIDE, TYPE, SIGN>::dft(out + k0, out + k1, out + k2, out + k3, out + k4, out + k5, out + k6, out + k7, out + k8, out + k9, out + k10, out + k11, out + k12, out + k13, out + k14, out + k15, out + k16, out + k17, out + k18,
				out[k0], z[l1] * out[k1], z[l2] * out[k2], z[l3] * out[k3], z[l4] * out[k4], z[l5] * out[k5], z[l6] * out[k6], z[l7] * out[k7], z[l8] * out[k8], z[l9] * out[k9], z[l10] * out[k10], z[l11] * out[k11], z[l12] * out[k12], z[l13] * out[k13], z[l14] * out[k14], z[l15] * out[k15], z[l16] * out[k16], z[l17] * out[k17], z[l18] * out[k18]);
	}
};

template <int STRIDE, typename TYPE>
struct Dit<23, 23, STRIDE, TYPE, -1>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE *out2, TYPE *out3, TYPE *out4, TYPE *out5, TYPE *out6, TYPE *out7, TYPE *out8, TYPE *out9, TYPE *out10, TYPE *out11, TYPE *out12, TYPE *out13, TYPE *out14, TYPE *out15, TYPE *out16, TYPE *out17, TYPE *out18, TYPE *out19, TYPE *out20, TYPE *out21, TYPE *out22,
			TYPE in0, TYPE in1, TYPE in2, TYPE in3, TYPE in4, TYPE in5, TYPE in6, TYPE in7, TYPE in8, TYPE in9, TYPE in10, TYPE in11, TYPE in12, TYPE in13, TYPE in14, TYPE in15, TYPE in16, TYPE in17, TYPE in18, TYPE in19, TYPE in20, TYPE in21, TYPE in22)
	{
		Dit<23, 23, STRIDE, TYPE, 1>::dft(out0, out22, out21, out20, out19, out18, out17, out16, out15, out14, out13, out12, out11, out10, out9, out8, out7, out6, out5, out4, out3, out2, out1,
			in0, in1, in2, in3, in4, in5, in6, in7, in8, in9, in10, in11, in12, in13, in14, in15, in16, in17, in18, in19, in20, in21, in22);
	}
	static inline void dit(TYPE *out, const TYPE *in, const TYPE *)
	{
		dft(out, out + 1, out + 2, out + 3, out + 4, out + 5, out + 6, out + 7, out + 8, out + 9, out + 10, out + 11, out + 12, out + 13, out + 14, out + 15, out + 16, out + 17, out + 18, out + 19, out + 20, out + 21, out + 22,
			in[0], in[STRIDE], in[2 * STRIDE], in[3 * STRIDE], in[4 * STRIDE], in[5 * STRIDE], in[6 * STRIDE], in[7 * STRIDE], in[8 * STRIDE], in[9 * STRIDE], in[10 * STRIDE], in[11 * STRIDE], in[12 * STRIDE], in[13 * STRIDE], in[14 * STRIDE], in[15 * STRIDE], in[16 * STRIDE], in[17 * STRIDE], in[18 * STRIDE], in[19 * STRIDE], in[20 * STRIDE], in[21 * STRIDE], in[22 * STRIDE]);
	}
};

template <int STRIDE, typename TYPE>
struct Dit<23, 23, STRIDE, TYPE, 1>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE *out2, TYPE *out3, TYPE *out4, TYPE *out5, TYPE *out6, TYPE *out7, TYPE *out8, TYPE *out9, TYPE *out10, TYPE *out11, TYPE *out12, TYPE *out13, TYPE *out14, TYPE *out15, TYPE *out16, TYPE *out17, TYPE *out18, TYPE *out19, TYPE *out20, TYPE *out21, TYPE *out22,
			TYPE in0, TYPE in1, TYPE in2, TYPE in3, TYPE in4, TYPE in5, TYPE in6, TYPE in7, TYPE in8, TYPE in9, TYPE in10, TYPE in11, TYPE in12, TYPE in13, TYPE in14, TYPE in15, TYPE in16, TYPE in17, TYPE in18, TYPE in19, TYPE in20, TYPE in21, TYPE in22)
	{
		TYPE a1(in1 + in22), t1(twiddle(in1, in22));
		TYPE a2(in2 + in21), t2(twiddle(in2, in21));
		TYPE a3(in3 + in20), t3(twiddle(in3, in20));
		TYPE a4(in4 + in19), t4(twiddle(in4, in19));
		TYPE a5(in5 + in18), t5(twiddle(in5, in18));
		TYPE a6(in6 + in17), t6(twiddle(in6, in17));
		TYPE a7(in7 + in16), t7(twiddle(in7, in16));
		TYPE a8(in8 + in15), t8(twiddle(in8, in15));
		TYPE a9(in9 + in14), t9(twiddle(in9, in14));
		TYPE a10(in10 + in13), t10(twiddle(in10, in13));
		TYPE a11(in11 + in12), t11(twiddle(in11, in12));
		TYPE c1(cx<1,23>(a1) + cx<2,23>(a2) + cx<3,23>(a3) + cx<4,23>(a4) + cx<5,23>(a5) + cx<6,23>(a6) + cx<7,23>(a7) + cx<8,23>(a8) + cx<9,23>(a9) + cx<10,23>(a10) + cx<11,23>(a11));
		TYPE s1(sx<1,23>(t1) + sx<2,23>(t2) + sx<3,23>(t3) + sx<4,23>(t4) + sx<5,23>(t5) + sx<6,23>(t6) + sx<7,23>(t7) + sx<8,23>(t8) + sx<9,23>(t9) + sx<10,23>(t10) + sx<11,23>(t11));
		TYPE c2(cx<2,23>(a1) + cx<4,23>(a2) + cx<6,23>(a3) + cx<8,23>(a4) + cx<10,23>(a5) + cx<11,23>(a6) + cx<9,23>(a7) + cx<7,23>(a8) + cx<5,23>(a9) + cx<3,23>(a10) + cx<1,23>(a11));
		TYPE s2(sx<2,23>(t1) + sx<4,23>(t2) + sx<6,23>(t3) + sx<8,23>(t4) + sx<10,23>(t5) - sx<11,23>(t6) - sx<9,23>(t7) - sx<7,23>(t8) - sx<5,23>(t9) - sx<3,23>(t10) - sx<1,23>(t11));
		TYPE c3(cx<3,23>(a1) + cx<6,23>(a2) + cx<9,23>(a3) + cx<11,23>(a4) + cx<8,23>(a5) + cx<5,23>(a6) + cx<2,23>(a7) + cx<1,23>(a8) + cx<4,23>(a9) + cx<7,23>(a10) + cx<10,23>(a11));
		TYPE s3(sx<3,23>(t1) + sx<6,23>(t2) + sx<9,23>(t3) - sx<11,23>(t4) - sx<8,23>(t5) - sx<5,23>(t6) - sx<2,23>(t7) + sx<1,23>(t8) + sx<4,23>(t9) + sx<7,23>(t10) + sx<10,23>(t11));
		TYPE c4(cx<4,23>(a1) + cx<8,23>(a2) + cx<11,23>(a3) + cx<7,23>(a4) + cx<3,23>(a5) + cx<1,23>(a6) + cx<5,23>(a7) + cx<9,23>(a8) + cx<10,23>(a9) + cx<6,23>(a10) + cx<2,23>(a11));
		TYPE s4(sx<4,23>(t1) + sx<8,23>(t2) - sx<11,23>(t3) - sx<7,23>(t4) - sx<3,23>(t5) + sx<1,23>(t6) + sx<5,23>(t7) + sx<9,23>(t8) - sx<10,23>(t9) - sx<6,23>(t10) - sx<2,23>(t11));
		TYPE c5(cx<5,23>(a1) + cx<10,23>(a2) + cx<8,23>(a3) + cx<3,23>(a4) + cx<2,23>(a5) + cx<7,23>(a6) + cx<11,23>(a7) + cx<6,23>(a8) + cx<1,23>(a9) + cx<4,23>(a10) + cx<9,23>(a11));
		TYPE s5(sx<5,23>(t1) + sx<10,23>(t2) - sx<8,23>(t3) - sx<3,23>(t4) + sx<2,23>(t5) + sx<7,23>(t6) - sx<11,23>(t7) - sx<6,23>(t8) - sx<1,23>(t9) + sx<4,23>(t10) + sx<9,23>(t11));
		TYPE c6(cx<6,23>(a1) + cx<11,23>(a2) + cx<5,23>(a3) + cx<1,23>(a4) + cx<7,23>(a5) + cx<10,23>(a6) + cx<4,23>(a7) + cx<2,23>(a8) + cx<8,23>(a9) + cx<9,23>(a10) + cx<3,23>(a11));
		TYPE s6(sx<6,23>(t1) - sx<11,23>(t2) - sx<5,23>(t3) + sx<1,23>(t4) + sx<7,23>(t5) - sx<10,23>(t6) - sx<4,23>(t7) + sx<2,23>(t8) + sx<8,23>(t9) - sx<9,23>(t10) - sx<3,23>(t11));
		TYPE c7(cx<7,23>(a1) + cx<9,23>(a2) + cx<2,23>(a3) + cx<5,23>(a4) + cx<11,23>(a5) + cx<4,23>(a6) + cx<3,23>(a7) + cx<10,23>(a8) + cx<6,23>(a9) + cx<1,23>(a10) + cx<8,23>(a11));
		TYPE s7(sx<7,23>(t1) - sx<9,23>(t2) - sx<2,23>(t3) + sx<5,23>(t4) - sx<11,23>(t5) - sx<4,23>(t6) + sx<3,23>(t7) + sx<10,23>(t8) - sx<6,23>(t9) + sx<1,23>(t10) + sx<8,23>(t11));
		TYPE c8(cx<8,23>(a1) + cx<7,23>(a2) + cx<1,23>(a3) + cx<9,23>(a4) + cx<6,23>(a5) + cx<2,23>(a6) + cx<10,23>(a7) + cx<5,23>(a8) + cx<3,23>(a9) + cx<11,23>(a10) + cx<4,23>(a11));
		TYPE s8(sx<8,23>(t1) - sx<7,23>(t2) + sx<1,23>(t3) + sx<9,23>(t4) - sx<6,23>(t5) + sx<2,23>(t6) + sx<10,23>(t7) - sx<5,23>(t8) + sx<3,23>(t9) + sx<11,23>(t10) - sx<4,23>(t11));
		TYPE c9(cx<9,23>(a1) + cx<5,23>(a2) + cx<4,23>(a3) + cx<10,23>(a4) + cx<1,23>(a5) + cx<8,23>(a6) + cx<6,23>(a7) + cx<3,23>(a8) + cx<11,23>(a9) + cx<2,23>(a10) + cx<7,23>(a11));
		TYPE s9(sx<9,23>(t1) - sx<5,23>(t2) + sx<4,23>(t3) - sx<10,23>(t4) - sx<1,23>(t5) + sx<8,23>(t6) - sx<6,23>(t7) + sx<3,23>(t8) - sx<11,23>(t9) - sx<2,23>(t10) + sx<7,23>(t11));
		TYPE c10(cx<10,23>(a1) + cx<3,23>(a2) + cx<7,23>(a3) + cx<6,23>(a4) + cx<4,23>(a5) + cx<9,23>(a6) + cx<1,23>(a7) + cx<11,23>(a8) + cx<2,23>(a9) + cx<8,23>(a10) + cx<5,23>(a11));
		TYPE s10(sx<10,23>(t1) - sx<3,23>(t2) + sx<7,23>(t3) - sx<6,23>(t4) + sx<4,23>(t5) - sx<9,23>(t6) + sx<1,23>(t7) + sx<11,23>(t8) - sx<2,23>(t9) + sx<8,23>(t10) - sx<5,23>(t11));
		TYPE c11(cx<11,23>(a1) + cx<1,23>(a2) + cx<10,23>(a3) + cx<2,23>(a4) + cx<9,23>(a5) + cx<3,23>(a6) + cx<8,23>(a7) + cx<4,23>(a8) + cx<7,23>(a9) + cx<5,23>(a10) + cx<6,23>(a11));
		TYPE s11(sx<11,23>(t1) - sx<1,23>(t2) + sx<10,23>(t3) - sx<2,23>(t4) + sx<9,23>(t5) - sx<3,23>(t6) + sx<8,23>(t7) - sx<4,23>(t8) + sx<7,23>(t9) - sx<5,23>(t10) + sx<6,23>(t11));
		*out0 = in0 + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11;
		*out1 = in0 + c1 - s1;
		*out2 = in0 + c2 - s2;
		*out3 = in0 + c3 - s3;
		*out4 = in0 + c4 - s4;
		*out5 = in0 + c5 - s5;
		*out6 = in0 + c6 - s6;
		*out7 = in0 + c7 - s7;
		*out8 = in0 + c8 - s8;
		*out9 = in0 + c9 - s9;
		*out10 = in0 + c10 - s10;
		*out11 = in0 + c11 - s11;
		*out12 = in0 + c11 + s11;
		*out13 = in0 + c10 + s10;
		*out14 = in0 + c9 + s9;
		*out15 = in0 + c8 + s8;
		*out16 = in0 + c7 + s7;
		*out17 = in0 + c6 + s6;
		*out18 = in0 + c5 + s5;
		*out19 = in0 + c4 + s4;
		*out20 = in0 + c3 + s3;
		*out21 = in0 + c2 + s2;
		*out22 = in0 + c1 + s1;
	}
	static inline void dit(TYPE *out, const TYPE *in, const TYPE *)
	{
		dft(out, out + 1, out + 2, out + 3, out + 4, out + 5, out + 6, out + 7, out + 8, out + 9, out + 10, out + 11, out + 12, out + 13, out + 14, out + 15, out + 16, out + 17, out + 18, out + 19, out + 20, out + 21, out + 22,
			in[0], in[STRIDE], in[2 * STRIDE], in[3 * STRIDE], in[4 * STRIDE], in[5 * STRIDE], in[6 * STRIDE], in[7 * STRIDE], in[8 * STRIDE], in[9 * STRIDE], in[10 * STRIDE], in[11 * STRIDE], in[12 * STRIDE], in[13 * STRIDE], in[14 * STRIDE], in[15 * STRIDE], in[16 * STRIDE], in[17 * STRIDE], in[18 * STRIDE], in[19 * STRIDE], in[20 * STRIDE], in[21 * STRIDE], in[22 * STRIDE]);
	}
};

template <int BINS, int STRIDE, typename TYPE, int SIGN>
struct Dit<23, BINS, STRIDE, TYPE, SIGN>
{
	static const int RADIX = 23;
	static const int QUOTIENT = BINS / RADIX;
	static void dit(TYPE *out, const TYPE *in, const TYPE *z)
	{
		for (int o = 0, i = 0; o < BINS; o += QUOTIENT, i += STRIDE)
			Dit<split(QUOTIENT), QUOTIENT, RADIX * STRIDE, TYPE, SIGN>::dit(out + o, in + i, z);
		for (int k0 = 0, k1 = QUOTIENT, k2 = 2 * QUOTIENT, k3 = 3 * QUOTIENT, k4 = 4 * QUOTIENT, k5 = 5 * QUOTIENT, k6 = 6 * QUOTIENT, k7 = 7 * QUOTIENT, k8 = 8 * QUOTIENT, k9 = 9 * QUOTIENT, k10 = 10 * QUOTIENT, k11 = 11 * QUOTIENT, k12 = 12 * QUOTIENT, k13 = 13 * QUOTIENT, k14 = 14 * QUOTIENT, k15 = 15 * QUOTIENT, k16 = 16 * QUOTIENT, k17 = 17 * QUOTIENT, k18 = 18 * QUOTIENT, k19 = 19 * QUOTIENT, k20 = 20 * QUOTIENT, k21 = 21 * QUOTIENT, k22 = 22 * QUOTIENT,
				l1 = 0, l2 = 0, l3 = 0, l4 = 0, l5 = 0, l6 = 0, l7 = 0, l8 = 0, l9 = 0, l10 = 0, l11 = 0, l12 = 0, l13 = 0, l14 = 0, l15 = 0, l16 = 0, l17 = 0, l18 = 0, l19 = 0, l20 = 0, l21 = 0, l22 = 0;
				k0 < QUOTIENT;
				++k0, ++k1, ++k2, ++k3, ++k4, ++k5, ++k6, ++k7, ++k8, ++k9, ++k10, ++k11, ++k12, ++k13, ++k14, ++k15, ++k16, ++k17, ++k18, ++k19, ++k20, ++k21, ++k22,
				l1 += STRIDE, l2 += 2 * STRIDE, l3 += 3 * STRIDE, l4 += 4 * STRIDE, l5 += 5 * STRIDE, l6 += 6 * STRIDE, l7 += 7 * STRIDE, l8 += 8 * STRIDE, l9 += 9 * STRIDE, l10 += 10 * STRIDE, l11 += 11 * STRIDE, l12 += 12 * STRIDE, l13 += 13 * STRIDE, l14 += 14 * STRIDE, l15 += 15 * STRIDE, l16 += 16 * STRIDE, l17 += 17 * STRIDE, l18 += 18 * STRIDE, l19 += 19 * STRIDE, l20 += 20 * STRIDE, l21 += 21 * STRIDE, l22 += 22 * STRIDE)
			Dit<RADIX, RADIX, STRIDE, TYPE, SIGN>::dft(out + k0, out + k1, out + k2, out + k3, out + k4, out + k5, out + k6, out + k7, out + k8, out + k9, out + k10, out + k11, out + k12, out + k13, out + k14, out + k15, out + k16, out + k17, out + k18, out + k19, out + k20, out + k21, out + k22,
				out[k0], z[l1] * out[k1], z[l2] * out[k2], z[l3] * out[k3], z[l4] * out[k4], z[l5] * out[k5], z[l6] * out[k6], z[l7] * out[k7], z[l8] * out[k8], z[l9] * out[k9], z[l10] * out[k10], z[l11] * out[k11], z[l12] * out[k12], z[l13] * out[k13], z[l14] * out[k14], z[l15] * out[k15], z[l16] * out[k16], z[l17] * out[k17], z[l18] * out[k18], z[l19] * out[k19], z[l20] * out[k20], z[l21] * out[k21], z[l22] * out[k22]);
	}
};

template <int STRIDE, typename TYPE>
struct Dit<29, 29, STRIDE, TYPE, -1>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE *out2, TYPE *out3, TYPE *out4, TYPE *out5, TYPE *out6, TYPE *out7, TYPE *out8, TYPE *out9, TYPE *out10, TYPE *out11, TYPE *out12, TYPE *out13, TYPE *out14, TYPE *out15, TYPE *out16, TYPE *out17, TYPE *out18, TYPE *out19, TYPE *out20, TYPE *out21, TYPE *out22, TYPE *out23, TYPE *out24, TYPE *out25, TYPE *out26, TYPE *out27, TYPE *out28,
			TYPE in0, TYPE in1, TYPE in2, TYPE in3, TYPE in4, TYPE in5, TYPE in6, TYPE in7, TYPE in8, TYPE in9, TYPE in10, TYPE in11, TYPE in12, TYPE in13, TYPE in14, TYPE in15, TYPE in16, TYPE in17, TYPE in18, TYPE in19, TYPE in20, TYPE in21, TYPE in22, TYPE in23, TYPE in24, TYPE in25, TYPE in26, TYPE in27, TYPE in28)
	{
		Dit<29, 29, STRIDE, TYPE, 1>::dft(out0, out28, out27, out26, out25, out24, out23, out22, out21, out20, out19, out18, out17, out16, out15, out14, out13, out12, out11, out10, out9, out8, out7, out6, out5, out4, out3, out2, out1,
			in0, in1, in2, in3, in4, in5, in6, in7, in8, in9, in10, in11, in12, in13, in14, in15, in16, in17, in18, in19, in20, in21, in22, in23, in24, in25, in26, in27, in28);
	}
	static inline void dit(TYPE *out, const TYPE *in, const TYPE *)
	{
		dft(out, out + 1, out + 2, out + 3, out + 4, out + 5, out + 6, out + 7, out + 8, out + 9, out + 10, out + 11, out + 12, out + 13, out + 14, out + 15, out + 16, out + 17, out + 18, out + 19, out + 20, out + 21, out + 22, out + 23, out + 24, out + 25, out + 26, out + 27, out + 28,
			in[0], in[STRIDE], in[2 * STRIDE], in[3 * STRIDE], in[4 * STRIDE], in[5 * STRIDE], in[6 * STRIDE], in[7 * STRIDE], in[8 * STRIDE], in[9 * STRIDE], in[10 * STRIDE], in[11 * STRIDE], in[12 * STRIDE], in[13 * STRIDE], in[14 * STRIDE], in[15 * STRIDE], in[16 * STRIDE], in[17 * STRIDE], in[18 * STRIDE], in[19 * STRIDE], in[20 * STRIDE], in[21 * STRIDE], in[22 * STRIDE], in[23 * STRIDE], in[24 * STRIDE], in[25 * STRIDE], in[26 * STRIDE], in[27 * STRIDE], in[28 * STRIDE]);
	}
};

template <int STRIDE, typename TYPE>
struct Dit<29, 29, STRIDE, TYPE, 1>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE *out2, TYPE *out3, TYPE *out4, TYPE *out5, TYPE *out6, TYPE *out7, TYPE *out8, TYPE *out9, TYPE *out10, TYPE *out11, TYPE *out12, TYPE *out13, TYPE *out14, TYPE *out15, TYPE *out16, TYPE *out17, TYPE *out18, TYPE *out19, TYPE *out20, TYPE *out21, TYPE *out22, TYPE *out23, TYPE *out24, TYPE *out25, TYPE *out26, TYPE *out27, TYPE *out28,
			TYPE in0, TYPE in1, TYPE in2, TYPE in3, TYPE in4, TYPE in5, TYPE in6, TYPE in7, TYPE in8, TYPE in9, TYPE in10, TYPE in11, TYPE in12, TYPE in13, TYPE in14, TYPE in15, TYPE in16, TYPE in17, TYPE in18, TYPE in19, TYPE in20, TYPE in21, TYPE in22, TYPE in23, TYPE in24, TYPE in25, TYPE in26, TYPE in27, TYPE in28)
	{
		TYPE a1(in1 + in28), t1(twiddle(in1, in28));
		TYPE a2(in2 + in27), t2(twiddle(in2, in27));
		TYPE a3(in3 + in26), t3(twiddle(in3, in26));
		TYPE a4(in4 + in25), t4(twiddle(in4, in25));
		TYPE a5(in5 + in24), t5(twiddle(in5, in24));
		TYPE a6(in6 + in23), t6(twiddle(in6, in23));
		TYPE a7(in7 + in22), t7(twiddle(in7, in22));
		TYPE a8(in8 + in21), t8(twiddle(in8, in21));
		TYPE a9(in9 + in20), t9(twiddle(in9, in20));
		TYPE a10(in10 + in19), t10(twiddle(in10, in19));
		TYPE a11(in11 + in18), t11(twiddle(in11, in18));
		TYPE a12(in12 + in17), t12(twiddle(in12, in17));
		TYPE a13(in13 + in16), t13(twiddle(in13, in16));
		TYPE a14(in14 + in15), t14(twiddle(in14, in15));
		TYPE c1(cx<1,29>(a1) + cx<2,29>(a2) + cx<3,29>(a3) + cx<4,29>(a4) + cx<5,29>(a5) + cx<6,29>(a6) + cx<7,29>(a7) + cx<8,29>(a8) + cx<9,29>(a9) + cx<10,29>(a10) + cx<11,29>(a11) + cx<12,29>(a12) + cx<13,29>(a13) + cx<14,29>(a14));
		TYPE s1(sx<1,29>(t1) + sx<2,29>(t2) + sx<3,29>(t3) + sx<4,29>(t4) + sx<5,29>(t5) + sx<6,29>(t6) + sx<7,29>(t7) + sx<8,29>(t8) + sx<9,29>(t9) + sx<10,29>(t10) + sx<11,29>(t11) + sx<12,29>(t12) + sx<13,29>(t13) + sx<14,29>(t14));
		TYPE c2(cx<2,29>(a1) + cx<4,29>(a2) + cx<6,29>(a3) + cx<8,29>(a4) + cx<10,29>(a5) + cx<12,29>(a6) + cx<14,29>(a7) + cx<13,29>(a8) + cx<11,29>(a9) + cx<9,29>(a10) + cx<7,29>(a11) + cx<5,29>(a12) + cx<3,29>(a13) + cx<1,29>(a14));
		TYPE s2(sx<2,29>(t1) + sx<4,29>(t2) + sx<6,29>(t3) + sx<8,29>(t4) + sx<10,29>(t5) + sx<12,29>(t6) + sx<14,29>(t7) - sx<13,29>(t8) - sx<11,29>(t9) - sx<9,29>(t10) - sx<7,29>(t11) - sx<5,29>(t12) - sx<3,29>(t13) - sx<1,29>(t14));
		TYPE c3(cx<3,29>(a1) + cx<6,29>(a2) + cx<9,29>(a3) + cx<12,29>(a4) + cx<14,29>(a5) + cx<11,29>(a6) + cx<8,29>(a7) + cx<5,29>(a8) + cx<2,29>(a9) + cx<1,29>(a10) + cx<4,29>(a11) + cx<7,29>(a12) + cx<10,29>(a13) + cx<13,29>(a14));
		TYPE s3(sx<3,29>(t1) + sx<6,29>(t2) + sx<9,29>(t3) + sx<12,29>(t4) - sx<14,29>(t5) - sx<11,29>(t6) - sx<8,29>(t7) - sx<5,29>(t8) - sx<2,29>(t9) + sx<1,29>(t10) + sx<4,29>(t11) + sx<7,29>(t12) + sx<10,29>(t13) + sx<13,29>(t14));
		TYPE c4(cx<4,29>(a1) + cx<8,29>(a2) + cx<12,29>(a3) + cx<13,29>(a4) + cx<9,29>(a5) + cx<5,29>(a6) + cx<1,29>(a7) + cx<3,29>(a8) + cx<7,29>(a9) + cx<11,29>(a10) + cx<14,29>(a11) + cx<10,29>(a12) + cx<6,29>(a13) + cx<2,29>(a14));
		TYPE s4(sx<4,29>(t1) + sx<8,29>(t2) + sx<12,29>(t3) - sx<13,29>(t4) - sx<9,29>(t5) - sx<5,29>(t6) - sx<1,29>(t7) + sx<3,29>(t8) + sx<7,29>(t9) + sx<11,29>(t10) - sx<14,29>(t11) - sx<10,29>(t12) - sx<6,29>(t13) - sx<2,29>(t14));
		TYPE c5(cx<5,29>(a1) + cx<10,29>(a2) + cx<14,29>(a3) + cx<9,29>(a4) + cx<4,29>(a5) + cx<1,29>(a6) + cx<6,29>(a7) + cx<11,29>(a8) + cx<13,29>(a9) + cx<8,29>(a10) + cx<3,29>(a11) + cx<2,29>(a12) + cx<7,29>(a13) + cx<12,29>(a14));
		TYPE s5(sx<5,29>(t1) + sx<10,29>(t2) - sx<14,29>(t3) - sx<9,29>(t4) - sx<4,29>(t5) + sx<1,29>(t6) + sx<6,29>(t7) + sx<11,29>(t8) - sx<13,29>(t9) - sx<8,29>(t10) - sx<3,29>(t11) + sx<2,29>(t12) + sx<7,29>(t13) + sx<12,29>(t14));
		TYPE c6(cx<6,29>(a1) + cx<12,29>(a2) + cx<11,29>(a3) + cx<5,29>(a4) + cx<1,29>(a5) + cx<7,29>(a6) + cx<13,29>(a7) + cx<10,29>(a8) + cx<4,29>(a9) + cx<2,29>(a10) + cx<8,29>(a11) + cx<14,29>(a12) + cx<9,29>(a13) + cx<3,29>(a14));
		TYPE s6(sx<6,29>(t1) + sx<12,29>(t2) - sx<11,29>(t3) - sx<5,29>(t4) + sx<1,29>(t5) + sx<7,29>(t6) + sx<13,29>(t7) - sx<10,29>(t8) - sx<4,29>(t9) + sx<2,29>(t10) + sx<8,29>(t11) + sx<14,29>(t12) - sx<9,29>(t13) - sx<3,29>(t14));
		TYPE c7(cx<7,29>(a1) + cx<14,29>(a2) + cx<8,29>(a3) + cx<1,29>(a4) + cx<6,29>(a5) + cx<13,29>(a6) + cx<9,29>(a7) + cx<2,29>(a8) + cx<5,29>(a9) + cx<12,29>(a10) + cx<10,29>(a11) + cx<3,29>(a12) + cx<4,29>(a13) + cx<11,29>(a14));
		TYPE s7(sx<7,29>(t1) + sx<14,29>(t2) - sx<8,29>(t3) - sx<1,29>(t4) + sx<6,29>(t5) + sx<13,29>(t6) - sx<9,29>(t7) - sx<2,29>(t8) + sx<5,29>(t9) + sx<12,29>(t10) - sx<10,29>(t11) - sx<3,29>(t12) + sx<4,29>(t13) + sx<11,29>(t14));
		TYPE c8(cx<8,29>(a1) + cx<13,29>(a2) + cx<5,29>(a3) + cx<3,29>(a4) + cx<11,29>(a5) + cx<10,29>(a6) + cx<2,29>(a7) + cx<6,29>(a8) + cx<14,29>(a9) + cx<7,29>(a10) + cx<1,29>(a11) + cx<9,29>(a12) + cx<12,29>(a13) + cx<4,29>(a14));
		TYPE s8(sx<8,29>(t1) - sx<13,29>(t2) - sx<5,29>(t3) + sx<3,29>(t4) + sx<11,29>(t5) - sx<10,29>(t6) - sx<2,29>(t7) + sx<6,29>(t8) + sx<14,29>(t9) - sx<7,29>(t10) + sx<1,29>(t11) + sx<9,29>(t12) - sx<12,29>(t13) - sx<4,29>(t14));
		TYPE c9(cx<9,29>(a1) + cx<11,29>(a2) + cx<2,29>(a3) + cx<7,29>(a4) + cx<13,29>(a5) + cx<4,29>(a6) + cx<5,29>(a7) + cx<14,29>(a8) + cx<6,29>(a9) + cx<3,29>(a10) + cx<12,29>(a11) + cx<8,29>(a12) + cx<1,29>(a13) + cx<10,29>(a14));
		TYPE s9(sx<9,29>(t1) - sx<11,29>(t2) - sx<2,29>(t3) + sx<7,29>(t4) - sx<13,29>(t5) - sx<4,29>(t6) + sx<5,29>(t7) + sx<14,29>(t8) - sx<6,29>(t9) + sx<3,29>(t10) + sx<12,29>(t11) - sx<8,29>(t12) + sx<1,29>(t13) + sx<10,29>(t14));
		TYPE c10(cx<10,29>(a1) + cx<9,29>(a2) + cx<1,29>(a3) + cx<11,29>(a4) + cx<8,29>(a5) + cx<2,29>(a6) + cx<12,29>(a7) + cx<7,29>(a8) + cx<3,29>(a9) + cx<13,29>(a10) + cx<6,29>(a11) + cx<4,29>(a12) + cx<14,29>(a13) + cx<5,29>(a14));
		TYPE s10(sx<10,29>(t1) - sx<9,29>(t2) + sx<1,29>(t3) + sx<11,29>(t4) - sx<8,29>(t5) + sx<2,29>(t6) + sx<12,29>(t7) - sx<7,29>(t8) + sx<3,29>(t9) + sx<13,29>(t10) - sx<6,29>(t11) + sx<4,29>(t12) + sx<14,29>(t13) - sx<5,29>(t14));
		TYPE c11(cx<11,29>(a1) + cx<7,29>(a2) + cx<4,29>(a3) + cx<14,29>(a4) + cx<3,29>(a5) + cx<8,29>(a6) + cx<10,29>(a7) + cx<1,29>(a8) + cx<12,29>(a9) + cx<6,29>(a10) + cx<5,29>(a11) + cx<13,29>(a12) + cx<2,29>(a13) + cx<9,29>(a14));
		TYPE s11(sx<11,29>(t1) - sx<7,29>(t2) + sx<4,29>(t3) - sx<14,29>(t4) - sx<3,29>(t5) + sx<8,29>(t6) - sx<10,29>(t7) + sx<1,29>(t8) + sx<12,29>(t9) - sx<6,29>(t10) + sx<5,29>(t11) - sx<13,29>(t12) - sx<2,29>(t13) + sx<9,29>(t14));
		TYPE c12(cx<12,29>(a1) + cx<5,29>(a2) + cx<7,29>(a3) + cx<10,29>(a4) + cx<2,29>(a5) + cx<14,29>(a6) + cx<3,29>(a7) + cx<9,29>(a8) + cx<8,29>(a9) + cx<4,29>(a10) + cx<13,29>(a11) + cx<1,29>(a12) + cx<11,29>(a13) + cx<6,29>(a14));
		TYPE s12(sx<12,29>(t1) - sx<5,29>(t2) + sx<7,29>(t3) - sx<10,29>(t4) + sx<2,29>(t5) + sx<14,29>(t6) - sx<3,29>(t7) + sx<9,29>(t8) - sx<8,29>(t9) + sx<4,29>(t10) - sx<13,29>(t11) - sx<1,29>(t12) + sx<11,29>(t13) - sx<6,29>(t14));
		TYPE c13(cx<13,29>(a1) + cx<3,29>(a2) + cx<10,29>(a3) + cx<6,29>(a4) + cx<7,29>(a5) + cx<9,29>(a6) + cx<4,29>(a7) + cx<12,29>(a8) + cx<1,29>(a9) + cx<14,29>(a10) + cx<2,29>(a11) + cx<11,29>(a12) + cx<5,29>(a13) + cx<8,29>(a14));
		TYPE s13(sx<13,29>(t1) - sx<3,29>(t2) + sx<10,29>(t3) - sx<6,29>(t4) + sx<7,29>(t5) - sx<9,29>(t6) + sx<4,29>(t7) - sx<12,29>(t8) + sx<1,29>(t9) + sx<14,29>(t10) - sx<2,29>(t11) + sx<11,29>(t12) - sx<5,29>(t13) + sx<8,29>(t14));
		TYPE c14(cx<14,29>(a1) + cx<1,29>(a2) + cx<13,29>(a3) + cx<2,29>(a4) + cx<12,29>(a5) + cx<3,29>(a6) + cx<11,29>(a7) + cx<4,29>(a8) + cx<10,29>(a9) + cx<5,29>(a10) + cx<9,29>(a11) + cx<6,29>(a12) + cx<8,29>(a13) + cx<7,29>(a14));
		TYPE s14(sx<14,29>(t1) - sx<1,29>(t2) + sx<13,29>(t3) - sx<2,29>(t4) + sx<12,29>(t5) - sx<3,29>(t6) + sx<11,29>(t7) - sx<4,29>(t8) + sx<10,29>(t9) - sx<5,29>(t10) + sx<9,29>(t11) - sx<6,29>(t12) + sx<8,29>(t13) - sx<7,29>(t14));
		*out0 = in0 + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12 + a13 + a14;
		*out1 = in0 + c1 - s1;
		*out2 = in0 + c2 - s2;
		*out3 = in0 + c3 - s3;
		*out4 = in0 + c4 - s4;
		*out5 = in0 + c5 - s5;
		*out6 = in0 + c6 - s6;
		*out7 = in0 + c7 - s7;
		*out8 = in0 + c8 - s8;
		*out9 = in0 + c9 - s9;
		*out10 = in0 + c10 - s10;
		*out11 = in0 + c11 - s11;
		*out12 = in0 + c12 - s12;
		*out13 = in0 + c13 - s13;
		*out14 = in0 + c14 - s14;
		*out15 = in0 + c14 + s14;
		*out16 = in0 + c13 + s13;
		*out17 = in0 + c12 + s12;
		*out18 = in0 + c11 + s11;
		*out19 = in0 + c10 + s10;
		*out20 = in0 + c9 + s9;
		*out21 = in0 + c8 + s8;
		*out22 = in0 + c7 + s7;
		*out23 = in0 + c6 + s6;
		*out24 = in0 + c5 + s5;
		*out25 = in0 + c4 + s4;
		*out26 = in0 + c3 + s3;
		*out27 = in0 + c2 + s2;
		*out28 = in0 + c1 + s1;
	}
	static inline void dit(TYPE *out, const TYPE *in, const TYPE *)
	{
		dft(out, out + 1, out + 2, out + 3, out + 4, out + 5, out + 6, out + 7, out + 8, out + 9, out + 10, out + 11, out + 12, out + 13, out + 14, out + 15, out + 16, out + 17, out + 18, out + 19, out + 20, out + 21, out + 22, out + 23, out + 24, out + 25, out + 26, out + 27, out + 28,
			in[0], in[STRIDE], in[2 * STRIDE], in[3 * STRIDE], in[4 * STRIDE], in[5 * STRIDE], in[6 * STRIDE], in[7 * STRIDE], in[8 * STRIDE], in[9 * STRIDE], in[10 * STRIDE], in[11 * STRIDE], in[12 * STRIDE], in[13 * STRIDE], in[14 * STRIDE], in[15 * STRIDE], in[16 * STRIDE], in[17 * STRIDE], in[18 * STRIDE], in[19 * STRIDE], in[20 * STRIDE], in[21 * STRIDE], in[22 * STRIDE], in[23 * STRIDE], in[24 * STRIDE], in[25 * STRIDE], in[26 * STRIDE], in[27 * STRIDE], in[28 * STRIDE]);
	}
};

template <int BINS, int STRIDE, typename TYPE, int SIGN>
struct Dit<29, BINS, STRIDE, TYPE, SIGN>
{
	static const int RADIX = 29;
	static const int QUOTIENT = BINS / RADIX;
	static void dit(TYPE *out, const TYPE *in, const TYPE *z)
	{
		for (int o = 0, i = 0; o < BINS; o += QUOTIENT, i += STRIDE)
			Dit<split(QUOTIENT), QUOTIENT, RADIX * STRIDE, TYPE, SIGN>::dit(out + o, in + i, z);
		for (int k0 = 0, k1 = QUOTIENT, k2 = 2 * QUOTIENT, k3 = 3 * QUOTIENT, k4 = 4 * QUOTIENT, k5 = 5 * QUOTIENT, k6 = 6 * QUOTIENT, k7 = 7 * QUOTIENT, k8 = 8 * QUOTIENT, k9 = 9 * QUOTIENT, k10 = 10 * QUOTIENT, k11 = 11 * QUOTIENT, k12 = 12 * QUOTIENT, k13 = 13 * QUOTIENT, k14 = 14 * QUOTIENT, k15 = 15 * QUOTIENT, k16 = 16 * QUOTIENT, k17 = 17 * QUOTIENT, k18 = 18 * QUOTIENT, k19 = 19 * QUOTIENT, k20 = 20 * QUOTIENT, k21 = 21 * QUOTIENT, k22 = 22 * QUOTIENT, k23 = 23 * QUOTIENT, k24 = 24 * QUOTIENT, k25 = 25 * QUOTIENT, k26 = 26 * QUOTIENT, k27 = 27 * QUOTIENT, k28 = 28 * QUOTIENT,
				l1 = 0, l2 = 0, l3 = 0, l4 = 0, l5 = 0, l6 = 0, l7 = 0, l8 = 0, l9 = 0, l10 = 0, l11 = 0, l12 = 0, l13 = 0, l14 = 0, l15 = 0, l16 = 0, l17 = 0, l18 = 0, l19 = 0, l20 = 0, l21 = 0, l22 = 0, l23 = 0, l24 = 0, l25 = 0, l26 = 0, l27 = 0, l28 = 0;
				k0 < QUOTIENT;
				++k0, ++k1, ++k2, ++k3, ++k4, ++k5, ++k6, ++k7, ++k8, ++k9, ++k10, ++k11, ++k12, ++k13, ++k14, ++k15, ++k16, ++k17, ++k18, ++k19, ++k20, ++k21, ++k22, ++k23, ++k24, ++k25, ++k26, ++k27, ++k28,
				l1 += STRIDE, l2 += 2 * STRIDE, l3 += 3 * STRIDE, l4 += 4 * STRIDE, l5 += 5 * STRIDE, l6 += 6 * STRIDE, l7 += 7 * STRIDE, l8 += 8 * STRIDE, l9 += 9 * STRIDE, l10 += 10 * STRIDE, l11 += 11 * STRIDE, l12 += 12 * STRIDE, l13 += 13 * STRIDE, l14 += 14 * STRIDE, l15 += 15 * STRIDE, l16 += 16 * STRIDE, l17 += 17 * STRIDE, l18 += 18 * STRIDE, l19 += 19 * STRIDE, l20 += 20 * STRIDE, l21 += 21 * STRIDE, l22 += 22 * STRIDE, l23 += 23 * STRIDE, l24 += 24 * STRIDE, l25 += 25 * STRIDE, l26 += 26 * STRIDE, l27 += 27 * STRIDE, l28 += 28 * STRIDE)
			Dit<RADIX, RADIX, STRIDE, TYPE, SIGN>::dft(out + k0, out + k1, out + k2, out + k3, out + k4, out + k5, out + k6, out + k7, out + k8, out + k9, out + k10, out + k11, out + k12, out + k13, out + k14, out + k15, out + k16, out + k17, out + k18, out + k19, out + k20, out + k21, out + k22, out + k23, out + k24, out + k25, out + k26, out + k27, out + k28,
				out[k0], z[l1] * out[k1], z[l2] * out[k2], z[l3] * out[k3], z[l4] * out[k4], z[l5] * out[k5], z[l6] * out[k6], z[l7] * out[k7], z[l8] * out[k8], z[l9] * out[k9], z[l10] * out[k10], z[l11] * out[k11], z[l12] * out[k12], z[l13] * out[k13], z[l14] * out[k14], z[l15] * out[k15], z[l16] * out[k16], z[l17] * out[k17], z[l18] * out[k18], z[l19] * out[k19], z[l20] * out[k20], z[l21] * out[k21], z[l22] * out[k22], z[l23] * out[k23], z[l24] * out[k24], z[l25] * out[k25], z[l26] * out[k26], z[l27] * out[k27], z[l28] * out[k28]);
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
