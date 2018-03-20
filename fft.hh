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
	return (!(N%13)) ? 13 : (!(N%11)) ? 11 : (!(N%7)) ? 7 : (!(N%5)) ? 5 : (!(N%3)) ? 3 : (!(N%8)&&pow8(N)) ? 8 : (!(N%4)&&pow4(N)) ? 4 : (!(N%2)) ? 2 : 1;
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
		TYPE c1a(cx<1,5>(a)), c1b(cx<1,5>(b)), s1c(sx<1,5>(c)), s1d(sx<1,5>(d));
		TYPE c2a(cx<2,5>(a)), c2b(cx<2,5>(b)), s2c(sx<2,5>(c)), s2d(sx<2,5>(d));
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
		TYPE c1a(cx<1,5>(a)), c1b(cx<1,5>(b)), s1c(sx<1,5>(c)), s1d(sx<1,5>(d));
		TYPE c2a(cx<2,5>(a)), c2b(cx<2,5>(b)), s2c(sx<2,5>(c)), s2d(sx<2,5>(d));
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
		TYPE c1a(cx<1,7>(a)), c1b(cx<1,7>(b)), c1c(cx<1,7>(c)), s1d(sx<1,7>(d)), s1e(sx<1,7>(e)), s1f(sx<1,7>(f));
		TYPE c2a(cx<2,7>(a)), c2b(cx<2,7>(b)), c2c(cx<2,7>(c)), s2d(sx<2,7>(d)), s2e(sx<2,7>(e)), s2f(sx<2,7>(f));
		TYPE c3a(cx<3,7>(a)), c3b(cx<3,7>(b)), c3c(cx<3,7>(c)), s3d(sx<3,7>(d)), s3e(sx<3,7>(e)), s3f(sx<3,7>(f));

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
		TYPE c1a(cx<1,7>(a)), c1b(cx<1,7>(b)), c1c(cx<1,7>(c)), s1d(sx<1,7>(d)), s1e(sx<1,7>(e)), s1f(sx<1,7>(f));
		TYPE c2a(cx<2,7>(a)), c2b(cx<2,7>(b)), c2c(cx<2,7>(c)), s2d(sx<2,7>(d)), s2e(sx<2,7>(e)), s2f(sx<2,7>(f));
		TYPE c3a(cx<3,7>(a)), c3b(cx<3,7>(b)), c3c(cx<3,7>(c)), s3d(sx<3,7>(d)), s3e(sx<3,7>(e)), s3f(sx<3,7>(f));

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

template <int STRIDE, typename TYPE>
struct Dit<11, 11, STRIDE, TYPE, -1>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE *out2, TYPE *out3, TYPE *out4, TYPE *out5, TYPE *out6, TYPE *out7, TYPE *out8, TYPE *out9, TYPE *out10,
			TYPE in0, TYPE in1, TYPE in2, TYPE in3, TYPE in4, TYPE in5, TYPE in6, TYPE in7, TYPE in8, TYPE in9, TYPE in10)
	{
		TYPE a(in1 + in10), b(in2 + in9), c(in3 + in8), d(in4 + in7), e(in5 + in6), f(twiddle(in1, in10)), g(twiddle(in2, in9)), h(twiddle(in3, in8)), i(twiddle(in4, in7)), j(twiddle(in5, in6));
		TYPE c1a(cx<1,11>(a)), c1b(cx<1,11>(b)), c1c(cx<1,11>(c)), c1d(cx<1,11>(d)), c1e(cx<1,11>(e)), s1f(sx<1,11>(f)), s1g(sx<1,11>(g)), s1h(sx<1,11>(h)), s1i(sx<1,11>(i)), s1j(sx<1,11>(j));
		TYPE c2a(cx<2,11>(a)), c2b(cx<2,11>(b)), c2c(cx<2,11>(c)), c2d(cx<2,11>(d)), c2e(cx<2,11>(e)), s2f(sx<2,11>(f)), s2g(sx<2,11>(g)), s2h(sx<2,11>(h)), s2i(sx<2,11>(i)), s2j(sx<2,11>(j));
		TYPE c3a(cx<3,11>(a)), c3b(cx<3,11>(b)), c3c(cx<3,11>(c)), c3d(cx<3,11>(d)), c3e(cx<3,11>(e)), s3f(sx<3,11>(f)), s3g(sx<3,11>(g)), s3h(sx<3,11>(h)), s3i(sx<3,11>(i)), s3j(sx<3,11>(j));
		TYPE c4a(cx<4,11>(a)), c4b(cx<4,11>(b)), c4c(cx<4,11>(c)), c4d(cx<4,11>(d)), c4e(cx<4,11>(e)), s4f(sx<4,11>(f)), s4g(sx<4,11>(g)), s4h(sx<4,11>(h)), s4i(sx<4,11>(i)), s4j(sx<4,11>(j));
		TYPE c5a(cx<5,11>(a)), c5b(cx<5,11>(b)), c5c(cx<5,11>(c)), c5d(cx<5,11>(d)), c5e(cx<5,11>(e)), s5f(sx<5,11>(f)), s5g(sx<5,11>(g)), s5h(sx<5,11>(h)), s5i(sx<5,11>(i)), s5j(sx<5,11>(j));

		*out0 = in0 + a + b + c + d + e;
		*out1 = in0 + c1a + c2b + c3c + c4d + c5e + s1f + s2g + s3h + s4i + s5j;
		*out2 = in0 + c2a + c4b + c5c + c3d + c1e + s2f + s4g - s5h - s3i - s1j;
		*out3 = in0 + c3a + c5b + c2c + c1d + c4e + s3f - s5g - s2h + s1i + s4j;
		*out4 = in0 + c4a + c3b + c1c + c5d + c2e + s4f - s3g + s1h + s5i - s2j;
		*out5 = in0 + c5a + c1b + c4c + c2d + c3e + s5f - s1g + s4h - s2i + s3j;
		*out6 = in0 + c5a + c1b + c4c + c2d + c3e - s5f + s1g - s4h + s2i - s3j;
		*out7 = in0 + c4a + c3b + c1c + c5d + c2e - s4f + s3g - s1h - s5i + s2j;
		*out8 = in0 + c3a + c5b + c2c + c1d + c4e - s3f + s5g + s2h - s1i - s4j;
		*out9 = in0 + c2a + c4b + c5c + c3d + c1e - s2f - s4g + s5h + s3i + s1j;
		*out10= in0 + c1a + c2b + c3c + c4d + c5e - s1f - s2g - s3h - s4i - s5j;
	}
	static inline void dit(TYPE *out, const TYPE *in, const TYPE *)
	{
		dft(out, out + 1, out + 2, out + 3, out + 4, out + 5, out + 6, out + 7, out + 8, out + 9, out + 10, in[0], in[STRIDE], in[2 * STRIDE], in[3 * STRIDE], in[4 * STRIDE], in[5 * STRIDE], in[6 * STRIDE], in[7 * STRIDE], in[8 * STRIDE], in[9 * STRIDE], in[10 * STRIDE]);
	}
};

template <int STRIDE, typename TYPE>
struct Dit<11, 11, STRIDE, TYPE, 1>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE *out2, TYPE *out3, TYPE *out4, TYPE *out5, TYPE *out6, TYPE *out7, TYPE *out8, TYPE *out9, TYPE *out10,
			TYPE in0, TYPE in1, TYPE in2, TYPE in3, TYPE in4, TYPE in5, TYPE in6, TYPE in7, TYPE in8, TYPE in9, TYPE in10)
	{
		TYPE a(in1 + in10), b(in2 + in9), c(in3 + in8), d(in4 + in7), e(in5 + in6), f(twiddle(in1, in10)), g(twiddle(in2, in9)), h(twiddle(in3, in8)), i(twiddle(in4, in7)), j(twiddle(in5, in6));
		TYPE c1a(cx<1,11>(a)), c1b(cx<1,11>(b)), c1c(cx<1,11>(c)), c1d(cx<1,11>(d)), c1e(cx<1,11>(e)), s1f(sx<1,11>(f)), s1g(sx<1,11>(g)), s1h(sx<1,11>(h)), s1i(sx<1,11>(i)), s1j(sx<1,11>(j));
		TYPE c2a(cx<2,11>(a)), c2b(cx<2,11>(b)), c2c(cx<2,11>(c)), c2d(cx<2,11>(d)), c2e(cx<2,11>(e)), s2f(sx<2,11>(f)), s2g(sx<2,11>(g)), s2h(sx<2,11>(h)), s2i(sx<2,11>(i)), s2j(sx<2,11>(j));
		TYPE c3a(cx<3,11>(a)), c3b(cx<3,11>(b)), c3c(cx<3,11>(c)), c3d(cx<3,11>(d)), c3e(cx<3,11>(e)), s3f(sx<3,11>(f)), s3g(sx<3,11>(g)), s3h(sx<3,11>(h)), s3i(sx<3,11>(i)), s3j(sx<3,11>(j));
		TYPE c4a(cx<4,11>(a)), c4b(cx<4,11>(b)), c4c(cx<4,11>(c)), c4d(cx<4,11>(d)), c4e(cx<4,11>(e)), s4f(sx<4,11>(f)), s4g(sx<4,11>(g)), s4h(sx<4,11>(h)), s4i(sx<4,11>(i)), s4j(sx<4,11>(j));
		TYPE c5a(cx<5,11>(a)), c5b(cx<5,11>(b)), c5c(cx<5,11>(c)), c5d(cx<5,11>(d)), c5e(cx<5,11>(e)), s5f(sx<5,11>(f)), s5g(sx<5,11>(g)), s5h(sx<5,11>(h)), s5i(sx<5,11>(i)), s5j(sx<5,11>(j));

		*out0 = in0 + a + b + c + d + e;
		*out1 = in0 + c1a + c2b + c3c + c4d + c5e - s1f - s2g - s3h - s4i - s5j;
		*out2 = in0 + c2a + c4b + c5c + c3d + c1e - s2f - s4g + s5h + s3i + s1j;
		*out3 = in0 + c3a + c5b + c2c + c1d + c4e - s3f + s5g + s2h - s1i - s4j;
		*out4 = in0 + c4a + c3b + c1c + c5d + c2e - s4f + s3g - s1h - s5i + s2j;
		*out5 = in0 + c5a + c1b + c4c + c2d + c3e - s5f + s1g - s4h + s2i - s3j;
		*out6 = in0 + c5a + c1b + c4c + c2d + c3e + s5f - s1g + s4h - s2i + s3j;
		*out7 = in0 + c4a + c3b + c1c + c5d + c2e + s4f - s3g + s1h + s5i - s2j;
		*out8 = in0 + c3a + c5b + c2c + c1d + c4e + s3f - s5g - s2h + s1i + s4j;
		*out9 = in0 + c2a + c4b + c5c + c3d + c1e + s2f + s4g - s5h - s3i - s1j;
		*out10= in0 + c1a + c2b + c3c + c4d + c5e + s1f + s2g + s3h + s4i + s5j;
	}
	static inline void dit(TYPE *out, const TYPE *in, const TYPE *)
	{
		dft(out, out + 1, out + 2, out + 3, out + 4, out + 5, out + 6, out + 7, out + 8, out + 9, out + 10, in[0], in[STRIDE], in[2 * STRIDE], in[3 * STRIDE], in[4 * STRIDE], in[5 * STRIDE], in[6 * STRIDE], in[7 * STRIDE], in[8 * STRIDE], in[9 * STRIDE], in[10 * STRIDE]);
	}
};

template <int STRIDE, typename TYPE>
struct Dit<13, 13, STRIDE, TYPE, -1>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE *out2, TYPE *out3, TYPE *out4, TYPE *out5, TYPE *out6, TYPE *out7, TYPE *out8, TYPE *out9, TYPE *out10, TYPE *out11, TYPE *out12,
			TYPE in0, TYPE in1, TYPE in2, TYPE in3, TYPE in4, TYPE in5, TYPE in6, TYPE in7, TYPE in8, TYPE in9, TYPE in10, TYPE in11, TYPE in12)
	{
		*out0 = ex<-0,13>(in0) + ex<-0,13>(in1) + ex<-0,13>(in2) + ex<-0,13>(in3) + ex<-0,13>(in4) + ex<-0,13>(in5) + ex<-0,13>(in6) + ex<-0,13>(in7) + ex<-0,13>(in8) + ex<-0,13>(in9) + ex<-0,13>(in10) + ex<-0,13>(in11) + ex<-0,13>(in12);
		*out1 = ex<-0,13>(in0) + ex<-1,13>(in1) + ex<-2,13>(in2) + ex<-3,13>(in3) + ex<-4,13>(in4) + ex<-5,13>(in5) + ex<-6,13>(in6) + ex<-7,13>(in7) + ex<-8,13>(in8) + ex<-9,13>(in9) + ex<-10,13>(in10) + ex<-11,13>(in11) + ex<-12,13>(in12);
		*out2 = ex<-0,13>(in0) + ex<-2,13>(in1) + ex<-4,13>(in2) + ex<-6,13>(in3) + ex<-8,13>(in4) + ex<-10,13>(in5) + ex<-12,13>(in6) + ex<-1,13>(in7) + ex<-3,13>(in8) + ex<-5,13>(in9) + ex<-7,13>(in10) + ex<-9,13>(in11) + ex<-11,13>(in12);
		*out3 = ex<-0,13>(in0) + ex<-3,13>(in1) + ex<-6,13>(in2) + ex<-9,13>(in3) + ex<-12,13>(in4) + ex<-2,13>(in5) + ex<-5,13>(in6) + ex<-8,13>(in7) + ex<-11,13>(in8) + ex<-1,13>(in9) + ex<-4,13>(in10) + ex<-7,13>(in11) + ex<-10,13>(in12);
		*out4 = ex<-0,13>(in0) + ex<-4,13>(in1) + ex<-8,13>(in2) + ex<-12,13>(in3) + ex<-3,13>(in4) + ex<-7,13>(in5) + ex<-11,13>(in6) + ex<-2,13>(in7) + ex<-6,13>(in8) + ex<-10,13>(in9) + ex<-1,13>(in10) + ex<-5,13>(in11) + ex<-9,13>(in12);
		*out5 = ex<-0,13>(in0) + ex<-5,13>(in1) + ex<-10,13>(in2) + ex<-2,13>(in3) + ex<-7,13>(in4) + ex<-12,13>(in5) + ex<-4,13>(in6) + ex<-9,13>(in7) + ex<-1,13>(in8) + ex<-6,13>(in9) + ex<-11,13>(in10) + ex<-3,13>(in11) + ex<-8,13>(in12);
		*out6 = ex<-0,13>(in0) + ex<-6,13>(in1) + ex<-12,13>(in2) + ex<-5,13>(in3) + ex<-11,13>(in4) + ex<-4,13>(in5) + ex<-10,13>(in6) + ex<-3,13>(in7) + ex<-9,13>(in8) + ex<-2,13>(in9) + ex<-8,13>(in10) + ex<-1,13>(in11) + ex<-7,13>(in12);
		*out7 = ex<-0,13>(in0) + ex<-7,13>(in1) + ex<-1,13>(in2) + ex<-8,13>(in3) + ex<-2,13>(in4) + ex<-9,13>(in5) + ex<-3,13>(in6) + ex<-10,13>(in7) + ex<-4,13>(in8) + ex<-11,13>(in9) + ex<-5,13>(in10) + ex<-12,13>(in11) + ex<-6,13>(in12);
		*out8 = ex<-0,13>(in0) + ex<-8,13>(in1) + ex<-3,13>(in2) + ex<-11,13>(in3) + ex<-6,13>(in4) + ex<-1,13>(in5) + ex<-9,13>(in6) + ex<-4,13>(in7) + ex<-12,13>(in8) + ex<-7,13>(in9) + ex<-2,13>(in10) + ex<-10,13>(in11) + ex<-5,13>(in12);
		*out9 = ex<-0,13>(in0) + ex<-9,13>(in1) + ex<-5,13>(in2) + ex<-1,13>(in3) + ex<-10,13>(in4) + ex<-6,13>(in5) + ex<-2,13>(in6) + ex<-11,13>(in7) + ex<-7,13>(in8) + ex<-3,13>(in9) + ex<-12,13>(in10) + ex<-8,13>(in11) + ex<-4,13>(in12);
		*out10= ex<-0,13>(in0) + ex<-10,13>(in1) + ex<-7,13>(in2) + ex<-4,13>(in3) + ex<-1,13>(in4) + ex<-11,13>(in5) + ex<-8,13>(in6) + ex<-5,13>(in7) + ex<-2,13>(in8) + ex<-12,13>(in9) + ex<-9,13>(in10) + ex<-6,13>(in11) + ex<-3,13>(in12);
		*out11= ex<-0,13>(in0) + ex<-11,13>(in1) + ex<-9,13>(in2) + ex<-7,13>(in3) + ex<-5,13>(in4) + ex<-3,13>(in5) + ex<-1,13>(in6) + ex<-12,13>(in7) + ex<-10,13>(in8) + ex<-8,13>(in9) + ex<-6,13>(in10) + ex<-4,13>(in11) + ex<-2,13>(in12);
		*out12= ex<-0,13>(in0) + ex<-12,13>(in1) + ex<-11,13>(in2) + ex<-10,13>(in3) + ex<-9,13>(in4) + ex<-8,13>(in5) + ex<-7,13>(in6) + ex<-6,13>(in7) + ex<-5,13>(in8) + ex<-4,13>(in9) + ex<-3,13>(in10) + ex<-2,13>(in11) + ex<-1,13>(in12);
	}
	static inline void dit(TYPE *out, const TYPE *in, const TYPE *)
	{
		dft(out, out + 1, out + 2, out + 3, out + 4, out + 5, out + 6, out + 7, out + 8, out + 9, out + 10, out + 11, out + 12, in[0], in[STRIDE], in[2 * STRIDE], in[3 * STRIDE], in[4 * STRIDE], in[5 * STRIDE], in[6 * STRIDE], in[7 * STRIDE], in[8 * STRIDE], in[9 * STRIDE], in[10 * STRIDE], in[11 * STRIDE], in[12 * STRIDE]);
	}
};

template <int STRIDE, typename TYPE>
struct Dit<13, 13, STRIDE, TYPE, 1>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE *out2, TYPE *out3, TYPE *out4, TYPE *out5, TYPE *out6, TYPE *out7, TYPE *out8, TYPE *out9, TYPE *out10, TYPE *out11, TYPE *out12,
			TYPE in0, TYPE in1, TYPE in2, TYPE in3, TYPE in4, TYPE in5, TYPE in6, TYPE in7, TYPE in8, TYPE in9, TYPE in10, TYPE in11, TYPE in12)
	{
		*out0 = ex<+0,13>(in0) + ex<+0,13>(in1) + ex<+0,13>(in2) + ex<+0,13>(in3) + ex<+0,13>(in4) + ex<+0,13>(in5) + ex<+0,13>(in6) + ex<+0,13>(in7) + ex<+0,13>(in8) + ex<+0,13>(in9) + ex<+0,13>(in10) + ex<+0,13>(in11) + ex<+0,13>(in12);
		*out1 = ex<+0,13>(in0) + ex<+1,13>(in1) + ex<+2,13>(in2) + ex<+3,13>(in3) + ex<+4,13>(in4) + ex<+5,13>(in5) + ex<+6,13>(in6) + ex<+7,13>(in7) + ex<+8,13>(in8) + ex<+9,13>(in9) + ex<+10,13>(in10) + ex<+11,13>(in11) + ex<+12,13>(in12);
		*out2 = ex<+0,13>(in0) + ex<+2,13>(in1) + ex<+4,13>(in2) + ex<+6,13>(in3) + ex<+8,13>(in4) + ex<+10,13>(in5) + ex<+12,13>(in6) + ex<+1,13>(in7) + ex<+3,13>(in8) + ex<+5,13>(in9) + ex<+7,13>(in10) + ex<+9,13>(in11) + ex<+11,13>(in12);
		*out3 = ex<+0,13>(in0) + ex<+3,13>(in1) + ex<+6,13>(in2) + ex<+9,13>(in3) + ex<+12,13>(in4) + ex<+2,13>(in5) + ex<+5,13>(in6) + ex<+8,13>(in7) + ex<+11,13>(in8) + ex<+1,13>(in9) + ex<+4,13>(in10) + ex<+7,13>(in11) + ex<+10,13>(in12);
		*out4 = ex<+0,13>(in0) + ex<+4,13>(in1) + ex<+8,13>(in2) + ex<+12,13>(in3) + ex<+3,13>(in4) + ex<+7,13>(in5) + ex<+11,13>(in6) + ex<+2,13>(in7) + ex<+6,13>(in8) + ex<+10,13>(in9) + ex<+1,13>(in10) + ex<+5,13>(in11) + ex<+9,13>(in12);
		*out5 = ex<+0,13>(in0) + ex<+5,13>(in1) + ex<+10,13>(in2) + ex<+2,13>(in3) + ex<+7,13>(in4) + ex<+12,13>(in5) + ex<+4,13>(in6) + ex<+9,13>(in7) + ex<+1,13>(in8) + ex<+6,13>(in9) + ex<+11,13>(in10) + ex<+3,13>(in11) + ex<+8,13>(in12);
		*out6 = ex<+0,13>(in0) + ex<+6,13>(in1) + ex<+12,13>(in2) + ex<+5,13>(in3) + ex<+11,13>(in4) + ex<+4,13>(in5) + ex<+10,13>(in6) + ex<+3,13>(in7) + ex<+9,13>(in8) + ex<+2,13>(in9) + ex<+8,13>(in10) + ex<+1,13>(in11) + ex<+7,13>(in12);
		*out7 = ex<+0,13>(in0) + ex<+7,13>(in1) + ex<+1,13>(in2) + ex<+8,13>(in3) + ex<+2,13>(in4) + ex<+9,13>(in5) + ex<+3,13>(in6) + ex<+10,13>(in7) + ex<+4,13>(in8) + ex<+11,13>(in9) + ex<+5,13>(in10) + ex<+12,13>(in11) + ex<+6,13>(in12);
		*out8 = ex<+0,13>(in0) + ex<+8,13>(in1) + ex<+3,13>(in2) + ex<+11,13>(in3) + ex<+6,13>(in4) + ex<+1,13>(in5) + ex<+9,13>(in6) + ex<+4,13>(in7) + ex<+12,13>(in8) + ex<+7,13>(in9) + ex<+2,13>(in10) + ex<+10,13>(in11) + ex<+5,13>(in12);
		*out9 = ex<+0,13>(in0) + ex<+9,13>(in1) + ex<+5,13>(in2) + ex<+1,13>(in3) + ex<+10,13>(in4) + ex<+6,13>(in5) + ex<+2,13>(in6) + ex<+11,13>(in7) + ex<+7,13>(in8) + ex<+3,13>(in9) + ex<+12,13>(in10) + ex<+8,13>(in11) + ex<+4,13>(in12);
		*out10= ex<+0,13>(in0) + ex<+10,13>(in1) + ex<+7,13>(in2) + ex<+4,13>(in3) + ex<+1,13>(in4) + ex<+11,13>(in5) + ex<+8,13>(in6) + ex<+5,13>(in7) + ex<+2,13>(in8) + ex<+12,13>(in9) + ex<+9,13>(in10) + ex<+6,13>(in11) + ex<+3,13>(in12);
		*out11= ex<+0,13>(in0) + ex<+11,13>(in1) + ex<+9,13>(in2) + ex<+7,13>(in3) + ex<+5,13>(in4) + ex<+3,13>(in5) + ex<+1,13>(in6) + ex<+12,13>(in7) + ex<+10,13>(in8) + ex<+8,13>(in9) + ex<+6,13>(in10) + ex<+4,13>(in11) + ex<+2,13>(in12);
		*out12= ex<+0,13>(in0) + ex<+12,13>(in1) + ex<+11,13>(in2) + ex<+10,13>(in3) + ex<+9,13>(in4) + ex<+8,13>(in5) + ex<+7,13>(in6) + ex<+6,13>(in7) + ex<+5,13>(in8) + ex<+4,13>(in9) + ex<+3,13>(in10) + ex<+2,13>(in11) + ex<+1,13>(in12);
	}
	static inline void dit(TYPE *out, const TYPE *in, const TYPE *)
	{
		dft(out, out + 1, out + 2, out + 3, out + 4, out + 5, out + 6, out + 7, out + 8, out + 9, out + 10, out + 11, out + 12, in[0], in[STRIDE], in[2 * STRIDE], in[3 * STRIDE], in[4 * STRIDE], in[5 * STRIDE], in[6 * STRIDE], in[7 * STRIDE], in[8 * STRIDE], in[9 * STRIDE], in[10 * STRIDE], in[11 * STRIDE], in[12 * STRIDE]);
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
