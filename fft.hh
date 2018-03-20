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
	return TYPE(a.real() + a.imag() - b.real() + b.imag(),
		a.imag() - a.real() - b.real() - b.imag());
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
		TYPE ac(in1 + in4), as(twiddle(in1, in4));
		TYPE bc(in2 + in3), bs(twiddle(in2, in3));
		TYPE c1(cx<1,5>(ac) + cx<2,5>(bc));
		TYPE s1(sx<1,5>(as) + sx<2,5>(bs));
		TYPE c2(cx<2,5>(ac) + cx<1,5>(bc));
		TYPE s2(sx<2,5>(as) - sx<1,5>(bs));
		*out0 = in0 + ac + bc;
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
		TYPE ac(in1 + in6), as(twiddle(in1, in6));
		TYPE bc(in2 + in5), bs(twiddle(in2, in5));
		TYPE cc(in3 + in4), cs(twiddle(in3, in4));
		TYPE c1(cx<1,7>(ac) + cx<2,7>(bc) + cx<3,7>(cc));
		TYPE s1(sx<1,7>(as) + sx<2,7>(bs) + sx<3,7>(cs));
		TYPE c2(cx<2,7>(ac) + cx<3,7>(bc) + cx<1,7>(cc));
		TYPE s2(sx<2,7>(as) - sx<3,7>(bs) - sx<1,7>(cs));
		TYPE c3(cx<3,7>(ac) + cx<1,7>(bc) + cx<2,7>(cc));
		TYPE s3(sx<3,7>(as) - sx<1,7>(bs) + sx<2,7>(cs));
		*out0 = in0 + ac + bc + cc;
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
		TYPE ac(in1 + in10), as(twiddle(in1, in10));
		TYPE bc(in2 + in9), bs(twiddle(in2, in9));
		TYPE cc(in3 + in8), cs(twiddle(in3, in8));
		TYPE dc(in4 + in7), ds(twiddle(in4, in7));
		TYPE ec(in5 + in6), es(twiddle(in5, in6));
		TYPE c1(cx<1,11>(ac) + cx<2,11>(bc) + cx<3,11>(cc) + cx<4,11>(dc) + cx<5,11>(ec));
		TYPE s1(sx<1,11>(as) + sx<2,11>(bs) + sx<3,11>(cs) + sx<4,11>(ds) + sx<5,11>(es));
		TYPE c2(cx<2,11>(ac) + cx<4,11>(bc) + cx<5,11>(cc) + cx<3,11>(dc) + cx<1,11>(ec));
		TYPE s2(sx<2,11>(as) + sx<4,11>(bs) - sx<5,11>(cs) - sx<3,11>(ds) - sx<1,11>(es));
		TYPE c3(cx<3,11>(ac) + cx<5,11>(bc) + cx<2,11>(cc) + cx<1,11>(dc) + cx<4,11>(ec));
		TYPE s3(sx<3,11>(as) - sx<5,11>(bs) - sx<2,11>(cs) + sx<1,11>(ds) + sx<4,11>(es));
		TYPE c4(cx<4,11>(ac) + cx<3,11>(bc) + cx<1,11>(cc) + cx<5,11>(dc) + cx<2,11>(ec));
		TYPE s4(sx<4,11>(as) - sx<3,11>(bs) + sx<1,11>(cs) + sx<5,11>(ds) - sx<2,11>(es));
		TYPE c5(cx<5,11>(ac) + cx<1,11>(bc) + cx<4,11>(cc) + cx<2,11>(dc) + cx<3,11>(ec));
		TYPE s5(sx<5,11>(as) - sx<1,11>(bs) + sx<4,11>(cs) - sx<2,11>(ds) + sx<3,11>(es));
		*out0 = in0 + ac + bc + cc + dc + ec;
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
		TYPE ac(in1 + in12), as(twiddle(in1, in12));
		TYPE bc(in2 + in11), bs(twiddle(in2, in11));
		TYPE cc(in3 + in10), cs(twiddle(in3, in10));
		TYPE dc(in4 + in9), ds(twiddle(in4, in9));
		TYPE ec(in5 + in8), es(twiddle(in5, in8));
		TYPE fc(in6 + in7), fs(twiddle(in6, in7));
		TYPE c1(cx<1,13>(ac) + cx<2,13>(bc) + cx<3,13>(cc) + cx<4,13>(dc) + cx<5,13>(ec) + cx<6,13>(fc));
		TYPE s1(sx<1,13>(as) + sx<2,13>(bs) + sx<3,13>(cs) + sx<4,13>(ds) + sx<5,13>(es) + sx<6,13>(fs));
		TYPE c2(cx<2,13>(ac) + cx<4,13>(bc) + cx<6,13>(cc) + cx<5,13>(dc) + cx<3,13>(ec) + cx<1,13>(fc));
		TYPE s2(sx<2,13>(as) + sx<4,13>(bs) + sx<6,13>(cs) - sx<5,13>(ds) - sx<3,13>(es) - sx<1,13>(fs));
		TYPE c3(cx<3,13>(ac) + cx<6,13>(bc) + cx<4,13>(cc) + cx<1,13>(dc) + cx<2,13>(ec) + cx<5,13>(fc));
		TYPE s3(sx<3,13>(as) + sx<6,13>(bs) - sx<4,13>(cs) - sx<1,13>(ds) + sx<2,13>(es) + sx<5,13>(fs));
		TYPE c4(cx<4,13>(ac) + cx<5,13>(bc) + cx<1,13>(cc) + cx<3,13>(dc) + cx<6,13>(ec) + cx<2,13>(fc));
		TYPE s4(sx<4,13>(as) - sx<5,13>(bs) - sx<1,13>(cs) + sx<3,13>(ds) - sx<6,13>(es) - sx<2,13>(fs));
		TYPE c5(cx<5,13>(ac) + cx<3,13>(bc) + cx<2,13>(cc) + cx<6,13>(dc) + cx<1,13>(ec) + cx<4,13>(fc));
		TYPE s5(sx<5,13>(as) - sx<3,13>(bs) + sx<2,13>(cs) - sx<6,13>(ds) - sx<1,13>(es) + sx<4,13>(fs));
		TYPE c6(cx<6,13>(ac) + cx<1,13>(bc) + cx<5,13>(cc) + cx<2,13>(dc) + cx<4,13>(ec) + cx<3,13>(fc));
		TYPE s6(sx<6,13>(as) - sx<1,13>(bs) + sx<5,13>(cs) - sx<2,13>(ds) + sx<4,13>(es) - sx<3,13>(fs));
		*out0 = in0 + ac + bc + cc + dc + ec + fc;
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
		TYPE ac(in1 + in16), as(twiddle(in1, in16));
		TYPE bc(in2 + in15), bs(twiddle(in2, in15));
		TYPE cc(in3 + in14), cs(twiddle(in3, in14));
		TYPE dc(in4 + in13), ds(twiddle(in4, in13));
		TYPE ec(in5 + in12), es(twiddle(in5, in12));
		TYPE fc(in6 + in11), fs(twiddle(in6, in11));
		TYPE gc(in7 + in10), gs(twiddle(in7, in10));
		TYPE hc(in8 + in9), hs(twiddle(in8, in9));
		TYPE c1(cx<1,17>(ac) + cx<2,17>(bc) + cx<3,17>(cc) + cx<4,17>(dc) + cx<5,17>(ec) + cx<6,17>(fc) + cx<7,17>(gc) + cx<8,17>(hc));
		TYPE s1(sx<1,17>(as) + sx<2,17>(bs) + sx<3,17>(cs) + sx<4,17>(ds) + sx<5,17>(es) + sx<6,17>(fs) + sx<7,17>(gs) + sx<8,17>(hs));
		TYPE c2(cx<2,17>(ac) + cx<4,17>(bc) + cx<6,17>(cc) + cx<8,17>(dc) + cx<7,17>(ec) + cx<5,17>(fc) + cx<3,17>(gc) + cx<1,17>(hc));
		TYPE s2(sx<2,17>(as) + sx<4,17>(bs) + sx<6,17>(cs) + sx<8,17>(ds) - sx<7,17>(es) - sx<5,17>(fs) - sx<3,17>(gs) - sx<1,17>(hs));
		TYPE c3(cx<3,17>(ac) + cx<6,17>(bc) + cx<8,17>(cc) + cx<5,17>(dc) + cx<2,17>(ec) + cx<1,17>(fc) + cx<4,17>(gc) + cx<7,17>(hc));
		TYPE s3(sx<3,17>(as) + sx<6,17>(bs) - sx<8,17>(cs) - sx<5,17>(ds) - sx<2,17>(es) + sx<1,17>(fs) + sx<4,17>(gs) + sx<7,17>(hs));
		TYPE c4(cx<4,17>(ac) + cx<8,17>(bc) + cx<5,17>(cc) + cx<1,17>(dc) + cx<3,17>(ec) + cx<7,17>(fc) + cx<6,17>(gc) + cx<2,17>(hc));
		TYPE s4(sx<4,17>(as) + sx<8,17>(bs) - sx<5,17>(cs) - sx<1,17>(ds) + sx<3,17>(es) + sx<7,17>(fs) - sx<6,17>(gs) - sx<2,17>(hs));
		TYPE c5(cx<5,17>(ac) + cx<7,17>(bc) + cx<2,17>(cc) + cx<3,17>(dc) + cx<8,17>(ec) + cx<4,17>(fc) + cx<1,17>(gc) + cx<6,17>(hc));
		TYPE s5(sx<5,17>(as) - sx<7,17>(bs) - sx<2,17>(cs) + sx<3,17>(ds) + sx<8,17>(es) - sx<4,17>(fs) + sx<1,17>(gs) + sx<6,17>(hs));
		TYPE c6(cx<6,17>(ac) + cx<5,17>(bc) + cx<1,17>(cc) + cx<7,17>(dc) + cx<4,17>(ec) + cx<2,17>(fc) + cx<8,17>(gc) + cx<3,17>(hc));
		TYPE s6(sx<6,17>(as) - sx<5,17>(bs) + sx<1,17>(cs) + sx<7,17>(ds) - sx<4,17>(es) + sx<2,17>(fs) + sx<8,17>(gs) - sx<3,17>(hs));
		TYPE c7(cx<7,17>(ac) + cx<3,17>(bc) + cx<4,17>(cc) + cx<6,17>(dc) + cx<1,17>(ec) + cx<8,17>(fc) + cx<2,17>(gc) + cx<5,17>(hc));
		TYPE s7(sx<7,17>(as) - sx<3,17>(bs) + sx<4,17>(cs) - sx<6,17>(ds) + sx<1,17>(es) + sx<8,17>(fs) - sx<2,17>(gs) + sx<5,17>(hs));
		TYPE c8(cx<8,17>(ac) + cx<1,17>(bc) + cx<7,17>(cc) + cx<2,17>(dc) + cx<6,17>(ec) + cx<3,17>(fc) + cx<5,17>(gc) + cx<4,17>(hc));
		TYPE s8(sx<8,17>(as) - sx<1,17>(bs) + sx<7,17>(cs) - sx<2,17>(ds) + sx<6,17>(es) - sx<3,17>(fs) + sx<5,17>(gs) - sx<4,17>(hs));
		*out0 = in0 + ac + bc + cc + dc + ec + fc + gc + hc;
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
		TYPE ac(in1 + in18), as(twiddle(in1, in18));
		TYPE bc(in2 + in17), bs(twiddle(in2, in17));
		TYPE cc(in3 + in16), cs(twiddle(in3, in16));
		TYPE dc(in4 + in15), ds(twiddle(in4, in15));
		TYPE ec(in5 + in14), es(twiddle(in5, in14));
		TYPE fc(in6 + in13), fs(twiddle(in6, in13));
		TYPE gc(in7 + in12), gs(twiddle(in7, in12));
		TYPE hc(in8 + in11), hs(twiddle(in8, in11));
		TYPE ic(in9 + in10), is(twiddle(in9, in10));
		TYPE c1(cx<1,19>(ac) + cx<2,19>(bc) + cx<3,19>(cc) + cx<4,19>(dc) + cx<5,19>(ec) + cx<6,19>(fc) + cx<7,19>(gc) + cx<8,19>(hc) + cx<9,19>(ic));
		TYPE s1(sx<1,19>(as) + sx<2,19>(bs) + sx<3,19>(cs) + sx<4,19>(ds) + sx<5,19>(es) + sx<6,19>(fs) + sx<7,19>(gs) + sx<8,19>(hs) + sx<9,19>(is));
		TYPE c2(cx<2,19>(ac) + cx<4,19>(bc) + cx<6,19>(cc) + cx<8,19>(dc) + cx<9,19>(ec) + cx<7,19>(fc) + cx<5,19>(gc) + cx<3,19>(hc) + cx<1,19>(ic));
		TYPE s2(sx<2,19>(as) + sx<4,19>(bs) + sx<6,19>(cs) + sx<8,19>(ds) - sx<9,19>(es) - sx<7,19>(fs) - sx<5,19>(gs) - sx<3,19>(hs) - sx<1,19>(is));
		TYPE c3(cx<3,19>(ac) + cx<6,19>(bc) + cx<9,19>(cc) + cx<7,19>(dc) + cx<4,19>(ec) + cx<1,19>(fc) + cx<2,19>(gc) + cx<5,19>(hc) + cx<8,19>(ic));
		TYPE s3(sx<3,19>(as) + sx<6,19>(bs) + sx<9,19>(cs) - sx<7,19>(ds) - sx<4,19>(es) - sx<1,19>(fs) + sx<2,19>(gs) + sx<5,19>(hs) + sx<8,19>(is));
		TYPE c4(cx<4,19>(ac) + cx<8,19>(bc) + cx<7,19>(cc) + cx<3,19>(dc) + cx<1,19>(ec) + cx<5,19>(fc) + cx<9,19>(gc) + cx<6,19>(hc) + cx<2,19>(ic));
		TYPE s4(sx<4,19>(as) + sx<8,19>(bs) - sx<7,19>(cs) - sx<3,19>(ds) + sx<1,19>(es) + sx<5,19>(fs) + sx<9,19>(gs) - sx<6,19>(hs) - sx<2,19>(is));
		TYPE c5(cx<5,19>(ac) + cx<9,19>(bc) + cx<4,19>(cc) + cx<1,19>(dc) + cx<6,19>(ec) + cx<8,19>(fc) + cx<3,19>(gc) + cx<2,19>(hc) + cx<7,19>(ic));
		TYPE s5(sx<5,19>(as) - sx<9,19>(bs) - sx<4,19>(cs) + sx<1,19>(ds) + sx<6,19>(es) - sx<8,19>(fs) - sx<3,19>(gs) + sx<2,19>(hs) + sx<7,19>(is));
		TYPE c6(cx<6,19>(ac) + cx<7,19>(bc) + cx<1,19>(cc) + cx<5,19>(dc) + cx<8,19>(ec) + cx<2,19>(fc) + cx<4,19>(gc) + cx<9,19>(hc) + cx<3,19>(ic));
		TYPE s6(sx<6,19>(as) - sx<7,19>(bs) - sx<1,19>(cs) + sx<5,19>(ds) - sx<8,19>(es) - sx<2,19>(fs) + sx<4,19>(gs) - sx<9,19>(hs) - sx<3,19>(is));
		TYPE c7(cx<7,19>(ac) + cx<5,19>(bc) + cx<2,19>(cc) + cx<9,19>(dc) + cx<3,19>(ec) + cx<4,19>(fc) + cx<8,19>(gc) + cx<1,19>(hc) + cx<6,19>(ic));
		TYPE s7(sx<7,19>(as) - sx<5,19>(bs) + sx<2,19>(cs) + sx<9,19>(ds) - sx<3,19>(es) + sx<4,19>(fs) - sx<8,19>(gs) - sx<1,19>(hs) + sx<6,19>(is));
		TYPE c8(cx<8,19>(ac) + cx<3,19>(bc) + cx<5,19>(cc) + cx<6,19>(dc) + cx<2,19>(ec) + cx<9,19>(fc) + cx<1,19>(gc) + cx<7,19>(hc) + cx<4,19>(ic));
		TYPE s8(sx<8,19>(as) - sx<3,19>(bs) + sx<5,19>(cs) - sx<6,19>(ds) + sx<2,19>(es) - sx<9,19>(fs) - sx<1,19>(gs) + sx<7,19>(hs) - sx<4,19>(is));
		TYPE c9(cx<9,19>(ac) + cx<1,19>(bc) + cx<8,19>(cc) + cx<2,19>(dc) + cx<7,19>(ec) + cx<3,19>(fc) + cx<6,19>(gc) + cx<4,19>(hc) + cx<5,19>(ic));
		TYPE s9(sx<9,19>(as) - sx<1,19>(bs) + sx<8,19>(cs) - sx<2,19>(ds) + sx<7,19>(es) - sx<3,19>(fs) + sx<6,19>(gs) - sx<4,19>(hs) + sx<5,19>(is));
		*out0 = in0 + ac + bc + cc + dc + ec + fc + gc + hc + ic;
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
