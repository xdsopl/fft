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
static inline TYPE cos2pi5(TYPE a) { return std::cos(typename TYPE::value_type(2 * M_PI / 5)) * a; }
template <typename TYPE>
static inline TYPE sin2pi5(TYPE a) { return std::sin(typename TYPE::value_type(2 * M_PI / 5)) * a; }
template <typename TYPE>
static inline TYPE cos4pi5(TYPE a) { return std::cos(typename TYPE::value_type(4 * M_PI / 5)) * a; }
template <typename TYPE>
static inline TYPE sin4pi5(TYPE a) { return std::sin(typename TYPE::value_type(4 * M_PI / 5)) * a; }

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

static constexpr int split(int N)
{
	return (!(N%5)) ? 5 : (!(N%3)) ? 3 : (!(N%2)&&pow4(N)) ? 4 : (!(N%2)) ? 2 : 1;
}

template <int RADIX, int N, int S, typename TYPE, int SIGN>
struct Dit {};

template <int S, typename TYPE, int SIGN>
struct Dit<1, 1, S, TYPE, SIGN>
{
	static inline void dit(TYPE *out, TYPE *in, TYPE *)
	{
		*out = *in;
	}
};

template <int S, typename TYPE, int SIGN>
struct Dit<2, 2, S, TYPE, SIGN>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE in0, TYPE in1)
	{
		*out0 = in0 + in1;
		*out1 = in0 - in1;
	}
	static inline void dit(TYPE *out, TYPE *in, TYPE *)
	{
		dft(out, out + 1, in[0], in[S]);
	}
};

template <int S, typename TYPE>
struct Dit<3, 3, S, TYPE, -1>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE *out2, TYPE in0, TYPE in1, TYPE in2)
	{
		TYPE a(in1 + in2), b(sqrt3(twiddle(in2, in1)));
		*out0 = in0 + a;
		*out1 = in0 - half(a + b);
		*out2 = in0 - half(a - b);
	}
	static inline void dit(TYPE *out, TYPE *in, TYPE *)
	{
		dft(out, out + 1, out + 2, in[0], in[S], in[2 * S]);
	}
};

template <int S, typename TYPE>
struct Dit<3, 3, S, TYPE, 1>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE *out2, TYPE in0, TYPE in1, TYPE in2)
	{
		TYPE a(in1 + in2), b(sqrt3(twiddle(in2, in1)));
		*out0 = in0 + a;
		*out1 = in0 - half(a - b);
		*out2 = in0 - half(a + b);
	}
	static inline void dit(TYPE *out, TYPE *in, TYPE *)
	{
		dft(out, out + 1, out + 2, in[0], in[S], in[2 * S]);
	}
};

template <int S, typename TYPE>
struct Dit<4, 4, S, TYPE, -1>
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
	static inline void dit(TYPE *out, TYPE *in, TYPE *)
	{
		dft(out, out + 1, out + 2, out + 3, in[0], in[S], in[2 * S], in[3 * S]);
	}
};

template <int S, typename TYPE>
struct Dit<4, 4, S, TYPE, 1>
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
	static inline void dit(TYPE *out, TYPE *in, TYPE *)
	{
		dft(out, out + 1, out + 2, out + 3, in[0], in[S], in[2 * S], in[3 * S]);
	}
};

template <int S, typename TYPE>
struct Dit<5, 5, S, TYPE, -1>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE *out2, TYPE *out3, TYPE *out4,
			TYPE in0, TYPE in1, TYPE in2, TYPE in3, TYPE in4)
	{
		TYPE a(in1 + in4), b(in2 + in3), c(twiddle(in1, in4)), d(twiddle(in2, in3));
		TYPE c1a(cos2pi5(a)), c2b(cos4pi5(b)), s1c(sin2pi5(c)), s2d(sin4pi5(d));
		TYPE c2a(cos4pi5(a)), c1b(cos2pi5(b)), s2c(sin4pi5(c)), s1d(sin2pi5(d));
		*out0 = in0 + a + b;
		*out1 = in0 + c1a + c2b + s1c + s2d;
		*out2 = in0 + c2a + c1b + s2c - s1d;
		*out3 = in0 + c2a + c1b - s2c + s1d;
		*out4 = in0 + c1a + c2b - s1c - s2d;
	}
	static inline void dit(TYPE *out, TYPE *in, TYPE *)
	{
		dft(out, out + 1, out + 2, out + 3, out + 4, in[0], in[S], in[2 * S], in[3 * S], in[4 * S]);
	}
};

template <int S, typename TYPE>
struct Dit<5, 5, S, TYPE, 1>
{
	static inline void dft(TYPE *out0, TYPE *out1, TYPE *out2, TYPE *out3, TYPE *out4,
			TYPE in0, TYPE in1, TYPE in2, TYPE in3, TYPE in4)
	{
		TYPE a(in1 + in4), b(in2 + in3), c(twiddle(in1, in4)), d(twiddle(in2, in3));
		TYPE c1a(cos2pi5(a)), c2b(cos4pi5(b)), s1c(sin2pi5(c)), s2d(sin4pi5(d));
		TYPE c2a(cos4pi5(a)), c1b(cos2pi5(b)), s2c(sin4pi5(c)), s1d(sin2pi5(d));
		*out0 = in0 + a + b;
		*out1 = in0 + c1a + c2b - s1c - s2d;
		*out2 = in0 + c2a + c1b - s2c + s1d;
		*out3 = in0 + c2a + c1b + s2c - s1d;
		*out4 = in0 + c1a + c2b + s1c + s2d;
	}
	static inline void dit(TYPE *out, TYPE *in, TYPE *)
	{
		dft(out, out + 1, out + 2, out + 3, out + 4, in[0], in[S], in[2 * S], in[3 * S], in[4 * S]);
	}
};

template <int N, int S, typename TYPE, int SIGN>
struct Dit<2, N, S, TYPE, SIGN>
{
	static void dit(TYPE *out, TYPE *in, TYPE *z)
	{
		for (int o = 0, i = 0; o < N; o += N / 2, i += S)
			Dit<split(N / 2), N / 2, 2 * S, TYPE, SIGN>::ditd(out + o, in + i, z);
		for (int k0 = 0, k1 = N / 2; k0 < N / 2; ++k0, ++k1)
			Dit<2, 2, S, TYPE, SIGN>::dft(out + k0, out + k1, out[k0], z[k0 * S] * out[k1]);
	}
};

template <int N, int S, typename TYPE, int SIGN>
struct Dit<3, N, S, TYPE, SIGN>
{
	static void dit(TYPE *out, TYPE *in, TYPE *z)
	{
		for (int o = 0, i = 0; o < N; o += N / 3, i += S)
			Dit<split(N / 3), N / 3, 3 * S, TYPE, SIGN>::dit(out + o, in + i, z);
		for (int k0 = 0, k1 = N / 3, k2 = 2 * N / 3; k0 < N / 3; ++k0, ++k1, ++k2)
			Dit<3, 3, S, TYPE, SIGN>::dft(out + k0, out + k1, out + k2,
				out[k0], z[k0 * S] * out[k1], z[2 * k0 * S] * out[k2]);
	}
};

template <int N, int S, typename TYPE, int SIGN>
struct Dit<4, N, S, TYPE, SIGN>
{
	static void dit(TYPE *out, TYPE *in, TYPE *z)
	{
		for (int o = 0, i = 0; o < N; o += N / 4, i += S)
			Dit<split(N / 4), N / 4, 4 * S, TYPE, SING>::dit(out + o, in + i, z);
		for (int k0 = 0, k1 = N / 4, k2 = 2 * N / 4, k3 = 3 * N / 4;
				k0 < N / 4; ++k0, ++k1, ++k2, ++k3)
			Dit<4, 4, S, TYPE, SIGN>::dft(out + k0, out + k1, out + k2, out + k3,
				out[k0], z[k0 * S] * out[k1], z[2 * k0 * S] * out[k2],
				z[3 * k0 * S] * out[k3]);
	}
};

template <int N, int S, typename TYPE, int SIGN>
struct Dit<5, N, S, TYPE, SIGN>
{
	static void dit(TYPE *out, TYPE *in, TYPE *z)
	{
		for (int o = 0, i = 0; o < N; o += N / 5, i += S)
			Dit<split(N / 5), N / 5, 5 * S, TYPE, SIGN>::dit(out + o, in + i, z);
		for (int k0 = 0, k1 = N / 5, k2 = 2 * N / 5, k3 = 3 * N / 5, k4 = 4 * N / 5;
				k0 < N / 5; ++k0, ++k1, ++k2, ++k3, ++k4)
			Dit<5, 5, S, TYPE, SIGN>::dft(out + k0, out + k1, out + k2, out + k3, out + k4,
				out[k0], z[k0 * S] * out[k1], z[2 * k0 * S] * out[k2],
				z[3 * k0 * S] * out[k3], z[4 * k0 * S] * out[k4]);
	}
};

template <int BINS, typename TYPE, int SIGN>
struct Factors
{
	typedef typename TYPE::value_type value_type;
	TYPE z[BINS];
	Factors()
	{
		for (int n = 0; n < BINS; ++n)
			z[n] = exp(TYPE(0, value_type(SIGN * 2 * M_PI) * value_type(n) / value_type(BINS)));
	}
};

template <int BINS, typename TYPE>
class Forward
{
	Factors<BINS, TYPE, -1> factors;
public:
	typedef typename TYPE::value_type value_type;
	inline void operator ()(TYPE *out, TYPE *in)
	{
		Dit<split(BINS), BINS, 1, TYPE, -1>::dit(out, in, factors.z);
	}
};

template <int BINS, typename TYPE>
class Backward
{
	Factors<BINS, TYPE, 1> factors;
public:
	typedef typename TYPE::value_type value_type;
	inline void operator ()(TYPE *out, TYPE *in)
	{
		Dit<split(BINS), BINS, 1, TYPE, 1>::dit(out, in, factors.z);
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
