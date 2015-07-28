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
static inline TYPE twiddle(TYPE a, TYPE b)
{
	return TYPE(a.imag() - b.imag(), b.real() - a.real());
}

template <typename TYPE>
static inline void fwd2(TYPE *out0, TYPE *out1, TYPE in0, TYPE in1)
{
	*out0 = in0 + in1;
	*out1 = in0 - in1;
}

template <typename TYPE>
static inline void bwd2(TYPE *out0, TYPE *out1, TYPE in0, TYPE in1)
{
	*out0 = in0 + in1;
	*out1 = in0 - in1;
}

template <typename TYPE>
static inline void fwd3(TYPE *out0, TYPE *out1, TYPE *out2, TYPE in0, TYPE in1, TYPE in2)
{
	typedef typename TYPE::value_type value_type;
	static const value_type half(value_type(1) / value_type(2)), sqrt3(std::sqrt(value_type(3)));
	TYPE a(in1 + in2), b(sqrt3 * twiddle(in2, in1));
	*out0 = in0 + a;
	*out1 = in0 - half * (a + b);
	*out2 = in0 - half * (a - b);
}

template <typename TYPE>
static inline void bwd3(TYPE *out0, TYPE *out1, TYPE *out2, TYPE in0, TYPE in1, TYPE in2)
{
	typedef typename TYPE::value_type value_type;
	static const value_type half(value_type(1) / value_type(2)), sqrt3(std::sqrt(value_type(3)));
	TYPE a(in1 + in2), b(sqrt3 * twiddle(in2, in1));
	*out0 = in0 + a;
	*out1 = in0 - half * (a - b);
	*out2 = in0 - half * (a + b);
}

template <typename TYPE>
static inline void fwd4(TYPE *out0, TYPE *out1, TYPE *out2, TYPE *out3,
		TYPE in0, TYPE in1, TYPE in2, TYPE in3)
{
	TYPE a(in0 + in2), b(in0 - in2), c(in1 + in3), d(twiddle(in1, in3));
	*out0 = a + c;
	*out1 = b + d;
	*out2 = a - c;
	*out3 = b - d;
}

template <typename TYPE>
static inline void bwd4(TYPE *out0, TYPE *out1, TYPE *out2, TYPE *out3,
		TYPE in0, TYPE in1, TYPE in2, TYPE in3)
{
	TYPE a(in0 + in2), b(in0 - in2), c(in1 + in3), d(twiddle(in1, in3));
	*out0 = a + c;
	*out1 = b - d;
	*out2 = a - c;
	*out3 = b + d;
}

template <typename TYPE>
static inline void fwd5(TYPE *out0, TYPE *out1, TYPE *out2, TYPE *out3, TYPE *out4,
		TYPE in0, TYPE in1, TYPE in2, TYPE in3, TYPE in4)
{
	typedef typename TYPE::value_type value_type;
	static const value_type r1(2 * M_PI / 5), c1(std::cos(r1)), s1(std::sin(r1));
	static const value_type r2(4 * M_PI / 5), c2(std::cos(r2)), s2(std::sin(r2));
	TYPE a(in1 + in4), b(in2 + in3), c(twiddle(in1, in4)), d(twiddle(in2, in3));
	TYPE c1a(c1 * a), c2b(c2 * b), s1c(s1 * c), s2d(s2 * d);
	TYPE c2a(c2 * a), c1b(c1 * b), s2c(s2 * c), s1d(s1 * d);
	*out0 = in0 + a + b;
	*out1 = in0 + c1a + c2b + s1c + s2d;
	*out2 = in0 + c2a + c1b + s2c - s1d;
	*out3 = in0 + c2a + c1b - s2c + s1d;
	*out4 = in0 + c1a + c2b - s1c - s2d;
}

template <typename TYPE>
static inline void bwd5(TYPE *out0, TYPE *out1, TYPE *out2, TYPE *out3, TYPE *out4,
		TYPE in0, TYPE in1, TYPE in2, TYPE in3, TYPE in4)
{
	typedef typename TYPE::value_type value_type;
	static const value_type r1(2 * M_PI / 5), c1(std::cos(r1)), s1(std::sin(r1));
	static const value_type r2(4 * M_PI / 5), c2(std::cos(r2)), s2(std::sin(r2));
	TYPE a(in1 + in4), b(in2 + in3), c(twiddle(in1, in4)), d(twiddle(in2, in3));
	TYPE c1a(c1 * a), c2b(c2 * b), s1c(s1 * c), s2d(s2 * d);
	TYPE c2a(c2 * a), c1b(c1 * b), s2c(s2 * c), s1d(s1 * d);
	*out0 = in0 + a + b;
	*out1 = in0 + c1a + c2b - s1c - s2d;
	*out2 = in0 + c2a + c1b - s2c + s1d;
	*out3 = in0 + c2a + c1b + s2c - s1d;
	*out4 = in0 + c1a + c2b + s1c + s2d;
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

template <int RADIX, int N, int S, typename TYPE>
struct DitFwd {};

template <int S, typename TYPE>
struct DitFwd<1, 1, S, TYPE>
{
	static inline void fwd(TYPE *out, TYPE *in, TYPE *)
	{
		*out = *in;
	}
};

template <int S, typename TYPE>
struct DitFwd<2, 2, S, TYPE>
{
	static inline void fwd(TYPE *out, TYPE *in, TYPE *)
	{
		fwd2(out, out + 1, in[0], in[S]);
	}
};

template <int S, typename TYPE>
struct DitFwd<3, 3, S, TYPE>
{
	static inline void fwd(TYPE *out, TYPE *in, TYPE *)
	{
		fwd3(out, out + 1, out + 2, in[0], in[S], in[2 * S]);
	}
};

template <int S, typename TYPE>
struct DitFwd<4, 4, S, TYPE>
{
	static inline void fwd(TYPE *out, TYPE *in, TYPE *)
	{
		fwd4(out, out + 1, out + 2, out + 3, in[0], in[S], in[2 * S], in[3 * S]);
	}
};

template <int S, typename TYPE>
struct DitFwd<5, 5, S, TYPE>
{
	static inline void fwd(TYPE *out, TYPE *in, TYPE *)
	{
		fwd5(out, out + 1, out + 2, out + 3, out + 4, in[0], in[S], in[2 * S], in[3 * S], in[4 * S]);
	}
};

template <int N, int S, typename TYPE>
struct DitFwd<2, N, S, TYPE>
{
	static void fwd(TYPE *out, TYPE *in, TYPE *z)
	{
		typedef DitFwd<split(N / 2), N / 2, 2 * S, TYPE> dit;
		dit::fwd(out, in, z);
		dit::fwd(out + N / 2, in + S, z);
		for (int k0 = 0, k1 = N / 2; k0 < N / 2; ++k0, ++k1)
			fwd2(out + k0, out + k1, out[k0], z[k0 * S] * out[k1]);
	}
};

template <int N, int S, typename TYPE>
struct DitFwd<3, N, S, TYPE>
{
	static void fwd(TYPE *out, TYPE *in, TYPE *z)
	{
		typedef DitFwd<split(N / 3), N / 3, 3 * S, TYPE> dit;
		dit::fwd(out, in, z);
		dit::fwd(out + N / 3, in + S, z);
		dit::fwd(out + 2 * N / 3, in + 2 * S, z);
		for (int k0 = 0, k1 = N / 3, k2 = 2 * N / 3; k0 < N / 3; ++k0, ++k1, ++k2)
			fwd3(out + k0, out + k1, out + k2,
				out[k0], z[k0 * S] * out[k1], z[2 * k0 * S] * out[k2]);
	}
};

template <int N, int S, typename TYPE>
struct DitFwd<4, N, S, TYPE>
{
	static void fwd(TYPE *out, TYPE *in, TYPE *z)
	{
		typedef DitFwd<split(N / 4), N / 4, 4 * S, TYPE> dit;
		dit::fwd(out, in, z);
		dit::fwd(out + N / 4, in + S, z);
		dit::fwd(out + 2 * N / 4, in + 2 * S, z);
		dit::fwd(out + 3 * N / 4, in + 3 * S, z);
		for (int k0 = 0, k1 = N / 4, k2 = 2 * N / 4, k3 = 3 * N / 4;
				k0 < N / 4; ++k0, ++k1, ++k2, ++k3)
			fwd4(out + k0, out + k1, out + k2, out + k3,
				out[k0], z[k0 * S] * out[k1], z[2 * k0 * S] * out[k2],
				z[3 * k0 * S] * out[k3]);
	}
};

template <int N, int S, typename TYPE>
struct DitFwd<5, N, S, TYPE>
{
	static void fwd(TYPE *out, TYPE *in, TYPE *z)
	{
		typedef DitFwd<split(N / 5), N / 5, 5 * S, TYPE> dit;
		dit::fwd(out, in, z);
		dit::fwd(out + N / 5, in + S, z);
		dit::fwd(out + 2 * N / 5, in + 2 * S, z);
		dit::fwd(out + 3 * N / 5, in + 3 * S, z);
		dit::fwd(out + 4 * N / 5, in + 4 * S, z);
		for (int k0 = 0, k1 = N / 5, k2 = 2 * N / 5, k3 = 3 * N / 5, k4 = 4 * N / 5;
				k0 < N / 5; ++k0, ++k1, ++k2, ++k3, ++k4)
			fwd5(out + k0, out + k1, out + k2, out + k3, out + k4,
				out[k0], z[k0 * S] * out[k1], z[2 * k0 * S] * out[k2],
				z[3 * k0 * S] * out[k3], z[4 * k0 * S] * out[k4]);
	}
};

template <int RADIX, int N, int S, typename TYPE>
struct DitBwd {};

template <int S, typename TYPE>
struct DitBwd<1, 1, S, TYPE>
{
	static inline void bwd(TYPE *out, TYPE *in, TYPE *)
	{
		*out = *in;
	}
};

template <int S, typename TYPE>
struct DitBwd<2, 2, S, TYPE>
{
	static inline void bwd(TYPE *out, TYPE *in, TYPE *) {
		bwd2(out, out + 1, in[0], in[S]);
	}
};

template <int S, typename TYPE>
struct DitBwd<3, 3, S, TYPE>
{
	static inline void bwd(TYPE *out, TYPE *in, TYPE *) {
		bwd3(out, out + 1, out + 2, in[0], in[S], in[2 * S]);
	}
};

template <int S, typename TYPE>
struct DitBwd<4, 4, S, TYPE>
{
	static inline void bwd(TYPE *out, TYPE *in, TYPE *) {
		bwd4(out, out + 1, out + 2, out + 3, in[0], in[S], in[2 * S], in[3 * S]);
	}
};

template <int S, typename TYPE>
struct DitBwd<5, 5, S, TYPE>
{
	static inline void bwd(TYPE *out, TYPE *in, TYPE *) {
		bwd5(out, out + 1, out + 2, out + 3, out + 4, in[0], in[S], in[2 * S], in[3 * S], in[4 * S]);
	}
};

template <int N, int S, typename TYPE>
struct DitBwd<2, N, S, TYPE>
{
	static void bwd(TYPE *out, TYPE *in, TYPE *z)
	{
		typedef DitBwd<split(N / 2), N / 2, 2 * S, TYPE> dit;
		dit::bwd(out, in, z);
		dit::bwd(out + N / 2, in + S, z);
		for (int k0 = 0, k1 = N / 2; k0 < N / 2; ++k0, ++k1)
			bwd2(out + k0, out + k1, out[k0], z[k0 * S] * out[k1]);
	}
};

template <int N, int S, typename TYPE>
struct DitBwd<3, N, S, TYPE>
{
	static void bwd(TYPE *out, TYPE *in, TYPE *z)
	{
		typedef DitBwd<split(N / 3), N / 3, 3 * S, TYPE> dit;
		dit::bwd(out, in, z);
		dit::bwd(out + N / 3, in + S, z);
		dit::bwd(out + 2 * N / 3, in + 2 * S, z);
		for (int k0 = 0, k1 = N / 3, k2 = 2 * N / 3; k0 < N / 3; ++k0, ++k1, ++k2)
			bwd3(out + k0, out + k1, out + k2,
				out[k0], z[k0 * S] * out[k1], z[2 * k0 * S] * out[k2]);
	}
};

template <int N, int S, typename TYPE>
struct DitBwd<4, N, S, TYPE>
{
	static void bwd(TYPE *out, TYPE *in, TYPE *z)
	{
		typedef DitBwd<split(N / 4), N / 4, 4 * S, TYPE> dit;
		dit::bwd(out, in, z);
		dit::bwd(out + N / 4, in + S, z);
		dit::bwd(out + 2 * N / 4, in + 2 * S, z);
		dit::bwd(out + 3 * N / 4, in + 3 * S, z);
		for (int k0 = 0, k1 = N / 4, k2 = 2 * N / 4, k3 = 3 * N / 4;
				k0 < N / 4; ++k0, ++k1, ++k2, ++k3)
			bwd4(out + k0, out + k1, out + k2, out + k3,
				out[k0], z[k0 * S] * out[k1], z[2 * k0 * S] * out[k2],
				z[3 * k0 * S] * out[k3]);
	}
};

template <int N, int S, typename TYPE>
struct DitBwd<5, N, S, TYPE>
{
	static void bwd(TYPE *out, TYPE *in, TYPE *z)
	{
		typedef DitBwd<split(N / 5), N / 5, 5 * S, TYPE> dit;
		dit::bwd(out, in, z);
		dit::bwd(out + N / 5, in + S, z);
		dit::bwd(out + 2 * N / 5, in + 2 * S, z);
		dit::bwd(out + 3 * N / 5, in + 3 * S, z);
		dit::bwd(out + 4 * N / 5, in + 4 * S, z);
		for (int k0 = 0, k1 = N / 5, k2 = 2 * N / 5, k3 = 3 * N / 5, k4 = 4 * N / 5;
				k0 < N / 5; ++k0, ++k1, ++k2, ++k3, ++k4)
			bwd5(out + k0, out + k1, out + k2, out + k3, out + k4,
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
		DitFwd<split(BINS), BINS, 1, TYPE>::fwd(out, in, factors.z);
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
		DitBwd<split(BINS), BINS, 1, TYPE>::bwd(out, in, factors.z);
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
