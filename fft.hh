/*
fft - mixed radix fft
Written in 2015 by <Ahmet Inan> <xdsopl@googlemail.com>
To the extent possible under law, the author(s) have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide. This software is distributed without any warranty.
You should have received a copy of the CC0 Public Domain Dedication along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/

#ifndef FFT_HH
#define FFT_HH

#include <complex>
#include <cmath>

namespace FFT {

template <typename TYPE>
static inline std::complex<TYPE> twiddle(std::complex<TYPE> a, std::complex<TYPE> b)
{
	return std::complex<TYPE>(a.imag() - b.imag(), b.real() - a.real());
}

template <typename TYPE>
static inline void fwd2(std::complex<TYPE> *out0, std::complex<TYPE> *out1,
		std::complex<TYPE> in0, std::complex<TYPE> in1)
{
	*out0 = in0 + in1;
	*out1 = in0 - in1;
}

template <typename TYPE>
static inline void bwd2(std::complex<TYPE> *out0, std::complex<TYPE> *out1,
		std::complex<TYPE> in0, std::complex<TYPE> in1)
{
	*out0 = in0 + in1;
	*out1 = in0 - in1;
}

template <typename TYPE>
static inline void fwd3(std::complex<TYPE> *out0, std::complex<TYPE> *out1, std::complex<TYPE> *out2,
		std::complex<TYPE> in0, std::complex<TYPE> in1, std::complex<TYPE> in2)
{
	static const TYPE half(TYPE(1) / TYPE(2)), sqrt3(std::sqrt(TYPE(3)));
	std::complex<TYPE> a(in1 + in2), b(sqrt3 * twiddle(in2, in1));
	*out0 = in0 + a;
	*out1 = in0 - half * (a + b);
	*out2 = in0 - half * (a - b);
}

template <typename TYPE>
static inline void bwd3(std::complex<TYPE> *out0, std::complex<TYPE> *out1, std::complex<TYPE> *out2,
		std::complex<TYPE> in0, std::complex<TYPE> in1, std::complex<TYPE> in2)
{
	static const TYPE half(TYPE(1) / TYPE(2)), sqrt3(std::sqrt(TYPE(3)));
	std::complex<TYPE> a(in1 + in2), b(sqrt3 * twiddle(in2, in1));
	*out0 = in0 + a;
	*out1 = in0 - half * (a - b);
	*out2 = in0 - half * (a + b);
}

template <typename TYPE>
static inline void fwd4(std::complex<TYPE> *out0, std::complex<TYPE> *out1, std::complex<TYPE> *out2, std::complex<TYPE> *out3,
		std::complex<TYPE> in0, std::complex<TYPE> in1, std::complex<TYPE> in2, std::complex<TYPE> in3)
{
	std::complex<TYPE> a(in0 + in2), b(in0 - in2), c(in1 + in3), d(twiddle(in1, in3));
	*out0 = a + c;
	*out1 = b + d;
	*out2 = a - c;
	*out3 = b - d;
}

template <typename TYPE>
static inline void bwd4(std::complex<TYPE> *out0, std::complex<TYPE> *out1, std::complex<TYPE> *out2, std::complex<TYPE> *out3,
		std::complex<TYPE> in0, std::complex<TYPE> in1, std::complex<TYPE> in2, std::complex<TYPE> in3)
{
	std::complex<TYPE> a(in0 + in2), b(in0 - in2), c(in1 + in3), d(twiddle(in1, in3));
	*out0 = a + c;
	*out1 = b - d;
	*out2 = a - c;
	*out3 = b + d;
}

template <typename TYPE>
static inline void fwd5(std::complex<TYPE> *out0, std::complex<TYPE> *out1, std::complex<TYPE> *out2, std::complex<TYPE> *out3, std::complex<TYPE> *out4,
		std::complex<TYPE> in0, std::complex<TYPE> in1, std::complex<TYPE> in2, std::complex<TYPE> in3, std::complex<TYPE> in4)
{
	static const TYPE r1(2 * M_PI / 5), c1(std::cos(r1)), s1(std::sin(r1));
	static const TYPE r2(4 * M_PI / 5), c2(std::cos(r2)), s2(std::sin(r2));
	std::complex<TYPE> a(in1 + in4), b(in2 + in3), c(twiddle(in1, in4)), d(twiddle(in2, in3));
	std::complex<TYPE> c1a(c1 * a), c2b(c2 * b), s1c(s1 * c), s2d(s2 * d);
	std::complex<TYPE> c2a(c2 * a), c1b(c1 * b), s2c(s2 * c), s1d(s1 * d);
	*out0 = in0 + a + b;
	*out1 = in0 + c1a + c2b + s1c + s2d;
	*out2 = in0 + c2a + c1b + s2c - s1d;
	*out3 = in0 + c2a + c1b - s2c + s1d;
	*out4 = in0 + c1a + c2b - s1c - s2d;
}

template <typename TYPE>
static inline void bwd5(std::complex<TYPE> *out0, std::complex<TYPE> *out1, std::complex<TYPE> *out2, std::complex<TYPE> *out3, std::complex<TYPE> *out4,
		std::complex<TYPE> in0, std::complex<TYPE> in1, std::complex<TYPE> in2, std::complex<TYPE> in3, std::complex<TYPE> in4)
{
	static const TYPE r1(2 * M_PI / 5), c1(std::cos(r1)), s1(std::sin(r1));
	static const TYPE r2(4 * M_PI / 5), c2(std::cos(r2)), s2(std::sin(r2));
	std::complex<TYPE> a(in1 + in4), b(in2 + in3), c(twiddle(in1, in4)), d(twiddle(in2, in3));
	std::complex<TYPE> c1a(c1 * a), c2b(c2 * b), s1c(s1 * c), s2d(s2 * d);
	std::complex<TYPE> c2a(c2 * a), c1b(c1 * b), s2c(s2 * c), s1d(s1 * d);
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

template <int RADIX, int N, int S, typename FACTORS, typename TYPE>
struct DitFwd {};

template <int S, typename FACTORS, typename TYPE>
struct DitFwd<1, 1, S, FACTORS, TYPE>
{
	static inline void fwd(std::complex<TYPE> *out, std::complex<TYPE> *in)
	{
		*out = *in;
	}
};

template <int S, typename FACTORS, typename TYPE>
struct DitFwd<2, 2, S, FACTORS, TYPE>
{
	static inline void fwd(std::complex<TYPE> *out, std::complex<TYPE> *in)
	{
		fwd2(out, out + 1, in[0], in[S]);
	}
};

template <int S, typename FACTORS, typename TYPE>
struct DitFwd<3, 3, S, FACTORS, TYPE>
{
	static inline void fwd(std::complex<TYPE> *out, std::complex<TYPE> *in)
	{
		fwd3(out, out + 1, out + 2, in[0], in[S], in[2 * S]);
	}
};

template <int S, typename FACTORS, typename TYPE>
struct DitFwd<4, 4, S, FACTORS, TYPE>
{
	static inline void fwd(std::complex<TYPE> *out, std::complex<TYPE> *in)
	{
		fwd4(out, out + 1, out + 2, out + 3, in[0], in[S], in[2 * S], in[3 * S]);
	}
};

template <int S, typename FACTORS, typename TYPE>
struct DitFwd<5, 5, S, FACTORS, TYPE>
{
	static inline void fwd(std::complex<TYPE> *out, std::complex<TYPE> *in)
	{
		fwd5(out, out + 1, out + 2, out + 3, out + 4, in[0], in[S], in[2 * S], in[3 * S], in[4 * S]);
	}
};

template <int N, int S, typename FACTORS, typename TYPE>
struct DitFwd<2, N, S, FACTORS, TYPE>
{
	static void fwd(std::complex<TYPE> *out, std::complex<TYPE> *in)
	{
		typedef DitFwd<split(N / 2), N / 2, 2 * S, FACTORS, TYPE> dit;
		FACTORS z;
		dit::fwd(out, in);
		dit::fwd(out + N / 2, in + S);
		for (int k0 = 0, k1 = N / 2; k0 < N / 2; ++k0, ++k1)
			fwd2(out + k0, out + k1, out[k0], z[k0 * S] * out[k1]);
	}
};

template <int N, int S, typename FACTORS, typename TYPE>
struct DitFwd<3, N, S, FACTORS, TYPE>
{
	static void fwd(std::complex<TYPE> *out, std::complex<TYPE> *in)
	{
		typedef DitFwd<split(N / 3), N / 3, 3 * S, FACTORS, TYPE> dit;
		FACTORS z;
		dit::fwd(out, in);
		dit::fwd(out + N / 3, in + S);
		dit::fwd(out + 2 * N / 3, in + 2 * S);
		for (int k0 = 0, k1 = N / 3, k2 = 2 * N / 3; k0 < N / 3; ++k0, ++k1, ++k2)
			fwd3(out + k0, out + k1, out + k2,
				out[k0], z[k0 * S] * out[k1], z[2 * k0 * S] * out[k2]);
	}
};

template <int N, int S, typename FACTORS, typename TYPE>
struct DitFwd<4, N, S, FACTORS, TYPE>
{
	static void fwd(std::complex<TYPE> *out, std::complex<TYPE> *in)
	{
		typedef DitFwd<split(N / 4), N / 4, 4 * S, FACTORS, TYPE> dit;
		FACTORS z;
		dit::fwd(out, in);
		dit::fwd(out + N / 4, in + S);
		dit::fwd(out + 2 * N / 4, in + 2 * S);
		dit::fwd(out + 3 * N / 4, in + 3 * S);
		for (int k0 = 0, k1 = N / 4, k2 = 2 * N / 4, k3 = 3 * N / 4;
				k0 < N / 4; ++k0, ++k1, ++k2, ++k3)
			fwd4(out + k0, out + k1, out + k2, out + k3,
				out[k0], z[k0 * S] * out[k1], z[2 * k0 * S] * out[k2],
				z[3 * k0 * S] * out[k3]);
	}
};

template <int N, int S, typename FACTORS, typename TYPE>
struct DitFwd<5, N, S, FACTORS, TYPE>
{
	static void fwd(std::complex<TYPE> *out, std::complex<TYPE> *in)
	{
		typedef DitFwd<split(N / 5), N / 5, 5 * S, FACTORS, TYPE> dit;
		FACTORS z;
		dit::fwd(out, in);
		dit::fwd(out + N / 5, in + S);
		dit::fwd(out + 2 * N / 5, in + 2 * S);
		dit::fwd(out + 3 * N / 5, in + 3 * S);
		dit::fwd(out + 4 * N / 5, in + 4 * S);
		for (int k0 = 0, k1 = N / 5, k2 = 2 * N / 5, k3 = 3 * N / 5, k4 = 4 * N / 5;
				k0 < N / 5; ++k0, ++k1, ++k2, ++k3, ++k4)
			fwd5(out + k0, out + k1, out + k2, out + k3, out + k4,
				out[k0], z[k0 * S] * out[k1], z[2 * k0 * S] * out[k2],
				z[3 * k0 * S] * out[k3], z[4 * k0 * S] * out[k4]);
	}
};

template <int RADIX, int N, int S, typename FACTORS, typename TYPE>
struct DitBwd {};

template <int S, typename FACTORS, typename TYPE>
struct DitBwd<1, 1, S, FACTORS, TYPE>
{
	static inline void bwd(std::complex<TYPE> *out, std::complex<TYPE> *in)
	{
		*out = *in;
	}
};

template <int S, typename FACTORS, typename TYPE>
struct DitBwd<2, 2, S, FACTORS, TYPE>
{
	static inline void bwd(std::complex<TYPE> *out, std::complex<TYPE> *in) {
		bwd2(out, out + 1, in[0], in[S]);
	}
};

template <int S, typename FACTORS, typename TYPE>
struct DitBwd<3, 3, S, FACTORS, TYPE>
{
	static inline void bwd(std::complex<TYPE> *out, std::complex<TYPE> *in) {
		bwd3(out, out + 1, out + 2, in[0], in[S], in[2 * S]);
	}
};

template <int S, typename FACTORS, typename TYPE>
struct DitBwd<4, 4, S, FACTORS, TYPE>
{
	static inline void bwd(std::complex<TYPE> *out, std::complex<TYPE> *in) {
		bwd4(out, out + 1, out + 2, out + 3, in[0], in[S], in[2 * S], in[3 * S]);
	}
};

template <int S, typename FACTORS, typename TYPE>
struct DitBwd<5, 5, S, FACTORS, TYPE>
{
	static inline void bwd(std::complex<TYPE> *out, std::complex<TYPE> *in) {
		bwd5(out, out + 1, out + 2, out + 3, out + 4, in[0], in[S], in[2 * S], in[3 * S], in[4 * S]);
	}
};

template <int N, int S, typename FACTORS, typename TYPE>
struct DitBwd<2, N, S, FACTORS, TYPE>
{
	static void bwd(std::complex<TYPE> *out, std::complex<TYPE> *in)
	{
		typedef DitBwd<split(N / 2), N / 2, 2 * S, FACTORS, TYPE> dit;
		FACTORS z;
		dit::bwd(out, in);
		dit::bwd(out + N / 2, in + S);
		for (int k0 = 0, k1 = N / 2; k0 < N / 2; ++k0, ++k1)
			bwd2(out + k0, out + k1, out[k0], z[k0 * S] * out[k1]);
	}
};

template <int N, int S, typename FACTORS, typename TYPE>
struct DitBwd<3, N, S, FACTORS, TYPE>
{
	static void bwd(std::complex<TYPE> *out, std::complex<TYPE> *in)
	{
		typedef DitBwd<split(N / 3), N / 3, 3 * S, FACTORS, TYPE> dit;
		FACTORS z;
		dit::bwd(out, in);
		dit::bwd(out + N / 3, in + S);
		dit::bwd(out + 2 * N / 3, in + 2 * S);
		for (int k0 = 0, k1 = N / 3, k2 = 2 * N / 3; k0 < N / 3; ++k0, ++k1, ++k2)
			bwd3(out + k0, out + k1, out + k2,
				out[k0], z[k0 * S] * out[k1], z[2 * k0 * S] * out[k2]);
	}
};

template <int N, int S, typename FACTORS, typename TYPE>
struct DitBwd<4, N, S, FACTORS, TYPE>
{
	static void bwd(std::complex<TYPE> *out, std::complex<TYPE> *in)
	{
		typedef DitBwd<split(N / 4), N / 4, 4 * S, FACTORS, TYPE> dit;
		FACTORS z;
		dit::bwd(out, in);
		dit::bwd(out + N / 4, in + S);
		dit::bwd(out + 2 * N / 4, in + 2 * S);
		dit::bwd(out + 3 * N / 4, in + 3 * S);
		for (int k0 = 0, k1 = N / 4, k2 = 2 * N / 4, k3 = 3 * N / 4;
				k0 < N / 4; ++k0, ++k1, ++k2, ++k3)
			bwd4(out + k0, out + k1, out + k2, out + k3,
				out[k0], z[k0 * S] * out[k1], z[2 * k0 * S] * out[k2],
				z[3 * k0 * S] * out[k3]);
	}
};

template <int N, int S, typename FACTORS, typename TYPE>
struct DitBwd<5, N, S, FACTORS, TYPE>
{
	static void bwd(std::complex<TYPE> *out, std::complex<TYPE> *in)
	{
		typedef DitBwd<split(N / 5), N / 5, 5 * S, FACTORS, TYPE> dit;
		FACTORS z;
		dit::bwd(out, in);
		dit::bwd(out + N / 5, in + S);
		dit::bwd(out + 2 * N / 5, in + 2 * S);
		dit::bwd(out + 3 * N / 5, in + 3 * S);
		dit::bwd(out + 4 * N / 5, in + 4 * S);
		for (int k0 = 0, k1 = N / 5, k2 = 2 * N / 5, k3 = 3 * N / 5, k4 = 4 * N / 5;
				k0 < N / 5; ++k0, ++k1, ++k2, ++k3, ++k4)
			bwd5(out + k0, out + k1, out + k2, out + k3, out + k4,
				out[k0], z[k0 * S] * out[k1], z[2 * k0 * S] * out[k2],
				z[3 * k0 * S] * out[k3], z[4 * k0 * S] * out[k4]);
	}
};

template <int BINS, typename TYPE, int SIGN>
class Factors
{
	static const struct StaticFactors
	{
		std::complex<TYPE> z[BINS];
		StaticFactors()
		{
			for (int n = 0; n < BINS; ++n)
				z[n] = exp(std::complex<TYPE>(0, TYPE(SIGN * 2 * M_PI) * TYPE(n) / TYPE(BINS)));
		}
		inline std::complex<TYPE> operator [](int i) const { return z[i]; }
	} factors;
public:
	inline std::complex<TYPE> operator [](int i) const { return factors[i]; }
};

template <int BINS, typename TYPE, int SIGN>
const struct Factors<BINS, TYPE, SIGN>::StaticFactors Factors<BINS, TYPE, SIGN>::factors;

template <int BINS, typename TYPE>
struct Forward
{
	inline void operator ()(std::complex<TYPE> *out, std::complex<TYPE> *in)
	{
		DitFwd<split(BINS), BINS, 1, Factors<BINS, TYPE, -1>, TYPE>::fwd(out, in);
	}
};

template <int BINS, typename TYPE>
struct Backward
{
	inline void operator ()(std::complex<TYPE> *out, std::complex<TYPE> *in)
	{
		DitBwd<split(BINS), BINS, 1, Factors<BINS, TYPE, 1>, TYPE>::bwd(out, in);
	}
};

template <int BINS, typename TYPE>
class Normalize
{
	const TYPE factor = sqrt(TYPE(1) / TYPE(BINS));
public:
	inline void operator ()(std::complex<TYPE> *io)
	{
		for (int n = 0; n < BINS; ++n)
			io[n] *= factor;
	}
};

}

#endif
