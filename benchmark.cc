/*
fft - mixed radix fft
Written in 2015 by <Ahmet Inan> <xdsopl@googlemail.com>
To the extent possible under law, the author(s) have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide. This software is distributed without any warranty.
You should have received a copy of the CC0 Public Domain Dedication along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/

#include <iostream>
#include <iomanip>
#include <random>
#include <complex>
#include "complex.hh"
#include "fft.hh"

template <int BINS, typename TYPE>
static void test()
{
	typedef typename TYPE::value_type value_type;
	std::random_device rd;
	std::default_random_engine generator(rd());
	std::uniform_real_distribution<value_type> noise_distribution(-1.0f, 1.0f);
	auto noise = std::bind(noise_distribution, generator);

	TYPE a[BINS], b[BINS], c[BINS];
	for (int i = 0; i < BINS; ++i)
		a[i] = TYPE(noise(), noise());

	FFT::Normalize<BINS, TYPE> norm;
	FFT::Backward<BINS, TYPE> bwd;
	FFT::Forward<BINS, TYPE> fwd;

	fwd(b, a);
	norm(b);
	bwd(c, b);
	norm(c);

	value_type max_error = 0;
	for (int i = 0; i < BINS; ++i)
		max_error = std::max(max_error, abs(a[i] - c[i]));

	int ffts = ~1 & (int)(100000000 / BINS / (log2(BINS) + 1));
	auto start = std::chrono::system_clock::now();
	for (int i = 0; i < ffts; i += 2) {
		fwd(b, c);
		norm(b);
		bwd(c, b);
		norm(c);
	}
	auto end = std::chrono::system_clock::now();
	auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

	value_type max_error_growth = 0;
	for (int i = 0; i < BINS; ++i)
		max_error_growth = std::max(max_error_growth, abs(a[i] - c[i]));

	std::cerr << "FFT size: " << std::setw(5) << BINS;
	std::cerr << " max error: " << std::setw(11) << max_error;
	std::cerr << " max error growth after " << std::setw(9) << ffts << " ffts: " << std::setw(11) << max_error_growth;
	std::cerr << " speed: " << std::setw(12) << (ffts * 1000LL) / std::max(1LL, msec.count()) << " ffts / sec" << std::endl;
}

int main()
{
	typedef double value_type;
#if 0
	typedef std::complex<value_type> complex_type;
#else
	typedef Complex<value_type> complex_type;
#endif
	test<16, complex_type>();
	test<64, complex_type>();
	test<256, complex_type>();
	test<1024, complex_type>();
	test<4096, complex_type>();
	test<16384, complex_type>();
}

