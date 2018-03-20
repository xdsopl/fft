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

	std::cerr << "FFT size: " << std::setw(4) << BINS;
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

	test<1, complex_type>();
	test<2, complex_type>();
	test<3, complex_type>();
	test<4, complex_type>();
	test<5, complex_type>();
	test<6, complex_type>();
	test<7, complex_type>();
	test<8, complex_type>();
	test<9, complex_type>();
	test<10, complex_type>();
	test<11, complex_type>();
	test<12, complex_type>();
	test<13, complex_type>();
	test<14, complex_type>();
	test<15, complex_type>();
	test<16, complex_type>();
	test<18, complex_type>();
	test<20, complex_type>();
	test<21, complex_type>();
	test<22, complex_type>();
	test<24, complex_type>();
	test<25, complex_type>();
	test<26, complex_type>();
	test<27, complex_type>();
	test<28, complex_type>();
	test<30, complex_type>();
	test<32, complex_type>();
	test<33, complex_type>();
	test<35, complex_type>();
	test<36, complex_type>();
	test<39, complex_type>();
	test<40, complex_type>();
	test<42, complex_type>();
	test<44, complex_type>();
	test<45, complex_type>();
	test<48, complex_type>();
	test<49, complex_type>();
	test<50, complex_type>();
	test<52, complex_type>();
	test<54, complex_type>();
	test<55, complex_type>();
	test<56, complex_type>();
	test<60, complex_type>();
	test<63, complex_type>();
	test<64, complex_type>();
	test<65, complex_type>();
	test<66, complex_type>();
	test<70, complex_type>();
	test<72, complex_type>();
	test<75, complex_type>();
	test<77, complex_type>();
	test<78, complex_type>();
	test<80, complex_type>();
	test<81, complex_type>();
	test<84, complex_type>();
	test<88, complex_type>();
	test<90, complex_type>();
	test<91, complex_type>();
	test<96, complex_type>();
	test<98, complex_type>();
	test<99, complex_type>();
	test<100, complex_type>();
	test<104, complex_type>();
	test<105, complex_type>();
	test<108, complex_type>();
	test<110, complex_type>();
	test<112, complex_type>();
	test<117, complex_type>();
	test<120, complex_type>();
	test<121, complex_type>();
	test<125, complex_type>();
	test<126, complex_type>();
	test<128, complex_type>();
	test<256, complex_type>();
	test<480, complex_type>();
	test<512, complex_type>();
	test<640, complex_type>();
	test<720, complex_type>();
	test<882, complex_type>();
	test<1024, complex_type>();
	test<1080, complex_type>();
	test<1280, complex_type>();
	test<1920, complex_type>();
	test<4096, complex_type>();
}

