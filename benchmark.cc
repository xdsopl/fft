/*
fft - mixed radix fft
Written in 2015 by <Ahmet Inan> <xdsopl@googlemail.com>
To the extent possible under law, the author(s) have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide. This software is distributed without any warranty.
You should have received a copy of the CC0 Public Domain Dedication along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/

#include <iostream>
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

	std::cerr << "FFT size: " << BINS << " max error: " << max_error << std::endl;
}

int main()
{
	typedef double type;
#if 0
	typedef std::complex<type> complex;
#else
	typedef Complex<type> complex;
#endif

	test<1, complex>();
	test<2, complex>();
	test<3, complex>();
	test<4, complex>();
	test<5, complex>();
	test<6, complex>();
	test<8, complex>();
	test<9, complex>();
	test<10, complex>();
	test<12, complex>();
	test<15, complex>();
	test<16, complex>();
	test<18, complex>();
	test<20, complex>();
	test<24, complex>();
	test<25, complex>();
	test<27, complex>();
	test<30, complex>();
	test<32, complex>();
	test<36, complex>();
	test<40, complex>();
	test<45, complex>();
	test<48, complex>();
	test<50, complex>();
	test<54, complex>();
	test<60, complex>();
	test<64, complex>();

	std::random_device rd;
	std::default_random_engine generator(rd());
	std::uniform_real_distribution<type> noise_distribution(-1.0f, 1.0f);
	auto noise = std::bind(noise_distribution, generator);

	const int bins = 2700;
	const int ffts = 10000;

	complex a[bins], b[bins], c[bins];
	for (int i = 0; i < bins; ++i)
		a[i] = complex(noise(), noise());

	FFT::Normalize<bins, complex> norm;
	FFT::Backward<bins, complex> bwd;
	FFT::Forward<bins, complex> fwd;

	auto start = std::chrono::system_clock::now();
	fwd(b, a);
	norm(b);
	for (int i = 2; i < ffts; i += 2) {
		bwd(c, b);
		norm(c);
		fwd(b, c);
		norm(b);
	}
	bwd(c, b);
	norm(c);
	auto end = std::chrono::system_clock::now();
	auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

	type max_error = 0;
	for (int i = 0; i < bins; ++i)
		max_error = std::max(max_error, abs(a[i] - c[i]));

	std::cerr << "FFT size: " << bins << " max error after " << ffts << " ffts: " << max_error;
	std::cerr << " speed: " << (ffts * 1000LL) / std::max(1LL, msec.count()) << " ffts / sec" << std::endl;
	//for (int i = 0; i < bins; ++i)
	//	std::cout << a[i].real() << " " << a[i].imag() << " " << c[i].real() << " " << c[i].imag() << std::endl;
}

