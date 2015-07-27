/*
fft - mixed radix fft
Written in 2015 by <Ahmet Inan> <xdsopl@googlemail.com>
To the extent possible under law, the author(s) have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide. This software is distributed without any warranty.
You should have received a copy of the CC0 Public Domain Dedication along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/

#include <iostream>
#include <random>
#include "fft.hh"

int main()
{
	const int bins = 2700;
	const int loops = 1000;
	const int ffts = 4 * loops;
	typedef float type;

	std::random_device rd;
	std::default_random_engine generator(rd());
	std::uniform_real_distribution<type> noise_distribution(-1.0f, 1.0f);
	auto noise = std::bind(noise_distribution, generator);

	std::complex<type> a[bins], b[bins], c[bins];
	for (int i = 0; i < bins; ++i)
		a[i] = noise();

	FFT::Normalize<bins, type> norm;
	FFT::Backward<bins, type> bwd;
	FFT::Forward<bins, type> fwd;

	auto start = std::chrono::system_clock::now();
	type max_error = 0;
	for (int i = 0; i < loops; ++i) {
		fwd(b, a);
		norm(b);
		bwd(c, b);
		norm(c);
		for (int i = 0; i < bins; ++i)
			max_error = std::max(max_error, std::abs(a[i] - c[i]));
		fwd(b, c);
		norm(b);
		bwd(a, b);
		norm(a);
		for (int i = 0; i < bins; ++i)
			max_error = std::max(max_error, std::abs(a[i] - c[i]));
	}
	auto end = std::chrono::system_clock::now();
	auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

	std::cerr << "FFT size: " << bins << std::endl;
	std::cerr << "max error: " << max_error << std::endl;
	std::cerr << "speed: " << (ffts * 1000) / std::max(1LL, msec.count()) << " ffts / sec" << std::endl;
	//for (int i = 0; i < bins; ++i)
	//	std::cout << a[i].real() << " " << a[i].imag() << " " << c[i].real() << " " << c[i].imag() << std::endl;
}

