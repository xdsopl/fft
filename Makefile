
CXXFLAGS = -stdlib=libc++ -std=c++11 -W -Wall -O3 -march=native
CXX = clang++

benchmark: benchmark.cc fft.hh complex.hh
	$(CXX) $(CXXFLAGS) $< -o $@

test: benchmark
	./benchmark > /dev/null

.PHONY: clean test

clean:
	rm -f benchmark

