#! /bin/bash

# Code generator for prime number FFT
# Copyright 2018 Ahmet Inan <xdsopl@gmail.com>

RADIX=$1

struct_header() {
	cat << EOF
template <int STRIDE, typename TYPE>
struct Dit<$RADIX, $RADIX, STRIDE, TYPE, $1>
{
EOF
}

struct_footer() {
	echo "	static inline void dit(TYPE *out, const TYPE *in, const TYPE *)"
	echo "	{"
	echo -n "		dft(out,"
	for ((x = 1; x < RADIX; x++)) ; do
		echo -n " out + $x,"
	done
	echo
	echo -n "			in[0], in[STRIDE]"
	for ((x = 2; x < RADIX; x++)) ; do
		echo -n ", in[$x * STRIDE]"
	done
	echo ");"

	echo "	}"
	echo "};"
}

dft_header() {
	echo -n "	static inline void dft(TYPE *out0,"
	for ((x = 1; x < RADIX; x++)) ; do
		echo -n " TYPE *out$x,"
	done
	echo
	echo -n "		"
	for ((x = 0; x < RADIX-1; x++)) ; do
		echo -n "TYPE in$x, "
	done
	echo "TYPE in$((RADIX-1)))"
	echo "	{"
}

echo
struct_header -1
dft_header
echo -n "		Dit<$RADIX, $RADIX, STRIDE, TYPE, 1>::dft(out0,"
for ((x = RADIX-1; x > 0; x--)) ; do
	echo -n " out$x,"
done
echo
echo -n "			"
for ((x = 0; x < RADIX-1; x++)) ; do
	echo -n "in$x, "
done
echo "in$((RADIX-1)));"
echo "	}"
struct_footer

echo

struct_header 1
dft_header

for ((y = 1; y <= RADIX/2; y++)) ; do
	echo "		TYPE a$y(in$y + in$((RADIX-y))), t$y(twiddle(in$y, in$((RADIX - y))));"
done

for ((x = 1; x <= RADIX/2; x++)) ; do
	echo -n "		TYPE c$x(cx<$x,$RADIX>(a1)"
	for ((y = 2; y <= RADIX/2; y++)) ; do
		xy=$(((x*y)%RADIX))
		NUM=$((xy > RADIX/2 ? RADIX - xy : xy))
		echo -n " + cx<$NUM,$RADIX>(a$y)"
	done
	echo ");"
	echo -n "		TYPE s$x(sx<$x,$RADIX>(t1)"
	for ((y = 2; y <= RADIX/2; y++)) ; do
		xy=$(((x*y)%RADIX))
		NUM=$((xy > RADIX/2 ? RADIX - xy : xy))
		SIGN=+
		((xy > RADIX/2)) && SIGN=-
		echo -n " $SIGN sx<$NUM,$RADIX>(t$y)"
	done
	echo ");"
done

echo -n "		*out0 = in0"
for ((y = 1; y <= RADIX/2; y++)) ; do
	echo -n " + a$y"
done
echo ";"

for ((x = 1; x <= RADIX/2; x++)) ; do
	echo "		*out$x = in0 + c$x - s$x;"
done

for ((x = RADIX/2+1; x < RADIX; x++)) ; do
	echo "		*out$x = in0 + c$((RADIX-x)) + s$((RADIX-x));"
done

echo "	}"

struct_footer

echo

cat << EOF
template <int BINS, int STRIDE, typename TYPE, int SIGN>
struct Dit<$RADIX, BINS, STRIDE, TYPE, SIGN>
{
	static const int RADIX = $RADIX;
	static const int QUOTIENT = BINS / RADIX;
	static void dit(TYPE *out, const TYPE *in, const TYPE *z)
	{
		for (int o = 0, i = 0; o < BINS; o += QUOTIENT, i += STRIDE)
			Dit<split(QUOTIENT), QUOTIENT, RADIX * STRIDE, TYPE, SIGN>::dit(out + o, in + i, z);
EOF

echo -n "		for (int k0 = 0, k1 = QUOTIENT,"
for ((x = 2; x < RADIX; x++)) ; do
	echo -n " k$x = $x * QUOTIENT,"
done
echo
echo -n "				l1 = 0"
for ((x = 2; x < RADIX; x++)) ; do
	echo -n ", l$x = 0"
done
echo ";"
echo "				k0 < QUOTIENT;"
echo -n "				++k0,"
for ((x = 1; x < RADIX; x++)) ; do
	echo -n " ++k$x,"
done
echo
echo -n "				l1 += STRIDE"
for ((x = 2; x < RADIX; x++)) ; do
	echo -n ", l$x += $x * STRIDE"
done
echo ")"
echo -n "			Dit<RADIX, RADIX, STRIDE, TYPE, SIGN>::dft(out + k0,"
for ((x = 1; x < RADIX; x++)) ; do
	echo -n " out + k$x,"
done
echo
echo -n "				out[k0]"
for ((x = 1; x < RADIX; x++)) ; do
	echo -n ", z[l$x] * out[k$x]"
done
echo ");"
echo "	}"
echo "};"

