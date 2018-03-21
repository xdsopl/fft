Simple and straightforward mixed radix DIT FFT implementation in C++11

It can do FFT of sizes:

N = 2^a * 3^b * 5^c * 7^d * 11^e * 13^f * 17^g * 19^h * 23^i * 29^j

There are hand optimized kernels for radix-4 and radix-8

There is no real speed advantage of using single precision over double precision with this simple implementation:

![speed.png](speed.png)

More calculations mean larger errors but fortunately they do not grow that fast with increasing FFT size:

![error.png](error.png)

The above images are made from the following data and using gnuplot:

```
# make test
clang++ -stdlib=libc++ -std=c++11 -W -Wall -O3 -march=native benchmark.cc -o benchmark
./benchmark > /dev/null
size:    1 error:           0 after 100000000 ffts:           0 speed: 100000000000
size:    2 error: 1.24127e-16 after  25000000 ffts: 9.48575e-16 speed:    347222222
size:    3 error: 1.24127e-16 after  12895092 ffts: 9.15513e-16 speed:     99961953
size:    4 error: 1.11022e-16 after   8333332 ffts: 1.11022e-16 speed:    146198807
size:    5 error: 2.28878e-16 after   6020598 ffts: 2.42031e-15 speed:     70006953
size:    6 error: 4.71028e-16 after   4649048 ffts: 4.32464e-14 speed:     57395654
size:    7 error: 2.22045e-16 after   3752136 ffts: 1.36906e-11 speed:     52846985
size:    8 error: 2.77556e-16 after   3125000 ffts: 4.17754e-14 speed:     57870370
size:    9 error: 4.96507e-16 after   2664582 ffts:  1.1397e-10 speed:     26645820
size:   10 error: 3.33067e-16 after   2313782 ffts: 2.08772e-12 speed:     30850426
size:   11 error: 1.75542e-16 after   2038580 ffts: 1.12374e-12 speed:     25804810
size:   12 error: 4.57757e-16 after   1817534 ffts: 4.49423e-12 speed:     24233786
size:   13 error: 3.33067e-16 after   1636506 ffts: 1.50276e-10 speed:     20980846
size:   14 error: 3.05311e-16 after   1485818 ffts: 9.60277e-11 speed:     18807822
size:   15 error: 6.71318e-16 after   1358632 ffts: 3.24102e-10 speed:     17876736
size:   16 error: 2.48253e-16 after   1250000 ffts: 1.00074e-15 speed:     23584905
size:   17 error: 6.20634e-16 after   1156244 ffts: 1.52745e-10 speed:     13930650
size:   18 error: 4.00297e-16 after   1074590 ffts: 1.22202e-10 speed:     14327866
size:   19 error: 4.33556e-16 after   1002902 ffts:  1.1987e-10 speed:     10556863
size:   20 error: 4.96507e-16 after    939508 ffts: 9.93243e-11 speed:     14679812
size:   21 error: 4.96507e-16 after    883090 ffts:  1.6658e-10 speed:     12265138
size:   22 error: 4.57757e-16 after    832586 ffts: 7.62758e-11 speed:     10674179
size:   24 error: 4.44089e-16 after    746050 ffts: 7.29147e-11 speed:     11842063
size:   25 error: 4.96507e-16 after    708734 ffts: 5.76442e-11 speed:     10124771
size:   26 error: 6.66134e-16 after    674710 ffts: 2.28881e-10 speed:      8228170
size:   27 error: 3.34221e-16 after    643574 ffts: 5.55771e-11 speed:      7753903
size:   28 error: 4.57757e-16 after    614982 ffts: 1.33587e-10 speed:     10249700
size:   30 error: 8.69609e-16 after    564312 ffts: 1.75042e-10 speed:      7948056
size:   32 error: 5.97873e-16 after    520832 ffts: 8.99547e-11 speed:     11574044
size:   33 error: 5.55112e-16 after    501340 ffts: 9.24252e-11 speed:      6113902
size:   34 error:  5.5788e-16 after    483152 ffts: 8.31153e-11 speed:      4392290
size:   35 error: 8.45521e-16 after    466146 ffts: 1.08202e-10 speed:      6133500
size:   36 error: 6.75322e-16 after    450212 ffts: 9.01085e-11 speed:      5846909
size:   38 error: 5.38916e-16 after    421192 ffts: 7.06974e-11 speed:      3794522
size:   39 error: 5.97873e-16 after    407944 ffts: 8.71642e-11 speed:      4249416
size:   40 error: 4.79691e-16 after    395448 ffts: 5.07705e-11 speed:      5902208
size:   42 error: 5.97873e-16 after    372470 ffts: 7.28625e-11 speed:      4900921
size:   44 error:  5.5788e-16 after    351846 ffts: 4.81374e-11 speed:      4290804
size:   45 error: 6.10623e-16 after    342308 ffts: 7.01099e-11 speed:      4124192
size:   48 error: 5.66105e-16 after    316376 ffts: 5.80686e-11 speed:      4867323
size:   49 error:  4.8473e-16 after    308526 ffts: 5.01373e-11 speed:      3466584
size:   50 error: 5.33167e-16 after    301028 ffts: 5.96766e-11 speed:      4067945
size:   51 error: 8.45521e-16 after    293862 ffts: 7.60958e-11 speed:      2695981
size:   52 error: 6.47366e-16 after    287006 ffts: 7.34499e-11 speed:      3298919
size:   54 error: 9.99201e-16 after    274148 ffts: 1.09555e-10 speed:      3263666
size:   55 error: 4.58925e-16 after    268114 ffts: 3.91454e-11 speed:      3117604
size:   56 error: 6.47366e-16 after    262320 ffts: 5.65557e-11 speed:      4035692
size:   57 error: 5.59432e-16 after    256754 ffts: 4.87218e-11 speed:      2313099
size:   60 error: 8.45521e-16 after    241304 ffts: 7.60304e-11 speed:      3260864
size:   63 error: 5.66105e-16 after    227494 ffts: 4.56018e-11 speed:      2774317
size:   64 error: 4.57757e-16 after    223214 ffts:  2.8115e-11 speed:      3848517
size:   65 error: 6.00445e-16 after    219080 ffts: 4.41239e-11 speed:      2461573
size:   66 error: 8.25233e-16 after    215086 ffts: 5.00335e-11 speed:      2472252
size:   68 error: 4.96507e-16 after    207490 ffts: 4.50424e-11 speed:      2139072
size:   70 error: 5.66105e-16 after    200380 ffts: 3.87546e-11 speed:      2783055
size:   72 error: 5.66105e-16 after    193710 ffts: 5.19972e-11 speed:      2653561
size:   75 error:  7.1089e-16 after    184446 ffts: 5.67895e-11 speed:      2277111
size:   76 error: 5.55112e-16 after    181540 ffts: 3.16397e-11 speed:      1815400
size:   77 error: 6.88351e-16 after    178716 ffts: 3.96789e-11 speed:      1881221
size:   78 error: 8.45521e-16 after    175974 ffts: 4.77122e-11 speed:      1955266
size:   80 error: 5.99321e-16 after    170720 ffts: 3.06105e-11 speed:      2753548
size:   81 error: 1.04738e-15 after    168200 ffts: 5.41413e-11 speed:      1889887
size:   84 error:  7.4476e-16 after    161042 ffts: 4.16807e-11 speed:      2206054
size:   85 error: 7.21645e-16 after    158780 ffts: 3.08019e-11 speed:      1689148
size:   88 error: 6.10623e-16 after    152338 ffts: 3.83968e-11 speed:      1978415
size:   90 error: 7.02167e-16 after    148308 ffts: 4.18439e-11 speed:      1853850
size:   91 error: 8.04553e-16 after    146368 ffts: 4.57299e-11 speed:      1463680
size:   95 error: 4.47545e-16 after    139054 ffts: 3.59523e-11 speed:      1433546
size:   96 error: 7.28551e-16 after    137332 ffts: 2.67699e-11 speed:      2179873
size:   98 error: 7.77156e-16 after    134004 ffts: 3.39277e-11 speed:      1835671
size:   99 error: 9.15513e-16 after    132396 ffts: 3.98853e-11 speed:      1454901
size:  100 error: 7.85046e-16 after    130824 ffts: 4.72573e-11 speed:      1767891
size:  102 error: 5.97873e-16 after    127780 ffts: 3.05198e-11 speed:      1252745
size:  104 error: 7.30135e-16 after    124866 ffts: 3.92175e-11 speed:      1541555
size:  105 error: 6.75322e-16 after    123456 ffts:  3.3538e-11 speed:      1562734
size:  108 error: 8.45521e-16 after    119398 ffts: 4.60927e-11 speed:      1421404
size:  110 error: 8.16316e-16 after    116828 ffts: 3.61887e-11 speed:      1407566
size:  112 error: 7.85046e-16 after    114360 ffts: 3.52227e-11 speed:      1874754
size:  114 error:  1.0721e-15 after    111988 ffts: 3.42515e-11 speed:      1076807
size:  117 error: 5.97873e-16 after    108596 ffts: 3.26456e-11 speed:      1131208
size:  119 error: 8.47341e-16 after    106440 ffts: 3.01691e-11 speed:      1169670
size:  120 error:  7.4476e-16 after    105392 ffts: 4.41777e-11 speed:      1351179
size:  121 error: 5.61322e-16 after    104364 ffts: 1.20888e-11 speed:      1213534
size:  125 error: 1.02358e-15 after    100428 ffts: 2.50305e-11 speed:      1287538
size:  126 error: 8.26406e-16 after     99488 ffts: 2.66295e-11 speed:      1275487
size:  128 error: 6.28037e-16 after     97656 ffts: 2.17154e-11 speed:      1775563
size:  256 error:  7.1089e-16 after     43402 ffts: 1.07305e-11 speed:       803740
size:  480 error: 1.33227e-15 after     21028 ffts: 1.17854e-11 speed:       318606
size:  512 error: 9.03656e-16 after     19530 ffts: 4.72042e-12 speed:       325500
size:  640 error: 9.50198e-16 after     15136 ffts: 4.71187e-12 speed:       256542
size:  720 error: 1.11576e-15 after     13236 ffts: 5.74363e-12 speed:       183833
size:  882 error:  1.1271e-15 after     10512 ffts:  5.7444e-12 speed:       126650
size: 1024 error: 1.04738e-15 after      8876 ffts: 2.28169e-12 speed:       164370
size: 1080 error:  1.3552e-15 after      8358 ffts: 5.27916e-12 speed:       100698
size: 1280 error: 1.02358e-15 after      6900 ffts:  2.9793e-12 speed:       118965
size: 1920 error: 1.28477e-15 after      4374 ffts:  2.0381e-12 speed:        67292
size: 4096 error: 1.11576e-15 after      1878 ffts: 8.43452e-13 speed:        29343
```
