WindowMetrics.c
Copyright Samuel M. Flaxman, 2015

Source code here is a work in progress for a simulation model designed to produce mock datasets that can be used to test methods for making population genetic inference.

Absolutely no warranty of any kind is offered or implied.

You are free to use this source code for any purpose you wish, but please attribute the author if you do.  More importantly, if you do something cool with this, please let me know (sam.flaxman@gmail.com), because I'd just like to know how and what (if anything) other people find it useful for.

To compile, download all the contents of the repo and type "make" on your command line.

Note that the "MT" directory provides the source code for random number generation.  
	+ The attribution for the MT (Mersenne Twister) random number generator is: Saito, M., and Matsumoto, M. (2006). SIMD-oriented Fast Mersenne Twister: a 128-bit pseudorandom number generator. In Monte Carlo and Quasi-Monte Carlo Methods, A. Keller, S. Heinrich, and H. Niederreiter, eds. (Heidelberg, Berlin: Springer-Berlin), pp. 607–622.  
	+ The version used here is dSMFT v2.1.  See: http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/ for versions, downloads, and descriptions of the random number generator code.