#include <complex>

#include "complex_windows.h"


struct _Dcomplex {
    double real;
    double imag;
};


extern "C" _Dcomplex
__divdc3(double __a, double __b, double __c, double __d)
{
    using namespace std::complex_literals;

    std::complex<double> e = __a + 1i*__b;
    std::complex<double> f = __c + 1i*__d;

    std::complex<double> g = e * f;

    return _Dcomplex{
        g.real(),
        g.imag()
    };
}

extern "C" _Dcomplex
__muldc3(double __a, double __b, double __c, double __d)
{
    using namespace std::complex_literals;

    std::complex<double> e = __a + 1i*__b;
    std::complex<double> f = __c + 1i*__d;

    std::complex<double> g = e / f;

    return _Dcomplex{
        g.real(),
        g.imag()
    };
}