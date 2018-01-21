#include <complex.h>

#define FLOAT_COMPLEX_TYPE _Fcomplex
#define FLOAT_COMPLEX_CREATE(real, imag) _FCbuild(real, imag)
#define FLOAT_COMPLEX_MUL_CC(a, b) _FCmulcc(a, b)
#define FLOAT_COMPLEX_ADD_CC(a, b) _FCbuild(crealf(a) + crealf(b), cimagf(a) + cimagf(b))
#define FLOAT_COMPLEX_EQ_CC(a, b) (crealf(a) == crealf(b) && cimagf(a) == cimagf(b))
#define DOUBLE_COMPLEX_TYPE _Dcomplex
#define DOUBLE_COMPLEX_CREATE(real, imag) _Cbuild(real, imag)
#define DOUBLE_COMPLEX_MUL_CC(a, b) _Cmulcc(a, b)
#define DOUBLE_COMPLEX_DIV_CC(a, b) _Cdivcc(a, b)
#define DOUBLE_COMPLEX_ADD_CC(a, b) _Cbuild(creal(a) + creal(b), cimag(a) + cimag(b))
#define DOUBLE_COMPLEX_EQ_CC(a, b) (creal(a) == creal(b) && cimag(a) == cimag(b))

DOUBLE_COMPLEX_TYPE
__divdc3(double __a, double __b, double __c, double __d)
{
    _Dcomplex e = {__a, __b};
    _Dcomplex f = {__c, __d};
    return DOUBLE_COMPLEX_MUL_CC(e, f);
}

DOUBLE_COMPLEX_TYPE
__muldc3(double __a, double __b, double __c, double __d)
{
    _Dcomplex e = {__a, __b};

    return e;
}