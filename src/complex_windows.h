#if defined(__clang__) && defined(_MSC_VER)
/* clang-cl */
#define complex _Complex
#define _Complex_I  (0.0F +  1.0iF)
#define I _Complex_I
#endif