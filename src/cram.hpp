#ifdef PYNE_CRAM_SOLVE_C
#error "Both cram.h (the C header for CRAM) and cram.hpp (C++) have been included, only one is allowed!"
#endif

#ifndef PYNE_CRAM_SOLVE_CPP
#define PYNE_CRAM_SOLVE_CPP


typedef struct pyne_cram_transmute_info_tag {
  int n;
  int nnz;
  int* i;
  int* j;
  char** nucs;
  int* nucids;
  double* decay_matrix;
} pyne_cram_transmute_info_t;

extern pyne_cram_transmute_info_t pyne_cram_transmute_info;

int pyne_cram_transmute_ij(int i, int j);

int pyne_cram_transmute_nucid_to_i(int nucid);


void pyne_cram_solve_double(double* A, double* b, double* x);
void pyne_cram_diag_add_double(double* A, double alpha);
void pyne_cram_dot_double(double* A, double* x, double* y);
void pyne_cram_scalar_times_vector_double(double, double*);

// double complex types are not allowed in the C++ type convention,
// which extern "C" forced us into.

void pyne_cram_expm_multiply6(double* A, double* b, double* x);
void pyne_cram_expm_multiply8(double* A, double* b, double* x);
void pyne_cram_expm_multiply10(double* A, double* b, double* x);
void pyne_cram_expm_multiply12(double* A, double* b, double* x);
void pyne_cram_expm_multiply14(double* A, double* b, double* x);
void pyne_cram_expm_multiply16(double* A, double* b, double* x);
void pyne_cram_expm_multiply18(double* A, double* b, double* x);
#endif // PYNE_CRAM_SOLVE_CPP
