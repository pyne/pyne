/* This file was generated automatically with transmutagen version 1.0.1. */
/* The command used to generate this file was: python -m transmutagen.gensolve --py-solve --namespace=pyne_cram --outfile=cram.c*/
#ifndef PYNE_CRAM_SOLVE_C
#define PYNE_CRAM_SOLVE_C


#include <complex.h>

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

void pyne_cram_solve_complex(double complex* A, double complex* b, double complex* x);
void pyne_cram_diag_add_complex(double complex* A, double complex alpha);
void pyne_cram_dot_complex(double complex* A, double complex* x, double complex* y);
void pyne_cram_scalar_times_vector_complex(double complex, double complex*);

void pyne_cram_expm_multiply6(double* A, double* b, double* x);
void pyne_cram_expm_multiply8(double* A, double* b, double* x);
void pyne_cram_expm_multiply10(double* A, double* b, double* x);
void pyne_cram_expm_multiply12(double* A, double* b, double* x);
void pyne_cram_expm_multiply14(double* A, double* b, double* x);
void pyne_cram_expm_multiply16(double* A, double* b, double* x);
void pyne_cram_expm_multiply18(double* A, double* b, double* x);
#endif
