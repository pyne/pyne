cdef extern from "cram.h":

    cdef struct pyne_cram_transmute_info_tag:
        int n
        int nnz
        int* i
        int* j
        char** nucs
        int* nucids
        double* decay_matrix

    ctypedef pyne_cram_transmute_info_tag pyne_cram_transmute_info_t
    cdef pyne_cram_transmute_info_t pyne_cram_transmute_info

    void pyne_cram_solve_double(double*, double*, double*)
    void pyne_cram_diag_add_double(double*, double)
    void pyne_cram_dot_double(double*, double*, double*)
    void pyne_cram_scalar_times_vector_double(double, double*)

    void pyne_cram_solve_complex(double complex*, double complex*, double complex*)
    void pyne_cram_diag_add_complex(double complex*, double complex)
    void pyne_cram_dot_complex(double complex*, double complex*, double complex*)
    void pyne_cram_scalar_times_vector_complex(double complex, double complex*)

    void pyne_cram_expm_multiply6(double*, double*, double*)
    void pyne_cram_expm_multiply8(double*, double*, double*)
    void pyne_cram_expm_multiply10(double*, double*, double*)
    void pyne_cram_expm_multiply12(double*, double*, double*)
    void pyne_cram_expm_multiply14(double*, double*, double*)
    void pyne_cram_expm_multiply16(double*, double*, double*)
    void pyne_cram_expm_multiply18(double*, double*, double*)
