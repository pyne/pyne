// MCNP5/dagmc/meshtal_funcs.h

#ifndef DAGMC_MESHTAL_IFACE_H
#define DAGMC_MESHTAL_IFACE_H

#ifdef __cplusplus
extern "C" {
#endif

/*************************************************************************
 * This file declares the functions forming the C/Fortran bridge for the
 * advanced DAGMC Tally features.  Current MOAB-based mesh tallies that
 * are available include unstructured track length and KDE options.
 *************************************************************************/

/**
 * Functions from fmesh_mod are implemented in src/fmesh_mod.F90 
 * and should only be called from within meshtal_funcs.cpp
 */

/** FORT_FUNC:
 * Macro to access symbol of fortran function 'func' in module 'mod' 
 */
#ifndef FORT_FUNC

/* gcc/gfortran 4.3 and above: name mangling is '__module_MOD_function' */
#if __GNUC__ > 4 || ( __GNUC__ == 4  && __GNUC_MINOR__ >= 3 )
#define FORT_FUNC( mod, func ) __##mod##_MOD_##func

/* gcc/gfortran < 4.3: name mangling is '__module__function' */
#elif __GNUC__ == 4 
#define FORT_FUNC( mod, func ) __##mod##__##func

/* Something we haven't encountered yet */
#else
/* Comment out this error to force compile to proceed; it may or may not work */
#error "DagMC: unknown compiler with unknown fortran name mangling scheme."
#define FORT_FUNC( mod, func ) __##mod##__##func

#endif
#endif /* FORT_FUNC */

#define FMESH_FUNC( func ) FORT_FUNC( fmesh_mod, func )

/* Make a valid Fortran pointer to a C array */
extern void FMESH_FUNC(dagmc_make_fortran_pointer)(void* fort_ref, double* array, int* size);

/**
 * The dagmc_fmesh_*_ functions are called from fortran to drive our advanced mesh tallies,
 * mostly from fmesh_mod.F90.  They should probably not be called from C++ code.
 * Per-function documentation is found in meshtal_funcs.cpp
 */
void dagmc_fmesh_setup_mesh_(int* fm_ipt, int* id, int* fmesh_idx,
                             double* energy_mesh, int* n_energy_mesh, int* tot_energy_bin, 
                             char* comment, int* n_comment_lines, int* is_collision_tally);

void dagmc_fmesh_end_history_();

void dagmc_fmesh_score_(int* ipt,
                        double *x, double *y, double *z,
                        double *u, double *v, double *w, 
                        double *erg,double *wgt,
                        double *d, int* icl);

void dagmc_fmesh_print_(double* sp_norm);

void dagmc_collision_score_(int* ipt,
                            double* x,   double* y, double* z, 
                            double* erg, double* wgt,
                            double* ple, int* icl);

void dagmc_update_multiplier_(int* fmesh_idx, double* value);

void dagmc_fmesh_get_tally_data_(int* tally_id, void* fortran_data_pointer);
void dagmc_fmesh_get_error_data_(int* tally_id, void* fortran_data_pointer);
void dagmc_fmesh_get_scratch_data_(int* tally_id, void* fortran_data_pointer);
void dagmc_fmesh_clear_data_();
void dagmc_fmesh_add_scratch_to_tally_(int* tally_id);
void dagmc_fmesh_add_scratch_to_error_(int* tally_id);

#ifdef __cplusplus
} /* extern "C" */
#endif 

#endif /* DAGMC_MESHTAL_IFACE_H */

// end of MCNP5/dagmc/meshtal_funcs.h
