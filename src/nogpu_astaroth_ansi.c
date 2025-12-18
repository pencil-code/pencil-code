/*                            nogpu_astaroth_ansi.c
                              ---------------------
  Dummy version of gpu_astaroth_ansi.c for use with start

*/
#include "headers_c.h"

/* ------------------------------------------------------------------- */
void FTNIZE(initialize_gpu_c)(REAL* f, FINT* comm_fint, double* t, int* nt, 
				FINT* lreads_all_vars_from_device_,
				FINT* lcpu_timestep_on_gpu_
				)
{
}
/* ------------------------------------------------------------------- */
void FTNIZE(register_gpu_c)()
{
}
/* ------------------------------------------------------------------- */
void FTNIZE(finalize_gpu_c)()
{
}
/* ------------------------------------------------------------------- */
void FTNIZE(get_farray_ptr_gpu_c)(REAL** p_f_in)
{
}
/* ------------------------------------------------------------------- */
void FTNIZE(rhs_gpu_c)(FINT *isubstep, double* t)
{
}
/* ------------------------------------------------------------------- */
void FTNIZE(before_boundary_gpu_c)(FINT* lrmv, FINT *isubstep, double *t)
{
}
/* ------------------------------------------------------------------- */
void FTNIZE(after_timestep_gpu_c)()
{
}
/* ---------------------------------------------------------------------- */
void FTNIZE(copy_farray_c)(REAL* f)
{
}
/* ------------------------------------------------------------------- */
void FTNIZE(load_farray_c)()
{
}
/* ------------------------------------------------------------------- */
void FTNIZE(reload_gpu_config_c)()
{
}
/* ------------------------------------------------------------------- */
void FTNIZE(update_on_gpu_scal_by_ind_c)(int *index, REAL* value)
{
}
/* ------------------------------------------------------------------- */
void FTNIZE(update_on_gpu_arr_by_ind_c)(int *index)
{
}
/* ------------------------------------------------------------------- */
int FTNIZE(update_on_gpu_scal_by_name_c)(char *varname, REAL* value)
{
  return 0;
}
/* ------------------------------------------------------------------- */
int FTNIZE(update_on_gpu_arr_by_name_c)(char *varname)
{
  return 0;
}
/* ------------------------------------------------------------------- */
void FTNIZE(test_rhs_c)(REAL* f_in, REAL* df_truth)
{
}
/* ------------------------------------------------------------------- */
void FTNIZE(gpu_set_dt_c)(double* t)
{
}
/* ------------------------------------------------------------------- */
void FTNIZE(torchinfer_c)(FINT* itsub)
{
}
/* ------------------------------------------------------------------- */
void FTNIZE(torchtrain_c)(REAL *loss, FINT* itsub)
{
}
/* ------------------------------------------------------------------- */
void FTNIZE(radtransfer_gpu_c)()
{
}
/* ------------------------------------------------------------------- */
void FTNIZE(get_gpu_reduced_vars_c)(REAL* dst)
{
}
/* ------------------------------------------------------------------- */
void FTNIZE(test_bcs_c)(void)
{
}
/* ------------------------------------------------------------------- */
void FTNIZE(update_after_substep_gpu_c)(void)
{
}
/* ------------------------------------------------------------------- */
void FTNIZE(split_update_gpu_c)(void)
{
}
/* ------------------------------------------------------------------- */

