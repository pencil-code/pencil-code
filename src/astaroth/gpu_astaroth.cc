/*
   Date:   8-Feb-2017
   Author: M. Rheinhardt & j. Pekkilae
   Description:
           ANSI C and standard library callable function wrappers for ASTAROTH-nucleus functions to be called from Fortran.
  Comments: 
*/
// General headers.
#include <sstream>
#include <math.h>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>
#include <sys/resource.h>
#include <fstream>
#include <dlfcn.h>
#include <thread>
#include <stack>

static std::thread GW_thread{};
static int device_id = 0;

//TP: defined here since mpi.h can have its own definition of DOUBLE_PRECISION
//    and we don't want to conflict it with it. This is at least true on my laptop
#if AC_DOUBLE_PRECISION
#define DOUBLE_PRECISION 1
#endif

#define CUDA_ERRCHK(X)
int nt = 0;
int counter = 0;
int train_counter = 0;
int done = 0;
int my_rank;
bool calculated_coeff_scales = false;

// Astaroth headers
#include "astaroth.h"
#include "fft.h"

AcTaskGraph* GW_timestep_graph  =  NULL;
AcTaskGraph* boundary_z_halo_exchange_graph =  NULL;
//TP: logs performance metrics of Astaroth
const bool performance_logs = false;

#define real AcReal
#include "math_utils.h"
#define EXTERN
#define FINT int

#if DOUBLE_PRECISION
  #define REAL_MAX DBL_MAX
#else
  #define REAL_MAX FLT_MAX
#endif
#include "../cparam_c.h"

//TP: these are ugly but for the moment we live with these
#if TRANSPILATION
  #define thread_block_loop_factors thread_block_loop_factors__mod__gpu
  #define luses_aa_pot2_top luses_aa_pot2_top__mod__cdata 
  #define luses_aa_pwd_top luses_aa_pwd_top__mod__cdata 
  #define lsplit_gw_rhs_from_rest_on_gpu lsplit_gw_rhs_from_rest_on_gpu__mod__gravitational_waves_htxk
  #define limplicit_diffusion_with_fft limplicit_diffusion_with_fft__mod__implicitdiffusion
  #define limplicit_diffusion_with_cg  limplicit_diffusion_with_cg__mod__implicitdiffusion
  #define limplicit_resistivity limplicit_resistivity__mod__magnetic
  #define limplicit_viscosity   limplicit_viscosity__mod__viscosity
  #define lcumulative_df_on_gpu lcumulative_df_on_gpu__mod__gpu
  #define lbidiagonal_derij lbidiagonal_derij__mod__cdata
  #define nu nu__mod__viscosity
  #define nu_hyper2 nu_hyper2__mod__viscosity
  #define nu_hyper3 nu_hyper3__mod__viscosity
  #define n_odevars n_odevars__mod__cdata 
  #define AC_f_ode AC_f_ode__mod__cdata 
  #define gamma gamma__mod__equationofstate

  #define eta eta__mod__magnetic
  #define eta_hyper2 eta_hyper2__mod__magnetic
  #define eta_hyper3 eta_hyper3__mod__magnetic

  #define chi chi__mod__energy
  #define chi_hyper2 chi_hyper2__mod__energy
  #define chi_hyper3 chi_hyper3__mod__energy

  #define lmorton_curve lmorton_curve__mod__cdata
  #define ltest_bcs     ltest_bcs__mod__gpu
  #define num_substeps  num_substeps__mod__cdata
  #define maux_vtxbuf_index maux_vtxbuf_index__mod__cdata
  #define read_vtxbuf_from_gpu read_vtxbuf_from_gpu__mod__cdata 
  #define ldt ldt__mod__cdata
  #define dt dt__mod__cdata
  #define it it__mod__cdata
  #define dx dx__mod__cdata
  #define dy dy__mod__cdata
  #define dz dz__mod__cdata
  #define mu0 mu0__mod__cdata
  #define ldensity_nolog ldensity_nolog__mod__cdata 
  #define cdt            cdt__mod__cdata 
  #define itorder        itorder__mod__cdata
  #define dtinc          dtinc__mod__cdata
  #define dtdec          dtdec__mod__cdata

  #define AC_ldensity_nolog AC_ldensity_nolog__mod__cdata
  #define AC_mu0  AC_mu0__mod__cdata

  #define cdtv  cdtv__mod__cdata
  #define cdtv2 cdtv2__mod__cdata
  #define cdtv3 cdtv3__mod__cdata

  #define lcylindrical_coords lcylindrical_coords__mod__cdata
  #define lspherical_coords   lspherical_coords__mod__cdata
  #define lcartesian_coords   lcartesian_coords__mod__cdata

  #define lequidist lequidist__mod__cdata

  #define dx_1 dx_1__mod__cdata
  #define dy_1 dy_1__mod__cdata
  #define dz_1 dz_1__mod__cdata

  #define dx_tilde dx_tilde__mod__cdata
  #define dy_tilde dy_tilde__mod__cdata
  #define dz_tilde dz_tilde__mod__cdata
  
  #define rcyl_mn1 rcyl_mn1__mod__cdata
  #define r1_mn    r1_mn__mod__cdata
  #define sin1th   sin1th__mod__cdata
  #define cotth    cotth__mod__cdata

  #define ltrained ltrained__mod__training
  #define lcpu_timestep_on_gpu   lcpu_timestep_on_gpu__mod__gpu
  #define lsingle_precision_timestep lsingle_precision_timestep__mod__gpu
  #define lperi                  lperi__mod__cdata
  #define lxyz                   lxyz__mod__cdata
  #define ldebug ldebug__mod__cdata
  
  #define deltay  deltay__mod__cdata
  #define eps_rkf eps_rkf__mod__cdata
  #define itauxx itauxx__mod__training
  #define itauxy itauxy__mod__training
  #define itauxz itauxz__mod__training
  #define itauyy itauyy__mod__training
  #define itauyz itauyz__mod__training
  #define itauzz itauzz__mod__training
  
  #define input_channels  input_channels__mod__training
  #define output_channels output_channels__mod__training

  #define lsecond_force lsecond_force__mod__forcing
  #define lforce_helical lforce_helical__mod__forcing

  #define lrmv lrmv__mod__cdata
  #define lconserve_total_mass lconserve_total_mass__mod__density
  #define tstart_selfgrav tstart_selfgrav__mod__selfgravity
#endif

static bool lread_all_vars_from_device = false;
static bool lcpu_timestep_on_gpu       = false;
AcReal dt1_interface;
static int rank;
static MPI_Comm comm_pencil = MPI_COMM_NULL;

// Astaroth objects instantiation.
static AcMesh mesh = acInitMesh();
//static AcMesh test_mesh;

bool torch_train_CAPI(int sub_dims[3], AcReal* input, AcReal* label, AcReal* loss_val,
		     const int input_fields, const int output_fields, const char* model_name);
bool torch_infer_CAPI(int sub_dims[3], AcReal* input, AcReal* label, 
		     const int input_fields, const int output_fields, const char* model_name, bool subsample);
bool torch_create_model_CAPI(const char* name, const char* config_fname, int device);
bool torch_create_distributed_model_CAPI(const char* name, const char* config_fname, MPI_Comm mpi_comm, int device);
bool torch_load_CAPI(const char* name, const char* fname);
bool torch_load_checkpoint_CAPI(const char* name, const char* checkpoint_dir, int64_t* step_train, int64_t* step_inference);
bool torch_save_model_CAPI(const char* name, const char* fname);
bool torch_save_checkpoint_CAPI(const char* name, const char* checkpoint_dir);

void scaling();
extern "C" void print_debug();
float MSE();

extern "C" void copyFarray(AcReal* f);    // forward declaration
extern "C" void loadFarray(); // forward declaration
void denormalize(std::string filename, AcRealSymmetricTensor &tau_means, AcRealSymmetricTensor &tau_stds);

/***********************************************************************************************/
AcReal cpu_pow(AcReal const val, AcReal exponent)
{
// Masks hip GPU power function.
    return std::pow(val, exponent);
}
/***********************************************************************************************/
AcReal sign(const AcReal a, const AcReal b)
{
	return (b < (AcReal)0.0) ? -abs(a) : abs(a);
}
/***********************************************************************************************/
bool has_nans(AcMesh mesh_in)
{
  bool res = false;
  AcMeshDims dims = acGetMeshDims(acGridGetLocalMeshInfo());
  const auto DEVICE_VTXBUF_IDX = [&](const int x, const int y, const int z)
  				{
					return acVertexBufferIdx(x,y,z,mesh.info);
  				};
  for (size_t i = dims.n0.x; i < dims.n1.x; i++)
  {
    for (size_t j = dims.n0.y; j < dims.n1.y; j++)
    {
      for (size_t k = dims.n0.z; k < dims.n1.z; k++)
      {
        for (size_t ivar = 0; ivar < acGetNumFields(); ivar++)
        {
          if (mesh_in.vertex_buffer[ivar] == NULL) continue;
          if (isnan(mesh_in.vertex_buffer[ivar][DEVICE_VTXBUF_IDX(i, j, k)]))
          {
            res = true;
            acLogFromRootProc(rank,"nan at %d,%d,%d\n", i, j, k);
            acLogFromRootProc(rank,"field = %d", ivar);
          }
        }
      }
    }
  }
  return res;
}
/***********************************************************************************************/
AcReal to_real(void* param, const char* name)
{
	if (param == NULL)
	{
		fprintf(stderr,"Passed NULL to pushparsc, expected real: %s!!\n",name);
		abort();
	}
	return *((AcReal*)param);
}
/***********************************************************************************************/
int to_int(void* param, const char* name)
{
	if (param == NULL)
	{
		fprintf(stderr,"Passed NULL to pushparsc, expected int: %s!!\n",name);
		abort();
	}
	return *((int*)param);
}
/***********************************************************************************************/
bool to_bool(void* param, const char* name)
{
	if (param == NULL)
	{
		fprintf(stderr,"Passed NULL to pushparsc, expected bool: %s!!\n",name);
		abort();
	}
	return *((bool*)param);
}
/***********************************************************************************************/
int3 to_int3(void* param, const char* name)
{
	if (param == NULL)
	{
		fprintf(stderr,"Passed NULL to pushparsc, expected int3: %s!!\n",name);
		abort();
	}
        int* arr = (int*)param;
        return (int3){arr[0],arr[1],arr[2]};
}
/***********************************************************************************************/
AcReal3 to_real3(void* param, const char* name)
{
	if (param == NULL)
	{
		fprintf(stderr,"Passed NULL to pushparsc, expected real3: %s!!\n",name);
		abort();
	}
        AcReal* arr = (AcReal*)param;
        return (AcReal3){arr[0],arr[1],arr[2]};
}
/***********************************************************************************************/
AcBool3 to_bool3(void* param, const char* name)
{
	if (param == NULL)
	{
		fprintf(stderr,"Passed NULL to pushparsc, expected bool3: %s!!\n",name);
		abort();
	}
        bool* arr = (bool*)param;
        return (AcBool3){arr[0],arr[1],arr[2]};
}
/***********************************************************************************************/
// PC interface headers.
#include "PC_moduleflags.h"

#if TRANSPILATION
#if LFORCING
torus_rect to_torus_rect(void* param, const char* name)
{
       if (param == NULL)
       {
               fprintf(stderr,"Passed NULL to pushparsc: %s!!\n",name);
               abort();
       }
       //TP: placeholder for now before testing is torus_react being POD sufficient for save access here
       return (torus_rect){};
}
#endif
#endif
/***********************************************************************************************/
#include "../sub_c.h"           // provides set_dt
#include "../boundcond_c.h"     // provides boundconds[xyz] etc.
#include "PC_module_parfuncs.h" // provides stuff from physics modules

#if LFORCING
  #include "../forcing_c.h"     // provides forcing_pars_hel
  #include "forcing.h"
#endif

/***********************************************************************************************/
template <typename T>
void
pushpars_special(const char* func_name, T p_pars)
{
 	void* handle = dlopen("src/special.so",RTLD_NOW | RTLD_LOCAL);
	auto sub = reinterpret_cast<void (*)(T)>(dlsym(handle,func_name));
	sub(p_pars);
	dlclose(handle);
}
/***********************************************************************************************/
const char*
get_cwd()
{
  static char cwd[9000];
  cwd[0] = '\0';
  const char* err = getcwd(cwd, sizeof(cwd));
  if (err == NULL) 
  {
	  fprintf(stderr,"Was not able to get cwd!\n");
	  exit(EXIT_FAILURE);
  }
  return cwd;
}
/***********************************************************************************************/
#define PCLoad acPushToConfig
/***********************************************************************************************/
int same_path(const char *p1, const char *p2) {
    char r1[PATH_MAX];
    char r2[PATH_MAX];

    if (!realpath(p1, r1) || !realpath(p2, r2)) {
        return 0; // failed to resolve
    }

    return strcmp(r1, r2) == 0;
}
static int lac_sparse_autotuning = 0;
/***********************************************************************************************/
void setupConfig(AcMeshInfo& config)
{
  config = acInitInfo();
 
  #include "PC_modulepars.h"

  //CUDA_AWARE_MPI is true by default, to turn it off set CUDA_AWARE_MPI=OFF in your config
  PCLoad(config, AC_use_cuda_aware_mpi,CUDA_AWARE_MPI);
  PCLoad(config, AC_bidiagonal_derij,lbidiagonal_derij);

  PCLoad(config, AC_x,x__mod__cdata);
  PCLoad(config, AC_y,y__mod__cdata);
  PCLoad(config, AC_z,z__mod__cdata);

  PCLoad(config, AC_thread_block_loop_factors,thread_block_loop_factors);

  PCLoad(config, AC_inv_cyl_r,rcyl_mn1);
  PCLoad(config, AC_inv_r,r1_mn);
  PCLoad(config, AC_inv_sin_theta,sin1th);
  PCLoad(config, AC_cot_theta,cotth);

//  loads for non-equidistant grids
  PCLoad(config,AC_nonequidistant_grid, (AcBool3){!lequidist.x,!lequidist.y,!lequidist.z});
  PCLoad(config,AC_inv_mapping_func_derivative_x,dx_1);
  PCLoad(config,AC_inv_mapping_func_derivative_y,dy_1);
  PCLoad(config,AC_inv_mapping_func_derivative_z,dz_1);

  PCLoad(config,AC_mapping_func_tilde_x,dx_tilde);
  PCLoad(config,AC_mapping_func_tilde_y,dy_tilde);
  PCLoad(config,AC_mapping_func_tilde_z,dz_tilde);

//  loads for slope-limited-diffusion related arrays
  PCLoad(config, AC_x12,x12__mod__cdata);
  PCLoad(config, AC_y12,y12__mod__cdata);
  PCLoad(config, AC_sinth12,sinth12__mod__cdata);
  PCLoad(config, AC_z12,z12__mod__cdata);

  PCLoad(config, AC_rk_order, itorder);
  PCLoad(config, AC_shear,lshear);
  PCLoad(config, AC_rk_cumulative_df,lcumulative_df_on_gpu);

  if (lcartesian_coords)
  {
          PCLoad(config, AC_coordinate_system, AC_CARTESIAN_COORDINATES);
  }
  if (lspherical_coords)
  {
          PCLoad(config, AC_coordinate_system, AC_SPHERICAL_COORDINATES);
  }
  if (lcylindrical_coords)
  {
          PCLoad(config, AC_coordinate_system, AC_CYLINDRICAL_COORDINATES);
  }
  //TP: this is needed to not run out of memory for spherical-gdisk-planet-thermo and backreact_infl on norlx51
  PCLoad(config,AC_only_default_stream_for_taskgraphs, lonly_default_stream_for_taskgraphs__mod__gpu);
#if LGRAVITATIONAL_WAVES_HTXK
  PCLoad(config,AC_GW_Fourier_precision,lsingle_precision_ffts_for_gw_update__mod__gravitational_waves_htxk ? AC_SINGLE_PRECISION : AC_REAL_PRECISION);
#endif
  PCLoad(config,AC_error_buffer_precision,lsingle_precision_timestep ? AC_SINGLE_PRECISION : AC_REAL_PRECISION);
  PCLoad(config, AC_domain_decomposition, (int3) {nprocx,nprocy,nprocz});
  PCLoad(config, AC_ngrid, (int3){nxgrid,nygrid,nzgrid});
  PCLoad(config, AC_skip_single_gpu_optim, true);

  PCLoad(config,AC_decompose_strategy,AC_DECOMPOSE_STRATEGY_EXTERNAL);
  PCLoad(config,AC_proc_mapping_strategy,lmorton_curve ? AC_PROC_MAPPING_STRATEGY_MORTON : AC_PROC_MAPPING_STRATEGY_LINEAR);

  PCLoad(config, AC_include_3d_halo_corners, ltraining);
  PCLoad(config,AC_MPI_comm_strategy,AC_MPI_COMM_STRATEGY_DUP_USER);
  config.comm->handle = comm_pencil;

// Grid and geometry related parameters

  PCLoad(config,AC_ds,(AcReal3){dx,dy,dz});
  PCLoad(config,AC_periodic_grid,lperi);
//  Overwrites Astaroth's default formula for AC_len in case it does not cover everything like non-equidistant grids
  PCLoad(config,AC_len,lxyz);

  #if LTRAINING
  AcRealSymmetricTensor tau_hydro_means{};
  AcRealSymmetricTensor tau_hydro_stds{};
	
  //The statistics for doing the inverse are only needed when using the trained tau
  if(ltrained) denormalize("data/training/normalizer.bin", tau_hydro_means, tau_hydro_stds);

  PCLoad(config,AC_tau_hydro_means,tau_hydro_means);
  PCLoad(config,AC_tau_hydro_stds ,tau_hydro_stds);
  #endif
  PCLoad(config,AC_sparse_autotuning,lac_sparse_autotuning);

  //Calculate derived variable values based on the assignments in the initializations
  acHostUpdateParams(&config); 

  //Redirect the runtime compilation log to a log file
  if (!ltraining) config.runtime_compilation_log_dst = "ac_compilation.log";
  char build_path[18000];
  sprintf(build_path,"%s/src/astaroth/submodule/build",get_cwd());
  //Tell Astaroth that we want the runtime build to end up under this folder
  config.runtime_compilation_build_path = strdup(build_path);
  //We find the DSL compiler under the sample which did the initial compilation
  //So in case of pc_newrun -s the new sample would use the compiled compiler under the old sample.
  //We do this because changing where the compiler path would kick in unnecessary recompilations.
  config.acc_compiler_path = ACC_COMPILER_PATH;

  sprintf(build_path,"%s/src/astaroth",get_cwd());
  //This is here to skip the unnecessary make in case we run pc_newrun -s
  if(!same_path(build_path,ORIGINAL_BUILD_PATH))
  {
	  config.runtime_compilation_skip_make_if_nothing_has_changed = true;
  }

}
/***********************************************************************************************/
//TP: x,y and z macros are too general
#undef x
#undef y
#undef z
/***********************************************************************************************/
AcReal calc_dt1_courant(const AcReal t)
{
      const AcReal gpu_max_dt1 = lsingle_precision_timestep 
	      				? (AcReal)acDeviceGetOutput(acGridGetDevice(),AC_dt1_max_single_precision)
					: acDeviceGetOutput(acGridGetDevice(),AC_dt1_max);
      if (lfractional_tstep_advance__mod__cdata)
      {
	      return std::max((double)gpu_max_dt1,1./(dt_incr__mod__cdata*t));
      }
      return gpu_max_dt1;
}
/***********************************************************************************************/
AcReal GpuCalcDt(const AcReal t)
{
	acGridSynchronizeStream(STREAM_ALL);
  	acDeviceSetInput(acGridGetDevice(), AC_step_num, (PC_SUB_STEP_NUMBER) 0);
	const auto graph = acGetOptimizedDSLTaskGraph(AC_calculate_timestep);
	acGridExecuteTaskGraph(graph,1);
	acGridSynchronizeStream(STREAM_ALL);
	acDeviceSwapBuffers(acGridGetDevice());
	return calc_dt1_courant(t);
}
/***********************************************************************************************/
extern "C" void radTransfer()
{
#if Lradiation_ray_MODULE
	for(int inu = 1; inu <= nnu__mod__radiation; ++inu)
	{
  		acDeviceSetInput(acGridGetDevice(), AC_frequency_bin,inu);

//TP: need some fields (like TT) on the boundaries to calculate Srad and kappa_rho
		acGridHaloExchange();
  		auto bcs = acGetOptimizedDSLTaskGraph(boundconds);	
		acGridExecuteTaskGraph(bcs,1);

		const auto info = acDeviceGetLocalConfig(acGridGetDevice());
		const Volume offsets = (Volume)
		{
			(size_t)nxgrid == 1 ? 0 : 1,
			(size_t)nygrid == 1 ? 0 : 1,
			(size_t)nzgrid == 1 ? 0 : 1
		};
		const Volume start = acGetMinNN(info)-offsets;
		const Volume end   = acGetMaxNN(info)+offsets;
		acGridExecuteTaskGraph(acGetOptimizedDSLTaskGraph(get_source_function_and_opacity,start,end),1);
        	const auto Qintrinsic_graph = acGetOptimizedDSLTaskGraph(Qintrinsic_steps,true,Qintrinsic_bcs);
		acGridExecuteTaskGraph(Qintrinsic_graph,1);
        	const auto Qextrinsic_graph = acGetOptimizedDSLTaskGraph(Qextrinsic_steps);
		acGridExecuteTaskGraph(bcs,1);
		acGridExecuteTaskGraph(Qextrinsic_graph,1);
	}
#endif
}
/***********************************************************************************************/
extern "C" void splitUpdate(const AcReal error_tolerance, const int max_steps)
{
#if LIMPLICIT_DIFFUSION
	if(limplicit_diffusion_with_fft)
	{
#if LMAGNETIC
		if(limplicit_resistivity && lbfield)
		{
			ac_fft_split_diffusion_update(acGetF_BX(),dt,eta,acGetSPLIT_DIFFUSION_UPDATE_BUFFER_REAL(),acGetSPLIT_DIFFUSION_UPDATE_BUFFER_IMAG(),acGetSPLIT_DIFFUSION_UPDATE_BUFFER_REAL(),acGetSPLIT_DIFFUSION_UPDATE_BUFFER_IMAG());
			ac_fft_split_diffusion_update(acGetF_BY(),dt,eta,acGetSPLIT_DIFFUSION_UPDATE_BUFFER_REAL(),acGetSPLIT_DIFFUSION_UPDATE_BUFFER_IMAG(),acGetSPLIT_DIFFUSION_UPDATE_BUFFER_REAL(),acGetSPLIT_DIFFUSION_UPDATE_BUFFER_IMAG());
			ac_fft_split_diffusion_update(acGetF_BZ(),dt,eta,acGetSPLIT_DIFFUSION_UPDATE_BUFFER_REAL(),acGetSPLIT_DIFFUSION_UPDATE_BUFFER_IMAG(),acGetSPLIT_DIFFUSION_UPDATE_BUFFER_REAL(),acGetSPLIT_DIFFUSION_UPDATE_BUFFER_IMAG());
		}
#endif
#if LVISCOSITY
		if(limplicit_viscosity)
		{
			ac_fft_split_diffusion_update(acGetUUX(),dt,nu,acGetSPLIT_DIFFUSION_UPDATE_BUFFER_REAL(),acGetSPLIT_DIFFUSION_UPDATE_BUFFER_IMAG(),acGetSPLIT_DIFFUSION_UPDATE_BUFFER_REAL(),acGetSPLIT_DIFFUSION_UPDATE_BUFFER_IMAG());
			ac_fft_split_diffusion_update(acGetUUY(),dt,nu,acGetSPLIT_DIFFUSION_UPDATE_BUFFER_REAL(),acGetSPLIT_DIFFUSION_UPDATE_BUFFER_IMAG(),acGetSPLIT_DIFFUSION_UPDATE_BUFFER_REAL(),acGetSPLIT_DIFFUSION_UPDATE_BUFFER_IMAG());
			ac_fft_split_diffusion_update(acGetUUZ(),dt,nu,acGetSPLIT_DIFFUSION_UPDATE_BUFFER_REAL(),acGetSPLIT_DIFFUSION_UPDATE_BUFFER_IMAG(),acGetSPLIT_DIFFUSION_UPDATE_BUFFER_REAL(),acGetSPLIT_DIFFUSION_UPDATE_BUFFER_IMAG());
		}
#endif
	}
	else if(limplicit_diffusion_with_cg)
	{
		acGridExecuteTaskGraph(acGetOptimizedDSLTaskGraph(implicit_diffusion_init_cg),1);
		acDeviceSetInput(acGridGetDevice(),AC_implicit_diffusion_coefficient,eta);
		for(int field = 0; field < 3; ++field)
		{
			int step_num = 0;
			acDeviceSetInput(acGridGetDevice(),AC_CG_FIELD,CG_FIELD(field));
			AcReal residual = 1e100;
			bool over_max_steps = false;
			acGridExecuteTaskGraph(acGetOptimizedDSLTaskGraph(implicit_diffusion_get_residual),1);
			residual = sqrt(acDeviceGetOutput(acGridGetDevice(),AC_implicit_diffusion_residual));
			while(residual > error_tolerance && !over_max_steps)
			{
				acGridExecuteTaskGraph(acGetOptimizedDSLTaskGraph(implicit_diffusion_cg_step),1);
				acGridExecuteTaskGraph(acGetOptimizedDSLTaskGraph(implicit_diffusion_get_residual),1);
				residual = sqrt(acDeviceGetOutput(acGridGetDevice(),AC_implicit_diffusion_residual));
				++step_num;
				over_max_steps |= (max_steps > 0 && step_num >= max_steps);
			}
		}
		acDeviceSetInput(acGridGetDevice(),AC_implicit_diffusion_coefficient,nu);
		for(int field = 3; field < 6; ++field)
		{
			acDeviceSetInput(acGridGetDevice(),AC_CG_FIELD,CG_FIELD(field));
			AcReal residual = 1e100;
			acGridExecuteTaskGraph(acGetOptimizedDSLTaskGraph(implicit_diffusion_get_residual),1);
			residual = sqrt(acDeviceGetOutput(acGridGetDevice(),AC_implicit_diffusion_residual));
			bool over_max_steps = false;
			int step_num = 0;
			while(residual > error_tolerance && !over_max_steps)
			{
				acGridExecuteTaskGraph(acGetOptimizedDSLTaskGraph(implicit_diffusion_cg_step),1);
				acGridExecuteTaskGraph(acGetOptimizedDSLTaskGraph(implicit_diffusion_get_residual),1);
				residual = sqrt(acDeviceGetOutput(acGridGetDevice(),AC_implicit_diffusion_residual));
				++step_num;
				over_max_steps |= (max_steps > 0 && step_num >= max_steps);
			}
		}
	}
	else
	{
		fprintf(stderr,"Unknown method for implicit diffusion update on GPU!\n");
		fflush(stderr);
		exit(EXIT_FAILURE);
	}

#endif
}
/***********************************************************************************************/
extern "C" void registerGPU()
{

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  mesh.info = acInitInfo();
#if AC_RUNTIME_COMPILATION
#else
  AcResult res = acCheckDeviceAvailability();

  if (res == AC_FAILURE) 
  {
	  fprintf(stderr,"No devices!\n");
	  exit(EXIT_FAILURE);
  }
#endif
}
/***********************************************************************************************/
void save_stats(){
#if LTRAINING
	#include "user_constants.h"
    if(rank==0){
        std::ofstream myFile;
        std::string fileString = "running_statistics.csv";	
        myFile.open(fileString);
        myFile << "count,acc_count,in_acc_sum_x,in_acc_sum_y,in_acc_sum_z,in_acc_sum_squared_x,in_acc_sum_squared_y,in_acc_sum_squared_z,out_acc_sum_xx,out_acc_sum_yy,out_acc_sum_zz,out_acc_sum_xy,out_acc_sum_yz,out_acc_sum_xz,out_acc_sum_squared_xx,out_acc_sum_squared_yy,out_acc_sum_squared_zz,out_acc_sum_squared_xy,out_acc_sum_squared_yz,out_acc_sum_squared_xz\n";
        
        auto in_acc_sum_x = acDeviceGetOutput(acGridGetDevice(), in_acc_sum[0]);
        auto in_acc_sum_y = acDeviceGetOutput(acGridGetDevice(), in_acc_sum[1]);
        auto in_acc_sum_z = acDeviceGetOutput(acGridGetDevice(), in_acc_sum[2]);

        auto in_acc_sum_squared_x = acDeviceGetOutput(acGridGetDevice(), in_acc_sum_squared[0]);
        auto in_acc_sum_squared_y = acDeviceGetOutput(acGridGetDevice(), in_acc_sum_squared[1]);
        auto in_acc_sum_squared_z = acDeviceGetOutput(acGridGetDevice(), in_acc_sum_squared[2]);

        auto out_acc_sum_xx = acDeviceGetOutput(acGridGetDevice(), out_acc_sum[0]);
        auto out_acc_sum_yy = acDeviceGetOutput(acGridGetDevice(), out_acc_sum[1]);
        auto out_acc_sum_zz = acDeviceGetOutput(acGridGetDevice(), out_acc_sum[2]);
        auto out_acc_sum_xy = acDeviceGetOutput(acGridGetDevice(), out_acc_sum[3]);
        auto out_acc_sum_yz = acDeviceGetOutput(acGridGetDevice(), out_acc_sum[4]);
        auto out_acc_sum_xz = acDeviceGetOutput(acGridGetDevice(), out_acc_sum[5]);

        auto out_acc_sum_squared_xx = acDeviceGetOutput(acGridGetDevice(), out_acc_sum_squared[0]);
        auto out_acc_sum_squared_yy = acDeviceGetOutput(acGridGetDevice(), out_acc_sum_squared[1]);
        auto out_acc_sum_squared_zz = acDeviceGetOutput(acGridGetDevice(), out_acc_sum_squared[2]);
        auto out_acc_sum_squared_xy = acDeviceGetOutput(acGridGetDevice(), out_acc_sum_squared[3]);
        auto out_acc_sum_squared_yz = acDeviceGetOutput(acGridGetDevice(), out_acc_sum_squared[4]);
        auto out_acc_sum_squared_xz = acDeviceGetOutput(acGridGetDevice(), out_acc_sum_squared[5]);

        myFile << train_counter << "," << nxgrid*nygrid*nzgrid << "," << in_acc_sum_x << "," << in_acc_sum_y << "," << in_acc_sum_z << "," << in_acc_sum_squared_x << "," << in_acc_sum_squared_y << "," << in_acc_sum_squared_z << "," <<  out_acc_sum_xx << "," <<  out_acc_sum_yy << "," <<  out_acc_sum_zz << "," <<  out_acc_sum_xy << "," << out_acc_sum_yz << "," << out_acc_sum_xz << "," << out_acc_sum_squared_xx << "," << out_acc_sum_squared_yy << "," << out_acc_sum_squared_zz << "," << out_acc_sum_squared_xy << "," << out_acc_sum_squared_yz << "," << out_acc_sum_squared_xz << "\n";

        myFile.close();
    }
#endif
}
/***********************************************************************************************/
void read_stats(){
#if LTRAINING
    if(rank==0){
        std::ifstream statistics("running_statistics.csv");
        int count;
        int domain;
        std::vector<AcReal> in_sum;
        std::vector<AcReal> in_sum_sq;
        std::vector<AcReal> out_sum;
        std::vector<AcReal> out_sum_sq;

		if (!statistics.is_open()){
			fprintf(stderr, "Could not open running_statistics file");
			fflush(stderr);
			exit(EXIT_FAILURE);
		}
        
        std::string line;
        std::getline(statistics, line);

        while(std::getline(statistics, line)){
            std::stringstream ss(line);
            std::string val;
            std::vector<std::string> row;
            
            while (std::getline(ss, val, ',')){
                row.push_back(val);
            }
            
            count = std::stoi(row[0]);
            domain = std::stoi(row[1]);

            in_sum.push_back(std::stod(row[2]));
            in_sum.push_back(std::stod(row[3]));
            in_sum.push_back(std::stod(row[4]));

            in_sum_sq.push_back(std::stod(row[5]));
            in_sum_sq.push_back(std::stod(row[6]));
            in_sum_sq.push_back(std::stod(row[7]));


            out_sum.push_back(std::stod(row[8]));
            out_sum.push_back(std::stod(row[9]));
            out_sum.push_back(std::stod(row[10]));
            out_sum.push_back(std::stod(row[11]));
            out_sum.push_back(std::stod(row[12]));
            out_sum.push_back(std::stod(row[13]));

            out_sum_sq.push_back(std::stod(row[14]));
            out_sum_sq.push_back(std::stod(row[15]));
            out_sum_sq.push_back(std::stod(row[16]));
            out_sum_sq.push_back(std::stod(row[17]));
            out_sum_sq.push_back(std::stod(row[18]));
            out_sum_sq.push_back(std::stod(row[19]));
        }
        
    }
#endif
}
/***********************************************************************************************/
extern "C" void tf_create_model_c_api(const char *model_name, const char* config_file_path, int comm_fint, bool ldist){
#if LTRAINING
	int ndevices = 0;
	
	if(acGetDeviceCount(&ndevices) != cudaSuccess)
	{
			fprintf(stderr, "initialize_training, acGetDeviceCount failed");
			fflush(stderr);
			exit(EXIT_FAILURE);
	}

  if(acSetDevice(rank % ndevices) != cudaSuccess)
  {
  		    fprintf(stderr,"Was not able to set device id!\n");
			fflush(stderr);
 			exit(EXIT_FAILURE);
  }

  acLogFromRootProc(rank,"CONFIG_FILE: %s\n", config_file_path);
  acLogFromRootProc(rank,"Model_NAME: %s\n", model_name);
	
	if(ldist){
  	    comm_pencil = MPI_Comm_f2c(comm_fint);
		bool success = torch_create_distributed_model_CAPI(model_name, config_file_path, comm_pencil, rank % ndevices);
		if (success != 0){
			acLogFromRootProc(rank, "Error when creating distributed model: %s from config file: %s\n", model_name, config_file_path);
		}
		else{
  		acLogFromRootProc(rank,"CREATED DISTRIBUTED TRAINING\n");
		}
	}

	else{
		bool success = torch_create_model_CAPI(model_name, config_file_path, 0);
		if (success != 0){
			acLogFromRootProc(rank, "Error when creating the model: %s from config file: %s\n", model_name, config_file_path);
		}
		else{
  		acLogFromRootProc(rank,"CREATED SINGLE TRAINING\n");
		}
	}
	fflush(stderr);
	fflush(stdout);
#endif
}
/***********************************************************************************************/
extern "C" void tf_load_model_c_api(const char* name, const char* fname){
#if LTRAINING
	bool success = torch_load_CAPI(name, fname);

	if (success != 0){
		acLogFromRootProc(rank, "Error when loading model %s\n", name);
	}
	else{
		acLogFromRootProc(rank, "Torchfort model %s loaded succesfully", name);
	}
	fflush(stdout);
	fflush(stderr);
#endif
}
/***********************************************************************************************/
extern "C" void tf_load_model_checkpoint_c_api(const char* name, const char* checkpoint_dir){
#if LTRAINING
	int64_t step_train=0;
	int64_t step_inference=0;
	bool success = torch_load_checkpoint_CAPI(name, checkpoint_dir, &step_train, &step_inference);
	if (success != 0){
		acLogFromRootProc(rank, "Error when loading model %s\n", name);
	}
	else{
		acLogFromRootProc(rank, "Torchfort model %s loaded succesfully\n", name);
	}
    

	fflush(stdout);
	fflush(stderr);
#endif
}
/***********************************************************************************************/
extern "C" void tf_save_model_c_api(const char* name, const char* fname){
#if LTRAINING
	bool success = torch_save_model_CAPI(name, fname);
	if(success != 0){
		acLogFromRootProc(rank, "save_model: Error when saving ML model: %s\n", name);
	}
	else{
        save_stats();
		acLogFromRootProc(rank, "save_model: Saving ML model to: %s\n", fname);
	}
	fflush(stdout);
	fflush(stderr);
#endif
}
/***********************************************************************************************/
extern "C" void tf_save_checkpoint_c_api(const char* name, const char* checkpoint_dir){
#if LTRAINING
	bool success = torch_save_checkpoint_CAPI(name, checkpoint_dir);
	if(success != 0){
		acLogFromRootProc(rank, "save_checkpoint: Error when checkpointing ML model: %s in directory: %s\n", name, checkpoint_dir);
	}
	else{
        save_stats();
		acLogFromRootProc(rank, "save_checkpoint: Checkpoint ML model to: %s\n", checkpoint_dir);
	}
	fflush(stdout);
	fflush(stderr);
#endif
}
/***********************************************************************************************/
// used as a flag to check if we need to calcualte the taus & uumean again
bool called_training = false;

// flags for printing if we are training or infering
bool calling_train = false;
bool calling_infer = false;

int randomNumber;
std::vector<float>val_loss;
std::vector<float>train_loss;
std::vector<double>train_time;
std::vector<double>val_time;
std::vector<int>train_nts;

bool loaded_stats = false;
bool snap_print = false;
/***********************************************************************************************/
void denormalize(std::string filename, AcRealSymmetricTensor &tau_means, AcRealSymmetricTensor &tau_stds)
{
	if(!loaded_stats){
		loaded_stats = true;
		std::ifstream f(filename, std::ios::binary);
		if (!f.is_open()){
			fprintf(stderr, "Could not open stats file");
			fflush(stderr);
			exit(EXIT_FAILURE);
		}

		float acc_count;
		float num_acc;
		std::vector<float> acc_sum;
		std::vector<float> acc_sum_squared;
		
		while(f.peek() != EOF){
			uint32_t name_len;
			f.read(reinterpret_cast<char*>(&name_len), sizeof(name_len));
			if (f.eof()) break;
			std::string name;
			name.resize(name_len);
			f.read(&name[0], name_len);

			std::vector<int> shape;
			std::vector<float> data;

			uint32_t ndim = 1;
			f.read(reinterpret_cast<char*>(&ndim), sizeof(ndim));
			shape.resize(ndim);
			f.read(reinterpret_cast<char*>(shape.data()), ndim * sizeof(uint32_t));

			size_t num_elem = 1;
			for (auto s : shape) num_elem *=s;
			data.resize(num_elem);
			f.read(reinterpret_cast<char*>(data.data()), num_elem * sizeof(float));

			if (name == "acc_count"){
				acc_count = data[0];
			}
			
			if (name == "num_acc"){
				num_acc = data[0];
			}

			if(name == "acc_sum"){
				acc_sum = data;
			}

			if (name == "acc_sum_squared"){
				acc_sum_squared = data;
			}
		}
	
	float safe_count = std::max(acc_count, 1.0f);
	std::vector<float> means(acc_sum.size());
	for (int index = 0; index < acc_sum.size(); index++){
		means[index] = acc_sum[index] / safe_count;
	}
	std::vector<float> stds(acc_sum_squared.size());
	for (int index = 0; index < acc_sum.size(); index++){
		float var = (acc_sum_squared[index]/ safe_count) - (means[index] * means[index]);
		stds[index] = std::max(std::sqrt(var), 1e-8f);
	}

	tau_means.xx = means[0];
	tau_means.yy = means[1];
	tau_means.zz = means[2];
	tau_means.xy = means[3];
	tau_means.yz = means[4];
	tau_means.xz = means[5];


	tau_stds.xx = stds[0];
	tau_stds.yy = stds[1];
	tau_stds.zz = stds[2];
	tau_stds.xy = stds[3];
	tau_stds.yz = stds[4];
	tau_stds.xz = stds[5];

	/*
	std::cout << "printing acc_count: " << acc_count <<  "\n" << std::flush;

	std::cout << "printing num_acc: " << num_acc <<  "\n" << std::flush;

	std::cout << "printing acc_sum:" << std::flush;
	for (auto &i : acc_sum){
		std::cout << i << "";
	}
	std::cout << "\n" << std::flush;


	std::cout << "printing acc_sum_squared:" << std::flush;
	for (auto &i : acc_sum_squared){
		std::cout << i << "";
	}
	std::cout << "\n" << std::flush;
	*/
	


	}
}
/***********************************************************************************************/
extern "C" void torch_infer_c_api(int itsub)
{	
#if LTRAINING
	#include "user_constants.h"
	if(!ltrained) return;
	if(!calling_infer){
		AcRealSymmetricTensor tau_means = mesh.info[AC_tau_hydro_means];
		AcRealSymmetricTensor tau_stds  = mesh.info[AC_tau_hydro_stds];
    acLogFromRootProc(rank,"Doing inference\n");
		fprintf(stderr,"means xx: %f, yy: %f, zz: %f, xy: %f, yz: %f, xz: %f\n", tau_means.xx, tau_means.yy, tau_means.zz, tau_means.xy, tau_means.yz, tau_means.xz);
		fprintf(stderr,"stds xx: %f, yy: %f, zz: %f, xy: %f, yz: %f, xz: %f\n", tau_stds.xx, tau_stds.yy, tau_stds.zz, tau_stds.xy, tau_stds.yz, tau_stds.xz);
		fflush(stderr);
		fflush(stdout);
	}
	calling_infer = true;

	auto calc_uumean_tau = acGetOptimizedDSLTaskGraph(initialize_uumean);
	acGridExecuteTaskGraph(calc_uumean_tau, 1);

        auto bcs = acGetOptimizedDSLTaskGraph(boundconds);	
	acGridExecuteTaskGraph(bcs,1);


	AcReal* out = NULL;

	AcReal* uumean_ptr = NULL;
	AcReal* tau_infer_ptr = NULL;
	
	/*
	if (randomNumber == 0) acDeviceGetVertexBufferPtrs(acGridGetDevice(), UUMEANBatch[0].x, &uumean_ptr, &out);
	if (randomNumber == 1) acDeviceGetVertexBufferPtrs(acGridGetDevice(), UUMEANBatch[1].x, &uumean_ptr, &out);
	if (randomNumber == 2) acDeviceGetVertexBufferPtrs(acGridGetDevice(), UUMEANBatch[2].x, &uumean_ptr, &out);
	if (randomNumber == 3) acDeviceGetVertexBufferPtrs(acGridGetDevice(), UUMEANBatch[3].x, &uumean_ptr, &out);
	if (randomNumber == 4) acDeviceGetVertexBufferPtrs(acGridGetDevice(), UUMEANBatch[4].x, &uumean_ptr, &out);
	if (randomNumber == 5) acDeviceGetVertexBufferPtrs(acGridGetDevice(), UUMEANBatch[5].x, &uumean_ptr, &out);
	*/
	

	//
	//check if trained on same grid size and do not smooth the averaged out vals
	//
	//
	acDeviceGetVertexBufferPtrs(acGridGetDevice(), uumean.x, &uumean_ptr, &out);
	acDeviceGetVertexBufferPtrs(acGridGetDevice(), TAU_HYDRO_INFERRED.xx, &tau_infer_ptr, &out);

	acGridHaloExchange();
	
	torch_infer_CAPI((int[]){mx,my,mz}, uumean_ptr, tau_infer_ptr,input_channels,output_channels,"stationary",calling_train);
	
	auto descale_inf = acGetOptimizedDSLTaskGraph(descale_inferred_taus);
	acGridExecuteTaskGraph(descale_inf, 1);

        bcs = acGetOptimizedDSLTaskGraph(boundconds);	
	acGridExecuteTaskGraph(bcs,1);
#endif
}
/***********************************************************************************************/
extern "C" void torch_train_c_api(AcReal *loss_val, int itsub, double t) {
#if LTRAINING
  #include "user_constants.h"
  #include <stdlib.h>
  if(itsub != 1) return;
  
  if(!calling_train){
    acLogFromRootProc(rank,"Doing training\n");
  	fflush(stderr);
  	fflush(stdout);
  }


  called_training = true;
  calling_train = true;


  acGridExecuteTaskGraph(acGetOptimizedDSLTaskGraph(get_taus),1);

  auto bcs = acGetOptimizedDSLTaskGraph(boundconds);	
  acGridExecuteTaskGraph(bcs,1);
  	
  acGridHaloExchange();

  train_counter++;
  acDeviceSetInput(acGridGetDevice(), AC_count, train_counter);

  acGridExecuteTaskGraph(acGetOptimizedDSLTaskGraph(normalize),1);
  acGridExecuteTaskGraph(acGetOptimizedDSLTaskGraph(boundconds),1);

  AcReal* out = NULL;
  AcReal* uumean_ptr = NULL;
  AcReal* TAU_ptr = NULL;
  *loss_val = 0.1;
  
  acDeviceGetVertexBufferPtrs(acGridGetDevice(), tau_hydro.xx, &TAU_ptr, &out);
  acDeviceGetVertexBufferPtrs(acGridGetDevice(), uumean.x, &uumean_ptr, &out);
 
  acGridHaloExchange();
  torch_train_CAPI((int[]){mx,my,mz}, uumean_ptr, TAU_ptr, loss_val,input_channels,output_channels,"stationary");
  train_loss.push_back(*loss_val);
  train_nts.push_back(it);
  //print_debug();
  if (it==nt){
  	std::ofstream myFile;
  	std::string fileString = "train_loss_" + std::to_string(my_rank)  + ".csv";	
  	myFile.open(fileString);
  	myFile << "epoch,train_loss\n";
  	for (int i=0;i<train_loss.size();i++){
  		myFile << i << "," << train_loss[i] << "\n";
  	}
  	myFile.close();
		
		std::string train_sample = "train_sample_" + std::to_string(my_rank) + ".csv";
  	myFile.open(train_sample);
  	myFile << "train_sample_nts\n";
  	for (int i=0;i<train_nts.size();i++){
  		myFile << train_nts[i] << "\n";
  	}
  	myFile.close();
  }
#endif
}
/***********************************************************************************************/
float MSE(){
#if LTRAINING
	#include "user_constants.h"

	auto calc_validation = acGetOptimizedDSLTaskGraph(calc_validation_loss);
	acGridExecuteTaskGraph(calc_validation, 1);

 	auto bcs = acGetOptimizedDSLTaskGraph(boundconds);	
	acGridExecuteTaskGraph(bcs,1);

	calculated_coeff_scales = true;
	return (acDeviceGetOutput(acGridGetDevice(), AC_l2_sum))/(6*nxgrid*nygrid*nzgrid);
#else
        return 0;
#endif
}
/***********************************************************************************************/
void load_f_ode()
{
	//This seems to be a sufficiently performant way to keep ODE variables in sync:
	//For backreact_infl 256^3 subdomain:
	//Loading ode variables took: 6.08199999987846e-04
	//BeforeBoundary took: 7.85426100003406e-03 to 1.9e-2 (TODO: why does it fluctuate so much?)
        //RHS TOOK:   1.451359e-02
	//So for large enough subdomains rhs dominates
        if (n_odevars > 0)
        {
                AcReal start = MPI_Wtime();
                acDeviceSynchronizeStream(acGridGetDevice(),STREAM_DEFAULT);
                acDeviceLoad(acGridGetDevice(), STREAM_DEFAULT, mesh.info, AC_f_ode);
                acDeviceSynchronizeStream(acGridGetDevice(),STREAM_DEFAULT);
	        if(performance_logs) acLogFromRootProc(rank,"Loading ode variables took: %.14e\n",MPI_Wtime()-start);
        }
}
/***********************************************************************************************/
void
reload_dynamically_changing_variables(bool lrmv, int isubstep, double t, bool lsubstepping_in_time)
{
// Loads the values of the ODE variables
// which are are advanced on the host in dspecial_dt_ode
	load_f_ode();
// Sets the values of input parameters to the kernels.
// Does not load them to the device but instead sets values of the Astaroth config
// from which input parameters to the kernels are read from
        acDeviceSetInput(acGridGetDevice(), AC_timestep_number,it);
        acDeviceSetInput(acGridGetDevice(), AC_step_num,(PC_SUB_STEP_NUMBER) (isubstep-1));
 	acDeviceSetInput(acGridGetDevice(), AC_lrmv, lrmv);
 	acDeviceSetInput(acGridGetDevice(), AC_t, AcReal(t));
 	acDeviceSetInput(acGridGetDevice(), AC_lsubstepping_in_time, lsubstepping_in_time);
}
/***********************************************************************************************/
extern "C" void beforeBoundaryGPU(bool lrmv, int isubstep, double t, bool lsubstepping_in_time)
{
	reload_dynamically_changing_variables(lrmv,isubstep,t,lsubstepping_in_time);
        
	AcReal start_time = MPI_Wtime();
// Execute all "before-boundary-actions", which do not update the halos, by separate task graph
	acGridExecuteTaskGraph(acGetOptimizedDSLTaskGraph(AC_before_boundary_steps),1);

//Some Fields are directly calculated on the halos like yH in ioncalc.
//Could reformulate the kernels in a way that the bc is simply the same kernel as the normal calculation
//But don't want to repeat calc too often so this is a somewhat easy way to do it
	const auto steps_updating_halos = acGetOptimizedDSLTaskGraph(AC_before_boundary_steps_including_halos,
					          (Volume){0,0,0},acGetLocalMM(acGridGetLocalMeshInfo()));
	if (!acGridTaskGraphIsEmpty(steps_updating_halos))
	{
		AcTaskGraph* bcs = acGetOptimizedDSLTaskGraph(boundconds);
		acGridExecuteTaskGraph(bcs,1);                            // apply boundconds
		//This is not anymore done since it is fine for the resulting inner halos to be incorrect.
		//If those inner halos are later needed they will be anyways communicated
		//acGridHaloExchange();                                     // halo communication
		acGridExecuteTaskGraph(steps_updating_halos,1);           // f-array update
	}
#if LSELFGRAVITY
	if (t>=tstart_selfgrav)
	{
		acGridExecuteTaskGraph(acGetOptimizedDSLTaskGraph(AC_calc_selfgravity_rhs),1);
		acDeviceFFTR2C(acGridGetDevice(),acGetRHS_POISSON(),RHS_POISSON_COMPLEX);                     // FFT of rhs of Poisson eq.
    		AcMeshDims dims = acGetMeshDims(acGridGetLocalMeshInfo());
  		acGridLaunchKernel(STREAM_DEFAULT, selfgravity_poisson_solve, dims.n0, dims.n1);
		acGridSynchronizeStream(STREAM_ALL);
		acDeviceFFTC2R(acGridGetDevice(),SELFGRAVITY_POTENTIAL_COMPLEX,acGetSELFGRAVITY_POTENTIAL()); // FFT backtransform. of solution
		//TP: A placeholder for iterative solvers, in the future choose the number of solving steps based on the norm of the residual
		/**
		for (int i = 0; i < 100; ++i)
		{
			acGridExecuteTaskGraph(acGetOptimizedDSLTaskGraph(AC_sor_step),1);
		}
		**/
		acGridExecuteTaskGraph(acGetOptimizedDSLTaskGraph(AC_calc_final_potential),1);
	}

#endif
#if LNEWTON_COOLING
	acGridExecuteTaskGraph(acGetOptimizedDSLTaskGraph(AC_integrate_tau),1);
#endif
	if(performance_logs) acLogFromRootProc(rank,"BeforeBoundary took: %.14e\n",MPI_Wtime()-start_time);
}
/***********************************************************************************************/
bool idx_init = false;
std::vector<size_t> idx_cache;
std::vector<double> buffer;
std::stack<int> snaps({20000, 15000, 10000, 5000, 1000, 500, 100, 50, 10});


extern "C" void print_debug() {
/*
	if (snaps.empty() || it< snaps.top()) return;

	while (!snaps.empty() && it >= snaps.top()) {
		snaps.pop();
	}
*/

#if LTRAINING
    #include "user_constants.h"
		
		std::string fname = "snapshots/snapshot_multi_normalized_"+ std::to_string(nxgrid*nygrid*nzgrid)  +"_rank_" + std::to_string(rank) + "_it_" + std::to_string(it) + ".bin";
		std::ifstream infile(fname, std::ios::binary);
		if (infile.good()) return;
		
		counter = it;
		
		/*
		MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
		acDeviceSetInput(acGridGetDevice(), AC_ranNum, randomNumber);
		*/

	  acGridHaloExchange();
          copyFarray(mesh.vertex_buffer[0]);

          AcMeshDims dims = acGetMeshDims(acGridGetLocalMeshInfo());

		/*
  	const auto DEVICE_VTXBUF_IDX = [&](const int x, const int y, const int z)
  				{
					return acVertexBufferIdx(x,y,z,mesh.info);
  				};
		*/


		int x_size = (dims.m1.x - dims.m0.x);
		int y_size = (dims.m1.y - dims.m0.y);
		int z_size = (dims.m1.z - dims.m0.z);

		const size_t n_points =  x_size * y_size * z_size;
		const int n_fields = 22;
		
		if(!idx_init){
			idx_cache.reserve(n_points);
			buffer.reserve(n_points * n_fields);

    	for (size_t i = dims.m0.x; i < dims.m1.x; i++) {
      	for (size_t j = dims.m0.y; j < dims.m1.y; j++) {
        	for (size_t k = dims.m0.z; k < dims.m1.z; k++) {
						idx_cache.push_back(acVertexBufferIdx(i,j,k,mesh.info));
					}
				}
			}
			idx_init = true;
		}
		
		buffer.clear();

    const AcReal* tau_xx_buf = mesh.vertex_buffer[tau_hydro.xx];
    const AcReal* tau_inferred_xx_buf = mesh.vertex_buffer[TAU_HYDRO_INFERRED.xx];
    const AcReal* tau_yy_buf = mesh.vertex_buffer[tau_hydro.yy];
    const AcReal* tau_inferred_yy_buf = mesh.vertex_buffer[TAU_HYDRO_INFERRED.yy];
    const AcReal* tau_zz_buf = mesh.vertex_buffer[tau_hydro.zz];
    const AcReal* tau_inferred_zz_buf = mesh.vertex_buffer[TAU_HYDRO_INFERRED.zz];
    const AcReal* tau_xy_buf = mesh.vertex_buffer[tau_hydro.xy];
    const AcReal* tau_inferred_xy_buf = mesh.vertex_buffer[TAU_HYDRO_INFERRED.xy];
    const AcReal* tau_yz_buf = mesh.vertex_buffer[tau_hydro.yz];
    const AcReal* tau_inferred_yz_buf = mesh.vertex_buffer[TAU_HYDRO_INFERRED.yz];
    const AcReal* tau_xz_buf = mesh.vertex_buffer[tau_hydro.xz];
    const AcReal* tau_inferred_xz_buf = mesh.vertex_buffer[TAU_HYDRO_INFERRED.xz];
    const AcReal* uumean_x_buf = mesh.vertex_buffer[uumean.x];
    const AcReal* uumean_y_buf = mesh.vertex_buffer[uumean.y];
    const AcReal* uumean_z_buf = mesh.vertex_buffer[uumean.z];
    const AcReal* uux_buf = mesh.vertex_buffer[UUX];
    const AcReal* uuy_buf = mesh.vertex_buffer[UUY];
    const AcReal* uuz_buf = mesh.vertex_buffer[UUZ];

    const double it_double = static_cast<double>(it);
    const int yz_size = y_size * z_size;

		for(size_t idx_i=0; idx_i < idx_cache.size(); ++idx_i){
			
			const size_t idx = idx_cache[idx_i];

			const int i = dims.m0.x + (idx_i / (yz_size));
			const int j = dims.m0.y + ((idx_i / z_size) % y_size);
			const int k = dims.m0.z + (idx_i % z_size);
			
			 buffer.push_back(it_double);

       buffer.push_back(tau_xx_buf[idx]);
       buffer.push_back(tau_inferred_xx_buf[idx]);

       buffer.push_back(tau_yy_buf[idx]);
       buffer.push_back(tau_inferred_yy_buf[idx]);

       buffer.push_back(tau_zz_buf[idx]);
       buffer.push_back(tau_inferred_zz_buf[idx]);

       buffer.push_back(tau_xy_buf[idx]);
       buffer.push_back(tau_inferred_xy_buf[idx]);

       buffer.push_back(tau_yz_buf[idx]);
       buffer.push_back(tau_inferred_yz_buf[idx]);

       buffer.push_back(tau_xz_buf[idx]);
       buffer.push_back(tau_inferred_xz_buf[idx]);


       buffer.push_back(uumean_x_buf[idx]);
       buffer.push_back(uumean_y_buf[idx]);
       buffer.push_back(uumean_z_buf[idx]);

       buffer.push_back(uux_buf[idx]);
       buffer.push_back(uuy_buf[idx]);
       buffer.push_back(uuz_buf[idx]);

       buffer.push_back(static_cast<double>(i));
       buffer.push_back(static_cast<double>(j));
       buffer.push_back(static_cast<double>(k));
		}
	

		//std::ostringstream fname;
    //fname << "snapshots/snapshot_38_rank_" + std::to_string(my_rank) + "_it_" + std::to_string(it) + ".bin";
    std::ofstream out(fname, std::ios::binary);
    out.write(reinterpret_cast<const char*>(buffer.data()), buffer.size() * sizeof(double));
    out.close();

	#endif
	if(called_training) called_training = false;
}
/***********************************************************************************************/
void
GW_update(const AcReal dt_gw)
{
#if LGRAVITATIONAL_WAVES_HTXK
	//We set the device id here since another thread than the master might be executing this
  	if(acSetDevice(device_id) != cudaSuccess)
  	{
  	        fprintf(stderr,"Was not able to set device id!\n");
  	        exit(EXIT_FAILURE);
  	}
	acDeviceFFTR2PlanarBatched(acGridGetDevice(), acGetF_STRESS_0(),acGetAC_tpq_re__mod__gravitational_waves_htxk_0(),acGetAC_tpq_im__mod__gravitational_waves_htxk_0(),6);
	//TP: do this if you want to test the performance of utilizing the conjugate symmetry
	//acDeviceFFTR2HermitianPlanarBatched(acGridGetDevice(), acGetF_STRESS_0(),acGetAC_tpq_re__mod__gravitational_waves_htxk_0(),acGetAC_tpq_im__mod__gravitational_waves_htxk_0(),6,STREAM_10);
        acDeviceSynchronizeStream(acGridGetDevice(),STREAM_10);
	acDeviceSetInput(acGridGetDevice(),AC_dt,dt_gw);
	acGridExecuteTaskGraph(GW_timestep_graph,1);
	acDeviceSetInput(acGridGetDevice(),AC_dt,dt);
#endif
}
/***********************************************************************************************/
extern "C" void afterSubStepGPU()
{
	static AcReal dt_gw=0.0;
#if LGRAVITATIONAL_WAVES_HTXK
	if (acDeviceGetInput(acGridGetDevice(), AC_step_num) == PC_FIRST_SUB_STEP)
	{
	   dt_gw += dt;
	   if((it+1) % ntimesteps_per_gw_step__mod__gravitational_waves_htxk == 0)
	   {
	     GW_update(dt_gw);
	     dt_gw = 0.0;
	   }
	   //Do this if you want to test the performance of overlapping GW FFTs with the normal RHS
	   //Usually this has given subpar performance
	   //GW_thread = std::thread(GW_update);
	}
#endif
	acGridExecuteTaskGraph(acGetOptimizedDSLTaskGraph(AC_after_timestep),1);
}
/***********************************************************************************************/
//Tested to work with solar-atmosphere-magnetic
void
fourier_boundary_conditions()
{
	if(luses_aa_pwd_top)
	{
		acKernelInputParams params{};
		params.bc_aa_pwd_kernel.boundary = BOUNDARY_Z_TOP;
		params.bc_aa_pwd_kernel.topbot  = AC_top;
		const AcKernel kernel = acGetOptimizedKernel(bc_aa_pwd_kernel,params);
		const size_t boundary_z = (size_t)mesh.info[AC_nlocal_max].z-1;
		for(int ghost = 1; ghost <= NGHOST; ++ghost)
		{
			acDeviceFFTR2PlanarXY(acGridGetDevice(), acGetAAX(), acGetAX_FOURIER_REAL(), acGetAX_FOURIER_IMAG(), boundary_z-ghost);
			acDeviceFFTR2PlanarXY(acGridGetDevice(), acGetAAY(), acGetAY_FOURIER_REAL(), acGetAY_FOURIER_IMAG(), boundary_z-ghost);
			acDeviceFFTR2PlanarXY(acGridGetDevice(), acGetAAZ(), acGetAZ_FOURIER_REAL(), acGetAZ_FOURIER_IMAG(), boundary_z-ghost);
		}
		acDeviceLaunchKernel(acGridGetDevice(), STREAM_DEFAULT, kernel, (Volume){NGHOST,NGHOST,boundary_z+1}, 
				(Volume){(size_t)mesh.info[AC_nlocal_max].x,(size_t)mesh.info[AC_nlocal_max].y, boundary_z+2}
				);
  		acDeviceSynchronizeStream(acGridGetDevice(),STREAM_DEFAULT);
		for(int ghost = 1; ghost <= NGHOST; ++ghost)
		{
			acDeviceFFTBackwardTransformPlanar2RXY(acGridGetDevice(),  acGetAX_FOURIER_REAL(), acGetAX_FOURIER_IMAG(), acGetAAX(), boundary_z+ghost);
			acDeviceFFTBackwardTransformPlanar2RXY(acGridGetDevice(),  acGetAY_FOURIER_REAL(), acGetAY_FOURIER_IMAG(), acGetAAY(), boundary_z+ghost);
			acDeviceFFTBackwardTransformPlanar2RXY(acGridGetDevice(),  acGetAZ_FOURIER_REAL(), acGetAZ_FOURIER_IMAG(), acGetAAZ(), boundary_z+ghost);
		}
		acGridExecuteTaskGraph(boundary_z_halo_exchange_graph,1);
	}
	if(luses_aa_pot2_top)
	{
		acKernelInputParams params{};
		params.bc_aa_pot_kernel.boundary = BOUNDARY_Z_TOP;
		params.bc_aa_pot_kernel.topbot  = AC_top;
		const AcKernel kernel = acGetOptimizedKernel(bc_aa_pot_kernel,params);
		const size_t boundary_z = (size_t)mesh.info[AC_nlocal_max].z-1;
		acDeviceFFTR2PlanarXY(acGridGetDevice(), acGetAAX(), acGetAX_FOURIER_REAL(), acGetAX_FOURIER_IMAG(), boundary_z);
		acDeviceFFTR2PlanarXY(acGridGetDevice(), acGetAAY(), acGetAY_FOURIER_REAL(), acGetAY_FOURIER_IMAG(), boundary_z);
		acDeviceFFTR2PlanarXY(acGridGetDevice(), acGetAAZ(), acGetAZ_FOURIER_REAL(), acGetAZ_FOURIER_IMAG(), boundary_z);
		acDeviceLaunchKernel(acGridGetDevice(), STREAM_DEFAULT, kernel, (Volume){NGHOST,NGHOST,boundary_z+1}, 
				(Volume){(size_t)mesh.info[AC_nlocal_max].x,(size_t)mesh.info[AC_nlocal_max].y, boundary_z+2}
				);
  		acDeviceSynchronizeStream(acGridGetDevice(),STREAM_DEFAULT);
		for(int ghost = 1; ghost <= NGHOST; ++ghost)
		{
			acDeviceFFTBackwardTransformPlanar2RXY(acGridGetDevice(),  acGetAX_FOURIER_REAL(), acGetAX_FOURIER_IMAG(), acGetAAX(), boundary_z+ghost);
			acDeviceFFTBackwardTransformPlanar2RXY(acGridGetDevice(),  acGetAY_FOURIER_REAL(), acGetAY_FOURIER_IMAG(), acGetAAY(), boundary_z+ghost);
			acDeviceFFTBackwardTransformPlanar2RXY(acGridGetDevice(),  acGetAZ_FOURIER_REAL(), acGetAZ_FOURIER_IMAG(), acGetAAZ(), boundary_z+ghost);
		}
		acGridExecuteTaskGraph(boundary_z_halo_exchange_graph,1);
	}
}
/***********************************************************************************************/
void
update_forcing(const int isubstep)
{
#if LFORCING
  //Update forcing params
   if (lsecond_force) 
   {
	   fprintf(stderr,"Second forcing force not yet implemented on GPU!\n");
	   exit(EXIT_FAILURE);
   }
   if (isubstep == num_substeps) 
	 {
	   forcing_params.Update();  // calculate on CPU and load into GPU
	 }
#endif
}
/***********************************************************************************************/
void
prepare_rhs(const int isubstep, double t)
{
  fourier_boundary_conditions();
  acDeviceSetInput(acGridGetDevice(), AC_step_num,(PC_SUB_STEP_NUMBER) (isubstep-1));
  if (lshear) 
  {
	  acDeviceSetInput(acGridGetDevice(), AC_shear_delta_y, deltay);
  }
  acDeviceSetInput(acGridGetDevice(), AC_t,(AcReal)t);
}
/***********************************************************************************************/
void
set_timestep(double t)
{
    static bool lfirst_timestep_calculated = false;
    //Lcpu_timestep_on_gpu enables the same timestep as PC when testing
    if (ldt && lcourant_dt && (!lfirst_timestep_calculated || lcpu_timestep_on_gpu)) 
    {
            dt1_interface = GpuCalcDt(AcReal(t));
    	lfirst_timestep_calculated = true;
    }
    //Done in this more complex manner to ensure the actually integrated time and the time reported by Pencil agree
    //if we call set_dt after the first timestep there would be slight shift in dt what Pencil sees and what is actually used for time integration
    if (ldt) 
    {
      set_dt(dt1_interface);
    }
    acDeviceSetInput(acGridGetDevice(), AC_dt,dt);
}
/***********************************************************************************************/
void
calc_timestep(double t)
{
  constexpr AcReal unit = 1.0;
  AcReal dt1_;
  if (!lcourant_dt)
  {
    const AcReal maximum_error = lsingle_precision_timestep 
	    				? ((AcReal)acDeviceGetOutput(acGridGetDevice(), AC_maximum_error_single_precision))/eps_rkf
	    				: ((AcReal)acDeviceGetOutput(acGridGetDevice(), AC_maximum_error))/eps_rkf;
    AcReal dt_;
    const AcReal dt_increase=-unit/(itorder+dtinc);
    const AcReal dt_decrease=-unit/(itorder-dtdec);
    constexpr AcReal safety=(AcReal)0.95;
    if (maximum_error > 1)
    {
    	// Step above error threshold so decrease the next time step
    	const AcReal dt_temp = safety*dt*pow(maximum_error,dt_decrease);
    	// Don't decrease the time step by more than a factor of ten
        constexpr AcReal decrease_factor = (AcReal)0.1;
    	dt_ = sign(max(abs(dt_temp), decrease_factor*abs(dt)), dt);
    } 
    else
    {
    	dt_ = dt*pow(maximum_error,dt_increase);
    }
    set_next_dt(dt_);
    dt1_ = unit/dt_;
  }
  else 
  {
    dt1_ = calc_dt1_courant(AcReal(t));
  }
  dt1_interface = dt1_;
}
/***********************************************************************************************/
extern "C" void substepGPU(int isubstep, double t)
//
//  Do the 'isubstep'th integration step on all GPUs on the node and handle boundaries.
//
{
  prepare_rhs(isubstep,t);
  if (isubstep == 1) 
  {
#if LGRAVITATIONAL_WAVES_HTXK
    if(GW_thread.joinable())
    {
            GW_thread.join();
    }
#endif
    set_timestep(t);
#if LGRAVITATIONAL_WAVES_HTXK
    if(lsplit_gw_rhs_from_rest_on_gpu)
    {
          acGridExecuteTaskGraph(acGetOptimizedDSLTaskGraph(AC_GW_rhs),1);
    }
#endif
  }
  //Important that forcing vectors are updated only after the timestep is calculated
  //since the timestep is used in the calculation
  update_forcing(isubstep);

  AcTaskGraph *rhs =  acGetOptimizedDSLTaskGraph(AC_rhs);
  auto start = MPI_Wtime();
  acGridExecuteTaskGraph(rhs, 1);
  auto end = MPI_Wtime();
  if (performance_logs) acLogFromRootProc(rank,"RHS TOOK: %14e\n",end-start);
  if (ldt && (   (isubstep == 5 && !lcourant_dt) 
              || (isubstep == 1 &&  lcourant_dt)
             )
     )
  {
	calc_timestep(t);
  }
#if LTIMEAVGS
  if(isubstep == num_substeps)
  {
	  acGridExecuteTaskGraph(acGetOptimizedDSLTaskGraph(AC_update_timeavgs));
  }
#endif
  return;
}
/***********************************************************************************************/
void copyFarray(AcReal* f)
{
  #include "user_constants.h"

  acGridSynchronizeStream(STREAM_ALL);
  //TP: for now only copy the advanced fields back
  //TODO: should auxiliaries be needed on the GPU like e.g. Shock be copied? They can always be recomputed on the host if needed
  //If doing training we read all since we might want TAU components to calculate e.g. validation error
  const int all_vtxbufs = acGetNumFields();
  const int end = ltraining ? acGetNumFields(): 
	  	  lread_all_vars_from_device ? mfarray : mvar;

  AcMesh* dst = &mesh;
  AcMesh tmp;
  if (nghost > 0 && (dimensionality <= 1 || (dimensionality == 2 && nygrid == 1)))
  {
  	acHostMeshCopy(mesh, &tmp);
	dst = &tmp;
  }
  for (int i = 0; i < all_vtxbufs; ++i)
  {
	  //Have to specialize for training for now since we are reading fields that do not exist in Fortran: TODO: allocate fields on the CPU in Fortran
	  const int index = ltraining ? i : 
		  	    i < mvar ? i : maux_vtxbuf_index[i];
	  bool read_var = (i < end);
	  read_var &= (index != -1);
	  read_var |= (read_vtxbuf_from_gpu[i]);
	  read_var |= ltraining;
          if (!read_var) continue;
	  acDeviceStoreVertexBuffer(acGridGetDevice(),STREAM_DEFAULT,VertexBufferHandle(index),dst);
  }
  acGridSynchronizeStream(STREAM_ALL);
  // Astaroth does not allocate ghost zones for inactive dimensions for 1d and 2d simulations, and unlike for xy  we cannot simply offset into the farray so have to manually copy the values
  //    This is fine since 1d simulations are anyways mainly for testing
  if (nghost > 0 && dimensionality == 1)
  {
  	for (int i = 0; i < end; ++i)
  	{
		if (i >= mvar && maux_vtxbuf_index[i] == -1) continue;
 		if (nxgrid != 1)
		{
			for (int x = 0; x < nx; ++x)
			{
				const size_t f_index = (NGHOST+x)+ mx*(NGHOST + my*(NGHOST));
				mesh.vertex_buffer[i][f_index] = dst->vertex_buffer[i][NGHOST+x];
			}
		}
 		if (nygrid != 1)
		{
			for (int y = 0; y < ny; ++y)
			{
				const size_t f_index = (NGHOST)+ mx*((y+NGHOST) + my*(NGHOST));
				mesh.vertex_buffer[i][f_index] = dst->vertex_buffer[i][NGHOST+y];
			}
		}
 		if (nzgrid != 1)
		{
			for (int z = 0; z < nz; ++z)
			{
				const size_t f_index = (NGHOST)+ mx*(NGHOST + my*(z+NGHOST));
				mesh.vertex_buffer[i][f_index] = dst->vertex_buffer[i][NGHOST+z];
			}
		}
	}
  	acHostMeshDestroy(&tmp);
  }
  if (nghost > 0 && dimensionality == 0)
  {
    	for (int i = 0; i < end; ++i)
  	{
		const int index = (i < mvar) ? i : maux_vtxbuf_index[i];
		if (index == -1) continue;
		const size_t f_index = (NGHOST) + mx*(NGHOST + my*(NGHOST));
		const size_t ac_index = 0;
 		mesh.vertex_buffer[i][f_index] = dst->vertex_buffer[i][ac_index];
	}
  }
  if (nghost > 0 && dimensionality == 2 && nygrid == 1)
  {
    	for (int i = 0; i < end; ++i)
  	{
		const int index = (i < mvar) ? i : maux_vtxbuf_index[i];
		if (index == -1) continue;
		for (int x = 0; x < nx; ++x)
		{
			for (int z = 0; z < nz; ++z)
			{
				const size_t f_index = (x+NGHOST) + mx*(NGHOST + my*(z+NGHOST));
				const size_t ac_index = (x+NGHOST) + mx*(z+NGHOST);
 				mesh.vertex_buffer[index][f_index] = dst->vertex_buffer[index][ac_index];
			}
		}
	}
  	acHostMeshDestroy(&tmp);
  }
}
/***********************************************************************************************/
extern "C" void getFArrayIn(AcReal **p_f_in)
{
  return;   //!!!
  AcReal* out = NULL;

  AcReal* uux_ptr = NULL;
  AcReal* uuy_ptr = NULL;
  AcReal* uuz_ptr = NULL;

  acDeviceGetVertexBufferPtrs(acGridGetDevice(),UUX,&uux_ptr,&out);
  acDeviceGetVertexBufferPtrs(acGridGetDevice(),UUY,&uuy_ptr,&out);
  acDeviceGetVertexBufferPtrs(acGridGetDevice(),UUZ,&uuz_ptr,&out);
  if (uux_ptr + mw != uuy_ptr) fprintf(stderr, "UU not contiguous\n");
  if (uuy_ptr + mw != uuz_ptr) fprintf(stderr, "UU not contiguous\n");
  acDeviceGetVertexBufferPtrs(acGridGetDevice(),VertexBufferHandle(0),p_f_in,&out);
}
/***********************************************************************************************/
extern "C" void copyVBApointers(AcReal **in, AcReal **out)
{
  AcReal* uux_ptr = NULL;
  AcReal* uuy_ptr = NULL;
  AcReal* uuz_ptr = NULL;

  acDeviceGetVertexBufferPtrs(acGridGetDevice(),UUX,&uux_ptr,out);
  acDeviceGetVertexBufferPtrs(acGridGetDevice(),UUY,&uuy_ptr,out);
  acDeviceGetVertexBufferPtrs(acGridGetDevice(),UUZ,&uuz_ptr,out);
  if (uux_ptr + mw != uuy_ptr) fprintf(stderr, "UU not contiguous\n");
  if (uuy_ptr + mw != uuz_ptr) fprintf(stderr, "UU not contiguous\n");

  acDeviceGetVertexBufferPtrs(acGridGetDevice(),VertexBufferHandle(0),in,out);
}
/***********************************************************************************************/
void autotune_all_integration_substeps()
{
  for (int i = 1; i <= num_substeps; ++i)
  {
  	acDeviceSetInput(acGridGetDevice(), AC_step_num,(PC_SUB_STEP_NUMBER)(i-1));
        if (rank==0 && ldebug) printf("memusage before GetOptimizedDSLTaskGraph= %f MBytes\n", acMemUsage()/1024.);
  	acDeviceSetInput(acGridGetDevice(), AC_lrmv,false);
	acGetOptimizedDSLTaskGraph(AC_rhs);
  	acDeviceSetInput(acGridGetDevice(), AC_lrmv,true);
	acGetOptimizedDSLTaskGraph(AC_rhs);
        if (rank==0 && ldebug) printf("memusage after GetOptimizedDSLTaskGraph= %f MBytes\n", acMemUsage()/1024.);
 	beforeBoundaryGPU(true,i,0.0,false);
        beforeBoundaryGPU(false,i,0.0,false);
  }
  radTransfer();
  //This is here mainly because GW_update will allocate buffers on the first call
  //and I don't those allocations to mess with timing measurements for short runs
  GW_update(0.0);
  splitUpdate(1e11,1);
}
/***********************************************************************************************/
AcMesh
get_f_src()
{
  if (nghost > 0 && dimensionality == 1)
  {
        AcMesh src;
  	acHostMeshCopy(mesh, &src);
    	for (int i = 0; i < mfarray; ++i)
  	{
		const int index = (i < mvar) ? i : maux_vtxbuf_index[i];
		if (index == -1) continue;
 		if (nxgrid != 1)
		{
			for (int x = 0; x < nx; ++x)
			{
				const size_t f_index = (NGHOST+x)+ mx*(NGHOST + my*(NGHOST));
				src.vertex_buffer[index][NGHOST+x]  = mesh.vertex_buffer[index][f_index];
			}
		}
		else if (nygrid != 1)
		{
			for (int y = 0; y < ny; ++y)
			{
				const size_t f_index = (NGHOST)+ mx*((NGHOST+y) + my*(NGHOST));
				src.vertex_buffer[index][NGHOST+y]  = mesh.vertex_buffer[index][f_index];
			}
		}
		else if (nzgrid != 1)
		{
			for (int z = 0; z < nz; ++z)
			{
				const size_t f_index = (NGHOST)+ mx*(NGHOST + my*(z+NGHOST));
				src.vertex_buffer[index][NGHOST+z]  = mesh.vertex_buffer[index][f_index];
			}
		}
	}
	return src;
  }
  if (nghost > 0 && dimensionality == 2 && nygrid == 1)
  {
        AcMesh src;
  	acHostMeshCopy(mesh, &src);
    	for (int i = 0; i < mfarray; ++i)
  	{
		const int index = (i < mvar) ? i : maux_vtxbuf_index[i];
		if (index == -1) continue;
		for (int x = 0; x < nx; ++x)
		{
			for (int z = 0; z < nz; ++z)
			{
				const size_t f_index = (x+NGHOST) + mx*(NGHOST + my*(z+NGHOST));
				const size_t ac_index = (x+NGHOST) + mx*(z+NGHOST);
				src.vertex_buffer[index][ac_index] = mesh.vertex_buffer[index][f_index];
			}
		}
	}
	return src;
  }
  if (nghost > 0 && dimensionality == 0)
  {
	AcMesh src;
  	acHostMeshCopy(mesh, &src);
    	for (int i = 0; i < mfarray; ++i)
  	{
		const int index = (i < mvar) ? i : maux_vtxbuf_index[i];
		if (index == -1) continue;
		const size_t f_index = (NGHOST) + mx*(NGHOST + my*(NGHOST));
		const size_t ac_index = 0;
		src.vertex_buffer[index][ac_index] = mesh.vertex_buffer[index][f_index];
	}
	return src;
  }
  return mesh;
}
/***********************************************************************************************/
extern "C" void loadFarray()
{
  AcMesh src = get_f_src();
  acGridSynchronizeStream(STREAM_ALL);

  for (int i = 0; i < mvar; ++i)
      acDeviceLoadVertexBuffer(acGridGetDevice(), STREAM_DEFAULT, src, VertexBufferHandle(i));

  int n_aux_on_gpu = 0;
  for (int i = 0; i < mfarray; ++i)
  {
    if (maux_vtxbuf_index[i] != -1)
    {
        n_aux_on_gpu++;
      	acDeviceLoadVertexBuffer(acGridGetDevice(), STREAM_DEFAULT, src, VertexBufferHandle(maux_vtxbuf_index[i]));
    }
  }
  for (int i = 0; i < mglobal; ++i)
  {
      acDeviceLoadVertexBuffer(acGridGetDevice(), STREAM_DEFAULT, src, VertexBufferHandle(mvar+n_aux_on_gpu+i));
  }

  acGridSynchronizeStream(STREAM_ALL);
  if (nghost > 0 && (dimensionality < 2 || (dimensionality == 2 && nygrid == 1))) acHostMeshDestroy(&src);
}
/***********************************************************************************************/
void generate_bcs()
{
	if(rank != 0) return;
	bool bc2ast_exists = (system("ls src/scripts/bc2ast > /dev/null 2>&1") == 0);
	if(!bc2ast_exists)
	{
		fprintf(stderr,"AC Warning: Did not find src/scripts/bc2ast so skipping possible bc generation!\n");
		return;
	}
	if(system("cd src && scripts/bc2ast 1> ../tmp_bcs 2> /dev/null && cd .."))
	{
		fprintf(stderr,"AC Error: Was not able to generate bcs!\n");
		exit(EXIT_FAILURE);
	}
	const int different = system("diff tmp_bcs src/astaroth/DSL/local/boundconds.h");
	if(different)
	{
		if(system("mv tmp_bcs src/astaroth/DSL/local/boundconds.h"))
		{
			fprintf(stderr,"AC Error: Was not able move bcs to DSL/local!\n");
			exit(EXIT_FAILURE);
		}
		fprintf(stderr,"BCs different: recompiling!\n");
	}
	else
	{
		if(system("rm tmp_bcs"))
		{
			fprintf(stderr,"AC Error: Was not able to remove tmp_bcs!\n");
			exit(EXIT_FAILURE);
		}
	}
}
/***********************************************************************************************/
extern "C" void testBCs();     // forward declaration
/***********************************************************************************************/
extern "C" void initializeGPU(AcReal *farr, int comm_fint, double t, int nt_,
				int lread_all_vars_from_device_,
				int lcpu_timestep_on_gpu_,
				int lac_sparse_autotuning_)  // MPI_Fint comm_fint
{
  lac_sparse_autotuning = lac_sparse_autotuning_;
  if(lread_all_vars_from_device_) lread_all_vars_from_device = true;
  if(lcpu_timestep_on_gpu_) lcpu_timestep_on_gpu = true;
  //Setup configurations used for initializing and running the GPU code
  nt = nt_;
  comm_pencil = MPI_Comm_f2c(comm_fint);
  setupConfig(mesh.info);

  if (rank==0 && ldebug) printf("memusage after pointer assign= %f MBytes\n", acMemUsage()/1024.);
#if AC_RUNTIME_COMPILATION
#include "cmake_options.h"

  //Not worth it to get this working inside the container
  const bool inside_container = ltraining;
  if(!inside_container) generate_bcs();
  MPI_Barrier(MPI_COMM_WORLD);
  acCompile(cmake_options,mesh.info);
  acLoadLibrary(rank == 0 ? stderr : NULL,mesh.info);
  acCheckDeviceAvailability();
  acLogFromRootProc(rank, "Done setupConfig && acCompile\n");
#else
  acLogFromRootProc(rank, "Done setupConfig\n");
#endif
  fflush(stdout);
  //Done after setupConfig and acCompile since we need maux_vtxbuf_index and acGetNumFields
  //This is an ugly way to do this but works for now
  {
    const size_t z_offset  = (dimensionality == 2 && nzgrid == 1) ? NGHOST*mx*my : 0;
    for (int i = 0; i < mvar; ++i)
    {
      mesh.vertex_buffer[VertexBufferHandle(i)] = &farr[mw*i+ z_offset];
    }

    int n_aux_on_gpu = 0;
    for (int i = 0; i < mfarray; ++i)
    {
      if (maux_vtxbuf_index[i] != -1)
      {
	++n_aux_on_gpu;
        mesh.vertex_buffer[maux_vtxbuf_index[i]] = &farr[mw*i + z_offset];
      }
    }
    for (int i = 0; i < mfarray-mvar-maux; ++i)
    {
        mesh.vertex_buffer[mvar+n_aux_on_gpu+i] = &farr[mw*(mvar+maux+i) + z_offset];
    }
    //TP: for now for training we have all slots filled since we might want to read TAU components to the host for calculating validation error
    if (ltraining)
    {
    	for (int i = 0; i < acGetNumFields(); ++i)
    	{
	   if (mesh.vertex_buffer[i] == NULL)
	   {
    	    	mesh.vertex_buffer[i] = (AcReal*)malloc(sizeof(AcReal)*mw);
	   }
    	}
    }
  }
  if (rank==0 && ldebug) printf("memusage grid_init= %f MBytes\n", acMemUsage()/1024.);
  acGridInit(mesh);
  if(acGetDevice(&device_id) != cudaSuccess)
  {
	  fprintf(stderr,"Was not able to get device id!\n");
	  exit(EXIT_FAILURE);
  }
  if (rank==0 && ldebug) printf("memusage after grid_init= %f MBytes\n", acMemUsage()/1024.);

  mesh.info = acGridDecomposeMeshInfo(mesh.info);
  //Important that this is before autotuning
  acDeviceSetInput(acGridGetDevice(), AC_step_num,(PC_SUB_STEP_NUMBER)0);
  acDeviceSetInput(acGridGetDevice(), AC_dt,dt);
  acDeviceSetInput(acGridGetDevice(), AC_t,AcReal(t));
  acDeviceSetInput(acGridGetDevice(), AC_shear_delta_y,deltay);
  GW_timestep_graph = acGetOptimizedDSLTaskGraph(AC_gravitational_waves_solve_and_stress);
  if(luses_aa_pwd_top || luses_aa_pot2_top)
  {
	Field AA_fields[3];
	AA_fields[0] = acGetAAX();
	AA_fields[1] = acGetAAY();
	AA_fields[2] = acGetAAZ();
	boundary_z_halo_exchange_graph = acGridBuildTaskGraph({
			acHaloExchangeBoundary(AA_fields,3,BOUNDARY_Z_TOP)
			});
  }

		
  if (ltest_bcs) testBCs();
  //This is for autotuning
  afterSubStepGPU();
  autotune_all_integration_substeps();
  if (rank==0 && ldebug) printf("memusage before store config= %f MBytes\n", acMemUsage()/1024.);
  acStoreConfig(acDeviceGetLocalConfig(acGridGetDevice()), "PC-AC.conf");
  if (rank==0 && ldebug) printf("memusage after store config= %f MBytes\n", acMemUsage()/1024.);
  acGridSynchronizeStream(STREAM_ALL);
  if (rank==0 && ldebug) printf("memusage after store synchronize stream= %f MBytes\n", acMemUsage()/1024.);
  acLogFromRootProc(rank, "DONE initializeGPU\n");
  fflush(stdout);

  constexpr AcReal unit = 1.0;
  dt1_interface = unit/dt;
}
/***********************************************************************************************/
extern "C" void reloadConfig()
{
  setupConfig(mesh.info);
#if AC_RUNTIME_COMPILATION
  //save the current values of the vtxbufs since the device arrays are freed by acGridQuit
  copyFarray(mesh.vertex_buffer[0]);
  acGridQuit();
  const AcResult closed_res = acCloseLibrary();
  if (closed_res != AC_SUCCESS)
  {
	  if (rank == 0) 
	  {
		  fprintf(stderr,"Was not successful in closing Astaroth lib!\n");
		  int err = system("touch ERROR");
		  if(err)
		    fprintf(stderr,"Unable to create ERROR file!\n");
	  }
  	  MPI_Barrier(MPI_COMM_WORLD);
	  exit(EXIT_FAILURE);
  }
#include "cmake_options.h"
  generate_bcs();
  MPI_Barrier(MPI_COMM_WORLD);
  if(acCompile(cmake_options,mesh.info) != AC_SUCCESS)
  {
	  if(rank == 0) 
	  {
		  fprintf(stderr,"Was not able to compile Astaroth!\n");
		  int err = system("touch ERROR");
		  if(err)
		    fprintf(stderr,"Unable to create ERROR file!\n");
	  }
  	  MPI_Barrier(MPI_COMM_WORLD);
	  exit(EXIT_FAILURE);
  }
  acLoadLibrary(rank == 0 ? stderr : NULL,mesh.info);
  acGridInit(mesh);
  acLogFromRootProc(rank, "Done setupConfig && acCompile\n");
  fflush(stdout);
  fflush(stderr);
  //TP: this is important that we don't overwrite the output buffer in middle of a timestep when the output buffer holds some meaning!
  autotune_all_integration_substeps();
  //TP: restore the vtxbuf values before quitting grid
  loadFarray();
#else
  acGridSynchronizeStream(STREAM_ALL);
  acDeviceUpdate(acGridGetDevice(), mesh.info);
  acGridSynchronizeStream(STREAM_ALL);
#endif
  acStoreConfig(acDeviceGetLocalConfig(acGridGetDevice()), "PC-AC.conf");
  acLogFromRootProc(rank, "DONE reloading on GPU\n");
  fflush(stdout);
  fflush(stderr);
}
/***********************************************************************************************/
extern "C" void updateInConfigArr(int index)
{
     if (mesh.info[(AcRealArrayParam)index] != nullptr)
        acDeviceLoadRealArray(acGridGetDevice(),STREAM_DEFAULT,mesh.info,static_cast<AcRealArrayParam>(index));
}
/***********************************************************************************************/
extern "C" void updateInConfigScal(int index, AcReal value)
{
     acDeviceLoadScalarUniform(acGridGetDevice(),STREAM_DEFAULT,static_cast<AcRealParam>(index),value);
     acDeviceSynchronizeStream(acGridGetDevice(),STREAM_DEFAULT);
}
/***********************************************************************************************/
extern "C" void updateInConfigVec(int index, AcReal* value_)
{
     const AcReal3 value = (AcReal3){value_[0],value_[1],value_[2]};
     acDeviceLoadVectorUniform(acGridGetDevice(),STREAM_DEFAULT,static_cast<AcReal3Param>(index),value);
     acDeviceSynchronizeStream(acGridGetDevice(),STREAM_DEFAULT);
}
/***********************************************************************************************/
extern "C" int updateInConfigArrName(char *name)
{
    int index = -1;
    for (int i=0; i<NUM_REAL_ARRAYS; i++){
       if (strcmp(get_array_info(static_cast<AcRealArrayParam>(i)).name,name)==0) index=i;
    }
    if (index>-1) updateInConfigArr(index);
    else
    {
       fprintf(stderr,"Astaroth WARNING: Did not entry named %s in config!!\n",name);
       fflush(stderr);
    }
    return index;
}
/***********************************************************************************************/
extern "C" int updateInConfigScalName(char *name, AcReal value)
{
    int index = -1;
    for (int i=0; i<NUM_REAL_PARAMS; i++){
       if (strcmp(realparam_names[i],name)==0) index=i;
    }
    if (index>-1) updateInConfigScal(index, value);
    else
    {
       fprintf(stderr,"Astaroth WARNING: Did not entry named %s in config!!\n",name);
       fflush(stderr);
    }
    return index;
}
/**********************************************************************************************/
extern "C" int updateInConfigVecName(char *name, AcReal* value_)
{
    int index = -1;
    for (int i=0; i<NUM_REAL3_PARAMS; i++){
       if (strcmp(real3param_names[i],name)==0) index=i;
    }
    if (index>-1) updateInConfigVec(index, value_);
    else
    {
       fprintf(stderr,"Astaroth WARNING: Did not entry named %s in config!!\n",name);
       fflush(stderr);
    }
    return index;
}
/**********************************************************************************************/
extern "C" void finalizeGPU()
{
  // Deallocate everything on the GPUs and reset
  AcResult res = acGridQuit();
  if(GW_thread.joinable())
  {
          GW_thread.join();
  }
	
	// write the loss values to a file

  std::ofstream myFile;
  std::string fileString = "train_loss_" + std::to_string(my_rank)  + ".csv";	
  myFile.open(fileString);
  myFile << "epoch,train_loss\n";
  for (int i=0;i<train_loss.size();i++){
    myFile << i << "," << train_loss[i] << "\n";
  }
  myFile.close();

  std::string train_sample = "train_sample_" + std::to_string(my_rank) + ".csv";
  myFile.open(train_sample);
  myFile << "train_sample_nts\n";
  for (int i=0;i<train_nts.size();i++){
  	myFile << train_nts[i] << "\n";
  }
  myFile.close();

}
/***********************************************************************************************/
extern "C" void random_initial_condition()
{
  return;
  //acGridSynchronizeStream(STREAM_ALL);
  //AcMeshDims dims = acGetMeshDims(acGridGetLocalMeshInfo());
  //acGridLaunchKernel(STREAM_DEFAULT, randomize, dims.n0, dims.n1);
  //acGridSynchronizeStream(STREAM_ALL);
}
/***********************************************************************************************/
void testBCs()
{

  //if (dimensionality != 3) return;
  // Set random seed for reproducibility
  srand(321654987);
  const auto DEVICE_VTXBUF_IDX = [&](const int x, const int y, const int z)
  				{
					return acVertexBufferIdx(x,y,z,acGridGetLocalMeshInfo());
  				};
  AcTaskGraph* rhs = acGetOptimizedDSLTaskGraph(AC_rhs);
  const bool all_periodic = acGridTaskGraphHasPeriodicBoundcondsX(rhs) && 
  			    acGridTaskGraphHasPeriodicBoundcondsY(rhs) && 
  			    acGridTaskGraphHasPeriodicBoundcondsZ(rhs); 
  if (all_periodic && !lshear) return;
  if (lshear) acLogFromRootProc(rank,"testBCS: deltay: %7e\n",deltay);
  auto bcs = acGetDSLTaskGraph(boundconds);

  AcMesh tmp_mesh_to_store;
  acHostMeshCopy(mesh, &tmp_mesh_to_store);
  
  acHostMeshRandomize(&mesh);

  AcMesh mesh_to_copy;
  acHostMeshCopy(mesh, &mesh_to_copy);

  acGridSynchronizeStream(STREAM_ALL);
  acDeviceLoadMesh(acGridGetDevice(), STREAM_DEFAULT, mesh);
  acGridSynchronizeStream(STREAM_ALL);

  int ivar1 = 1;
  int ivar2 = mvar;
  boundconds_x_c(mesh.vertex_buffer[0],&ivar1,&ivar2);
  boundconds_y_c(mesh.vertex_buffer[0],&ivar1,&ivar2);
  boundconds_z_c(mesh.vertex_buffer[0],&ivar1,&ivar2);

  fourier_boundary_conditions();
  acGridSynchronizeStream(STREAM_ALL);
  acGridExecuteTaskGraph(bcs,1);
  acGridSynchronizeStream(STREAM_ALL);
  acGridSynchronizeStream(STREAM_ALL);
  acDeviceStoreMesh(acGridGetDevice(), STREAM_DEFAULT, &mesh_to_copy);
  acGridSynchronizeStream(STREAM_ALL);

  AcMeshDims dims = acGetMeshDims(acGridGetLocalMeshInfo());
  bool passed = true;
  AcReal max_abs_not_passed_val=-1.0;
  AcReal true_pair{};
  AcReal max_abs_relative_difference =-1.0;
  AcReal max_abs_value = -1.0;
  AcReal min_abs_value = 1.0;
  AcReal gpu_val_for_largest_diff{};
  AcReal true_val_for_largest_diff{};
  AcReal epsilon = (AcReal)pow(10.0,-12.0);
  int4 largest_diff_point{};

  auto skip_x = nxgrid == 1 || acGridTaskGraphHasPeriodicBoundcondsX(rhs);
  auto skip_y = nygrid == 1 || acGridTaskGraphHasPeriodicBoundcondsY(rhs);
  auto skip_z = nzgrid == 1 || acGridTaskGraphHasPeriodicBoundcondsZ(rhs);

  auto start_x =  skip_x ? NGHOST : 0;
  auto start_y =  skip_y ? NGHOST : 0;
  auto start_z =  skip_z ? NGHOST : 0;

  auto end_x =  skip_x ? dims.n1.x : dims.m1.x;
  auto end_y =  skip_y ? dims.n1.y : dims.m1.y;
  auto end_z =  skip_z ? dims.n1.z : dims.m1.z;

  int num_of_points_where_different[NUM_VTXBUF_HANDLES]{};
  bool different_in[NUM_VTXBUF_HANDLES][27]{};
  //TP: stupid levels of safety
  memset(different_in,0,NUM_VTXBUF_HANDLES*27*sizeof(bool));
  int num_of_points = 0;

  for (size_t i = start_x; i < end_x; i++)
  {
    for (size_t j = start_y; j < end_y; j++)
    {
      for (size_t k = start_z; k < end_z; k++)
      {
	const int x_idx = i < NGHOST ? -1 : i >= dims.n1.x ? 1 : 0;
	const int y_idx = j < NGHOST ? -1 : j >= dims.n1.y ? 1 : 0;
	const int z_idx = k < NGHOST ? -1 : k >= dims.n1.z ? 1 : 0;

	const bool x_in = i >= NGHOST && i < dims.n1.x;
	const bool y_in = j >= NGHOST && j < dims.n1.y;
	const bool z_in = k >= NGHOST && k < dims.n1.z;
	if (
		x_in && y_in && z_in
	  ) continue;
	  //if (i < NGHOST | i > dims.n1.x  || j < NGHOST || j > dims.n1.y) continue;
	++num_of_points;
        for (int ivar = 0; ivar < mvar; ivar++)
        {

	  const int index = (x_idx+1) + 3*((y_idx+1) + 3*(z_idx+1));
          AcReal out_val = mesh_to_copy.vertex_buffer[ivar][DEVICE_VTXBUF_IDX(i, j, k)];
          AcReal true_val = mesh.vertex_buffer[ivar][DEVICE_VTXBUF_IDX(i, j, k)];
          AcReal abs_diff = fabs(out_val - true_val);
          if (fabs(true_val) > max_abs_value) max_abs_value = fabs(true_val);
          if (fabs(true_val) < min_abs_value) min_abs_value = fabs(true_val);
          if (fabs(abs_diff/true_val) > epsilon || (true_val == (AcReal)0.0 && fabs(out_val) > (AcReal)pow(0.1,13)) || (epsilon == (AcReal)0.0 && true_val != out_val))
          {
	    different_in[ivar][index] |= true;
            passed = false;
            num_of_points_where_different[ivar]++;
            if (max_abs_not_passed_val<abs(out_val)){
              max_abs_not_passed_val = abs(out_val);
              true_pair = true_val;
            }
            if (true_val != (AcReal)0.0){
              if (max_abs_relative_difference<(abs_diff/true_val)){
                max_abs_relative_difference=(abs_diff/true_val);
                gpu_val_for_largest_diff = out_val;
                true_val_for_largest_diff = true_val;
		largest_diff_point = (int4){(int)i,(int)j,(int)k,(int)ivar};
              }
            }  
          }
          if (isnan(out_val))
          {
            acLogFromRootProc(rank,"nan before at %d,%d,%d,%s!\n!",i,j,k,acGetFieldName(Field(ivar)));
            acLogFromRootProc(rank,"%.7e\n",(double)out_val);
          }
        }
      }
    }
  }


  passed &= !has_nans(mesh);
  if (!passed)
  {
  	for (int ivar=0;ivar<mvar;ivar++)
        {
    		acLogFromRootProc(0,"ratio of values wrong for field: %s\t %f\n",acGetFieldName(Field(ivar)),(double)num_of_points_where_different[ivar]/num_of_points);
		for (int x = -1; x <= 1; ++x)
		{
			for (int y = -1; y <= 1; ++y)
			{
				for (int z = -1; z <= 1; ++z)
				{
	  				const int index = (x+1) + 3*((y+1) + 3*(z+1));
					if (different_in[ivar][index]) acLogFromRootProc(0,"different (%d,%d,%d)\n",x,y,z);
				}
			}
		}
	}
  	acLogFromRootProc(0,"max abs not passed val: %.7e\t%.7e\n",(double)max_abs_not_passed_val, (double)fabs(true_pair));
  	acLogFromRootProc(0,"max abs relative difference val: %.7e\n",(double)max_abs_relative_difference);
  	acLogFromRootProc(0,"Point where biggest rel diff: %d,%d,%d in Field: %s, diff: %.14e\n",largest_diff_point.x,largest_diff_point.y,largest_diff_point.z,acGetFieldName(Field(largest_diff_point.w)),max_abs_relative_difference);
  	acLogFromRootProc(0,"largest difference: %.7e\t%.7e\n",(double)gpu_val_for_largest_diff, (double)true_val_for_largest_diff);
  	acLogFromRootProc(0,"abs range: %.7e-%7e\n",(double)min_abs_value,(double)max_abs_value);
    	acLogFromRootProc(0,"Did not pass BC test :(\n");
	fprintf(stderr,"Did not pass BC\n");

    	exit(EXIT_FAILURE);
  }
  else
  {
	  fprintf(stderr,"Passed BC test :)\n");
  }
  fflush(stdout);
  fflush(stderr);

  acHostMeshDestroy(&mesh_to_copy);
  acHostMeshCopyVertexBuffers(tmp_mesh_to_store,mesh);
  acHostMeshDestroy(&tmp_mesh_to_store);
  acGridDestroyTaskGraph(bcs);
  acHostMeshDestroy(&tmp_mesh_to_store);
  exit(EXIT_SUCCESS);
}
/***********************************************************************************************/
extern "C" void prepareForFirstSubstep(double t)
{
	acGridSynchronizeStream(STREAM_ALL);
 	acDeviceSetInput(acGridGetDevice(), AC_t,AcReal(t));
	beforeBoundaryGPU(false,0,t,false);
	acGridExecuteTaskGraph(acGetOptimizedDSLTaskGraph(AC_initialize_sums),1);
	if (!lcourant_dt || !ldt)
	{
		return;
	}
	//Not needed but for extra safety
  	acDeviceSetInput(acGridGetDevice(), AC_step_num, (PC_SUB_STEP_NUMBER) 0);
	const auto graph = acGetOptimizedDSLTaskGraph(AC_calculate_timestep);

	acGridExecuteTaskGraph(graph,1);
	acGridSynchronizeStream(STREAM_ALL);
	AcReal dt1_ = calc_dt1_courant(AcReal(t));
	set_dt(dt1_);
	dt1_interface = dt1_;
        acDeviceSwapBuffers(acGridGetDevice());
	//Not strictly needed but for extra safety
	loadFarray();
}
/***********************************************************************************************/
//Reads the results of reductions and copies them to the host.
//At the moment the only use case are sums used in ODE integrations.
extern "C" void getGPUReducedVars(AcReal* dst)
{
        acDeviceSynchronizeStream(acGridGetDevice(),STREAM_DEFAULT);
#if LAXIONU1BACK
	dst[0] = acDeviceGetOutput(acGridGetDevice(), AC_edotb_sum__mod__axionu1back);
	dst[1] = acDeviceGetOutput(acGridGetDevice(), AC_rhoe__mod__axionu1back);
	dst[2] = acDeviceGetOutput(acGridGetDevice(), AC_rhob__mod__axionu1back);
#endif
#if LAXIONSU2BACK
	dst[0] = acDeviceGetOutput(acGridGetDevice(), AC_grand_sum);
	dst[1] = acDeviceGetOutput(acGridGetDevice(), AC_dgrant_sum);
	dst[2] = acDeviceGetOutput(acGridGetDevice(), AC_trdoteff2km_sum);
	dst[3] = acDeviceGetOutput(acGridGetDevice(), AC_trdoteff2m_sum);
	dst[4] = acDeviceGetOutput(acGridGetDevice(), AC_treff2km_sum);
	dst[5] = acDeviceGetOutput(acGridGetDevice(), AC_treff2m_sum);
	dst[6] = acDeviceGetOutput(acGridGetDevice(), AC_tldoteff2km_sum);
	dst[7] = acDeviceGetOutput(acGridGetDevice(), AC_tldoteff2m_sum);
	dst[8] = acDeviceGetOutput(acGridGetDevice(), AC_tleff2km_sum);
	dst[9] = acDeviceGetOutput(acGridGetDevice(), AC_tleff2m_sum);
#endif
#if LBACKREACT_INFL
	dst[0] = acDeviceGetOutput(acGridGetDevice(), AC_a2rhom_all__mod__backreact_infl);
	dst[1] = acDeviceGetOutput(acGridGetDevice(), AC_a2rhopm_all__mod__backreact_infl);
	dst[2] = acDeviceGetOutput(acGridGetDevice(), AC_a2rhophim_all__mod__backreact_infl);
	dst[3] = acDeviceGetOutput(acGridGetDevice(), AC_a2rhogphim_all__mod__backreact_infl);
	dst[4] = acDeviceGetOutput(acGridGetDevice(), AC_sige1m_all_nonaver__mod__backreact_infl);
	dst[5] = acDeviceGetOutput(acGridGetDevice(), AC_sigb1m_all_nonaver__mod__backreact_infl);
	dst[6] = acDeviceGetOutput(acGridGetDevice(), AC_ddotam_all__mod__backreact_infl);
	dst[7] = acDeviceGetOutput(acGridGetDevice(), AC_e2m_all__mod__backreact_infl);
	dst[8] = acDeviceGetOutput(acGridGetDevice(), AC_b2m_all__mod__backreact_infl);
#endif
#if LKLEIN_GORDON
	dst[0] = acDeviceGetOutput(acGridGetDevice(), AC_a2rhom_all__mod__klein_gordon);
	dst[1] = acDeviceGetOutput(acGridGetDevice(), AC_a2rhopm_all__mod__klein_gordon);
	dst[2] = acDeviceGetOutput(acGridGetDevice(), AC_a2rhophim_all__mod__klein_gordon);
	dst[3] = acDeviceGetOutput(acGridGetDevice(), AC_a2rhogphim_all__mod__klein_gordon);
	dst[4] = acDeviceGetOutput(acGridGetDevice(), AC_sige1m_all_nonaver__mod__klein_gordon);
	dst[5] = acDeviceGetOutput(acGridGetDevice(), AC_sigb1m_all_nonaver__mod__klein_gordon);
	dst[6] = acDeviceGetOutput(acGridGetDevice(), AC_ddotam_all__mod__klein_gordon);
	dst[7] = acDeviceGetOutput(acGridGetDevice(), AC_e2m_all__mod__klein_gordon);
	dst[8] = acDeviceGetOutput(acGridGetDevice(), AC_b2m_all__mod__klein_gordon);
#endif
}
/***********************************************************************************************/
