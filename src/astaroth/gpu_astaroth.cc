/*
   Date:   8-Feb-2017
   Author: M. Rheinhardt & j. Pekkilae
   Description:
           ANSI C and standard library callable function wrappers for ASTAROTH-nucleus functions to be called from Fortran.
  Comments: 
*/
// General headers.
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

#define CUDA_ERRCHK(X)

int counter = 0;
int done = 0;

// Astaroth headers.
#include "astaroth.h"
bool
has_nans(AcMesh mesh_in);

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
  #define nu nu__mod__viscosity
  #define nu_hyper2 nu_hyper2__mod__viscosity
  #define nu_hyper3 nu_hyper3__mod__viscosity
 
  #define gamma gamma__mod__equationofstate

  #define eta eta__mod__magnetic
  #define eta_hyper2 eta_hyper2__mod__magnetic
  #define eta_hyper3 eta_hyper3__mod__magnetic

  #define chi chi__mod__energy
  #define chi_hyper2 chi_hyper2__mod__energy
  #define chi_hyper3 chi_hyper3__mod__energy

  #define lmorton_curve lmorton_curve__mod__cdata
  #define ltest_bcs     ltest_bcs__mod__cdata
  #define num_substeps  num_substeps__mod__cdata
  #define maux_vtxbuf_index maux_vtxbuf_index__mod__cdata
  #define ldt ldt__mod__cdata
  #define lcourant_dt lcourant_dt__mod__cdata
  #define dt dt__mod__cdata
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

  #define lcpu_timestep_on_gpu   lcpu_timestep_on_gpu__mod__cdata
  #define lperi                  lperi__mod__cdata
  #define lxyz                   lxyz__mod__cdata
  #define lac_sparse_autotuning  lac_sparse_autotuning__mod__cdata
  #define ldebug ldebug__mod__cdata
  
  #define deltay  deltay__mod__cdata
  #define eps_rkf eps_rkf__mod__cdata
  
  #define itauxx itauxx__mod__training
  #define itauxy itauxy__mod__training
  #define itauxz itauxz__mod__training
  #define itauyy itauyy__mod__training
  #define itauyz itauyz__mod__training
  #define itauzz itauzz__mod__training

  #define lread_all_vars_from_device lread_all_vars_from_device__mod__cdata
  #define lsecond_force lsecond_force__mod__forcing
  #define lforce_helical lforce_helical__mod__forcing
#endif

AcReal dt1_interface{};
static int rank;
static AcMesh mesh = acInitMesh();
extern "C" void copyFarray(AcReal* f);

extern "C" void torch_trainCAPI(float* input, float* label, float* loss_val);
extern "C" void torch_inferCAPI(float* input, float* label);
extern "C" void torch_createmodel(const char* name, const char* config_fname, MPI_Comm mpi_comm, int device);


void print_debug() {
    #if TRAINING
    #include "user_constants.h"
	/*	
    auto temp = acGetOptimizedDSLTaskGraph(descale);
    acGridSynchronizeStream(STREAM_ALL);
    acGridExecuteTaskGraph(temp, 1);
    acGridSynchronizeStream(STREAM_ALL);
*/

		std::ofstream myFile;
		std::string fileString = "slices/slices_" + std::to_string(counter) + ".csv";	
		myFile.open(fileString);

    myFile << "it,TAU_xx,TAU_xx_inferred,TAU_yy,TAU_yy_inferred,TAU_zz,TAU_zz_inferred,"
           << "TAU_xy,TAU_xy_inferred,TAU_yz,TAU_yz_inferred,TAU_xz,TAU_xz_inferred,UUMEAN_x,UUMEAN_y,UUMEAN_z,"
           << "UUX,UUY,UUZ\n";


    copyFarray(NULL);
    const auto DEVICE_VTXBUF_IDX = [&](const int x, const int y, const int z) {
        return acVertexBufferIdx(x, y, z, mesh.info);
    };

    AcMeshDims dims = acGetMeshDims(acGridGetLocalMeshInfo());
	

    for (size_t i = dims.m0.x; i < dims.m1.x; i++) {
        for (size_t j = dims.m0.y; j < dims.m1.y; j++) {
            for (size_t k = dims.m0.z; k < dims.m1.z; k++) {

                myFile << counter << ","
                       << mesh.vertex_buffer[TAU.xx][DEVICE_VTXBUF_IDX(i, j, k)] << ","
                       << mesh.vertex_buffer[TAU_INFERRED.xx][DEVICE_VTXBUF_IDX(i, j, k)] << ","
                       << mesh.vertex_buffer[TAU.yy][DEVICE_VTXBUF_IDX(i, j, k)] << ","
                       << mesh.vertex_buffer[TAU_INFERRED.yy][DEVICE_VTXBUF_IDX(i, j, k)] << ","
                       << mesh.vertex_buffer[TAU.zz][DEVICE_VTXBUF_IDX(i, j, k)] << ","
                       << mesh.vertex_buffer[TAU_INFERRED.zz][DEVICE_VTXBUF_IDX(i, j, k)] << ","
                       << mesh.vertex_buffer[TAU.xy][DEVICE_VTXBUF_IDX(i, j, k)] << ","
                       << mesh.vertex_buffer[TAU_INFERRED.xy][DEVICE_VTXBUF_IDX(i, j, k)] << ","
                       << mesh.vertex_buffer[TAU.yz][DEVICE_VTXBUF_IDX(i, j, k)] << ","
                       << mesh.vertex_buffer[TAU_INFERRED.yz][DEVICE_VTXBUF_IDX(i, j, k)] << ","
                       << mesh.vertex_buffer[TAU.xz][DEVICE_VTXBUF_IDX(i, j, k)] << ","
                       << mesh.vertex_buffer[TAU_INFERRED.xz][DEVICE_VTXBUF_IDX(i, j, k)] << ","
                       << mesh.vertex_buffer[UUMEAN.x][DEVICE_VTXBUF_IDX(i, j, k)] << ","
                       << mesh.vertex_buffer[UUMEAN.y][DEVICE_VTXBUF_IDX(i, j, k)] << ","
                       << mesh.vertex_buffer[UUMEAN.z][DEVICE_VTXBUF_IDX(i, j, k)] << ","
                       << mesh.vertex_buffer[UUX][DEVICE_VTXBUF_IDX(i, j, k)] << ","
                       << mesh.vertex_buffer[UUY][DEVICE_VTXBUF_IDX(i, j, k)] << ","
                       << mesh.vertex_buffer[UUZ][DEVICE_VTXBUF_IDX(i, j, k)] << "\n";
            }
        }
    }
    #endif
}



extern "C" void torch_train_c_api(AcReal *loss_val){
	
#if TRAINING
	#include "user_constants.h"

	std::cout << "Counter is: " << counter << "\n"; 
	

	auto temp = acGetOptimizedDSLTaskGraph(train_prepare);
	acGridSynchronizeStream(STREAM_ALL);
	acGridExecuteTaskGraph(temp,1);
	acGridSynchronizeStream(STREAM_ALL);

	AcReal* out = NULL;
	
	AcReal* uumean_ptr = NULL;
	AcReal* TAU_ptr = NULL;
	*loss_val = 0.1;

	acDeviceGetVertexBufferPtrs(acGridGetDevice(), TAU.xx, &TAU_ptr, &out);
	acDeviceGetVertexBufferPtrs(acGridGetDevice(), UUMEAN.x, &uumean_ptr, &out);

	
	
  auto bcs = acGetOptimizedDSLTaskGraph(boundconds);	
	acGridSynchronizeStream(STREAM_ALL);
	acGridExecuteTaskGraph(bcs,1);
	acGridSynchronizeStream(STREAM_ALL);
	float avgloss = 0;

	for(int batch = 0; batch<5; batch++){
		torch_trainCAPI(uumean_ptr, TAU_ptr, loss_val);
		avgloss = avgloss + *loss_val;
	}
	printf("Loss after training: %f\n", avgloss/5);
	counter++;
#endif
}

float MSE(){
#if TRAINING
	#include "user_constants.h"
	
	auto temp = acGetOptimizedDSLTaskGraph(calc_validation_loss);
	acGridSynchronizeStream(STREAM_ALL);
	acGridExecuteTaskGraph(temp,1);
	acGridSynchronizeStream(STREAM_ALL);


  auto bcs = acGetOptimizedDSLTaskGraph(boundconds);	
	acGridSynchronizeStream(STREAM_ALL);
	acGridExecuteTaskGraph(bcs,1);
	acGridSynchronizeStream(STREAM_ALL);

	const auto DEVICE_VTXBUF_IDX = [&](const int x, const int y, const int z)
					{
						return acVertexBufferIdx(x, y, z, mesh.info);
					};
	AcMeshDims dims = acGetMeshDims(acGridGetLocalMeshInfo());


	copyFarray(NULL);

	return (acDeviceGetOutput(acGridGetDevice(), AC_l2_sum))/(6*32*32*32);
#else
#endif
}




extern "C" void torch_infer_c_api(int flag){	
#if TRAINING
	#include "user_constants.h"
	
		auto temp = acGetOptimizedDSLTaskGraph(train_prepare);
		acGridSynchronizeStream(STREAM_ALL);
		acGridExecuteTaskGraph(temp,1);
		acGridSynchronizeStream(STREAM_ALL);


		AcReal* out = NULL;
	
		AcReal* uumean_ptr = NULL;
		AcReal* tau_infer_ptr = NULL;

		acDeviceGetVertexBufferPtrs(acGridGetDevice(), TAU_INFERRED.xx, &tau_infer_ptr, &out);
		acDeviceGetVertexBufferPtrs(acGridGetDevice(), UUMEAN.x, &uumean_ptr, &out);

	
	
  	auto bcs = acGetOptimizedDSLTaskGraph(boundconds);	
		acGridSynchronizeStream(STREAM_ALL);
		acGridExecuteTaskGraph(bcs,1);
		acGridSynchronizeStream(STREAM_ALL);

		torch_inferCAPI(uumean_ptr, tau_infer_ptr);
		

		float vloss = MSE();
		
		
  	bcs = acGetOptimizedDSLTaskGraph(boundconds);	
		acGridSynchronizeStream(STREAM_ALL);
		acGridExecuteTaskGraph(bcs,1);
		acGridSynchronizeStream(STREAM_ALL);

		printf("Validation error is: %f\n", vloss);
		print_debug();
	#endif
}



/***********************************************************************************************/
int memusage()
  {
    struct rusage usage;
    int res=getrusage(RUSAGE_SELF,&usage);

    return usage.ru_maxrss;
  }
/***********************************************************************************************/
extern "C" void testRHS(AcReal *farray_in, AcReal *dfarray_truth)
{
  const auto DEVICE_VTXBUF_IDX = [&](const int x, const int y, const int z)
  				{
					return acVertexBufferIdx(x,y,z,mesh.info);
  				};
  // __energy_MOD_pushpars2c(p_par_energy);
  // mesh.info.profiles[AC_dlnhcond_prof]=dlnhcond_prof; // [2-1] dlnhcond_prof real(nz)
  // for (int i=20;i<30;++i)
  //   acLogFromRootProc(rank,"C: dlnhcond_prof %d=%f\n",i,mesh.info.profiles[AC_dlnhcond_prof][i]);
  // fflush(stdout);
  // return;
  // make_tasks(diagnostics_func, reduction_func, finalize_func, write_func);
  // thread_pool.WaitAll();
  constexpr AcReal alpha[3] = {0.0, (AcReal)-(5.0 / 9.0), (AcReal)-(153.0 / 128.0)};
  constexpr AcReal beta[3] = {(AcReal)(1.0 / 3.0), (AcReal)(15.0 / 16.0), (AcReal)(8.0 / 15.0)};
  constexpr int num_of_steps = 100;

  AcMesh mesh_true;
  AcMesh mesh_test;
  const AcReal epsilon = (AcReal)pow(0.1,12);
  // constexpr AcReal epsilon = 0.0;
  constexpr AcReal local_dt = (AcReal)0.001;
  Device dev = acGridGetDevice();
  acDeviceSetInput(acGridGetDevice(),AC_dt,local_dt);
  acGridSynchronizeStream(STREAM_ALL);

  size_t offset = 0;
  for (int i = 0; i < NUM_VTXBUF_HANDLES; ++i)
  {
    mesh_test.vertex_buffer[VertexBufferHandle(i)] = &farray_in[offset];
    offset += mw;
  }
  offset = 0;
	for (int i = 0; i < NUM_VTXBUF_HANDLES; ++i)
  {
    mesh_true.vertex_buffer[VertexBufferHandle(i)] = &dfarray_truth[offset];
    offset += mw;
  }
  AcMeshDims dims = acGetMeshDims(acGridGetLocalMeshInfo());
  acLogFromRootProc(rank,"n0: %d,%d,%d\trank: %d\n", dims.n0.x, dims.n0.y, dims.n0.z, rank);
  acLogFromRootProc(rank,"n1: %d,%d,%d\trank: %d\n", dims.n1.x, dims.n1.y, dims.n1.z, rank);

  //dryrun
  // acGridLaunchKernel(STREAM_DEFAULT, twopass_solve_intermediate_step0, dims.n0,dims.n1);
  // acGridSynchronizeStream(STREAM_ALL);
  // acGridSwapBuffers();
  // acGridSynchronizeStream(STREAM_ALL);
  // acGridLaunchKernel(STREAM_DEFAULT, twopass_solve_final_step0, dims.n0,dims.n1);
  // acGridSynchronizeStream(STREAM_ALL);
  // acGridSwapBuffers();
  // acGridSynchronizeStream(STREAM_ALL);
  // acGridLaunchKernel(STREAM_DEFAULT, twopass_solve_intermediate_step1, dims.n0,dims.n1);
  // acGridSynchronizeStream(STREAM_ALL);
  // acGridSwapBuffers();
  // acGridSynchronizeStream(STREAM_ALL);

  // acGridExecuteTaskGraph(rhs_test_graph,1);
  AcTaskGraph *rhs =  acGetOptimizedDSLTaskGraph(AC_rhs);
  for (int i = 0; i < 3; ++i)
  {
  	acDeviceSetInput(acGridGetDevice(), AC_step_num, (PC_SUB_STEP_NUMBER) i);
  	acGridExecuteTaskGraph(rhs,1);
  	acGridSynchronizeStream(STREAM_ALL);
  }

  acGridSynchronizeStream(STREAM_ALL);
  acDeviceLoadMesh(acGridGetDevice(), STREAM_DEFAULT, mesh_test);
  acGridSynchronizeStream(STREAM_ALL);
  acDeviceStoreMesh(acGridGetDevice(), STREAM_DEFAULT, &mesh);
  acGridSynchronizeStream(STREAM_ALL);

  bool loaded_correct = true;
  for (size_t i = 0; i < mx; i++)
  {
    for (size_t j = 0; j < my; j++)
    {
      for (size_t k = 0; k < mz; k++)
      {
        for (size_t ivar = 0; ivar < mfarray; ivar++)
        {
          AcReal out_val = mesh.vertex_buffer[ivar][DEVICE_VTXBUF_IDX(i, j, k)];
          AcReal true_val = mesh_test.vertex_buffer[ivar][DEVICE_VTXBUF_IDX(i, j, k)];
          if (out_val != true_val)
          {
            loaded_correct = false;
            acLogFromRootProc(rank,"Loaded val wrong at %d,%d,%d\n", i, j, k);
            acLogFromRootProc(rank,"field = %d", ivar);
            acLogFromRootProc(rank,"Loaded val: %f\tTRUE val: %f\n", (double)out_val, (double)true_val);
          }
        }
      }
    }
  }
  //AcReal derx_ux;
  //AcReal derx_normal;
  //AcReal dery_uy;
  //AcReal derz_uz;
  //test calculating divu different ways
  //const int y = 42;
  //const int x = 42;
  //const int z = 42;
  //const AcReal AC_inv_dsx = 2.0;

    // [0][0][-3] = -AC_inv_dsx * DER1_3,
    // [0][0][-2] = -AC_inv_dsx * DER1_2,
    // [0][0][-1] = -AC_inv_dsx * DER1_1,
    // [0][0][1]  = AC_inv_dsx * DER1_1,
    // [0][0][2]  = AC_inv_dsx * DER1_2,
    // [0][0][3]  = AC_inv_dsx * DER1_3

  // derx_normal = -AC_inv_dsx*DER1_3*(mesh.vertex_buffer[0][DEVICE_VTXBUF_IDX(x-3,y,z)]);
  // derx_normal += -AC_inv_dsx*DER1_2*(mesh.vertex_buffer[0][DEVICE_VTXBUF_IDX(x-2,y,z)]);
  // derx_normal += -AC_inv_dsx*DER1_1*(mesh.vertex_buffer[0][DEVICE_VTXBUF_IDX(x-1,y,z)]);
  // derx_normal += AC_inv_dsx*DER1_1*(mesh.vertex_buffer[0][DEVICE_VTXBUF_IDX(x+1,y,z)]);
  // derx_normal += AC_inv_dsx*DER1_2*(mesh.vertex_buffer[0][DEVICE_VTXBUF_IDX(x+2,y,z)]);
  // derx_normal += AC_inv_dsx*DER1_3*(mesh.vertex_buffer[0][DEVICE_VTXBUF_IDX(x+3,y,z)]);
  // acLogFromRootProc(rank,"test derx_ux: %.7e\n",derx_ux);
  // acLogFromRootProc(rank,"normal derx_ux. %.7e\n",derx_normal);
  fflush(stdout);
  // return;
  if (loaded_correct)
  {
    acLogFromRootProc(rank,"loaded correct data\n");
  }
  else
  {
    acLogFromRootProc(rank,"loaded incorrect data :(\n");
  }

  //set output buffer to 0 since if we are reading from it we don't want NaNs
  // acGridExecuteTaskGraph(rhs_test_graph, 1);
  // acGridSynchronizeStream(STREAM_ALL);
  // acGridExecuteTaskGraph(rhs_test_rhs_1, 1);
  // acGridSynchronizeStream(STREAM_ALL);
  acGridLaunchKernel(STREAM_DEFAULT, AC_BUILTIN_RESET, dims.n0, dims.n1);
  acGridSynchronizeStream(STREAM_ALL);

  //actual run
  //TP: works for only rk3 but is anyone anyways using this? Probably not
  for (int i=0;i<num_of_steps;i++){
    for (int substep = 0; substep <3; ++substep)
    {
    	acDeviceSetInput(acGridGetDevice(), AC_step_num, (PC_SUB_STEP_NUMBER)substep);
    	acGridExecuteTaskGraph(acGetOptimizedDSLTaskGraph(AC_rhs),1);
    	acGridSynchronizeStream(STREAM_ALL);
    }
    acGridSynchronizeStream(STREAM_ALL);
    acDeviceStoreMesh(acGridGetDevice(), STREAM_DEFAULT, &mesh);
  //  acGridSynchronizeStream(STREAM_ALL);
  //  AcReal max_uux2 = 0.0;
  //  AcReal max_uuy2 = 0.0;
  //  AcReal max_uuz2 = 0.0;
  //for (int i = dims.n0.x; i < dims.n1.x; i++)
  //{
  //  for (int j = dims.n0.y; j < dims.n1.y; j++)
  //  {
  //    for (int k = dims.n0.z; k < dims.n1.z; k++)
  //    {
  //        max_uux2 = max(max_uux2,pow(mesh.vertex_buffer[0][DEVICE_VTXBUF_IDX(i,j,k)],2));
  //        max_uuy2 = max(max_uuy2,pow(mesh.vertex_buffer[1][DEVICE_VTXBUF_IDX(i,j,k)],2));
  //        max_uuz2 = max(max_uuz2,pow(mesh.vertex_buffer[2][DEVICE_VTXBUF_IDX(i,j,k)],2));
  //    }
  //  }
  //}
  //acLogFromRootProc(rank,"GPU: uumax: %.7e\n", pow(max_uux2+max_uuy2+max_uuz2,0.5));
  //  acGridSynchronizeStream(STREAM_ALL);
  }
    acGridSynchronizeStream(STREAM_ALL);
  // acGridLaunchKernel(STREAM_DEFAULT, twopass_solve_intermediate_step0, dims.n0,dims.n1);
  // acGridSynchronizeStream(STREAM_ALL);
  // acGridSwapBuffers();
  // acGridSynchronizeStream(STREAM_ALL);
  // acGridLaunchKernel(STREAM_DEFAULT, twopass_solve_final_step0, dims.n0,dims.n1);
  // acGridSynchronizeStream(STREAM_ALL);
  // acGridSwapBuffers();
  // acGridSynchronizeStream(STREAM_ALL);
  // acGridLaunchKernel(STREAM_DEFAULT, twopass_solve_intermediate_step1, dims.n0,dims.n1);
  // acGridSynchronizeStream(STREAM_ALL);
  // acGridSwapBuffers();
  // acGridSynchronizeStream(STREAM_ALL);

  // acGridLaunchKernel(STREAM_DEFAULT, twopass_solve_final_step1, dims.n0,dims.n1);
  // acGridSynchronizeStream(STREAM_ALL);
  // acGridSwapBuffers();
  // acGridSynchronizeStream(STREAM_ALL);

  // acGridLaunchKernel(STREAM_DEFAULT, twopass_solve_intermediate_step2, dims.n0,dims.n1);
  // acGridSynchronizeStream(STREAM_ALL);
  // acGridSwapBuffers();
  // acGridSynchronizeStream(STREAM_ALL);

  // acGridLaunchKernel(STREAM_DEFAULT, twopass_solve_final_step2, dims.n0,dims.n1);
  // acGridSynchronizeStream(STREAM_ALL);
  // acGridSwapBuffers();
  // acGridSynchronizeStream(STREAM_ALL);

  acGridSynchronizeStream(STREAM_ALL);
  acDeviceStoreMesh(acGridGetDevice(), STREAM_DEFAULT, &mesh);

  bool passed = true;
  AcReal max_abs_not_passed_val=-1.0;
  AcReal true_pair{};
  AcReal max_abs_relative_difference =-1.0;
  AcReal max_abs_value = -1.0;
  AcReal min_abs_value = 1.0;
  AcReal gpu_val_for_largest_diff{};
  AcReal true_val_for_largest_diff{};
  int num_of_points_where_different[NUM_VTXBUF_HANDLES] = {0};

  for (size_t i = dims.n0.x; i < dims.n1.x; i++)
  {
    for (size_t j = dims.n0.y; j < dims.n1.y; j++)
    {
      for (size_t k = dims.n0.z; k < dims.n1.z; k++)
      {
        for (size_t ivar = 0; ivar < NUM_VTXBUF_HANDLES; ivar++)
        {
          AcReal out_val = mesh.vertex_buffer[ivar][DEVICE_VTXBUF_IDX(i, j, k)];
          AcReal true_val = mesh_true.vertex_buffer[ivar][DEVICE_VTXBUF_IDX(i, j, k)];
          AcReal abs_diff = fabs(out_val - true_val);
          if (fabs(true_val) > max_abs_value) max_abs_value = fabs(true_val);
          if (fabs(true_val) < min_abs_value) min_abs_value = fabs(true_val);
          if ((AcReal)(abs_diff/true_val) > epsilon || (true_val == (AcReal)0.0 && (AcReal)fabs(out_val) > (AcReal)pow(0.1,13)) || (epsilon == (AcReal)0.0 && true_val != out_val))
          {
            passed = false;
            num_of_points_where_different[ivar]++;
            // acLogFromRootProc(rank,"rhs val wrong at %d,%d,%d\n", i, j, k);
            // acLogFromRootProc(rank,"field = %d", ivar);
            // acLogFromRootProc(rank,"GPU val: %.7e\tTRUE val: %.7e\n", out_val, true_val);
            // acLogFromRootProc(rank,"PID: %d\n", rank);
            if (max_abs_not_passed_val<abs(out_val)){
              max_abs_not_passed_val = abs(out_val);
              true_pair = true_val;
            }
            if (true_val != (AcReal)0.0){
              if (max_abs_relative_difference<(abs_diff/true_val)){
                max_abs_relative_difference=(abs_diff/true_val);
                gpu_val_for_largest_diff = out_val;
                true_val_for_largest_diff = true_val;
              }
            }  
          }
          if (isnan(out_val))
          {
            acLogFromRootProc(rank,"nan before at %d,%d,%d,%d!\n!",i,j,k,ivar);
            acLogFromRootProc(rank,"%.7e\n",(double)out_val);
          }
        }
      }
    }
  }
  auto volume_size = [](const auto& a)
  {
	  return a.x*a.y*a.z;
  };
  for (int ivar=0;ivar<NUM_VTXBUF_HANDLES;ivar++)
    acLogFromRootProc(rank,"ratio of values wrong for field: %d\t %f\n",ivar,(double)num_of_points_where_different[ivar]/volume_size(dims.n1-dims.n0));
  passed &= !has_nans(mesh);
  if (passed)
  {
    acLogFromRootProc(rank,"Passed GPU test :)\n");
  }
  else
  {
    acLogFromRootProc(rank,"Did not pass GPU test :(\n");
  }
  acLogFromRootProc(rank,"max abs not passed val: %.7e\t%.7e\n",(double)max_abs_not_passed_val, (double)fabs(true_pair));
  acLogFromRootProc(rank,"max abs relative difference val: %.7e\n",(double)max_abs_relative_difference);
  acLogFromRootProc(rank,"largest difference: %.7e\t%.7e\n",(double)gpu_val_for_largest_diff, (double)true_val_for_largest_diff);
  acLogFromRootProc(rank,"abs range: %.7e-%7e\n",(double)min_abs_value,(double)max_abs_value);
  fflush(stdout);
}

AcReal
to_real(void* param, const char* name)
{
	if(param == NULL)
	{
		fprintf(stderr,"Passed NULL to pushparsc: %s!!\n",name);
		abort();
	}
	return *((AcReal*)param);
}
int
to_int(void* param, const char* name)
{
	if(param == NULL)
	{
		fprintf(stderr,"Passed NULL to pushparsc: %s!!\n",name);
		abort();
	}
	return *((int*)param);
}
bool
to_bool(void* param, const char* name)
{
	if(param == NULL)
	{
		fprintf(stderr,"Passed NULL to pushparsc: %s!!\n",name);
		abort();
	}
	return *((bool*)param);
}
int3
to_int3(void* param, const char* name)
{
	if(param == NULL)
	{
		fprintf(stderr,"Passed NULL to pushparsc: %s!!\n",name);
		abort();
	}
        int* arr = (int*)param;
        return (int3){arr[0],arr[1],arr[2]};
}
AcReal3
to_real3(void* param, const char* name)
{
	if(param == NULL)
	{
		fprintf(stderr,"Passed NULL to pushparsc: %s!!\n",name);
		abort();
	}
        AcReal* arr = (AcReal*)param;
        return (AcReal3){arr[0],arr[1],arr[2]};
}
AcBool3
to_bool3(void* param, const char* name)
{
	if(param == NULL)
	{
		fprintf(stderr,"Passed NULL to pushparsc: %s!!\n",name);
		abort();
	}
        bool* arr = (bool*)param;
        return (AcBool3){arr[0],arr[1],arr[2]};
}


__thread int tp_int;
typedef void (*rangefunc)(const int a, const int b);

AcReal cpu_pow(AcReal const val, AcReal exponent)
{
// masks hip GPU power function.
  return std::pow(val, exponent);
}
// PC interface headers.
#include "PC_moduleflags.h"
//#include "../cdata_c.h"
#include "../sub_c.h"           // provides set_dt
#include "../boundcond_c.h"     // provides boundconds[xyz] etc.
//#include "../mpicomm_c.h"       // provides finalize_sendrcv_bdry
#include "PC_module_parfuncs.h" // provides stuff from physics modules
				//
				//
//TP: x,y and z macros are too general

#if PACKED_DATA_TRANSFERS
  #include "loadStore.h"
#endif
#if LFORCING
  #include "../forcing_c.h"     // provides forcing_pars_hel
  #include "forcing.h"
#endif
// Astaroth objects instantiation.
//static AcMesh test_mesh;
static AcTaskGraph *randomize_graph;
static AcTaskGraph *rhs_test_graph;
static AcTaskGraph *rhs_test_rhs_1;

// Other.
int halo_xz_size[2] = {0, 0}, halo_yz_size[2] = {0, 0};
//static AcReal *xtop_buffer, *xbot_buffer, *ytop_buffer, *ybot_buffer;

/***********************************************************************************************/
int DCONST(const AcIntParam param)
{
  return mesh.info[param];
}
/***********************************************************************************************/
int3 DCONST(const AcInt3Param param)
{
  return mesh.info[param];
}
/***********************************************************************************************/
AcReal DCONST(const AcRealParam param)
{
  return mesh.info[param];
}
/***********************************************************************************************/
AcReal3 DCONST(const AcReal3Param param)
{
  return mesh.info[param];
}
/***********************************************************************************************/
//TP: for testing not usually used
/***
static const AcReal DER1_3_E = 1.0;
static const AcReal DER1_2_E = -9.0;
static const AcReal DER1_1_E = 45.0;
static const AcReal DER1_E_DIV = 60.0;
static const AcReal AC_inv_dsx = 20.3718327157626;
static const AcReal AC_inv_dsy = 20.3718327157626;
static const AcReal AC_inv_dsz = 20.3718327157626;
static const bool symmetric_der = true;
AcReal
derxx(const int x, const int y, const int z, AcMesh mesh, const int field)
{
    AcReal inv = AC_inv_dsx*AC_inv_dsx;
    AcReal derxx = inv*DER2_0*mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x,y,z)];
    derxx += inv*DER2_1*(mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x-1,y,z)]+mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x+1,y,z)]);
    derxx += inv*DER2_2*(mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x-2,y,z)]+mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x+2,y,z)]);
    derxx += inv*DER2_3*(mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x-3,y,z)]+mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x+3,y,z)]);
    return derxx;
}
AcReal
derx_astaroth(const int x, const int y, const int z, AcMesh mesh, const int field)
{
  AcReal derx = (-AC_inv_dsx*DER1_3)*(mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x-3,y,z)]);
  derx += (-AC_inv_dsx*DER1_2)*(mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x-2,y,z)]);
  derx += (-AC_inv_dsx*DER1_1)*(mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x-1,y,z)]);

  derx += (AC_inv_dsx*DER1_1)*(mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x+1,y,z)]);
  derx += (AC_inv_dsx*DER1_2)*(mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x+2,y,z)]);
  derx += (AC_inv_dsx*DER1_3)*(mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x+3,y,z)]);

  return derx;
}
AcReal
derx_symmetric(const int x, const int y, const int z, AcMesh mesh, const int field)
{
  AcReal derx = AC_inv_dsx*DER1_3*(mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x+3,y,z)]-mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x-3,y,z)]);
  derx += AC_inv_dsx*DER1_2*(mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x+2,y,z)]-mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x-2,y,z)]);
  derx += AC_inv_dsx*DER1_1*(mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x+1,y,z)]-mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x-1,y,z)]);
  return derx;
}
AcReal
derx_exact(const int x, const int y, const int z, AcMesh mesh, const int field)
{
  const AcReal divider = AC_inv_dsx*(1/DER1_E_DIV);
  AcReal derx = DER1_1_E*(mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x+1,y,z)]-mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x-1,y,z)]);
  derx += DER1_2_E*(mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x+2,y,z)]-mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x-2,y,z)]);
  derx += DER1_3_E*(mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x+3,y,z)]-mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x-3,y,z)]);
  return divider*derx;
}
AcReal
derx(const int x, const int y, const int z, AcMesh mesh, const int field)
{
  if (symmetric_der) return derx_symmetric(x,y,z,mesh,field);
  return derx_astaroth(x,y,z,mesh,field);
}
AcReal
dery_exact(const int x, const int y, const int z, AcMesh mesh, const int field)
{
  const AcReal divider = AC_inv_dsy*(1/DER1_E_DIV);
  AcReal dery = DER1_1_E*(mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x,y+1,z)]-mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x,y-1,z)]);
  dery += DER1_2_E*(mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x,y+2,z)]-mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x,y-2,z)]);
  dery += DER1_3_E*(mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x,y+3,z)]-mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x,y-3,z)]);
  return divider*dery;
}
AcReal
dery_symmetric(const int x, const int y, const int z, AcMesh mesh, const int field)
{
  AcReal dery = AC_inv_dsy*DER1_3*(mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x,y+3,z)]-mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x,y-3,z)]);
  dery += AC_inv_dsy*DER1_2*(mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x,y+2,z)]-mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x,y-2,z)]);
  dery += AC_inv_dsy*DER1_1*(mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x,y+1,z)]-mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x,y-1,z)]);
  return dery;
}
AcReal
dery_astaroth(const int x, const int y, const int z, AcMesh mesh, const int field)
{
  AcReal dery = (-AC_inv_dsy*DER1_3)*(mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x,y-3,z)]);
  dery += (-AC_inv_dsy*DER1_2)*(mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x,y-2,z)]);
  dery += (-AC_inv_dsy*DER1_1)*(mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x,y-1,z)]);

  dery += (AC_inv_dsy*DER1_1)*(mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x,y+1,z)]);
  dery += (AC_inv_dsy*DER1_2)*(mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x,y+2,z)]);
  dery += (AC_inv_dsy*DER1_3)*(mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x,y+3,z)]);

  return dery;
}
AcReal
dery(const int x, const int y, const int z, AcMesh mesh, const int field)
{
  if (symmetric_der) return dery_symmetric(x,y,z,mesh,field);
  return dery_astaroth(x,y,z,mesh,field);
}
AcReal
derz_exact(const int x, const int y, const int z, AcMesh mesh, const int field)
{
  const AcReal divider = AC_inv_dsz*(1/DER1_E_DIV);

  AcReal derz = DER1_1_E*(mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x,y,z+1)]-mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x,y,z-1)]);
  derz += DER1_2_E*(mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x,y,z+2)]-mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x,y,z-2)]);
  derz += DER1_3_E*(mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x,y,z+3)]-mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x,y,z-3)]);
  return divider*derz;
}
AcReal
derz_symmetric(const int x, const int y, const int z, AcMesh mesh, const int field)
{
  AcReal derz = AC_inv_dsz*DER1_3*(mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x,y,z+3)]-mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x,y,z-3)]);
  derz += AC_inv_dsz*DER1_2*(mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x,y,z+2)]-mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x,y,z-2)]);
  derz += AC_inv_dsz*DER1_1*(mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x,y,z+1)]-mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x,y,z-1)]);
  return derz;
}
AcReal
derz_astaroth(const int x, const int y, const int z, AcMesh mesh, const int field)
{
  AcReal derz = (-AC_inv_dsz*DER1_3)*(mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x,y,z-3)]);
  derz += (-AC_inv_dsz*DER1_2)*(mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x,y,z-2)]);
  derz += (-AC_inv_dsz*DER1_1)*(mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x,y,z-1)]);

  derz += (AC_inv_dsz*DER1_1)*(mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x,y,z+1)]);
  derz += (AC_inv_dsz*DER1_2)*(mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x,y,z+2)]);
  derz += (AC_inv_dsz*DER1_3)*(mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x,y,z+3)]);

  return derz;
}
AcReal
derz(const int x, const int y, const int z, AcMesh mesh, const int field)
{
  if (symmetric_der) return derz_symmetric(x,y,z,mesh,field);
  return derz_astaroth(x,y,z,mesh,field);
}
AcReal
divergence_symmetric(const int x, const int y, const int z, AcMesh mesh, const int field)
{
  return derx_symmetric(x,y,z,mesh,field) + dery_symmetric(x,y,z,mesh,field+1) + derz_symmetric(x,y,z,mesh,field+2);
}
AcReal
divergence_exact(const int x, const int y, const int z, AcMesh mesh, const int field)
{
  return derx_exact(x,y,z,mesh,field) + dery_exact(x,y,z,mesh,field+1) + derz_exact(x,y,z,mesh,field+2);
}
AcReal
divergence_astaroth(const int x, const int y, const int z, AcMesh mesh, const int field)
{
  return derx_astaroth(x,y,z,mesh,field) + dery_astaroth(x,y,z,mesh,field+1) + derz_astaroth(x,y,z,mesh,field+2);
}
AcReal
divergence(const int x, const int y, const int z, AcMesh mesh, const int field)
{
  if (symmetric_der) return divergence_symmetric(x,y,z,mesh,field);
  return divergence_astaroth(x,y,z,mesh,field);
}
AcReal3
gradient_exact(const int x, const int y, const int z, AcMesh mesh, const int field)
{
  return {derx_exact(x,y,z,mesh,field), dery_exact(x,y,z,mesh,field), derz_exact(x,y,z,mesh,field)};
}
AcReal3
gradient_symmetric(const int x, const int y, const int z, AcMesh mesh, const int field)
{
  return {derx_symmetric(x,y,z,mesh,field), dery_symmetric(x,y,z,mesh,field), derz_symmetric(x,y,z,mesh,field)};
}
AcReal3
gradient_astaroth(const int x, const int y, const int z, AcMesh mesh, const int field)
{
  return {derx_astaroth(x,y,z,mesh,field), dery_astaroth(x,y,z,mesh,field), derz_astaroth(x,y,z,mesh,field)};
}
AcReal3
gradient(const int x, const int y, const int z, AcMesh mesh, const int field)
{
  if (symmetric_der) return gradient_symmetric(x,y,z,mesh,field);
  return gradient_astaroth(x,y,z,mesh,field);
}
AcMatrix
gradients_symmetric(const int x, const int y, const int z, AcMesh mesh, const int field)
{
  return AcMatrix(gradient_symmetric(x,y,z,mesh,field),gradient_symmetric(x,y,z,mesh,field+1),gradient_symmetric(x,y,z,mesh,field+2));
}
AcMatrix
gradients_astaroth(const int x, const int y, const int z, AcMesh mesh, const int field)
{
  return AcMatrix(gradient_astaroth(x,y,z,mesh,field),gradient_astaroth(x,y,z,mesh,field+1),gradient_astaroth(x,y,z,mesh,field+2));
}
AcReal3
vecvalue(const int x, const int y, const int z, AcMesh mesh, const int field)
{
  return{mesh.vertex_buffer[field][DEVICE_VTXBUF_IDX(x,y,z)],mesh.vertex_buffer[field+1][DEVICE_VTXBUF_IDX(x,y,z)],mesh.vertex_buffer[field+2][DEVICE_VTXBUF_IDX(x,y,z)]};
}
AcMatrix
gradients(const int x, const int y, const int z, AcMesh mesh, const int field)
{
  if (symmetric_der) return gradients_symmetric(x,y,z,mesh,field);
  return gradients_astaroth(x,y,z,mesh,field);
}
***/
/***********************************************************************************************/
/****
void print_diagnostics(const int pid, const int step, const AcReal dt_, const AcReal simulation_time,
                       FILE *diag_file, const AcReal sink_mass, const AcReal accreted_mass,
                       int *found_nan)
{
  AcReal buf_rms, buf_max, buf_min;
  const int max_name_width = 16;

#if LHYDRO
  // Calculate rms, min and max from the velocity vector field
  acGridReduceVec(STREAM_DEFAULT, RTYPE_MAX, UUX, UUY, UUZ, &buf_max);
  acGridReduceVec(STREAM_DEFAULT, RTYPE_MIN, UUX, UUY, UUZ, &buf_min);
  acGridReduceVec(STREAM_DEFAULT, RTYPE_RMS, UUX, UUY, UUZ, &buf_rms);
#endif

  if (pid == 0)
  {
    //(pid, "Step %d, t_step %.3e, dt %e s\n", step, double(simulation_time), double(dt_));
    //(pid, "  %*s: min %.3e,\trms %.3e,\tmax %.3e\n", max_name_width, "uu total",
    // double(buf_min), double(buf_rms), double(buf_max));
    fprintf(diag_file, "%d %e %e %e %e %e ", step, double(simulation_time), double(dt_),
            double(buf_min), double(buf_rms), double(buf_max));
  }

#if LBFIELD
  acGridReduceVec(STREAM_DEFAULT, RTYPE_MAX, BFIELDX, BFIELDY, BFIELDZ, &buf_max);
  acGridReduceVec(STREAM_DEFAULT, RTYPE_MIN, BFIELDX, BFIELDY, BFIELDZ, &buf_min);
  acGridReduceVec(STREAM_DEFAULT, RTYPE_RMS, BFIELDX, BFIELDY, BFIELDZ, &buf_rms);

  acLogFromRootProc(pid, "  %*s: min %.3e,\trms %.3e,\tmax %.3e\n", max_name_width, "bb total",
                    double(buf_min), double(buf_rms), double(buf_max));
  if (pid == 0) fprintf(diag_file, "%e %e %e ", double(buf_min), double(buf_rms), double(buf_max));

  acGridReduceVecScal(STREAM_DEFAULT, RTYPE_ALFVEN_MAX, BFIELDX, BFIELDY, BFIELDZ, RHO, &buf_max);
  acGridReduceVecScal(STREAM_DEFAULT, RTYPE_ALFVEN_MIN, BFIELDX, BFIELDY, BFIELDZ, RHO, &buf_min);
  acGridReduceVecScal(STREAM_DEFAULT, RTYPE_ALFVEN_RMS, BFIELDX, BFIELDY, BFIELDZ, RHO, &buf_rms);

  acLogFromRootProc(pid, "  %*s: min %.3e,\trms %.3e,\tmax %.3e\n", max_name_width, "vA total",
                    double(buf_min), double(buf_rms), double(buf_max));
  if (pid == 0) fprintf(diag_file, "%e %e %e ", double(buf_min), double(buf_rms), double(buf_max));
#endif

  // Calculate rms, min and max from the variables as scalars

  for (int i = 0; i < NUM_VTXBUF_HANDLES; ++i)
  {
    acGridReduceScal(STREAM_DEFAULT, RTYPE_MAX, VertexBufferHandle(i), &buf_max);
    acGridReduceScal(STREAM_DEFAULT, RTYPE_MIN, VertexBufferHandle(i), &buf_min);
    acGridReduceScal(STREAM_DEFAULT, RTYPE_RMS, VertexBufferHandle(i), &buf_rms);

    acLogFromRootProc(pid, "  %*s: min %.3e,\trms %.3e,\tmax %.3e\n", max_name_width,
                      vtxbuf_names[i], double(buf_min), double(buf_rms), double(buf_max));
    if (pid == 0)
    {
      fprintf(diag_file, "%e %e %e ", double(buf_min), double(buf_rms), double(buf_max));
    }

    if (isnan(buf_max) || isnan(buf_min) || isnan(buf_rms))
    {
      *found_nan = 1;
    }
  }

  if ((sink_mass >= AcReal(0.0)) || (accreted_mass >= AcReal(0.0)))
  {
    if (pid == 0)
    {
      fprintf(diag_file, "%e %e ", double(sink_mass), double(accreted_mass));
    }
  }

  if (pid == 0)
  {
    fprintf(diag_file, "\n");
  }

  fflush(diag_file);
  fflush(stdout);
}
***/
/***********************************************************************************************/
std::array<AcReal,3>
visc_get_max_diffus()
{
	return
#if LVISCOSITY
	  {nu,nu_hyper2,nu_hyper3};
#else
	  {0.0,0.0,0.0};
#endif
}
std::array<AcReal,3>
magnetic_get_max_diffus()
{
	return
#if LMAGNETIC
	  {eta,eta_hyper2,eta_hyper3};
#else
	  {0.0,0.0,0.0};
#endif
}
std::array<AcReal,3>
energy_get_max_diffus()
{
 	constexpr AcReal chi_hyper2=0.;
	return
#if LENTROPY
	  {gamma*chi,gamma*chi_hyper2,gamma*chi_hyper3};
#else
	  {0.0,0.0,0.0};
#endif
}
std::array<AcReal,3>
elem_wise_max(const std::array<AcReal,3>& a,const std::array<AcReal,3>& b,const std::array<AcReal,3>& c)
{
	return 
	{
	  std::max(std::max(a[0],b[0]),c[0]),
	  std::max(std::max(a[1],b[1]),c[1]),
	  std::max(std::max(a[2],b[2]),c[2])
	};
}

AcReal max_diffus(AcReal );
/***********************************************************************************************/
bool
has_nans(AcMesh mesh_in);
/***********************************************************************************************/
AcReal
sign(const AcReal a, const AcReal b)
{
	if (b < (AcReal)0.0)
		return -abs(a);
	else
		return abs(a);
}
/***********************************************************************************************/
AcReal
calc_dt1_courant()
{
#if TRANSPILATION
	return acDeviceGetOutput(acGridGetDevice(),AC_dt1_max);
#endif
	//TP: temporary measure should not be done but at the moment we want to compile TG without TRANSPILATION=on
#if TRAINING
	return acDeviceGetOutput(acGridGetDevice(),AC_dt1_max);
#endif
      AcReal maxadvec = 0.;
#if LHYDRO
      maxadvec = acDeviceGetOutput(acGridGetDevice(), AC_maxadvec)/cdt;
      //if (rank==0) printf("rank, maxadvec= %d %e \n", rank, maxadvec);
#endif
      AcReal maxchi_dyn = 0.;
#if LENTROPY
      maxchi_dyn = acDeviceGetOutput(acGridGetDevice(), AC_maxchi);
#endif
      //fprintf(stderr, "HMM MAX ADVEC, DIFFUS: %14e, %14e\n",maxadvec,max_diffus());
      return (AcReal)sqrt(pow(maxadvec, 2) + pow(max_diffus(maxchi_dyn), 2));
      //acDeviceSetInput(acGridGetDevice(),AC_dt,dt);
      //if (rank==0) printf("rank, maxadvec, maxdiffus, dt1_= %d %e %e %e \n", rank, maxadvec,max_diffus(maxchi_dyn), dt1_);
}
/***********************************************************************************************/
AcReal
GpuCalcDt()
{
	acGridSynchronizeStream(STREAM_ALL);
  	acDeviceSetInput(acGridGetDevice(), AC_step_num, (PC_SUB_STEP_NUMBER) 0);
	const auto graph = acGetOptimizedDSLTaskGraph(AC_calculate_timestep);
	acGridExecuteTaskGraph(graph,1);
	acGridSynchronizeStream(STREAM_ALL);
        acDeviceSwapBuffers(acGridGetDevice());
	return calc_dt1_courant();
}
/***********************************************************************************************/
extern "C" void substepGPU(int isubstep)
//
//  Do the 'isubstep'th integration step on all GPUs on the node and handle boundaries.
//
{
   //TP: with on does timestepping the way PC does it
   //TP: logs performance metrics of Astaroth
   const bool log = false;
#if LFORCING
  //Update forcing params
   if (lsecond_force) 
   {
	   fprintf(stderr,"Second forcing force not yet implemented on GPU!\n");
	   exit(EXIT_FAILURE);
   }
   if (isubstep == itorder) forcing_params.Update();  // calculate on CPU and load into GPU
#endif

  acDeviceSetInput(acGridGetDevice(), AC_step_num,(PC_SUB_STEP_NUMBER) (isubstep-1));
  if(lshear && isubstep == 1) acDeviceSetInput(acGridGetDevice(), AC_shear_delta_y,deltay);
  Device dev = acGridGetDevice();
  //TP: done in this more complex manner to ensure the actually integrated time and the time reported by Pencil agree
  //if we call set_dt after the first timestep there would be slight shift in dt what Pencil sees and what is actually used for time integration
  
  if (isubstep == 1) 
  {
	  //TP: done to have the same timestep as PC when testing
	  if (ldt && lcpu_timestep_on_gpu) dt1_interface = GpuCalcDt();
	  if (ldt) set_dt(dt1_interface);
	  acDeviceSetInput(acGridGetDevice(), AC_dt,dt);
  }
  //fprintf(stderr,"before acGridExecuteTaskGraph");
  AcTaskGraph *rhs =  acGetOptimizedDSLTaskGraph(AC_rhs);
  auto start = MPI_Wtime();
  acGridExecuteTaskGraph(rhs, 1);
  auto end = MPI_Wtime();
  if (log && !rank) fprintf(stderr,"RHS TOOK: %14e\n",end-start);
  if (ldt &&
        ((isubstep == 5 && !lcourant_dt) 
        || (isubstep == 1 && lcourant_dt)
     ))
  {
    constexpr AcReal unit = 1.0;
    AcReal dt1_{};
    if (!lcourant_dt)
    {
      const AcReal maximum_error = acDeviceGetOutput(acGridGetDevice(), AC_maximum_error)/eps_rkf;
      AcReal dt_{};
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
      dt1_ = unit/dt_;
    }
    else 
    {
      dt1_ = calc_dt1_courant();
    }
    dt1_interface = dt1_;
  }
  return;
}
/***********************************************************************************************/
/**
extern "C" void testBcKernel(AcReal *farray_in, AcReal *farray_truth)
{
  AcMesh mesh_true;
  AcMesh mesh_test;
  AcReal epsilon = 0.00001;

  size_t offset = 0;
  for (int i = 0; i < NUM_VTXBUF_HANDLES; ++i)
  {
    mesh_test.vertex_buffer[VertexBufferHandle(i)] = &farray_in[offset];
    offset += mw;
  }
  offset = 0;
  for (int i = 0; i < NUM_VTXBUF_HANDLES; ++i)
  {
    mesh_true.vertex_buffer[VertexBufferHandle(i)] = &farray_truth[offset];
    offset += mw;
  }

  // AcReal cv1 = 1.66666663;
  // lnrho0 = 0;
  // cp = 1.0;
  // AcReal gamma_m1 = 0.666666627;

  //Emulate the gpu code serially
  int3 dims = {mx, my, 1};
  for (int i = 0; i < dims.x; i++)
  {
    for (int j = 0; j < dims.y; j++)
    {
      for (int k = 0; k < dims.z; k++)
      {
        //Remember to convert to zero based index :))))))
        // AcReal ftopktop = 1.0;
        // AcReal rho_xy=exp(mesh.vertex_buffer[ilnrho-1][DEVICE_VTXBUF_IDX(i,j,n2-1)]);
        // AcReal cs2_xy = mesh.vertex_buffer[iss-1][DEVICE_VTXBUF_IDX(i,j,n2-1)];
        // cs2_xy=cs20*exp(gamma_m1*(mesh.vertex_buffer[ilnrho-1][DEVICE_VTXBUF_IDX(i,j,n2-1)]-lnrho0)+cv1*cs2_xy);
        // AcReal tmp_xy=ftopktop/cs2_xy;
        // for (int z=1;z<=3;z++){
        //     rho_xy = mesh.vertex_buffer[ilnrho-1][DEVICE_VTXBUF_IDX(i,j,n2+z-1)]-mesh.vertex_buffer[ilnrho-1][DEVICE_VTXBUF_IDX(i,j,n2-z-1)];
        //     mesh.vertex_buffer[iss-1][DEVICE_VTXBUF_IDX(i,j,n2+z-1)] = mesh.vertex_buffer[iss-1][DEVICE_VTXBUF_IDX(i,j,n2-z-1)] + cp*(cp-cv)*(-rho_xy-1.0*tmp_xy);
        // }
        //#include "res.cuh"
      }
    }
  }
  bool passed = true;
  for (int i = 0; i < mx; i++)
  {

    for (int j = 0; j < my; j++)
    {
      for (int k = 0; k < mz; k++)
      {
        for (int ivar = 0; ivar < mfarray; ivar++)
        {
          AcReal out_val = mesh_test.vertex_buffer[ivar][DEVICE_VTXBUF_IDX(i, j, k)];
          AcReal true_val = mesh_true.vertex_buffer[ivar][DEVICE_VTXBUF_IDX(i, j, k)];
          if (fabs(out_val - true_val) > epsilon)
          {
            passed = false;
            acLogFromRootProc(rank,"C val wrong at %d,%d,%d\n", i, j, k);
            acLogFromRootProc(rank,"field = %d", ivar);
            acLogFromRootProc(rank,"C val: %f\tF val: %f\n", out_val, true_val);
          }
        }
      }
    }
  }
  if (passed)
  {
    acLogFromRootProc(rank,"Passed C test :)\n");
  }
  else
  {
    acLogFromRootProc(rank,"Did not pass C test :(\n");
  }
  if (!passed) return;
  acLogFromRootProc(rank,"Starting GPUtest\n");
  fflush(stdout);
  acGridSynchronizeStream(STREAM_ALL);
  // acGridLoadMesh(STREAM_DEFAULT,mesh);
  // acLogFromRootProc(rank,"loaded mesh\n");
  // acGridTestBCKernel({mx,my,1});

  acGridSynchronizeStream(STREAM_ALL);
  acLogFromRootProc(rank,"after bc kernel\n");
  fflush(stdout);
  acGridStoreMesh(STREAM_DEFAULT, &mesh);
  acGridSynchronizeStream(STREAM_ALL);
  acLogFromRootProc(rank,"after store\n");
  fflush(stdout);
  AcReal max_abs_not_passed=-1.0;
  for (int i = 0; i < mx; i++)
  {
    for (int j = 0; j < my; j++)
    {
      for (int k = 0; k < mz; k++)
      {
        for (int ivar = 0; ivar < mfarray; ivar++)
        {
          AcReal out_val = mesh.vertex_buffer[ivar][DEVICE_VTXBUF_IDX(i, j, k)];
          AcReal true_val = mesh_true.vertex_buffer[ivar][DEVICE_VTXBUF_IDX(i, j, k)];
          if (fabs(out_val - true_val) > epsilon)
          {
            passed = false;
            acLogFromRootProc(rank,"GPU val wrong at %d,%d,%d\n", i, j, k);
            acLogFromRootProc(rank,"field = %d", ivar);
            acLogFromRootProc(rank,"GPU val: %f\tTRUE val: %f\tDIFF: %f\n", out_val, true_val, fabs(out_val - true_val));
            if (fabs(out_val)>max_abs_not_passed) max_abs_not_passed = out_val;
          }
        }
      }
    }
  }
  acLogFromRootProc(rank,"maximum abs incorrect val: %f\n",max_abs_not_passed);
  if (passed)
  {
    acLogFromRootProc(rank,"Passed GPU test :)\n");
  }
  else
  {
    acLogFromRootProc(rank,"Did not pass GPU test :(\n");
  }
  return;
}
**/
/***********************************************************************************************/
extern "C" void registerGPU()
{
  // AcReal* profile_x_host = (AcReal*)malloc(sizeof(AcReal)*mx);
  // for (int i=0;i<mx;i++){
  //     profile_x_host[i] = (AcReal)i;
  //     acLogFromRootProc(rank,"profile_x_host[%d]=%f\n",i,profile_x_host[i]);
  // }
  // mesh.profiles[PROFILE_X] = profile_x_host;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  mesh.info = acInitInfo();
#if AC_RUNTIME_COMPILATION
#else
  AcResult res = acCheckDeviceAvailability();

  if(res == AC_FAILURE) 
  {
	  fprintf(stderr,"No devices!\n");
	  exit(EXIT_FAILURE);
  }
#endif
}
/***********************************************************************************************/
extern "C" void initGPU()
{
  // Check whether there is (at least) one GPU available
  //TP: moved to initializeGPU since with runtime compilation should call only after Astaroth is loaded
  //AcResult res = acCheckDeviceAvailability();
}
/***********************************************************************************************/
#define PCLoad acPushToConfig
MPI_Comm comm_pencil = MPI_COMM_NULL;
void modulepars(AcMeshInfo& config){
  // Enter basic parameters in config.
  #include "PC_modulepars.h"
}
void setupConfig(AcMeshInfo& config)
{ 
  modulepars(config);
  //TP: loads for non-cartesian derivatives
#if TRANSPILATION
  PCLoad(config, AC_inv_cyl_r,rcyl_mn1);
  PCLoad(config, AC_inv_r,r1_mn);
  PCLoad(config, AC_inv_sin_theta,sin1th);
  PCLoad(config, AC_cot_theta,cotth);

  //TP: loads for non-equidistant grids
  PCLoad(config,AC_nonequidistant_grid, (AcBool3){!lequidist.x,!lequidist.y,!lequidist.z});
  PCLoad(config,AC_inv_mapping_func_derivative_x,dx_1);
  PCLoad(config,AC_inv_mapping_func_derivative_y,dy_1);
  PCLoad(config,AC_inv_mapping_func_derivative_z,dz_1);

  PCLoad(config,AC_mapping_func_tilde_x,dx_tilde);
  PCLoad(config,AC_mapping_func_tilde_y,dy_tilde);
  PCLoad(config,AC_mapping_func_tilde_z,dz_tilde);
#endif

  PCLoad(config, AC_rk_order, itorder);
  PCLoad(config, AC_shear,lshear);

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
  PCLoad(config, AC_domain_decomposition, (int3) {nprocx,nprocy,nprocz});
  PCLoad(config, AC_ngrid, (int3){nxgrid,nygrid,nzgrid});
  PCLoad(config, AC_skip_single_gpu_optim, true);

  PCLoad(config,AC_decompose_strategy,AC_DECOMPOSE_STRATEGY_EXTERNAL);
  if (lmorton_curve)
    PCLoad(config,AC_proc_mapping_strategy,AC_PROC_MAPPING_STRATEGY_MORTON);
  else
    PCLoad(config,AC_proc_mapping_strategy,AC_PROC_MAPPING_STRATEGY_LINEAR);

	
	PCLoad(config, AC_include_3d_halo_corners, ltraining);
  PCLoad(config,AC_MPI_comm_strategy,AC_MPI_COMM_STRATEGY_DUP_USER);
  config.comm->handle = comm_pencil;

// grid and geometry related parameters

  PCLoad(config,AC_ds,(AcReal3){dx,dy,dz});
  PCLoad(config,AC_periodic_grid,lperi);
  //TP: Astaroth has default formula for AC_len but we want to make sure AC_len = lxyz, in case the Astaroth formula does not cover everything like non-equidistant grids
  PCLoad(config,AC_len,lxyz);

  PCLoad(config,AC_sparse_autotuning,lac_sparse_autotuning);
  // Enter physics related parameters in config.

  #if LDENSITY
    PCLoad(config,AC_ldensity_nolog,ldensity_nolog);
    //printf("ldensity_nolog is %d \n",config[AC_ldensity_nolog]);//ldensity_nolog);
  #endif
  //Device dev = acGridGetDevice();
#if LHYDRO
  //dev->output.real_outputs[AC_maxadvec]=0.;
#endif
#if LENTROPY
  //dev->output.real_outputs[AC_maxchi]=0.;
#endif
  acHostUpdateParams(&config); 
  if(!ltraining) config.runtime_compilation_log_dst = "ac_compilation.log";
  char cwd[9000];
  cwd[0] = '\0';
  const char* err = getcwd(cwd, sizeof(cwd));
  if(err == NULL) 
  {
	  fprintf(stderr,"Was not able to get cwd!\n");
	  exit(EXIT_FAILURE);
  }
  char build_path[18000];
  sprintf(build_path,"%s/src/astaroth/submodule/build",cwd);
  config.runtime_compilation_build_path = strdup(build_path);
}

#undef x
#undef y
#undef z
/***********************************************************************************************/
void checkConfig(AcMeshInfo &config)
{
 acLogFromRootProc(rank,"Check that config is correct\n");
 acLogFromRootProc(rank,"d[xyz]: %.14f %.14f %.14f \n", dx, dy, dz);
// acLogFromRootProc(rank,"rank= %d: l1, l2, n1, n2, m1, m2= %d %d %d %d %d %d \n", rank, l1, l2, n1, n2, m1, m2);
// acLogFromRootProc(rank,"zlen= %.14f %.14f \n", config[AC_len].z, lxyz[2]);
 /*
#if LHYDRO
 acLogFromRootProc(rank,"lpressuregradientgas= %d %d \n", lpressuregradient_gas, config[AC_lpressuregradient_gas]);
#endif
#if LENTROPY
 acLogFromRootProc(rank,"chi= %f %f \n", chi, config[AC_chi]);
 acLogFromRootProc(rank,"nkramers= %f %f \n", nkramers, config[AC_nkramers]);
 acLogFromRootProc(rank,"hcond0_kramers= %f %f \n", hcond0_kramers, config[AC_hcond0_kramers]);
 acLogFromRootProc(rank,"hcond_Kconst= %f %f \n", hcond_Kconst, config[AC_hcond_Kconst]);
 //acLogFromRootProc(rank,"Fbot= %f %f \n", Fbot, config[AC_Fbot]);
 acLogFromRootProc(rank,"chi_t= %f %f \n", chi_t, config[AC_chi_t]);
#endif
#if LVISCOSITY
 acLogFromRootProc(rank,"nu= %f %f \n", nu, config[AC_nu]);
 acLogFromRootProc(rank,"zeta= %f %f \n", zeta, config[AC_zeta]);
#endif
#if LMAGNETIC
  acLogFromRootProc(rank,"eta= %f %f \n", eta, config[AC_eta]);
#endif
#if LEOS
  acLogFromRootProc(rank,"cs20= %f %f \n", cs20, config[AC_cs20]);
  //  acLogFromRootProc(rank,"gamma= %f %f \n", gamma, get_real_param(config,comp_infoAC_gamma));
  acLogFromRootProc(rank,"gamma_m1= %f %f \n", gamma_m1, config[AC_gamma_m1]);
  acLogFromRootProc(rank,"gamma1= %f %f \n", gamma1, config[AC_gamma1]);
  acLogFromRootProc(rank,"cv= %f %f \n", cv, config[AC_cv]);
  acLogFromRootProc(rank,"cp= %f %f \n", cp, config[AC_cp]);
  acLogFromRootProc(rank,"lnT0= %f %f \n", lnTT0, config[AC_lnTT0]);
  acLogFromRootProc(rank,"lnrho0= %f %f \n", lnrho0, config[AC_lnrho0]);
#endif
#if LFORCING
  acLogFromRootProc(rank,"iforcing_zsym= %f %f \n", iforcing_zsym, config[AC_iforcing_zsym]);
  acLogFromRootProc(rank,"k1_ff= %f %f \n", k1_ff, config[AC_k1_ff]);
  acLogFromRootProc(rank,"tforce_stop= %f %f \n", tforce_stop, config[AC_tforce_stop]);
  acLogFromRootProc(rank,"k1_ff,profx_ampl, val= %f %d %lf %lf\n", k1_ff, profx_ampl, profx_ampl[0], profx_ampl[nx-1]);
#endif
  acLogFromRootProc(rank,"mu0= %f %f \n", (double)mu0, (double)config[AC_mu0]);
*/
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
void
testBCs();     // forward declaration
/***********************************************************************************************/
extern "C" void loadFarray(); // forward declaration
/***********************************************************************************************/
extern "C" void reloadConfig(); // forward declaration

void
autotune_all_integration_substeps()
{
  for (int i = 0; i < num_substeps; ++i)
  {
  	acDeviceSetInput(acGridGetDevice(), AC_step_num,(PC_SUB_STEP_NUMBER)i);
if (rank==0 && ldebug) printf("memusage before GetOptimizedDSLTaskGraph= %f MBytes\n", memusage()/1024.);
	acGetOptimizedDSLTaskGraph(AC_rhs);
if (rank==0 && ldebug) printf("memusage after GetOptimizedDSLTaskGraph= %f MBytes\n", memusage()/1024.);
  }
}

/***********************************************************************************************/
extern "C" void initializeGPU(AcReal *farr, int comm_fint)
{
  //Setup configurations used for initializing and running the GPU code
#if PACKED_DATA_TRANSFERS
  //initLoadStore();
#endif
  comm_pencil = MPI_Comm_f2c(comm_fint);
  setupConfig(mesh.info);
#if TRAINING
  #include "user_constants.h"
  fprintf(stderr,"INDEX OF TAU.xx: %d\n",acGetTAU_Xx());
/**
	{
		#include "user_constants.h"
		if(itauxx-1 != TAU.xx)
		{
			fprintf(stderr,"Mismatch of indeces for tauxx : %d,%d!!\n",itauxx,TAU.xx);
			exit(EXIT_FAILURE);
		}
		if(itauxy-1 != TAU.xy)
		{
			fprintf(stderr,"Mismatch of indeces for tauxy !!\n");
			exit(EXIT_FAILURE);
		}
		if(itauxz-1 != TAU.xz)
		{
			fprintf(stderr,"Mismatch of indeces for tauxz !!\n");
			exit(EXIT_FAILURE);
		}
		if(itauyy-1 != TAU.yy)
		{
			fprintf(stderr,"Mismatch of indeces for tauyy!!\n");
			exit(EXIT_FAILURE);
		}
		if(itauyz-1 != TAU.yz)
		{
			fprintf(stderr,"Mismatch of indeces for tauyz !!\n");
			exit(EXIT_FAILURE);
		}
		if(itauzz-1 != TAU.zz)
		{
			fprintf(stderr,"Mismatch of indeces for tauzz !!\n");
			exit(EXIT_FAILURE);
		}
	}
	**/
#endif
  //TP: done after setupConfig since we need maux_vtxbuf_index
  //TP: this is an ugly way to do this but works for now
  {
    size_t offset = 0;
    for (int i = 0; i < mvar; ++i)
    {
      mesh.vertex_buffer[VertexBufferHandle(i)] = &farr[offset];
      offset += mw;
    }

    for (int i = 0; i < mfarray; ++i)
    {
      if (maux_vtxbuf_index[i])
      {
        mesh.vertex_buffer[maux_vtxbuf_index[i]] = &farr[mw*i];
      }
    }
    //TP: for now for training we have all slots filled since we might want to read TAU components to the host for calculating validation error
    if(ltraining)
    {
    	for(int i = 0; i < NUM_VTXBUF_HANDLES; ++i)
    	{
	   if(mesh.vertex_buffer[i] == NULL)
	   {
    	    	mesh.vertex_buffer[i] = (AcReal*)malloc(sizeof(AcReal)*mw);
	   }
    	}
    }
  }
if (rank==0 && ldebug) printf("memusage after pointer assign= %f MBytes\n", memusage()/1024.);
#if AC_RUNTIME_COMPILATION
#include "cmake_options.h"
  acCompile(cmake_options,mesh.info);
  acLoadLibrary(rank == 0 ? stderr : NULL,mesh.info);
  acCheckDeviceAvailability();
  acLogFromRootProc(rank, "Done setupConfig && acCompile\n");
  fflush(stdout);
#else
  acCheckDeviceAvailability();
  acLogFromRootProc(rank, "Done setupConfig\n");
  fflush(stdout);
#endif
  checkConfig(mesh.info);
if (rank==0 && ldebug) printf("memusage grid_init= %f MBytes\n", memusage()/1024.);
  acGridInit(mesh);
if (rank==0 && ldebug) printf("memusage after grid_init= %f MBytes\n", memusage()/1024.);

  mesh.info = acGridDecomposeMeshInfo(mesh.info);
  //TP: important to do before autotuning
  acDeviceSetInput(acGridGetDevice(), AC_step_num,(PC_SUB_STEP_NUMBER)0);
  acDeviceSetInput(acGridGetDevice(), AC_dt,dt);
  acDeviceSetInput(acGridGetDevice(), AC_shear_delta_y,deltay);
		
  if (ltest_bcs) testBCs();
  autotune_all_integration_substeps();
if (rank==0 && ldebug) printf("memusage before store config= %f MBytes\n", memusage()/1024.);
  acStoreConfig(acDeviceGetLocalConfig(acGridGetDevice()), "PC-AC.conf");
if (rank==0 && ldebug) printf("memusage after store config= %f MBytes\n", memusage()/1024.);
  acGridSynchronizeStream(STREAM_ALL);
if (rank==0 && ldebug) printf("memusage after store synchronize stream= %f MBytes\n", memusage()/1024.);
  acLogFromRootProc(rank, "DONE initializeGPU\n");
  fflush(stdout);
  constexpr AcReal unit = 1.0;
  dt1_interface = unit/dt;


}
/***********************************************************************************************/
extern "C" void copyFarray(AcReal* f)
{
  /*
  if (has_nans(mesh)){
    acLogFromRootProc(rank,"found nans while copying\n");
    exit(0);
  }
  */
  AcMesh mesh_to_copy;
  size_t offset = 0;
  for (int i = 0; i < NUM_VTXBUF_HANDLES; ++i)
  {
    mesh_to_copy.vertex_buffer[VertexBufferHandle(i)] = &f[offset];
//printf("mesh %d %p \n",i, mesh.vertex_buffer[VertexBufferHandle(i)]);
    offset += mw;
  }
  mesh_to_copy.info = mesh.info;

  acGridSynchronizeStream(STREAM_ALL);
  //acDeviceStoreMesh(acGridGetDevice(), STREAM_DEFAULT, &mesh_to_copy);
  //acDeviceStoreMesh(acGridGetDevice(), STREAM_DEFAULT, &mesh);
  //TP: for now only copy the advanced fields back
  //TODO: should auxiliaries needed on the GPU like e.g. Shock be copied? They can always be recomputed on the host if needed
  //If doing training we read all since we might want TAU components to calculate e.g. validation error
  const int end = ltraining ? NUM_VTXBUF_HANDLES : 
	  	  lread_all_vars_from_device ? mfarray :
		  mvar;
  for (int i = 0; i < end; ++i)
  {
	  acDeviceStoreVertexBuffer(acGridGetDevice(),STREAM_DEFAULT,VertexBufferHandle(i),&mesh);
  }
  acGridSynchronizeStream(STREAM_ALL);
}
/***********************************************************************************************/
void
reloadConfig()
{
  mesh.info.run_consts = acInitCompInfo();
  setupConfig(mesh.info);
  acGridSynchronizeStream(STREAM_ALL);
  acDeviceUpdate(acGridGetDevice(), mesh.info);
  acGridSynchronizeStream(STREAM_ALL);
#if AC_RUNTIME_COMPILATION
  //TP: if quitting grid have to safe the current values of the vtxbufs since the device arrays are freed!
  copyFarray(mesh.vertex_buffer[0]);
  acGridQuit();
  acCloseLibrary();
#include "cmake_options.h"
  acCompile(cmake_options,mesh.info);
  acLoadLibrary(rank == 0 ? stderr : NULL,mesh.info);
  acGridInit(mesh);
  acLogFromRootProc(rank, "Done setupConfig && acCompile\n");
  fflush(stdout);
  fflush(stderr);
  //TP: this is important that we don't overwrite the output buffer in middle of a timestep when the output buffer holds some meaning!
  autotune_all_integration_substeps();
  //TP: restore the vtxbuf values before quitting grid
  loadFarray();
#endif
  acLogFromRootProc(rank, "DONE initializeGPU\n");
  fflush(stdout);
  fflush(stderr);
}
/***********************************************************************************************/
void loadFarray()
{
  /**
  if (has_nans(mesh)){
    acLogFromRootProc(rank,"found nans while copying\n");
    exit(0);
  }
  **/
  acGridSynchronizeStream(STREAM_ALL);
  {
    for (int i = 0; i < mvar; ++i)
  	acDeviceLoadVertexBuffer(acGridGetDevice(), STREAM_DEFAULT, mesh, VertexBufferHandle(i));

    for (int i = 0; i < mfarray; ++i)
      if (maux_vtxbuf_index[i])
  		acDeviceLoadVertexBuffer(acGridGetDevice(), STREAM_DEFAULT, mesh, VertexBufferHandle(maux_vtxbuf_index[i]));
  }
  acGridSynchronizeStream(STREAM_ALL);
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
}
/***********************************************************************************************/
extern "C" int updateInConfigArrName(char *name)
{
    int index = -1;
    for (int i=0; i<NUM_REAL_ARRAYS; i++){
       if (strcmp(get_array_info(static_cast<AcRealArrayParam>(i)).name,name)==0) index=i;
    }
    if (index>-1) updateInConfigArr(index);
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
    return index;
}
/**********************************************************************************************/
extern "C" void finalizeGPU()
{
#if PACKED_DATA_TRANSFERS
  //acDeviceStoreMesh(acGridGetDevice(), STREAM_DEFAULT, &mesh);  // needed?
#endif
  // Deallocate everything on the GPUs and reset
  AcResult res = acGridQuit();
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
bool
has_nans(AcMesh mesh_in)
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
        for (size_t ivar = 0; ivar < NUM_VTXBUF_HANDLES; ivar++)
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
AcReal max_diffus(AcReal maxchi_dyn)
{
  AcReal3 dxyz_vals = get_dxyzs();
  auto max_diffusions = elem_wise_max(visc_get_max_diffus(), magnetic_get_max_diffus(), energy_get_max_diffus());
#if LENTROPY
  max_diffusions[0] = std::max(max_diffusions[0],maxchi_dyn);
#endif
  return max_diffusions[0]*dxyz_vals.x/cdtv + max_diffusions[1]*dxyz_vals.y/cdtv2 + max_diffusions[2]*dxyz_vals.z/cdtv3;
}
/***********************************************************************************************/
//TP: this is not written the the most optimally since it needs two extra copies of the mesh where at least the tmp
//could be circumvented by temporarily using the output buffers on the GPU to store the f-array and load back from there
//but if we truly hit the mem limit for now the user can of course simply test the bcs with a smaller mesh and skip the test with a larger mesh

void
sym_z(AcMesh mesh_in)
{
  const auto DEVICE_VTXBUF_IDX = [&](const int x, const int y, const int z)
  				{
					return acVertexBufferIdx(x,y,z,acGridGetLocalMeshInfo());
  				};
  AcMeshDims dims = acGetMeshDims(acGridGetLocalMeshInfo());
  for (size_t i = 0; i < dims.m1.x; i++)
  {
    for (size_t j = 0; j < dims.m1.y; j++)
    {
      for (size_t k = 0; k < dims.m1.z; k++)
      {
	if (
	   i >= NGHOST && i < dims.n1.x &&
	   j >= NGHOST && j < dims.n1.y &&
	   k >= NGHOST && k < dims.n1.z
	  )continue;
	 if (k >= NGHOST && k < dims.n1.z) continue;
        for (int ivar = 0; ivar < NUM_VTXBUF_HANDLES; ivar++)
        {
	 //BOT
	 if (k < NGHOST)
	 {
	 	const auto idx = DEVICE_VTXBUF_IDX(i,j,k);
		const auto offset = NGHOST-k;
	 	const auto domain_z= NGHOST+offset;
	 	const auto domain_idx = DEVICE_VTXBUF_IDX(i,j,domain_z);
		mesh_in.vertex_buffer[ivar][idx] = mesh_in.vertex_buffer[ivar][domain_idx];
		//mesh_in.vertex_buffer[ivar][idx] = mesh.vertex_buffer[ivar][idx];
	 }
	 //TOP:
	 else
	 {
	 	const auto idx = DEVICE_VTXBUF_IDX(i,j,k);
		const auto offset = k-mesh_in.info[AC_nlocal_max].z+1;
	 	const auto domain_z= mesh_in.info[AC_nlocal_max].z-offset;
	 	const auto domain_idx = DEVICE_VTXBUF_IDX(i,j,domain_z);
		mesh_in.vertex_buffer[ivar][idx] = mesh_in.vertex_buffer[ivar][domain_idx];
		//mesh_in.vertex_buffer[ivar][idx] = mesh.vertex_buffer[ivar][idx];
	 }
	}
      }
    }
  }
}
/***********************************************************************************************/
void
check_sym_z(AcMesh mesh_in)
{
  const auto DEVICE_VTXBUF_IDX = [&](const int x, const int y, const int z)
  				{
					return acVertexBufferIdx(x,y,z,acGridGetLocalMeshInfo());
  				};
  if (rank == 1)
  {
  	//printf("HMM: %14e\n",mesh_in.vertex_buffer[0][DEVICE_VTXBUF_IDX(14,15,2)]);
  	//printf("HMM: %14e\n",mesh_in.vertex_buffer[0][DEVICE_VTXBUF_IDX(14,15,2)]);
  	//printf("HMM: %14e\n",mesh_in.vertex_buffer[0][DEVICE_VTXBUF_IDX(26,14,1)]);
  	//printf("HMM: %14e\n",mesh_in.vertex_buffer[0][DEVICE_VTXBUF_IDX(26,14,5)]);
  }
  AcMeshDims dims = acGetMeshDims(acGridGetLocalMeshInfo());
  for (size_t i = 0; i < dims.m1.x; i++)
  {
    for (size_t j = 0; j < dims.m1.y; j++)
    {
      for (size_t k = 0; k < dims.m1.z; k++)
      {
	if (
	   i >= NGHOST && i < dims.n1.x &&
	   j >= NGHOST && j < dims.n1.y &&
	   k >= NGHOST && k < dims.n1.z
	  )continue;
	 if (k >= NGHOST && k < dims.n1.z) continue;
        for (int ivar = 0; ivar < NUM_VTXBUF_HANDLES; ivar++)
        {
	 //BOT
	 if (k < NGHOST)
	 {
	 	const auto idx = DEVICE_VTXBUF_IDX(i,j,k);
		const auto offset = NGHOST-k;
	 	const auto domain_z= NGHOST+offset;
	 	const auto domain_idx = DEVICE_VTXBUF_IDX(i,j,domain_z);
		if (mesh_in.vertex_buffer[ivar][idx] !=  mesh_in.vertex_buffer[ivar][domain_idx])
			printf("WRONG\n");
	 }
	 //TOP:
	 else
	 {
	 	const auto idx = DEVICE_VTXBUF_IDX(i,j,k);
		const auto offset = k-mesh_in.info[AC_nlocal_max].z+1;
	 	const auto domain_z= mesh_in.info[AC_nlocal_max].z-offset;
	 	const auto domain_idx = DEVICE_VTXBUF_IDX(i,j,domain_z);
		mesh_in.vertex_buffer[ivar][idx] = mesh.vertex_buffer[ivar][domain_idx];
		if (mesh_in.vertex_buffer[ivar][idx] !=  mesh_in.vertex_buffer[ivar][domain_idx])
			printf("WRONG\n");
	 }
	}
      }
    }
  }
}
/***********************************************************************************************/
void
check_sym_x(const AcMesh mesh_in)
{
  const auto DEVICE_VTXBUF_IDX = [&](const int x, const int y, const int z)
  				{
					return acVertexBufferIdx(x,y,z,acGridGetLocalMeshInfo());
  				};
  if (rank == 1)
  {
  	//printf("HMM: %14e\n",mesh_in.vertex_buffer[0][DEVICE_VTXBUF_IDX(14,15,2)]);
  	//printf("HMM: %14e\n",mesh_in.vertex_buffer[0][DEVICE_VTXBUF_IDX(14,15,2)]);
  	//printf("HMM: %14e\n",mesh_in.vertex_buffer[0][DEVICE_VTXBUF_IDX(26,14,1)]);
  	//printf("HMM: %14e\n",mesh_in.vertex_buffer[0][DEVICE_VTXBUF_IDX(26,14,5)]);
  }
  AcMeshDims dims = acGetMeshDims(acGridGetLocalMeshInfo());
  for (size_t i = 0; i < dims.m1.x; i++)
  {
    for (size_t j = 0; j < dims.m1.y; j++)
    {
      for (size_t k = 0; k < dims.m1.z; k++)
      {
	if (
	   i >= NGHOST && i < dims.n1.x &&
	   j >= NGHOST && j < dims.n1.y &&
	   k >= NGHOST && k < dims.n1.z
	  )continue;
	 if (i >= NGHOST && i < dims.n1.x) continue;
        for (int ivar = UUY; ivar <= UUY; ivar++)
        {
	 //BOT
	 if (i < NGHOST)
	 {
	 	const auto idx = DEVICE_VTXBUF_IDX(i,j,k);
		const auto offset = NGHOST-i;
	 	const auto domain_x= NGHOST+offset;
	 	const auto domain_idx = DEVICE_VTXBUF_IDX(domain_x,j,k);
		if (mesh_in.vertex_buffer[ivar][idx] !=  mesh_in.vertex_buffer[ivar][domain_idx])
		{
			fprintf(stderr,"WRONG\n");
			exit(EXIT_FAILURE);

		}
	 }
	 //TOP:
	 else
	 {
	 	const auto idx = DEVICE_VTXBUF_IDX(i,j,k);
		const auto offset = i-mesh_in.info[AC_nlocal_max].x+1;
	 	const auto domain_x= mesh_in.info[AC_nlocal_max].x-offset;
	 	const auto domain_idx = DEVICE_VTXBUF_IDX(domain_x,j,k);
		if (mesh_in.vertex_buffer[ivar][idx] !=  mesh_in.vertex_buffer[ivar][domain_idx])
		{
			fprintf(stderr,"WRONG\n");
			exit(EXIT_FAILURE);

		}
	 }
	}
      }
    }
  }
}
/***********************************************************************************************/
void
testBCs()
{
  // Set random seed for reproducibility
  srand(321654987);
  const auto DEVICE_VTXBUF_IDX = [&](const int x, const int y, const int z)
  				{
					return acVertexBufferIdx(x,y,z,acGridGetLocalMeshInfo());
  				};
  AcTaskGraph* rhs = acGetOptimizedDSLTaskGraph(AC_rhs);
  const bool all_periodic = 
  				acGridTaskGraphHasPeriodicBoundcondsX(rhs) && 
  				acGridTaskGraphHasPeriodicBoundcondsY(rhs) && 
  				acGridTaskGraphHasPeriodicBoundcondsZ(rhs) 
				;
  if (all_periodic && !lshear) return;
  if(lshear) acLogFromRootProc(rank,"testBCS: deltay: %7e\n",deltay);
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

  acGridSynchronizeStream(STREAM_ALL);
  acGridExecuteTaskGraph(bcs,1);
  acGridSynchronizeStream(STREAM_ALL);

  acGridSynchronizeStream(STREAM_ALL);
  acDeviceStoreMesh(acGridGetDevice(), STREAM_DEFAULT, &mesh_to_copy);
  acGridSynchronizeStream(STREAM_ALL);

  //check_sym_x(mesh);
  //check_sym_x(mesh_to_copy);

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
  int3 largest_diff_point{};
  //AcReal epsilon = 0.0;

  auto skip_x = (acGridTaskGraphHasPeriodicBoundcondsX(rhs) && !lshear) || false;
  auto skip_y = acGridTaskGraphHasPeriodicBoundcondsY(rhs) || false;
  auto skip_z = acGridTaskGraphHasPeriodicBoundcondsZ(rhs) || false;

  auto start_x =  skip_x ? NGHOST : 0;
  auto start_y =  skip_y ? NGHOST : 0;
  auto start_z =  skip_z ? NGHOST : 0;

  auto end_x =  skip_x ? dims.n1.x : dims.m1.x;
  auto end_y =  skip_y ? dims.n1.y : dims.m1.y;
  auto end_z =  skip_z ? dims.n1.z : dims.m1.z;

  int num_of_points_where_different[NUM_VTXBUF_HANDLES]{};
  bool different_in[NUM_VTXBUF_HANDLES][6]{};
  //TP: stupid levels of safety
  memset(different_in,0,NUM_VTXBUF_HANDLES*6*sizeof(bool));
  const int bot_x = 0;
  const int top_x = 1;

  const int bot_y = 2;
  const int top_y = 3;

  const int bot_z = 4;
  const int top_z = 5;
  int num_of_points = 0;

  for (size_t i = start_x; i < end_x; i++)
  {
    for (size_t j = start_y; j < end_y; j++)
    {
      for (size_t k = start_z; k < end_z; k++)
      {
	if (
	   i >= NGHOST && i < dims.n1.x &&
	   j >= NGHOST && j < dims.n1.y &&
	   k >= NGHOST && k < dims.n1.z
	  )continue;
	  //if (i < NGHOST | i > dims.n1.x  || j < NGHOST || j > dims.n1.y) continue;
	++num_of_points;
        for (int ivar = 0; ivar < mvar; ivar++)
        {
          AcReal out_val = mesh_to_copy.vertex_buffer[ivar][DEVICE_VTXBUF_IDX(i, j, k)];
          AcReal true_val = mesh.vertex_buffer[ivar][DEVICE_VTXBUF_IDX(i, j, k)];
          AcReal abs_diff = fabs(out_val - true_val);
          if (fabs(true_val) > max_abs_value) max_abs_value = fabs(true_val);
          if (fabs(true_val) < min_abs_value) min_abs_value = fabs(true_val);
          if ((abs_diff/true_val) > epsilon || (true_val == (AcReal)0.0 && fabs(out_val) > (AcReal)pow(0.1,13)) || (epsilon == (AcReal)0.0 && true_val != out_val))
          {
	    different_in[ivar][bot_x] |= i < NGHOST;
	    different_in[ivar][top_x] |= i >= dims.n1.x;

	    different_in[ivar][bot_y] |= j < NGHOST;
	    different_in[ivar][top_y] |= j >= dims.n1.y;

	    different_in[ivar][bot_z] |= k < NGHOST;
	    different_in[ivar][top_z] |= k >= dims.n1.z;

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
		largest_diff_point = (int3){(int)i,(int)j,(int)k};
              }
            }  
          }
          if (isnan(out_val))
          {
            acLogFromRootProc(rank,"nan before at %d,%d,%d,%d!\n!",i,j,k,ivar);
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
    		acLogFromRootProc(0,"ratio of values wrong for field: %s\t %f\n",field_names[ivar],(double)num_of_points_where_different[ivar]/num_of_points);
		if (different_in[ivar][bot_x]) acLogFromRootProc(0,"different in BOT_X\n");
		if (different_in[ivar][top_x]) acLogFromRootProc(0,"different in TOP_X\n");

		if (different_in[ivar][bot_y]) acLogFromRootProc(0,"different in BOT_Y\n");
		if (different_in[ivar][top_y]) acLogFromRootProc(0,"different in TOP_Y\n");

		if (different_in[ivar][bot_z]) acLogFromRootProc(0,"different in BOT_Z\n");
		if (different_in[ivar][top_z]) acLogFromRootProc(0,"different in TOP_Z\n");
	}
  	acLogFromRootProc(0,"max abs not passed val: %.7e\t%.7e\n",(double)max_abs_not_passed_val, (double)fabs(true_pair));
  	acLogFromRootProc(0,"max abs relative difference val: %.7e\n",(double)max_abs_relative_difference);
  	acLogFromRootProc(0,"Point where biggest rel diff: %d,%d,%d\n",largest_diff_point.x,largest_diff_point.y,largest_diff_point.z);
  	acLogFromRootProc(0,"largest difference: %.7e\t%.7e\n",(double)gpu_val_for_largest_diff, (double)true_val_for_largest_diff);
  	acLogFromRootProc(0,"abs range: %.7e-%7e\n",(double)min_abs_value,(double)max_abs_value);
    	acLogFromRootProc(0,"Did not pass BC test :(\n");
	fprintf(stderr,"Did not pass BC\n");

    	exit(EXIT_FAILURE);
  }
  fflush(stdout);

  acHostMeshDestroy(&mesh_to_copy);
  acHostMeshCopyVertexBuffers(tmp_mesh_to_store,mesh);
  acHostMeshDestroy(&tmp_mesh_to_store);
  acGridDestroyTaskGraph(bcs);
  acHostMeshDestroy(&tmp_mesh_to_store);
  acGridDestroyTaskGraph(bcs);
}
/***********************************************************************************************/
extern "C" void
gpuSetDt()
{
	acGridSynchronizeStream(STREAM_ALL);
	if(!lcourant_dt)
	{
		fprintf(stderr,"gpuSetDt works only for Courant timestep!!\n");
		exit(EXIT_FAILURE);
	}
	//TP: not needed but for extra safety
  	acDeviceSetInput(acGridGetDevice(), AC_step_num, (PC_SUB_STEP_NUMBER) 0);
	const auto graph = acGetOptimizedDSLTaskGraph(AC_calculate_timestep);

	acGridExecuteTaskGraph(graph,1);
	acGridSynchronizeStream(STREAM_ALL);
	AcReal dt1_ = calc_dt1_courant();
	set_dt(dt1_);
	dt1_interface = dt1_;
        acDeviceSwapBuffers(acGridGetDevice());
	loadFarray();
	//TP: not strictly needed but for extra safety
}
/***********************************************************************************************/
