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

#define CUDA_ERRCHK(X)

// Astaroth headers.
#include "astaroth.h"
#include "math_utils.h"
#define real AcReal
#define EXTERN
#define FINT int

#if DOUBLE_PRECISION
  #define REAL_MAX DBL_MAX
#else
  #define REAL_MAX FLT_MAX
#endif

__thread int tp_int;
typedef void (*rangefunc)(const int a, const int b);

AcReal cpu_pow(AcReal const val, AcReal exponent)
{
// masks hip GPU power function.
  return std::pow(val, exponent);
}
// PC interface headers.
#include "PC_moduleflags.h"
#include "../cparam_c.h"
#include "../cdata_c.h"
#if LFORCING
  #include "../forcing_c.h"     // provides forcing_pars_hel
#endif
#include "../sub_c.h"           // provides set_dt
#include "../boundcond_c.h"     // provides boundconds[xyz] etc.
#include "../mpicomm_c.h"       // provides finalize_sendrcv_bdry
#include "PC_module_parfuncs.h" // provides stuff from physics modules

#if PACKED_DATA_TRANSFERS
  #include "loadStore.h"
#endif
#if LFORCING
  #include "forcing.h"
#endif
// Astaroth objects instantiation.
static AcMesh mesh;
//static AcMesh test_mesh;
static AcTaskGraph *rhs;
static AcTaskGraph *randomize_graph;
static AcTaskGraph *rhs_test_graph;
static AcTaskGraph *rhs_test_rhs_1;

// Other.
static int rank;
static MPI_Comm comm_pencil;
int halo_xz_size[2] = {0, 0}, halo_yz_size[2] = {0, 0};
//static AcReal *xtop_buffer, *xbot_buffer, *ytop_buffer, *ybot_buffer;


/***********************************************************************************************/
int DCONST(const AcIntParam param)
{
  return mesh.info.int_params[param];
}
/***********************************************************************************************/
int3 DCONST(const AcInt3Param param)
{
  return mesh.info.int3_params[param];
}
/***********************************************************************************************/
AcReal DCONST(const AcRealParam param)
{
  return mesh.info.real_params[param];
}
/***********************************************************************************************/
AcReal3 DCONST(const AcReal3Param param)
{
  return mesh.info.real3_params[param];
}
/***********************************************************************************************/
int DEVICE_VTXBUF_IDX(const int x_coordinate, const int y_coordinate, const int z_coordinate)
{
  return x_coordinate + mx * y_coordinate + mx * my * z_coordinate;
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
AcReal max_advec()
{
#if LHYDRO
  AcReal umax = 0.;
  acGridReduceVec(STREAM_DEFAULT, RTYPE_MAX, UUX, UUY, UUZ, &umax);
  return umax/sqrt(get_dxyzs().x);
#else
  return 0.
#endif
  
}
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

AcReal max_diffus()
{
  AcReal3 dxyz_vals = get_dxyzs();
  auto max_diffusions = elem_wise_max(visc_get_max_diffus(), magnetic_get_max_diffus(), energy_get_max_diffus());
  return max_diffusions[0]*dxyz_vals.x/cdtv + max_diffusions[1]*dxyz_vals.y/cdtv2 + max_diffusions[2]*dxyz_vals.z/cdtv3;
}
/***********************************************************************************************/
int id_to_tag(int3 id)
{
  return ((3 + id.x) % 3) * 9 + ((3 + id.y) % 3) * 3 + (3 + id.z) % 3;
}
/***********************************************************************************************/
int3 tag_to_id(int _tag)
{
  int3 _id = (int3){(_tag) / 9, ((_tag) % 9) / 3, (_tag) % 3};
  _id.x = _id.x == 2 ? -1 : _id.x;
  _id.y = _id.y == 2 ? -1 : _id.y;
  _id.z = _id.z == 2 ? -1 : _id.z;
  ERRCHK_ALWAYS(id_to_tag(_id) == _tag);
  return _id;
}
/***********************************************************************************************/
bool
has_nans(AcMesh mesh_in)
{
  bool res = false;
  AcMeshDims dims = acGetMeshDims(acGridGetLocalMeshInfo());
  for (int i = dims.n0.x; i < dims.n1.x; i++)
  {
    for (int j = dims.n0.y; j < dims.n1.y; j++)
    {
      for (int k = dims.n0.z; k < dims.n1.z; k++)
      {
        for (int ivar = 0; ivar < NUM_VTXBUF_HANDLES; ivar++)
        {
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
extern "C" void substepGPU(int isubstep)
//
//  Do the 'isubstep'th integration step on all GPUs on the node and handle boundaries.
//
{
#if LFORCING
  //Update forcing params
   if (isubstep == itorder) forcing_params.Update();  // calculate on CPU and load into GPU
#endif

  //acGridSynchronizeStream(STREAM_DEFAULT);
  //Transfer the updated ghost zone to the device(s) in the node

  //if (full)
  //{
  //  if (has_nans(mesh))
  //  {
  //    acLogFromRootProc(rank,"had nans before starting GPU comp\n");
  //    exit(0);
  //  }
  //  acLogFromRootProc(rank,"doing full i.e. loading\n");
  //  acGridSynchronizeStream(STREAM_ALL);
  //  acDeviceLoadMesh(acGridGetDevice(), STREAM_DEFAULT, mesh);
  //  acGridSynchronizeStream(STREAM_ALL);
  //  //set output buffer to 0 since if we are reading from it we don't want NaNs
  //  AcMeshDims dims = acGetMeshDims(acGridGetLocalMeshInfo());
  //  acGridSynchronizeStream(STREAM_ALL);
  //}
  acDeviceSetInput(acGridGetDevice(), AC_step_num,isubstep-1);
  //acGridSynchronizeStream(STREAM_ALL);
  Device dev = acGridGetDevice();
  if (isubstep == 1)
    acDeviceSetInput(acGridGetDevice(), AC_dt,dt);
  acGridExecuteTaskGraph(rhs, 1);
  if (isubstep == 1 && ldt)
  {
      acGridSynchronizeStream(STREAM_ALL);
      acGridFinalizeReduceLocal(rhs);
#if LHYDRO
      AcReal maxadvec = acDeviceGetOutput(acGridGetDevice(), AC_maxadvec)/cdt;
#endif
#if LENTROPY
      //AcReal maxchi    = dev->output.real_outputs[AC_maxchi];
      //maxdiffus = maxdiffus+maxchi/pow(dx,2)
#endif
      AcReal dt1_ = sqrt(pow(maxadvec, 2) + pow(max_diffus(), 2));
//printf("maxadvec, maxchi,maxdiffus= %e %e %e \n", maxadvec, maxchi, maxdiffus);
//printf("maxadvec, maxdiffus= %e %e %e \n", maxadvec, max_diffus());
      set_dt(dt1_);
      acDeviceSetInput(acGridGetDevice(),AC_dt,dt);
  }
  //acGridSynchronizeStream(STREAM_ALL);
  //acGridSynchronizeStream(STREAM_ALL);
  // acLogFromRootProc(rank,"Done substep: %d\n",isubstep);
  // fflush(stdout);
  // int found_nan;
  // FILE* diag_file = fopen("astaroth_timeseries.ts", "a");
  // print_diagnostics(rank, 1, dt, dt,diag_file, 0.0001, 0.0001, &found_nan);
  return;
}
/***********************************************************************************************/
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
/***********************************************************************************************/
extern "C" void testRHS(AcReal *farray_in, AcReal *dfarray_truth)
{
  // __energy_MOD_pushpars2c(p_par_energy);
  // mesh.info.profiles[AC_dlnhcond_prof]=dlnhcond_prof; // [2-1] dlnhcond_prof real(nz)
  // for (int i=20;i<30;++i)
  //   acLogFromRootProc(rank,"C: dlnhcond_prof %d=%f\n",i,mesh.info.profiles[AC_dlnhcond_prof][i]);
  // fflush(stdout);
  // return;
  // make_tasks(diagnostics_func, reduction_func, finalize_func, write_func);
  // thread_pool.WaitAll();
  constexpr real alpha[3] = {0.0, -(5.0 / 9.0), -(153.0 / 128.0)};
  constexpr real beta[3] = {(1.0 / 3.0), (15.0 / 16.0), (8.0 / 15.0)};
  constexpr int num_of_steps = 100;

  AcMesh mesh_true;
  AcMesh mesh_test;
  const AcReal epsilon = pow(0.1,12);
  // constexpr AcReal epsilon = 0.0;
  constexpr AcReal local_dt = 0.001;
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
  for(int i = 0; i < 3; ++i)
  {
  	acDeviceSetInput(acGridGetDevice(), AC_step_num, i);
  	acGridExecuteTaskGraph(rhs,1);
  	acGridSynchronizeStream(STREAM_ALL);
  }

  acGridSynchronizeStream(STREAM_ALL);
  acDeviceLoadMesh(acGridGetDevice(), STREAM_DEFAULT, mesh_test);
  acGridSynchronizeStream(STREAM_ALL);
  acDeviceStoreMesh(acGridGetDevice(), STREAM_DEFAULT, &mesh);
  acGridSynchronizeStream(STREAM_ALL);

  bool loaded_correct = true;
  for (int i = 0; i < mx; i++)
  {
    for (int j = 0; j < my; j++)
    {
      for (int k = 0; k < mz; k++)
      {
        for (int ivar = 0; ivar < mfarray; ivar++)
        {
          AcReal out_val = mesh.vertex_buffer[ivar][DEVICE_VTXBUF_IDX(i, j, k)];
          AcReal true_val = mesh_test.vertex_buffer[ivar][DEVICE_VTXBUF_IDX(i, j, k)];
          if (out_val != true_val)
          {
            loaded_correct = false;
            acLogFromRootProc(rank,"Loaded val wrong at %d,%d,%d\n", i, j, k);
            acLogFromRootProc(rank,"field = %d", ivar);
            acLogFromRootProc(rank,"Loaded val: %f\tTRUE val: %f\n", out_val, true_val);
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
  acGridLaunchKernel(STREAM_DEFAULT, acGetKernelByName("AC_BUILTIN_RESET"), dims.n0, dims.n1);
  acGridSynchronizeStream(STREAM_ALL);

  //actual run
  for (int i=0;i<num_of_steps;i++){
    for(int i = 0; i < 3; ++i)
    {
    	acDeviceSetInput(acGridGetDevice(), AC_step_num, i);
    	acGridExecuteTaskGraph(rhs,1);
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
  AcReal true_pair;
  AcReal max_abs_relative_difference =-1.0;
  AcReal max_abs_value = -1.0;
  AcReal min_abs_value = 1.0;
  AcReal gpu_val_for_largest_diff;
  AcReal true_val_for_largest_diff;
  int num_of_points_where_different[NUM_VTXBUF_HANDLES] = {0};

  for (int i = dims.n0.x; i < dims.n1.x; i++)
  {
    for (int j = dims.n0.y; j < dims.n1.y; j++)
    {
      for (int k = dims.n0.z; k < dims.n1.z; k++)
      {
        for (int ivar = 0; ivar < NUM_VTXBUF_HANDLES; ivar++)
        {
          AcReal out_val = mesh.vertex_buffer[ivar][DEVICE_VTXBUF_IDX(i, j, k)];
          AcReal true_val = mesh_true.vertex_buffer[ivar][DEVICE_VTXBUF_IDX(i, j, k)];
          AcReal abs_diff = fabs(out_val - true_val);
          if (fabs(true_val) > max_abs_value) max_abs_value = fabs(true_val);
          if (fabs(true_val) < min_abs_value) min_abs_value = fabs(true_val);
          if ((abs_diff/true_val) > epsilon || (true_val == 0.0 && fabs(out_val) > pow(0.1,13)) || (epsilon == 0.0 && true_val != out_val))
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
            if (true_val != 0.0){
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
            acLogFromRootProc(rank,"%.7e\n",out_val);
          }
        }
      }
    }
  }
  auto volume_size = [](int3 a)
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
  acLogFromRootProc(rank,"max abs not passed val: %.7e\t%.7e\n",max_abs_not_passed_val, fabs(true_pair));
  acLogFromRootProc(rank,"max abs relative difference val: %.7e\n",max_abs_relative_difference);
  acLogFromRootProc(rank,"largest difference: %.7e\t%.7e\n",gpu_val_for_largest_diff, true_val_for_largest_diff);
  acLogFromRootProc(rank,"abs range: %.7e-%7e\n",min_abs_value,max_abs_value);
  fflush(stdout);
}
/***********************************************************************************************/
extern "C" void registerGPU(AcReal *farray)
{
  // AcReal* profile_x_host = (AcReal*)malloc(sizeof(AcReal)*mx);
  // for (int i=0;i<mx;i++){
  //     profile_x_host[i] = (AcReal)i;
  //     acLogFromRootProc(rank,"profile_x_host[%d]=%f\n",i,profile_x_host[i]);
  // }
  // mesh.profiles[PROFILE_X] = profile_x_host;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  size_t offset = 0;
  for (int i = 0; i < NUM_VTXBUF_HANDLES; ++i)
  {
    mesh.vertex_buffer[VertexBufferHandle(i)] = &farray[offset];
    //test_mesh.vertex_buffer[VertexBufferHandle(i)] = (AcReal *)malloc(sizeof(AcReal) * mw);
    offset += mw;
  }
}
/***********************************************************************************************/
extern "C" void initGPU()
{
  // Initialize GPUs in the node
  //AcResult res = acCheckDeviceAvailability();
}
/***********************************************************************************************/
#define PCLoad acPushToConfig
void setupConfig(AcMeshInfo& config, AcCompInfo& comp_info)
{ 
  // Enter basic parameters in config.

  PCLoad(config,comp_info, AC_domain_decomposition, (int3) {nprocx,nprocy,nprocz});
  PCLoad(config,comp_info, AC_nxgrid, nxgrid);
  PCLoad(config,comp_info, AC_nygrid, nygrid);
  PCLoad(config,comp_info, AC_nzgrid, nzgrid);
  PCLoad(config,comp_info,AC_decompose_strategy,(int)AcDecomposeStrategy::External);
  if (lmorton_curve)
    PCLoad(config,comp_info,AC_proc_mapping_strategy,(int)AcProcMappingStrategy::Morton);
  else
    PCLoad(config,comp_info,AC_proc_mapping_strategy,(int)AcProcMappingStrategy::Linear);
  PCLoad(config,comp_info,AC_MPI_comm_strategy,(int)AcMPICommStrategy::DuplicateUserComm);
  config.comm = comm_pencil;

// grid and geometry related parameters

  PCLoad(config,comp_info,AC_dsx,dx);
  PCLoad(config,comp_info,AC_dsy,dy);
  PCLoad(config,comp_info,AC_dsz,dz);

  PCLoad(config,comp_info,AC_xlen,  Lxyz[0]);
  PCLoad(config,comp_info,AC_ylen,  Lxyz[1]);
  PCLoad(config,comp_info,AC_zlen,  Lxyz[2]);
  PCLoad(config,comp_info,AC_xorig, xyz0[0]);
  PCLoad(config,comp_info,AC_yorig, xyz0[1]);
  PCLoad(config,comp_info,AC_zorig, xyz0[2]);
  PCLoad(config,comp_info,AC_x,     x);
  PCLoad(config,comp_info,AC_y,     y);
  PCLoad(config,comp_info,AC_z,     z);

// physics related parameters

  PCLoad(config,comp_info,AC_mu0,mu0);

// parameter arrays for boundary conditions

  PCLoad(config,comp_info,AC_fbcx  , fbcx);
  PCLoad(config,comp_info,AC_fbcx_2, fbcx_2);
  PCLoad(config,comp_info,AC_fbcy  , fbcy);
  PCLoad(config,comp_info,AC_fbcy_1, fbcy_1);
  PCLoad(config,comp_info,AC_fbcy_2, fbcy_2);
  PCLoad(config,comp_info,AC_fbcz  , fbcz);
  PCLoad(config,comp_info,AC_fbcz_1, fbcz_1);
  PCLoad(config,comp_info,AC_fbcz_2, fbcz_2);

  // Enter physics related parameters in config.

  #include "PC_modulepars.h"
  #if LDENSITY
    PCLoad(config,comp_info,AC_ldensity_nolog,ldensity_nolog);
    //printf("ldensity_nolog is %d \n",config.int_params[AC_ldensity_nolog]);//ldensity_nolog);
  #endif
  //Device dev = acGridGetDevice();
#if LHYDRO
  //dev->output.real_outputs[AC_maxadvec]=0.;
#endif
#if LENTROPY
  //dev->output.real_outputs[AC_maxchi]=0.;
#endif

}

int
get_int_param(const AcMeshInfo &config, AcCompInfo& comp_info, AcIntParam param)
{
	return config.int_params[(int)param];
}

int
get_int_param(const AcMeshInfo &config, AcCompInfo& comp_info, AcIntCompParam param)
{
	return comp_info.config.int_params[(int)param];
}

AcReal
get_real_param(const AcMeshInfo &config, AcCompInfo& comp_info, AcRealParam param)
{
	return config.real_params[(int)param];
}

AcReal
get_real_param(const AcMeshInfo &config, AcCompInfo& comp_info, AcRealCompParam param)
{
	return comp_info.config.real_params[(int)param];
}
/***********************************************************************************************/
void checkConfig(AcMeshInfo &config, AcCompInfo& comp_info)
{
 acLogFromRootProc(rank,"Check that config is correct\n");
 acLogFromRootProc(rank,"n[xyz]grid, d[xyz]: %d %d %d %.14f %.14f %.14f \n", nxgrid, nygrid, nzgrid, dx, dy, dz);
 acLogFromRootProc(rank,"rank= %d: l1, l2, n1, n2, m1, m2= %d %d %d %d %d %d \n", rank, l1, l2, n1, n2, m1, m2);
 acLogFromRootProc(rank,"rank= %d: zlen= %.14f %.14f \n", config.real_params[AC_zlen], Lxyz[2]);

#if LENTROPY
 acLogFromRootProc(rank,"lpressuregradientgas= %d %d \n", lpressuregradient_gas, get_int_param(config,comp_info,AC_lpressuregradient_gas));
 acLogFromRootProc(rank,"chi= %f %f \n", chi, get_real_param(config,comp_info,AC_chi));
 acLogFromRootProc(rank,"nkramers= %f %f \n", nkramers, get_real_param(config,comp_info,AC_nkramers));
 acLogFromRootProc(rank,"hcond0_kramers= %f %f \n", hcond0_kramers, get_real_param(config,comp_info,AC_hcond0_kramers));
 acLogFromRootProc(rank,"hcond_Kconst= %f %f \n", hcond_Kconst, get_real_param(config,comp_info,AC_hcond_Kconst));
#endif
#if LVISCOSITY
 acLogFromRootProc(rank,"nu= %f %f \n", nu, get_real_param(config,comp_info,AC_nu));
 acLogFromRootProc(rank,"zeta= %f %f \n", zeta, get_real_param(config,comp_info,AC_zeta));
#endif
#if LMAGNETIC
  acLogFromRootProc(rank,"eta= %f %f \n", eta, get_real_param(config,comp_info,AC_eta));
#endif
#if LEOS
  acLogFromRootProc(rank,"cs20= %f %f \n", cs20, get_real_param(config,comp_info,AC_cs20));
  //  acLogFromRootProc(rank,"gamma= %f %f \n", gamma, get_real_param(config,comp_infoAC_gamma));
  acLogFromRootProc(rank,"gamma_m1= %f %f \n", gamma_m1, get_real_param(config,comp_info,AC_gamma_m1));
  acLogFromRootProc(rank,"gamma1= %f %f \n", gamma1, get_real_param(config,comp_info,AC_gamma1));
  acLogFromRootProc(rank,"cv= %f %f \n", cv, get_real_param(config,comp_info,AC_cv));
  acLogFromRootProc(rank,"cp= %f %f \n", cp, get_real_param(config,comp_info,AC_cp));
  acLogFromRootProc(rank,"lnT0= %f %f \n", lnTT0, get_real_param(config,comp_info,AC_lnTT0));
  acLogFromRootProc(rank,"lnrho0= %f %f \n", lnrho0, get_real_param(config,comp_info,AC_lnrho0));
#endif
#if LFORCING
  acLogFromRootProc(rank,"iforcing_zsym= %f %f \n", iforcing_zsym, get_int_param(config,comp_info,AC_iforcing_zsym));
  acLogFromRootProc(rank,"k1_ff= %f %f \n", k1_ff, get_real_param(config,comp_info,AC_k1_ff));
  acLogFromRootProc(rank,"tforce_stop= %f %f \n", tforce_stop, get_real_param(config,comp_info,AC_tforce_stop));
  acLogFromRootProc(rank,"k1_ff,profx_ampl, val= %f %d %lf %lf\n", k1_ff, profx_ampl, profx_ampl[0], profx_ampl[nx-1]);
#endif
  acLogFromRootProc(rank,"mu0= %f %f \n", mu0, get_real_param(config,comp_info,AC_mu0));
}
/***********************************************************************************************/
extern "C" void getFArrayIn(AcReal **p_f_in)
{
  auto VBA = acGridGetVBA();
  *p_f_in = VBA.in[0];
}
/***********************************************************************************************/
extern "C" void copyVBApointers(AcReal **in, AcReal **out)
{
  Device device = acGridGetDevice();
  auto VBA = acGridGetVBA();
  *in =  VBA.in[0];
  *out = VBA.out[0];
}
/***********************************************************************************************/
extern "C" void initializeGPU(AcReal **farr_GPU_in, AcReal **farr_GPU_out, int comm_fint)
{
  //Setup configurations used for initializing and running the GPU code

#if PACKED_DATA_TRANSFERS
  //initLoadStore();
#endif

  comm_pencil = MPI_Comm_f2c(comm_fint);
  AcCompInfo comp_info = acInitCompInfo();
  setupConfig(mesh.info,comp_info);
#if AC_RUNTIME_COMPILATION
  if(rank == 0)
  {
#include "cmake_options.h"
	  acCompile(cmake_options,comp_info);
  }
  MPI_Barrier(comm_pencil);
  acLoadLibrary();
  acLogFromRootProc(rank, "Done setupConfig && acCompile\n");
  fflush(stdout);
#else
  acLogFromRootProc(rank, "Done setupConfig\n");
  fflush(stdout);
#endif
  checkConfig(mesh.info,comp_info);
  acCheckDeviceAvailability();
  acGridInit(mesh.info);
  rhs = acGetDSLTaskGraph(AC_rhs);
  acGridSynchronizeStream(STREAM_ALL);
  acLogFromRootProc(rank, "DONE initializeGPU\n");
  fflush(stdout);
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
    offset += mw;
  }
  mesh_to_copy.info = mesh.info;

  acGridSynchronizeStream(STREAM_ALL);
  //acDeviceStoreMesh(acGridGetDevice(), STREAM_DEFAULT, &mesh_to_copy);
  acDeviceStoreMesh(acGridGetDevice(), STREAM_DEFAULT, &mesh);
  acGridSynchronizeStream(STREAM_ALL);
}
/***********************************************************************************************/
extern "C" void reloadConfig()
{
  AcCompInfo comp_info = acInitCompInfo();
  setupConfig(mesh.info, comp_info);
  acGridSynchronizeStream(STREAM_ALL);
  acDeviceUpdate(acGridGetDevice(), mesh.info);
  acGridSynchronizeStream(STREAM_ALL);
}
/***********************************************************************************************/
extern "C" void loadFarray()
{
  /**
  if (has_nans(mesh)){
    acLogFromRootProc(rank,"found nans while copying\n");
    exit(0);
  }
  **/
  acGridSynchronizeStream(STREAM_ALL);
  acDeviceLoadMesh(acGridGetDevice(), STREAM_DEFAULT, mesh);
  acGridSynchronizeStream(STREAM_ALL);
}
/***********************************************************************************************/
extern "C" void updateInConfigArr(int index)
{
     if (mesh.info.real_arrays[index] != nullptr)
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
    int ind = -1;
    for (int i=0; i<NUM_REAL_ARRAYS; i++){
       if (strcmp(get_array_info(static_cast<AcRealArrayParam>(i)).name,name)==0) ind=i;
    }
    if (ind>-1) updateInConfigArr(ind);
    return ind;
}
/***********************************************************************************************/
extern "C" int updateInConfigScalName(char *name, AcReal value)
{
    int ind = -1;
    for (int i=0; i<NUM_REAL_PARAMS; i++){
       if (strcmp(realparam_names[i],name)==0) ind=i;
    }
    if (ind>-1) updateInConfigScal(ind, value);
    return ind;
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
