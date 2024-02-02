/*
                               --------------------

   Date:   8-Feb-2017
   Author: M. Rheinhardt & J. Pekkilae
   Description:
           ANSI C and standard library callable function wrappers for ASTAROTH-nucleus functions to be called from Fortran.
  Comments: 
*/

//General
//#include <cmath>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>

#define CUDA_ERRCHK(X)

#include "submodule/acc-runtime/api/errchk.h"
#include "submodule/acc-runtime/api/math_utils.h"
#include "submodule/include/astaroth.h"
#include "submodule/build/acc-runtime/api/user_defines.h"
#include "submodule/src/core/kernels/kernels.h"
#include "submodule/src/core/task.h"
#include "submodule/acc-runtime/api/math_utils.h"
#include "submodule/include/astaroth_utils.h"
#define real AcReal
#define EXTERN
#define FINT int
__thread int tp_int;
typedef void (*rangefunc)(const int a, const int b);
//PC interface
#include "PC_moduleflags.h"
#include "../cparam_c.h"
#include "../cdata_c.h"
#include "../sub_c.h"       // provides set_dt
#include "../boundcond_c.h" // provides boundconds[xyz] etc.
#include "../mpicomm_c.h"   // provides finalize_sendrcv_bdry
#include "./stencil_loader.h"
//#include "diagnostics/diagnostics.h"
#if PACKED_DATA_TRANSFERS
//#include "loadStore.h"
#endif
#if LFORCING
//   #include "forcing.h"
#endif

#include "PC_module_parfuncs.h"
static const AcReal DER1_3_E = 1.0;
static const AcReal DER1_2_E = -9.0;
static const AcReal DER1_1_E = 45.0;
static const AcReal DER1_E_DIV = 60.0;
static const AcReal AC_inv_dsx = 20.3718327157626;
static const AcReal AC_inv_dsy = 20.3718327157626;
static const AcReal AC_inv_dsz = 20.3718327157626;
static const bool symmetric_der = true;
static AcMesh mesh;
static AcMesh test_mesh;
static AcTaskGraph *graph_1;
static AcTaskGraph *graph_2;
static AcTaskGraph *graph_3;
static AcTaskGraph *randomize_graph;
static AcTaskGraph *rhs_test_graph;
static AcTaskGraph *rhs_test_graph_2;
static int pid;
static int rank;
DeviceConfiguration devConfig;
int halo_xz_size[2] = {0, 0}, halo_yz_size[2] = {0, 0};
static int l0 = 1;

static AcReal *xtop_buffer;
static AcReal *xbot_buffer;

static AcReal *ytop_buffer;
static AcReal *ybot_buffer;
#if LFORCING
// static ForcingParams forcing_params;
// #include "../forcing_c.h"
#endif
int DCONST(const AcIntParam param)
{
  return mesh.info.int_params[param];
}
int3 DCONST(const AcInt3Param param)
{
  return mesh.info.int3_params[param];
}
AcReal
DCONST(const AcRealParam param)
{
  return mesh.info.real_params[param];
}
AcReal3
DCONST(const AcReal3Param param)
{
  return mesh.info.real3_params[param];
}
int DEVICE_VTXBUF_IDX(const int x, const int y, const int z)
{
  return x + mx * y + mx * my * z;
}
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
  if(symmetric_der) return derx_symmetric(x,y,z,mesh,field);
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
  if(symmetric_der) return dery_symmetric(x,y,z,mesh,field);
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
  if(symmetric_der) return derz_symmetric(x,y,z,mesh,field);
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
  if(symmetric_der) return divergence_symmetric(x,y,z,mesh,field);
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
  if(symmetric_der) return gradient_symmetric(x,y,z,mesh,field);
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
  if(symmetric_der) return gradients_symmetric(x,y,z,mesh,field);
  return gradients_astaroth(x,y,z,mesh,field);
}
/***********************************************************************************************/
void print_diagnostics(const int pid, const int step, const AcReal dt, const AcReal simulation_time,
                       FILE *diag_file, const AcReal sink_mass, const AcReal accreted_mass,
                       int *found_nan)
{

  AcReal buf_rms, buf_max, buf_min;
  const int max_name_width = 16;

  // Calculate rms, min and max from the velocity vector field
  acGridReduceVec(STREAM_DEFAULT, RTYPE_MAX, VTXBUF_UUX, VTXBUF_UUY, VTXBUF_UUZ, &buf_max);
  acGridReduceVec(STREAM_DEFAULT, RTYPE_MIN, VTXBUF_UUX, VTXBUF_UUY, VTXBUF_UUZ, &buf_min);
  acGridReduceVec(STREAM_DEFAULT, RTYPE_RMS, VTXBUF_UUX, VTXBUF_UUY, VTXBUF_UUZ, &buf_rms);

  if (pid == 0)
    (pid, "Step %d, t_step %.3e, dt %e s\n", step, double(simulation_time),
     double(dt));
  if (pid == 0)
    (pid, "  %*s: min %.3e,\trms %.3e,\tmax %.3e\n", max_name_width, "uu total",
     double(buf_min), double(buf_rms), double(buf_max));
  if (pid == 0)
  {
    fprintf(diag_file, "%d %e %e %e %e %e ", step, double(simulation_time), double(dt),
            double(buf_min), double(buf_rms), double(buf_max));
  }

#if LBFIELD
  acGridReduceVec(STREAM_DEFAULT, RTYPE_MAX, BFIELDX, BFIELDY, BFIELDZ, &buf_max);
  acGridReduceVec(STREAM_DEFAULT, RTYPE_MIN, BFIELDX, BFIELDY, BFIELDZ, &buf_min);
  acGridReduceVec(STREAM_DEFAULT, RTYPE_RMS, BFIELDX, BFIELDY, BFIELDZ, &buf_rms);

  acLogFromRootProc(pid, "  %*s: min %.3e,\trms %.3e,\tmax %.3e\n", max_name_width, "bb total",
                    double(buf_min), double(buf_rms), double(buf_max));
  if (pid == 0)
  {
    fprintf(diag_file, "%e %e %e ", double(buf_min), double(buf_rms), double(buf_max));
  }

  acGridReduceVecScal(STREAM_DEFAULT, RTYPE_ALFVEN_MAX, BFIELDX, BFIELDY, BFIELDZ, VTXBUF_LNRHO,
                      &buf_max);
  acGridReduceVecScal(STREAM_DEFAULT, RTYPE_ALFVEN_MIN, BFIELDX, BFIELDY, BFIELDZ, VTXBUF_LNRHO,
                      &buf_min);
  acGridReduceVecScal(STREAM_DEFAULT, RTYPE_ALFVEN_RMS, BFIELDX, BFIELDY, BFIELDZ, VTXBUF_LNRHO,
                      &buf_rms);

  acLogFromRootProc(pid, "  %*s: min %.3e,\trms %.3e,\tmax %.3e\n", max_name_width, "vA total",
                    double(buf_min), double(buf_rms), double(buf_max));
  if (pid == 0)
  {
    fprintf(diag_file, "%e %e %e ", double(buf_min), double(buf_rms), double(buf_max));
  }
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

#if LSINK
  acLogFromRootProc(pid, "sink mass is: %.15e \n", double(sink_mass));
  acLogFromRootProc(pid, "accreted mass is: %.15e \n", double(accreted_mass));
#endif

  fflush(diag_file);
  fflush(stdout);
}
/***********************************************************************************************/
AcReal max_advec()
{
  AcReal maxadvec_ = 0.;
#if LHYDRO
  AcReal umax;
  acGridReduceVec(STREAM_DEFAULT, RTYPE_MAX, VTXBUF_UUX, VTXBUF_UUY, VTXBUF_UUZ, &umax);
#endif
  return maxadvec_;
}
/***********************************************************************************************/
AcReal max_diffus()
{
  AcReal dxyz_2_val = get_dxyz_2();
  AcReal maxdiffus_ = 0.;
#if LVISCOSITY
  maxdiffus_ = nu * dxyz_2_val;
#endif
#if LMAGNETIC
  maxdiffus_ = std::max(maxdiffus_, eta * dxyz_2_val);
        //maxdiffus_=nu*dxyz_2[nghost-1];
        //maxdiffus_=nu/(dx*dx);
#endif
#if LMAGNETIC
        //maxdiffus_=std::max(maxdiffus_,eta*dxyz_2[nghost-1]);
        //maxdiffus_=eta/(dx*dx);
#endif
  return maxdiffus_;
}
int id_to_tag(int3 id)
{
  return ((3 + id.x) % 3) * 9 + ((3 + id.y) % 3) * 3 + (3 + id.z) % 3;
}
int3 tag_to_id(int _tag)
{
  int3 _id = (int3){(_tag) / 9, ((_tag) % 9) / 3, (_tag) % 3};
  _id.x = _id.x == 2 ? -1 : _id.x;
  _id.y = _id.y == 2 ? -1 : _id.y;
  _id.z = _id.z == 2 ? -1 : _id.z;
  ERRCHK_ALWAYS(id_to_tag(_id) == _tag);
  return _id;
}

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
            printf("nan at %d,%d,%d\n", i, j, k);
            printf("field = %d", ivar);
          }
        }
      }
    }
  }
  return res;
}
/***********************************************************************************************/
//Do the 'isubstep'th integration step on all GPUs on the node and handle boundaries.
//
extern "C" void substepGPU(int isubstep, bool full = false, bool early_finalize = false)
{
#if LFORCING
  // //Update forcing params

  // if (isubstep == itorder)
  //      forcing_params.Update();  // calculate on CPU and load into GPU
#endif
  if (lfirst && ldt)
  {
    AcReal dt1_advec = max_advec() / cdt;
    AcReal dt1_diffus = max_diffus() / cdtv;
    AcReal dt1_ = sqrt(pow(dt1_advec, 2) + pow(dt1_diffus, 2));
    set_dt(dt1_);
  }
  acGridLoadScalarUniform(STREAM_DEFAULT, AC_dt, dt);
  acGridSynchronizeStream(STREAM_DEFAULT);
  //Transfer the updated ghost zone to the device(s) in the node

  if (full)
  {
    acGridSynchronizeStream(STREAM_DEFAULT);
    acDeviceLoadMesh(acGridGetDevice(), STREAM_DEFAULT, mesh);
    acGridSynchronizeStream(STREAM_ALL);
    //set output buffer to 0 since if we are reading from it we don't want NaNs
    AcMeshDims dims = acGetMeshDims(acGridGetLocalMeshInfo());
    acGridLaunchKernel(STREAM_DEFAULT, reset, dims.n0,dims.n1);
    acGridSynchronizeStream(STREAM_ALL);
  }
  acGridSynchronizeStream(STREAM_ALL);
  acGridLoadScalarUniform(STREAM_DEFAULT, AC_dt, dt);
  if (isubstep == 1)
  {
    acGridExecuteTaskGraph(graph_1, 1);
  }
  if (isubstep == 2)
  {
    acGridExecuteTaskGraph(graph_2, 1);
  }
  if (isubstep == 3)
  {
    acGridExecuteTaskGraph(graph_3, 1);
  }
  acGridSynchronizeStream(STREAM_ALL);
  // printf("Done substep: %d\n",isubstep);
  // fflush(stdout);
  // int found_nan;
  // FILE* diag_file = fopen("astaroth_timeseries.ts", "a");
  // int rank;
  // print_diagnostics(rank, 1, dt, dt,diag_file, 0.0001, 0.0001,
  //               &found_nan);
  return;
}
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
  //Run the gpu code serially
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
        // for(int z=1;z<=3;z++){
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
            printf("C val wrong at %d,%d,%d\n", i, j, k);
            printf("field = %d", ivar);
            printf("C val: %f\tF val: %f\n", out_val, true_val);
          }
        }
      }
    }
  }
  if (passed)
  {
    printf("Passed C test :)\n");
  }
  else
  {
    printf("Did not pass C test :(\n");
  }
  if (!passed)
    return;
  printf("Starting GPUtest\n");
  fflush(stdout);
  acGridSynchronizeStream(STREAM_ALL);
  // acGridLoadMesh(STREAM_DEFAULT,mesh);
  // printf("loaded mesh\n");
  // acGridTestBCKernel({mx,my,1});

  acGridSynchronizeStream(STREAM_ALL);
  printf("after bc kernel\n");
  fflush(stdout);
  acGridStoreMesh(STREAM_DEFAULT, &mesh);
  acGridSynchronizeStream(STREAM_ALL);
  printf("after store\n");
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
            printf("GPU val wrong at %d,%d,%d\n", i, j, k);
            printf("field = %d", ivar);
            printf("GPU val: %f\tTRUE val: %f\tDIFF: %f\n", out_val, true_val, fabs(out_val - true_val));
            if(fabs(out_val)>max_abs_not_passed)
              max_abs_not_passed = out_val;
          }
        }
      }
    }
  }
  printf("maximum abs incorrect val: %f\n",max_abs_not_passed);
  if (passed)
  {
    printf("Passed GPU test :)\n");
  }
  else
  {
    printf("Did not pass GPU test :(\n");
  }
  return;
}
extern "C" void testRHS(AcReal *farray_in, AcReal *dfarray_truth)
{
  // __energy_MOD_pushpars2c(p_par_energy);
  // mesh.info.profiles[AC_dlnhcond_prof]=dlnhcond_prof; // [2-1] dlnhcond_prof real(nz)
  // for(int i=20;i<30;++i)
  //   printf("C: dlnhcond_prof %d=%f\n",i,mesh.info.profiles[AC_dlnhcond_prof][i]);
  // fflush(stdout);
  // return;
  // make_tasks(diagnostics_func, reduction_func, finalize_func, write_func);
  // thread_pool.WaitAll();
  constexpr real alpha[3] = {0.0, -(5.0 / 9.0), -(153.0 / 128.0)};
  constexpr real beta[3] = {(1.0 / 3.0), (15.0 / 16.0), (8.0 / 15.0)};
  constexpr int num_of_steps = 1;
  printf("HI from testRHS\n");
  fflush(stdout);
  AcMesh mesh_true;
  AcMesh mesh_test;
  AcMesh df_mesh;
  AcMesh ds_mesh;
  AcMesh next_mesh;
  for (int i = 0; i < NUM_VTXBUF_HANDLES; ++i)
  {
    ds_mesh.vertex_buffer[VertexBufferHandle(i)] = (AcReal *)malloc(sizeof(AcReal) * mw);
    df_mesh.vertex_buffer[VertexBufferHandle(i)] = (AcReal *)malloc(sizeof(AcReal) * mw);
    next_mesh.vertex_buffer[VertexBufferHandle(i)] = (AcReal *)malloc(sizeof(AcReal) * mw);
  }
  const AcReal epsilon = pow(0.1,12);
  // constexpr AcReal epsilon = 0.0;
  constexpr AcReal local_dt = 0.001;
  acGridLoadScalarUniform(STREAM_DEFAULT, AC_dt, local_dt);
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
  printf("n0: %d,%d,%d\tpid: %d\n", dims.n0.x, dims.n0.y, dims.n0.z, pid);
  printf("n1: %d,%d,%d\tpid: %d\n", dims.n1.x, dims.n1.y, dims.n1.z, pid);

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

    acGridLaunchKernel(STREAM_DEFAULT, twopass_solve_final_step1, dims.n0,dims.n1);
    acGridSynchronizeStream(STREAM_ALL);
    acGridSwapBuffers();
    acGridSynchronizeStream(STREAM_ALL);

    acGridLaunchKernel(STREAM_DEFAULT, twopass_solve_intermediate_step2, dims.n0,dims.n1);
    acGridSynchronizeStream(STREAM_ALL);
    acGridSwapBuffers();
    acGridSynchronizeStream(STREAM_ALL);

    acGridLaunchKernel(STREAM_DEFAULT, twopass_solve_final_step2, dims.n0,dims.n1);
    acGridSynchronizeStream(STREAM_ALL);
    acGridSwapBuffers();
    acGridSynchronizeStream(STREAM_ALL);

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
            printf("Loaded val wrong at %d,%d,%d\n", i, j, k);
            printf("field = %d", ivar);
            printf("Loaded val: %f\tTRUE val: %f\n", out_val, true_val);
          }
        }
      }
    }
  }
  AcReal derx_ux;
  AcReal derx_normal;
  AcReal dery_uy;
  AcReal derz_uz;
  //test calculating divu different ways
  const int y = 42;
  const int x = 42;
  const int z = 42;
  const AcReal AC_inv_dsx = 2.0;

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
  // printf("test derx_ux: %.7e\n",derx_ux);
  // printf("normal derx_ux. %.7e\n",derx_normal);
  fflush(stdout);
  // return;
  if (loaded_correct)
  {
    printf("loaded correct data\n");
  }
  else
  {
    printf("loaded incorrect data :(\n");
  }


  //set output buffer to 0 since if we are reading from it we don't want NaNs
  acGridLaunchKernel(STREAM_DEFAULT, reset, dims.n0,dims.n1);
  acGridSynchronizeStream(STREAM_ALL);
    // acGridExecuteTaskGraph(rhs_test_graph, 1);

    // acGridSynchronizeStream(STREAM_ALL);

    // acGridExecuteTaskGraph(rhs_test_graph_2, 1);

    // acGridSynchronizeStream(STREAM_ALL);


  //actual run
  for(int i=0;i<num_of_steps;i++){
    acGridExecuteTaskGraph(rhs_test_graph,1);
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
  }
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
          if(fabs(true_val)>max_abs_value)
            max_abs_value = fabs(true_val);
          if(fabs(true_val) < min_abs_value)
            min_abs_value = fabs(true_val);
          if((abs_diff/true_val) > epsilon || true_val == 0.0 && fabs(out_val) > pow(0.1,13) || epsilon == 0.0 && true_val != out_val)
          {
            passed = false;
            num_of_points_where_different[ivar]++;
            // printf("rhs val wrong at %d,%d,%d\n", i, j, k);
            // printf("field = %d", ivar);
            // printf("GPU val: %.7e\tTRUE val: %.7e\n", out_val, true_val);
            // printf("PID: %d\n", pid);
            if(max_abs_not_passed_val<abs(out_val)){
              max_abs_not_passed_val = abs(out_val);
              true_pair = true_val;
            }
            if(true_val != 0.0){
              if(max_abs_relative_difference<(abs_diff/true_val)){
                max_abs_relative_difference=(abs_diff/true_val);
                gpu_val_for_largest_diff = out_val;
                true_val_for_largest_diff = true_val;
              }
            }  
          }
          if(isnan(out_val))
          {
            printf("TP: nan before at %d,%d,%d,%d!\n!",i,j,k,ivar);
            printf("%.7e\n",out_val);
          }
        }
      }
    }
  }
  for(int ivar=0;ivar<NUM_VTXBUF_HANDLES;ivar++)
    printf("ratio of values wrong for field: %d\t %f\n",ivar,(double)num_of_points_where_different[ivar]/volume_size(dims.n1-dims.n0));
  passed &= !has_nans(mesh);
  if (passed)
  {
    printf("Passed GPU test :)\n");
  }
  else
  {
    printf("Did not pass GPU test :(\n");
  }
  printf("max abs not passed val: %.7e\t%.7e\n",max_abs_not_passed_val, fabs(true_pair));
  printf("max abs relative difference val: %.7e\n",max_abs_relative_difference);
  printf("largest difference: %.7e\t%.7e\n",gpu_val_for_largest_diff, true_val_for_largest_diff);
  printf("abs range: %.7e-%7e\n",min_abs_value,max_abs_value);
  fflush(stdout);
}
/***********************************************************************************************/
extern "C" void registerGPU(AcReal *farray)
{

  // AcReal* profile_x_host = (AcReal*)malloc(sizeof(AcReal)*mx);
  // for(int i=0;i<mx;i++){
  //     profile_x_host[i] = (AcReal)i;
  //     printf("profile_x_host[%d]=%f\n",i,profile_x_host[i]);
  // }
  // mesh.profiles[PROFILE_X] = profile_x_host;

  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  size_t offset = 0;
  for (int i = 0; i < NUM_VTXBUF_HANDLES; ++i)
  {
    mesh.vertex_buffer[VertexBufferHandle(i)] = &farray[offset];
    test_mesh.vertex_buffer[VertexBufferHandle(i)] = (AcReal *)malloc(sizeof(AcReal) * mw);
    offset += mw;
  }
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
}
/***********************************************************************************************/
extern "C" void initGPU()
{
  // Initialize GPUs in the node
  AcResult res = acCheckDeviceAvailability();
}
/***********************************************************************************************/
void setupConfig(AcMeshInfo &config)
{

  printf("nx etc. %d %d %d %.14f %.14f %.14f \n", nxgrid, nygrid, nzgrid, dx, dy, dz);
  config.int_params[AC_nx] = nxgrid;
  config.int_params[AC_ny] = nygrid;
  config.int_params[AC_nz] = nzgrid;
  config.int_params[AC_mx] = mxgrid;
  config.int_params[AC_my] = mygrid;
  config.int_params[AC_mz] = mzgrid;
  config.int_params[AC_nx_min] = l1;
  config.int_params[AC_nx_max] = l2;
  config.int_params[AC_ny_min] = m1;
  config.int_params[AC_ny_max] = m2;
  config.int_params[AC_nz_min] = n1;
  config.int_params[AC_nz_max] = n2;
  config.real_params[AC_dsx] = dx;
  config.real_params[AC_dsy] = dy;
  config.real_params[AC_dsz] = dz;
  printf("%d: l1i etc. %d %d %d %d %d %d \n", pid, l1i, l2i, n1i, n2i, m1i, m2i);
  printf("%d: l1 etc. %d %d %d %d %d %d \n", pid, l1, l2, n1, n2, m1, m2);
  config.real_params[AC_dsmin] = std::min(dx, std::min(dy, dz));
  config.real_params[AC_xlen] = 6.28318548;
  config.real_params[AC_ylen] = 6.28318548;
  config.real_params[AC_zlen] = 6.28318548;
  config.real_params[AC_xorig] = -3.14159274;
  config.real_params[AC_yorig] = -3.14159274;
  config.real_params[AC_zorig] = -3.14159274;
  config.real_params[AC_mu0] = mu0;

#include "PC_modulepars.h"
  printf("Done setUpConfig\n");
  fflush(stdout);
}
/***********************************************************************************************/
void checkConfig(AcMeshInfo &config)
{
  printf("check that config is correct\n");
//printf("setupConfig:mesh.info.real_params[AC_k1_ff]= %f \n",config.real_params[AC_k1_ff]);
#if LENTROPY
  printf("lpressuregradientgas= %d %d \n", lpressuregradient_gas, config.int_params[AC_lpressuregradient_gas]);
#endif
#if LENTROPY
  printf("chi= %f %f \n", chi, config.real_params[AC_chi]);
#endif
#if LVISCOSITY
  printf("nu= %f %f \n", nu, config.real_params[AC_nu]);
  printf("zeta= %f %f \n", zeta, config.real_params[AC_zeta]);
#endif
#if LMAGNETIC
  printf("eta= %f %f \n", eta, config.real_params[AC_eta]);
#endif
#if LEOS
  printf("cs20= %f %f \n", cs20, config.real_params[AC_cs20]);
  //  printf("gamma= %f %f \n", gamma, config.real_params[AC_gamma]);
  printf("gamma_m1= %f %f \n", gamma_m1, config.real_params[AC_gamma_m1]);
  printf("gamma1= %f %f \n", gamma1, config.real_params[AC_gamma1]);
  printf("cv= %f %f \n", cv, config.real_params[AC_cv]);
  printf("cp= %f %f \n", cp, config.real_params[AC_cp]);
  printf("lnT0= %f %f \n", lnTT0, config.real_params[AC_lnTT0]);
  printf("lnrho0= %f %f \n", lnrho0, config.real_params[AC_lnrho0]);
#endif
#if LFORCING
  //  printf("iforcing_zsym= %f %f \n", iforcing_zsym, config.int_params[AC_iforcing_zsym]);
  //  printf("k1_ff= %f %f \n", k1_ff, config.real_params[AC_k1_ff]);
  //  printf("tforce_stop= %f %f \n", tforce_stop, config.real_params[AC_tforce_stop]);
  //printf("k1_ff,profx_ampl, val= %f %d %lf %lf\n", k1_ff, profx_ampl, profx_ampl[0], profx_ampl[nx-1]);
#endif
  printf("mu0= %f %f \n", mu0, config.real_params[AC_mu0]);
}
/***********************************************************************************************/

void loadProfiles(AcMeshInfo &config)
{
  // #if LFORCING
  //      PUT(profx_ampl,nx,0,0)
  //      PUT(profy_ampl,0,my,0)
  //      PUT(profz_ampl,0,0,mz)
  //      PUT(profx_hel,nx,0,0)
  //      PUT(profy_hel,0,my,0)
  //      PUT(profz_hel,0,0,mz)
  // #endif
  // acGridLoadProfile(STREAM_DEFAULT, PROFILE_X, config);
}
/***********************************************************************************************/
extern "C" void copyVBApointers(AcReal **in, AcReal **out)
{
  Device device = acGridGetDevice();
  *in = device->vba.in[0];
  *out = device->vba.out[0];
}

extern "C" void getFArrayIn(AcReal **p_f_in)
{
  Device device = acGridGetDevice();
  *p_f_in = device->vba.in[0];
}
extern "C" void initializeGPU(AcReal **farr_GPU_in, AcReal **farr_GPU_out)
{
  //Setup configurations used for initializing and running the GPU code
#if PACKED_DATA_TRANSFERS
  // initLoadStore();
#endif
  setupConfig(mesh.info);

  printf("Setting domain int3\n");
  int3 decomp = {nprocx, nprocy, nprocz};
  fflush(stdout);
  acCheckDeviceAvailability();
  printf("Setted domain int3\n");
  fflush(stdout);

  mesh.info.int3_params[AC_domain_decomposition] = decomp;
  checkConfig(mesh.info);
  fflush(stdout);
  //TODO make cleaner.
  //Needed since based on the sample used some profiles might not have been 
  //allocated
  acGridInit(mesh.info);

  VertexBufferHandle all_fields[NUM_VTXBUF_HANDLES];
  for (int i = 0; i < NUM_VTXBUF_HANDLES; i++)
  {
    all_fields[i] = (VertexBufferHandle)i;
  }

  randomize_graph = acGridBuildTaskGraph({

      acHaloExchange(all_fields),
      acBoundaryCondition(BOUNDARY_XYZ, BOUNDCOND_PERIODIC, all_fields),
      acCompute(KERNEL_randomize, all_fields),
  });
  rhs_test_graph = acGridBuildTaskGraph({
      acHaloExchange(all_fields),
      acBoundaryCondition(BOUNDARY_XYZ, BOUNDCOND_PERIODIC, all_fields),
      acCompute(KERNEL_twopass_solve_intermediate_step0, all_fields),
      acCompute(KERNEL_twopass_solve_final_step0, all_fields),

      acHaloExchange(all_fields),
      acBoundaryCondition(BOUNDARY_XYZ, BOUNDCOND_PERIODIC, all_fields),
      acCompute(KERNEL_twopass_solve_intermediate_step1, all_fields),
      acCompute(KERNEL_twopass_solve_final_step1, all_fields),

      acHaloExchange(all_fields),
      acBoundaryCondition(BOUNDARY_XYZ, BOUNDCOND_PERIODIC, all_fields),
      acCompute(KERNEL_twopass_solve_intermediate_step2, all_fields),
      acCompute(KERNEL_twopass_solve_final_step2, all_fields), //
  });
  rhs_test_graph_2 = acGridBuildTaskGraph({
      acCompute(KERNEL_twopass_solve_intermediate_step1, all_fields),
      acCompute(KERNEL_twopass_solve_final_step1, all_fields),
  });
  
  graph_1 = acGridBuildTaskGraph({
    acHaloExchange(all_fields),
    acBoundaryCondition(BOUNDARY_XYZ, BOUNDCOND_PERIODIC, all_fields),
    acCompute(KERNEL_twopass_solve_intermediate_step0, all_fields),
    acCompute(KERNEL_twopass_solve_final_step0, all_fields),
  });
  graph_2 = acGridBuildTaskGraph(
    {
        acHaloExchange(all_fields),
        acBoundaryCondition(BOUNDARY_XYZ, BOUNDCOND_PERIODIC, all_fields),
        acCompute(KERNEL_twopass_solve_intermediate_step1, all_fields),
        acCompute(KERNEL_twopass_solve_final_step1, all_fields),
    });
  graph_3 = acGridBuildTaskGraph(
    {
        acHaloExchange(all_fields),
        acBoundaryCondition(BOUNDARY_XYZ, BOUNDCOND_PERIODIC, all_fields),
        acCompute(KERNEL_twopass_solve_intermediate_step2, all_fields),
        acCompute(KERNEL_twopass_solve_final_step2, all_fields),
    });
  printf("BUILD graphs\n");
  acGridExecuteTaskGraph(rhs_test_graph,1);
  acGridExecuteTaskGraph(graph_1,1);
  acGridExecuteTaskGraph(graph_2,1);
  acGridExecuteTaskGraph(graph_3,1);
  load_stencil_from_config(mesh.info);
}
/***********************************************************************************************/
extern "C" void copyFarray()
{
  //    AcResult res=acGridStoreMesh(STREAM_DEFAULT, &mesh);
  // acGridGetDevice();
  // AcReal* vba_in = acGridGetVBApointerIn();
  // acGridSynchronizeStream(STREAM_ALL);
  // // hipMemcpy(&mesh.vertex_buffer[0],vba_in,mw*NUM_VTXBUF_HANDLES*sizeof(AcReal),hipMemcpyDeviceToHost);
  acGridSynchronizeStream(STREAM_ALL);
  acDeviceStoreMesh(acGridGetDevice(), STREAM_DEFAULT, &mesh);
  // acDeviceStoreMesh(acGridGetDevice(), STREAM_DEFAULT, &mesh);
  acGridSynchronizeStream(STREAM_ALL);
  // printf("store all %d \n",res); fflush(stdout);
}
/***********************************************************************************************/
extern "C" void finalizeGPU()
{
#if PACKED_DATA_TRANSFERS
  acDeviceStoreMesh(acGridGetDevice(), STREAM_DEFAULT, &mesh);
#endif
  // Deallocate everything on the GPUs and reset
  AcResult res = acGridQuit();
}
extern "C" void random_initial_condition()
{
  acGridSynchronizeStream(STREAM_ALL);
  acGridExecuteTaskGraph(randomize_graph, 1);
  acGridSynchronizeStream(STREAM_ALL);
}
/***********************************************************************************************/
