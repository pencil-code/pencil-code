#include <thrust/complex.h>
#include "rk3_cuda_generic.cuh"
#include "diff_cuda_generic.cuh"

#include "gpu/cuda/core/errorhandler_cuda.cuh"
#ifdef GPU_ASTAROTH
#include "common/PC_moduleflags.h"
#endif

typedef struct {
    real *s_lnrho, *s_uux, *s_uuy, *s_uuz;
    real *r_lnrho, *r_uux, *r_uuy, *r_uuz;
} HydroStencil;

typedef struct {
    real *s_Ax, *s_Ay, *s_Az;
    real *r_Ax, *r_Ay, *r_Az;
} InductionStencil;

static __device__ __inline__ void load_halos(int smem_idx_base, int grid_idx_base, real* __restrict__ s_scal, const real* __restrict__ d_scal)
{
    if (threadIdx.x < BOUND_SIZE) {
        //Load left        
        {
            const int smem_idx = smem_idx_base - BOUND_SIZE;
            const int grid_idx = grid_idx_base - BOUND_SIZE;
            s_scal[smem_idx] = d_scal[grid_idx];
        }
    
        //Load right
        {
            const int smem_idx = smem_idx_base + blockDim.x;
            const int grid_idx = grid_idx_base + blockDim.x;
            s_scal[smem_idx] = d_scal[grid_idx];
        }
    }

    if (threadIdx.y < BOUND_SIZE) {      
        //Load bottom        
        {
            const int smem_idx = smem_idx_base - BOUND_SIZE*SMEM_WIDTH;
            const int grid_idx = grid_idx_base - BOUND_SIZE*d_mx;
            s_scal[smem_idx] = d_scal[grid_idx];
        }
    
        //Load top
        {
            const int smem_idx = smem_idx_base + RK_THREADS_Y*SMEM_WIDTH;
            const int grid_idx = grid_idx_base + RK_THREADS_Y*d_mx;
            s_scal[smem_idx] = d_scal[grid_idx];
        }
    }


    if (threadIdx.x < BOUND_SIZE && threadIdx.y < BOUND_SIZE) {
        //Load bottom left   
        {
            const int smem_idx = smem_idx_base - BOUND_SIZE - BOUND_SIZE*SMEM_WIDTH;
            const int grid_idx = grid_idx_base - BOUND_SIZE - BOUND_SIZE*d_mx;
            s_scal[smem_idx] = d_scal[grid_idx];
        }
        //Load bottom right   
        {
            const int smem_idx = smem_idx_base + blockDim.x - BOUND_SIZE*SMEM_WIDTH;
            const int grid_idx = grid_idx_base + blockDim.x - BOUND_SIZE*d_mx;
            s_scal[smem_idx] = d_scal[grid_idx];
        }
        //Load top left   
        {
            const int smem_idx = smem_idx_base - BOUND_SIZE + RK_THREADS_Y*SMEM_WIDTH;
            const int grid_idx = grid_idx_base - BOUND_SIZE + RK_THREADS_Y*d_mx;
            s_scal[smem_idx] = d_scal[grid_idx];
        } 
        //Load top right   
        {
            const int smem_idx = smem_idx_base + blockDim.x + RK_THREADS_Y*SMEM_WIDTH;
            const int grid_idx = grid_idx_base + blockDim.x + RK_THREADS_Y*d_mx;
            s_scal[smem_idx] = d_scal[grid_idx];
        }
    }
}

static __device__ __inline__ real laplace(const int smem_idx, real* s_scal, real* r_scal)
{
    const real d2x_scal = der2_scalx(smem_idx, s_scal);
    const real d2y_scal = der2_scaly(smem_idx, s_scal);
    const real d2z_scal = der2_scalz(r_scal);
    return d2x_scal + d2y_scal + d2z_scal;
}

#ifdef HYDRO
static __device__ real continuity(const int smem_idx, HydroStencil& stncl)
{
    const real ddx_lnrho = der_scalx(smem_idx, stncl.s_lnrho);
    const real ddx_uux   = der_scalx(smem_idx, stncl.s_uux);

    const real ddy_lnrho = der_scaly(smem_idx, stncl.s_lnrho);
    const real ddy_uuy   = der_scaly(smem_idx, stncl.s_uuy);

    const real ddz_lnrho = der_scalz(stncl.r_lnrho);
    const real ddz_uuz   = der_scalz(stncl.r_uuz);

    //Continuity
    const real res = - stncl.r_uux[3] * ddx_lnrho 
                     - stncl.r_uuy[3] * ddy_lnrho
		             - stncl.r_uuz[3] * ddz_lnrho
		             - ddx_uux - ddy_uuy - ddz_uuz;
    return res;
}

typedef enum {X_AXIS, Y_AXIS, Z_AXIS} AXIS;
template <AXIS axis>
static __device__ real momentum(const int smem_idx, HydroStencil& stncl)
{
    const real ddx_lnrho = der_scalx(smem_idx, stncl.s_lnrho);
    const real ddx_uux   = der_scalx(smem_idx, stncl.s_uux);
    const real ddx_uuy   = der_scalx(smem_idx, stncl.s_uuy);
    const real ddx_uuz   = der_scalx(smem_idx, stncl.s_uuz);

    const real ddy_lnrho = der_scaly(smem_idx, stncl.s_lnrho);
    const real ddy_uux   = der_scaly(smem_idx, stncl.s_uux);
    const real ddy_uuy   = der_scaly(smem_idx, stncl.s_uuy);
    const real ddy_uuz   = der_scaly(smem_idx, stncl.s_uuz);

    const real ddz_lnrho = der_scalz(stncl.r_lnrho);
    const real ddz_uux   = der_scalz(stncl.r_uux); 
    const real ddz_uuy   = der_scalz(stncl.r_uuy); 
    const real ddz_uuz   = der_scalz(stncl.r_uuz); 

    //S_grad_lnrho  //Eq(.9)
    const real divu3 = (1.0/3.0)*(ddx_uux + ddy_uuy + ddz_uuz);
    //const real Sxx = (2.0/3.0)*ddx_uux - (1.0/3.0)*(ddy_uuy + ddz_uuz);
    const real Sxx = ddx_uux - divu3;
    const real Sxy = 0.5*(ddy_uux + ddx_uuy);
    const real Sxz = 0.5*(ddz_uux + ddx_uuz);
    //const real Syy = (2.0/3.0)*ddy_uuy - (1.0/3.0)*(ddx_uux + ddz_uuz);
    const real Syy = ddy_uuy - divu3;
    const real Syz = 0.5*(ddz_uuy + ddy_uuz);
    //const real Szz = (2.0/3.0)*ddz_uuz - (1.0/3.0)*(ddx_uux + ddy_uuy);
    const real Szz = ddz_uuz - divu3;

    if (axis == X_AXIS) {
        const real d2x_uux = der2_scalx(smem_idx, stncl.s_uux);
        const real nu_const_uux = laplace(smem_idx, stncl.s_uux, stncl.r_uux);

        const real d2xy_uuy = der2_scalxy(smem_idx, stncl.s_uuy);

        const real res =   - stncl.r_uux[3] * ddx_uux //vec_dot_nabla_scal
                           - stncl.r_uuy[3] * ddy_uux
                           - stncl.r_uuz[3] * ddz_uux
                           - d_CS20*ddx_lnrho //ddx part of grad lnrho
                           + d_NU * nu_const_uux //nu_const
                           + 2.0*d_NU*(Sxx*ddx_lnrho + Sxy*ddy_lnrho + Sxz*ddz_lnrho)
                           + d_NU*(1.0/3.0)*(d2x_uux + d2xy_uuy); //S_grad_lnrho
        return res;
    } else if (axis == Y_AXIS) {
        const real d2y_uuy = der2_scaly(smem_idx, stncl.s_uuy);
        const real nu_const_uuy =  laplace(smem_idx, stncl.s_uuy, stncl.r_uuy);

        const real d2xy_uux = der2_scalxy(smem_idx, stncl.s_uux);

        const real res =   - stncl.r_uux[3] * ddx_uuy //vec_dot_nabla_scal
                           - stncl.r_uuy[3] * ddy_uuy
                           - stncl.r_uuz[3] * ddz_uuy
                           - d_CS20*ddy_lnrho //ddx part of grad lnrho
                           + d_NU * nu_const_uuy //nu_const
                           + 2.0*d_NU*(Sxy*ddx_lnrho + Syy*ddy_lnrho + Syz*ddz_lnrho)
                           + d_NU*(1.0/3.0)*(d2xy_uux + d2y_uuy); //S_grad_lnrho 
        return res;
    } else {
        const real d2z_uuz = der2_scalz(stncl.r_uuz);
        const real nu_const_uuz =  laplace(smem_idx, stncl.s_uuz, stncl.r_uuz);

        const real res =   - stncl.r_uux[3] * ddx_uuz //vec_dot_nabla_scal
                           - stncl.r_uuy[3] * ddy_uuz
                           - stncl.r_uuz[3] * ddz_uuz
                           - d_CS20*ddz_lnrho //ddx part of grad lnrho
                           + d_NU * nu_const_uuz //nu_const
                           + 2.0*d_NU*(Sxz*ddx_lnrho + Syz*ddy_lnrho + Szz*ddz_lnrho)
                           + d_NU*(1.0/3.0)*(d2z_uuz);
        return res;
    }
}
#endif

#ifdef MAGNETIC
template <AXIS axis>
static __device__ real induction(const int smem_idx, InductionStencil& stncl, const real uux, const real uuy, const real uuz)
{
    const real ddx_Az = der_scalx(smem_idx, stncl.s_Az);
    const real ddx_Ay = der_scalx(smem_idx, stncl.s_Ay);
    const real ddy_Ax = der_scaly(smem_idx, stncl.s_Ax);
    const real ddy_Az = der_scaly(smem_idx, stncl.s_Az);
    const real ddz_Ay = der_scalz(stncl.r_Ay);
    const real ddz_Ax = der_scalz(stncl.r_Ax);

    const real Bx = ddy_Az - ddz_Ay;
    const real By = ddz_Ax - ddx_Az;
    const real Bz = ddx_Ay - ddy_Ax;
    
    if (axis == X_AXIS) {
        const real laplace_Ax =  laplace(smem_idx, stncl.s_Ax, stncl.r_Ax);
        const real d2x_Ax = der2_scalx(smem_idx, stncl.s_Ax);
        const real d2xy_Ay = der2_scalxy(smem_idx, stncl.s_Ay);
        const real part_grad_div_Ax = d2x_Ax + d2xy_Ay;// + d2xz_Az;

        const real res = uuy*Bz - uuz*By - d_ETA * (-laplace_Ax + part_grad_div_Ax);
        return res;
    } else if (axis == Y_AXIS) {
        const real laplace_Ay = laplace(smem_idx, stncl.s_Ay, stncl.r_Ay);

        const real d2xy_Ax = der2_scalxy(smem_idx, stncl.s_Ax);
        const real d2y_Ay = der2_scaly(smem_idx, stncl.s_Ay);
        const real part_grad_div_Ay = d2xy_Ax + d2y_Ay;

        const real res = uuz*Bx - uux*Bz - d_ETA * (-laplace_Ay + part_grad_div_Ay);
        return res;
    } else {
        const real laplace_Az = laplace(smem_idx, stncl.s_Az, stncl.r_Az);
        const real d2z_Az = der2_scalz(stncl.r_Az);
        const real part_grad_div_Az = d2z_Az;// + d2xz_Ax;

        const real res = uux*By - uuy*Bx - d_ETA * (-laplace_Az + part_grad_div_Az);
        return res;
    }
}
#endif

//Nonhelical forcing adapted from astaroth_legacy
#ifdef FORCING
#ifdef GPU_ASTAROTH
template <AXIS axis>
static __device__ __inline__ real forcing(const int tx, const int ty, const int tz)
{
	vec3r profyz_hel_coef2;
        profyz_hel_coef2.x = d_PROFY_HEL[ty]*d_PROFZ_HEL[tz]*d_FORCING_COEF2[0];
        profyz_hel_coef2.y = d_PROFY_HEL[ty]*d_PROFZ_HEL[tz]*d_FORCING_COEF2[1];
        profyz_hel_coef2.z = d_PROFY_HEL[ty]*d_PROFZ_HEL[tz]*d_FORCING_COEF2[2];

	thrust::complex<REAL> fxyz = thrust::complex<REAL>(d_FORCING_FX[tx][0],d_FORCING_FX[tx][1])
                                    *thrust::complex<REAL>(d_FORCING_FY[ty][0],d_FORCING_FY[ty][1])
                                    *thrust::complex<REAL>(d_FORCING_FZ[tz][0],d_FORCING_FZ[tz][1]), tmp;
	REAL force_ampl=d_PROFX_AMPL[tx]*d_PROFY_AMPL[ty]*d_PROFZ_AMPL[tz];

#undef real //!!!
    if (axis == X_AXIS){
        tmp=thrust::complex<REAL>(d_FORCING_COEF1[0],d_PROFX_HEL[tx]*profyz_hel_coef2.x)*fxyz;
        return force_ampl*tmp.real()*d_FORCING_FDA[0];
    }
    else if (axis == Y_AXIS){
        tmp=thrust::complex<REAL>(d_FORCING_COEF1[1],d_PROFX_HEL[tx]*profyz_hel_coef2.y)*fxyz;
        return force_ampl*tmp.real()*d_FORCING_FDA[1];
    }
    else if (axis == Z_AXIS){
        tmp=thrust::complex<REAL>(d_FORCING_COEF1[2],d_PROFX_HEL[tx]*profyz_hel_coef2.z)*fxyz;
        return force_ampl*tmp.real()*d_FORCING_FDA[2];
    }
    else
      return 0.;
#define real REAL
}
#else
template <AXIS axis>
static __device__ __inline__ real forcing(const int tx, const int ty, const int tz)
{
    const real k_dot_x = (d_DSX*tx + d_DSX_OFFSET - d_XORIG)*d_KK_VEC_X 
                       + (d_DSY*ty + d_DSY_OFFSET - d_YORIG)*d_KK_VEC_Y 
                       + (d_DSZ*tz + d_DSZ_OFFSET - d_ZORIG)*d_KK_VEC_Z;
    const real waves = cos(k_dot_x)*cos(d_PHI) - sin(k_dot_x)*sin(d_PHI);

    if (axis == X_AXIS)
        return d_FORCING_KK_PART_X*waves;
    else if (axis == Y_AXIS)
        return d_FORCING_KK_PART_Y*waves;
    else
        return d_FORCING_KK_PART_Z*waves;
}
#endif
#endif
// Front: the outer computational domain at z=ZBOUND_SIZE...(ZBOUND_SIZE+BOUND_SIZE)
// Mid: segment at z=(ZBOUND_SIZE+BOUND_SIZE)...(ZBOUND_SIZE + d_nz - BOUND_SIZE)
// Back: segment at z=(ZBOUND_SIZE + d_nz - BOUND_SIZE)...(ZBOUND_SIZE + d_nz)
typedef enum {SEGMENT_FRONT=0, SEGMENT_MID, SEGMENT_BACK, SEGMENT_FULL, NUM_SEGMENTS} SegmentType;

#ifdef HYDRO
template <int step_number>
__launch_bounds__(RK_THREADS_PER_BLOCK, 1)
static __global__ void hydro_step(const real* __restrict__ d_lnrho,  //SOURCE
                        const real* __restrict__ d_uux, 
                        const real* __restrict__ d_uuy, 
                        const real* __restrict__ d_uuz, 
                        real* __restrict__ d_lnrho_dst,              //DESTINATION
                        real* __restrict__ d_uux_dst, 
                        real* __restrict__ d_uuy_dst, 
                        real* __restrict__ d_uuz_dst, 
                        const real dt,
                        SegmentType segtype)
{
    int zstart;
    int zmax;
    switch (segtype) {
        case SEGMENT_FRONT://OK
            zstart = d_nz_min;
            zmax = d_nz_min + 2*BOUND_SIZE;
            break;
        case SEGMENT_MID://OK
            zstart = d_nz_min + BOUND_SIZE;
            zmax = d_nz_max;
            break;
        case SEGMENT_BACK:
            zstart = d_nz_max - BOUND_SIZE;
            zmax = d_nz_max+BOUND_SIZE;
            break;
        default: //SEGMENT_FULL otherwise
            zstart = d_nz_min;
            zmax = d_nz_max + BOUND_SIZE;
            break;
    }
#ifdef GPU_ASTAROTH
    #define alphas d_ALPHA_TS
    #define betas d_BETA_TS
    //if (threadIdx.x==0 && threadIdx.y && threadIdx.z==0 &&
    //    blockIdx.x==0 && blockIdx.y==0 && blockIdx.z==0) printf("d_CS20= %f \n",d_CS20);
    /*if (threadIdx.x==0 && threadIdx.y && threadIdx.z==0 &&
        blockIdx.x==0 && blockIdx.y==0 && blockIdx.z==0) printf("d_NU= %f \n",d_NU);
    if (threadIdx.x==0 && threadIdx.y && threadIdx.z==0 &&
        blockIdx.x==0 && blockIdx.y==0 && blockIdx.z==0) printf("step_number,ALPHA,BETA= %d %f %f \n",step_number,ALPHA,BETA);*/
#else
    //const real alphas[] = {0.0, -0.53125, -1.1851851851851851};
    //const real betas[]  = {0.25, 0.88888888888888884, 0.75, 0.0};
    // from Williamson (1980)
    const real alphas[] = {0.0 ,  -5/9., -153/128. };
    const real betas[]  = {1/3., 15/16.,     8/15. };
#endif    
    const real ALPHA = alphas[step_number];
    const real BETA = betas[step_number];
    const real INVBETAPREV = step_number==0 ? 0. : 1./betas[step_number-1];      //1.0 / betas[(4+step_number-1) % 4];

    const int tx = threadIdx.x + blockIdx.x*blockDim.x + XBOUND_SIZE;//Start within comp domain
    const int ty = threadIdx.y + blockIdx.y*blockDim.y + YBOUND_SIZE;//Start within comp domain
    const int tz = threadIdx.z + blockIdx.z*blockDim.z*RK_ELEMS_PER_THREAD + (zstart-BOUND_SIZE);//Start from bound zone
    const int grid_idx = tx + ty*d_mx + tz*d_mxy;

    //Registers/////////////////////////////////////////////////////////////////
    //Z pencil
    const int Z_PENCIL_LENGTH = 2*BOUND_SIZE + 1;
    register real r_lnrho[Z_PENCIL_LENGTH] = {NAN};
    register real r_uux[Z_PENCIL_LENGTH]   = {NAN};
    register real r_uuy[Z_PENCIL_LENGTH]   = {NAN};
    register real r_uuz[Z_PENCIL_LENGTH]   = {NAN};

    //Partial momentum
    const int PART_MOM_SIZE = 2*BOUND_SIZE+1;
    register real mom_x[PART_MOM_SIZE] = {NAN};
    register real mom_y[PART_MOM_SIZE] = {NAN};
    register real mom_z[PART_MOM_SIZE] = {NAN};    
    ////////////////////////////////////////////////////////////////////////////

    //Shared memory///////////////////////////////////////////////////////////// 
    const int SMEM_SIZE = SMEM_WIDTH * SMEM_HEIGHT * SMEM_DEPTH;
    __shared__ real s_lnrho[SMEM_SIZE];
    __shared__ real s_uux[SMEM_SIZE];
    __shared__ real s_uuy[SMEM_SIZE];
    __shared__ real s_uuz[SMEM_SIZE];
    const int smem_idx = threadIdx.x + BOUND_SIZE + (threadIdx.y+BOUND_SIZE)*SMEM_WIDTH;
    ////////////////////////////////////////////////////////////////////////////

    //Special case: initialize registers near initial boundary 
    #pragma unroll 
    for (int k=BOUND_SIZE; k < Z_PENCIL_LENGTH-1; ++k) {
        const int curr_idx = grid_idx + (k - BOUND_SIZE)*d_mx*d_my;
        r_lnrho[k] = d_lnrho[curr_idx];
        r_uux [k] = d_uux[curr_idx];    
        r_uuy [k] = d_uuy[curr_idx];    
        r_uuz [k] = d_uuz[curr_idx];       
    }

    for (int k=0; k < RK_ELEMS_PER_THREAD + 2*ZBOUND_SIZE; ++k) {
        if (tz + k >= zmax)
            break;//Continue or break, depends on whether we want to unroll this or not

        const int curr_idx = grid_idx + k*d_mxy;

        //Update the current smem slab
        __syncthreads();
        s_lnrho[smem_idx] = r_lnrho[3]; load_halos(smem_idx, curr_idx, s_lnrho, d_lnrho);
        s_uux[smem_idx]   = r_uux[3];   load_halos(smem_idx, curr_idx, s_uux, d_uux);
        s_uuy[smem_idx]   = r_uuy[3];   load_halos(smem_idx, curr_idx, s_uuy, d_uuy);
        s_uuz[smem_idx]   = r_uuz[3];   load_halos(smem_idx, curr_idx, s_uuz, d_uuz);
        __syncthreads();

        real preloaded_lnrho_dst, preloaded_uux_dst, preloaded_uuy_dst, preloaded_uuz_dst;
        if (k >= ZBOUND_SIZE && k < ZBOUND_SIZE+RK_ELEMS_PER_THREAD)
            preloaded_lnrho_dst = d_lnrho_dst[curr_idx];
        if (k >= 2*ZBOUND_SIZE) {   
            preloaded_uux_dst = d_uux_dst[curr_idx - BOUND_SIZE*d_mxy];
            preloaded_uuy_dst = d_uuy_dst[curr_idx - BOUND_SIZE*d_mxy];
            preloaded_uuz_dst = d_uuz_dst[curr_idx - BOUND_SIZE*d_mxy];
        }

        //Update the leading slab in registers
        if (k+BOUND_SIZE < RK_ELEMS_PER_THREAD + 2*ZBOUND_SIZE && tz + k + BOUND_SIZE < zmax) {
            const int next_idx = curr_idx + BOUND_SIZE*d_mxy; 
            assert(next_idx < d_mxy*d_mz);
            r_lnrho[6] = d_lnrho[next_idx];
            r_uux  [6] = d_uux[next_idx];
            r_uuy  [6] = d_uuy[next_idx];
            r_uuz  [6] = d_uuz[next_idx];
        }

        HydroStencil stncl = {s_lnrho, s_uux, s_uuy, s_uuz, r_lnrho, r_uux, r_uuy, r_uuz};

        //Solve partial divergence
        
        mom_x[0] -= d_NU*(1.0/3.0)*der2_scalxz<3>(smem_idx, s_uuz);
        mom_x[1] -= d_NU*(1.0/3.0)*der2_scalxz<2>(smem_idx, s_uuz);
        mom_x[2] -= d_NU*(1.0/3.0)*der2_scalxz<1>(smem_idx, s_uuz);
        mom_x[4] += d_NU*(1.0/3.0)*der2_scalxz<1>(smem_idx, s_uuz);
        mom_x[5] += d_NU*(1.0/3.0)*der2_scalxz<2>(smem_idx, s_uuz);
        mom_x[6]  = d_NU*(1.0/3.0)*der2_scalxz<3>(smem_idx, s_uuz);

        mom_y[0] -= d_NU*(1.0/3.0)*der2_scalyz<3>(smem_idx, s_uuz);
        mom_y[1] -= d_NU*(1.0/3.0)*der2_scalyz<2>(smem_idx, s_uuz);
        mom_y[2] -= d_NU*(1.0/3.0)*der2_scalyz<1>(smem_idx, s_uuz);
        mom_y[4] += d_NU*(1.0/3.0)*der2_scalyz<1>(smem_idx, s_uuz);
        mom_y[5] += d_NU*(1.0/3.0)*der2_scalyz<2>(smem_idx, s_uuz);
        mom_y[6]  = d_NU*(1.0/3.0)*der2_scalyz<3>(smem_idx, s_uuz);

        mom_z[0] -= d_NU*(1.0/3.0)*(der2_scalxz<3>(smem_idx, s_uux) + der2_scalyz<3>(smem_idx, s_uuy));
        mom_z[1] -= d_NU*(1.0/3.0)*(der2_scalxz<2>(smem_idx, s_uux) + der2_scalyz<2>(smem_idx, s_uuy));
        mom_z[2] -= d_NU*(1.0/3.0)*(der2_scalxz<1>(smem_idx, s_uux) + der2_scalyz<1>(smem_idx, s_uuy));
        mom_z[4] += d_NU*(1.0/3.0)*(der2_scalxz<1>(smem_idx, s_uux) + der2_scalyz<1>(smem_idx, s_uuy));
        mom_z[5] += d_NU*(1.0/3.0)*(der2_scalxz<2>(smem_idx, s_uux) + der2_scalyz<2>(smem_idx, s_uuy));
        mom_z[6]  = d_NU*(1.0/3.0)*(der2_scalxz<3>(smem_idx, s_uux) + der2_scalyz<3>(smem_idx, s_uuy));

        if (k >= ZBOUND_SIZE && k < ZBOUND_SIZE+RK_ELEMS_PER_THREAD && tz + k < zmax - BOUND_SIZE) {
            //if (threadIdx.x == threadIdx.y && threadIdx.x == 0 && blockIdx.x == blockIdx.y && blockIdx.y == 0)
                //printf("Solving cont at %d (seg %d). Comp domain at (%d, %d)\n", tz + k, (int)segtype, d_nz_min, d_nz_max);
            const real cont_res = continuity(smem_idx, stncl);
            if (!step_number)
                d_lnrho_dst[curr_idx] = r_lnrho[3] + BETA*dt*cont_res;
            else
                d_lnrho_dst[curr_idx] = r_lnrho[3] + BETA*(ALPHA*(r_lnrho[3] - preloaded_lnrho_dst) * INVBETAPREV + dt*cont_res);

            mom_x[3] += momentum<X_AXIS>(smem_idx, stncl);
            mom_y[3] += momentum<Y_AXIS>(smem_idx, stncl);
            mom_z[3] += momentum<Z_AXIS>(smem_idx, stncl);
        }
        if (k >= 2*ZBOUND_SIZE && tz + k < zmax) {
            const int write_idx = curr_idx - BOUND_SIZE*d_mxy;
            if (!step_number) {
                d_uux_dst[write_idx] = r_uux[0] + BETA*dt*mom_x[0];
                d_uuy_dst[write_idx] = r_uuy[0] + BETA*dt*mom_y[0];
                d_uuz_dst[write_idx] = r_uuz[0] + BETA*dt*mom_z[0];
            } else {
                real uux_res = r_uux[0] + BETA*(dt*mom_x[0] 
                                                + ALPHA*(r_uux[0] - preloaded_uux_dst) * INVBETAPREV);
                real uuy_res = r_uuy[0] + BETA*(dt*mom_y[0] 
                                                + ALPHA*(r_uuy[0] - preloaded_uuy_dst) * INVBETAPREV);
                real uuz_res = r_uuz[0] + BETA*(dt*mom_z[0] 
                                                + ALPHA*(r_uuz[0] - preloaded_uuz_dst) * INVBETAPREV);

                #ifdef FORCING
                if (step_number == 2) {
                    const int tz_offset = tz + k - ZBOUND_SIZE;
                    uux_res += forcing<X_AXIS>(tx, ty, tz_offset);
                    uuy_res += forcing<Y_AXIS>(tx, ty, tz_offset);
                    uuz_res += forcing<Z_AXIS>(tx, ty, tz_offset);
                }
                #endif                

                d_uux_dst[write_idx] = uux_res;
                d_uuy_dst[write_idx] = uuy_res;
                d_uuz_dst[write_idx] = uuz_res;
            }
        }

        #pragma unroll
        for (int i=0; i < Z_PENCIL_LENGTH-1; ++i) {
            r_lnrho[i] = r_lnrho[i+1];
            r_uux [i] = r_uux [i+1];    
            r_uuy [i] = r_uuy [i+1];    
            r_uuz [i] = r_uuz [i+1];   
        }

        #pragma unroll
        for (int i=0; i < PART_MOM_SIZE-1; ++i) {
            mom_x[i] = mom_x[i+1];
            mom_y[i] = mom_y[i+1];
            mom_z[i] = mom_z[i+1];
        }
    }
}
#endif

#ifdef MAGNETIC
template <int step_number>
__launch_bounds__(RK_THREADS_PER_BLOCK, 1)
static __global__ void induction_step(  const real* __restrict__ d_Ax, 
                                        const real* __restrict__ d_Ay, 
                                        const real* __restrict__ d_Az,
                                        const real* __restrict__ d_uux, 
                                        const real* __restrict__ d_uuy, 
                                        const real* __restrict__ d_uuz,
                                        real* __restrict__ d_Ax_dst, 
                                        real* __restrict__ d_Ay_dst, 
                                        real* __restrict__ d_Az_dst, 
                                        const real dt,
                                        SegmentType segtype)
{
    int zstart;
    int zmax;
    switch (segtype) {
        case SEGMENT_FRONT://OK
            zstart = d_nz_min;
            zmax = d_nz_min + 2*BOUND_SIZE;
            break;
        case SEGMENT_MID://OK
            zstart = d_nz_min + BOUND_SIZE;
            zmax = d_nz_max;
            break;
        case SEGMENT_BACK:
            zstart = d_nz_max - BOUND_SIZE;
            zmax = d_nz_max+BOUND_SIZE;
            break;
        default: //SEGMENT_FULL otherwise
            zstart = d_nz_min;
            zmax = d_nz_max + BOUND_SIZE;
            break;
        
    }
#ifdef GPU_ASTAROTH
    #define alphas d_ALPHA_TS
    #define betas d_BETA_TS
#else
    //const real alphas[] = {0.0, -0.53125, -1.1851851851851851};
    //const real betas[]  = {0.25, 0.88888888888888884, 0.75, 0.0};
    const real alphas[] = {0.0 ,  -5/9., -153/128. };
    const real betas[]  = {1/3., 15/16.,     8/15. };
#endif
    const real ALPHA = alphas[step_number];
    const real BETA = betas[step_number];
    const real INVBETAPREV = step_number==0 ? 0. : 1./betas[step_number-1];      //1.0 / betas[(4+step_number-1) % 4];

    const int tx = threadIdx.x + blockIdx.x*blockDim.x + XBOUND_SIZE;//Start within comp domain
    const int ty = threadIdx.y + blockIdx.y*blockDim.y + YBOUND_SIZE;//Start within comp domain
    const int tz = threadIdx.z + blockIdx.z*blockDim.z*RK_ELEMS_PER_THREAD + (zstart-BOUND_SIZE);//Start from bound zone
    const int grid_idx = tx + ty*d_mx + tz*d_mxy;


    //Registers/////////////////////////////////////////////////////////////////
    //Z pencil
    const int Z_PENCIL_LENGTH = 2*BOUND_SIZE + 1;
    register real r_Ax[Z_PENCIL_LENGTH]   = {NAN};
    register real r_Ay[Z_PENCIL_LENGTH]   = {NAN};
    register real r_Az[Z_PENCIL_LENGTH]   = {NAN};

    //Partial magnetic vector potential
    const int PART_A_SIZE = 2*BOUND_SIZE + 1;
    register real part_Ax[PART_A_SIZE] = {NAN};
    register real part_Ay[PART_A_SIZE] = {NAN};
    register real part_Az[PART_A_SIZE] = {NAN};    
    ////////////////////////////////////////////////////////////////////////////

    //Shared memory///////////////////////////////////////////////////////////// 
    const int SMEM_SIZE = SMEM_WIDTH * SMEM_HEIGHT * SMEM_DEPTH;
    __shared__ real s_Ax[SMEM_SIZE];
    __shared__ real s_Ay[SMEM_SIZE];
    __shared__ real s_Az[SMEM_SIZE];
    const int smem_idx = threadIdx.x + BOUND_SIZE + (threadIdx.y+BOUND_SIZE)*SMEM_WIDTH;
    ////////////////////////////////////////////////////////////////////////////

    //Special case: initialize registers near initial boundary 
    #pragma unroll 
    for (int k=BOUND_SIZE; k < Z_PENCIL_LENGTH-1; ++k) {
        const int curr_idx = grid_idx + (k - BOUND_SIZE)*d_mx*d_my;
        r_Ax [k] = d_Ax[curr_idx];    
        r_Ay [k] = d_Ay[curr_idx];    
        r_Az [k] = d_Az[curr_idx];       
    }


    for (int k=0; k < RK_ELEMS_PER_THREAD + 2*ZBOUND_SIZE; ++k) {
        if (tz + k >= zmax)
            break;//Continue or break, depends on whether we want to unroll this or not

        const int curr_idx = grid_idx + k*d_mxy;

        //Update the current smem slab
        __syncthreads();
        s_Ax[smem_idx]   = r_Ax[3];   load_halos(smem_idx, curr_idx, s_Ax, d_Ax);
        s_Ay[smem_idx]   = r_Ay[3];   load_halos(smem_idx, curr_idx, s_Ay, d_Ay);
        s_Az[smem_idx]   = r_Az[3];   load_halos(smem_idx, curr_idx, s_Az, d_Az);
        __syncthreads();

        //Load local uu
        const real uux = d_uux[curr_idx];
        const real uuy = d_uuy[curr_idx];
        const real uuz = d_uuz[curr_idx]; 
        //

        //Update the leading slab in registers
        if (k+BOUND_SIZE < RK_ELEMS_PER_THREAD + 2*ZBOUND_SIZE && tz + k + BOUND_SIZE < zmax) {
            const int next_idx = curr_idx + BOUND_SIZE*d_mxy; 
            assert(next_idx < d_mxy*d_mz);
            r_Ax  [6] = d_Ax[next_idx];
            r_Ay  [6] = d_Ay[next_idx];
            r_Az  [6] = d_Az[next_idx];
        }

        InductionStencil stncl = {s_Ax, s_Ay, s_Az, r_Ax, r_Ay, r_Az};

        //Solve partial divergence
        part_Ax[0] -= - d_ETA * der2_scalxz<3>(smem_idx, s_Az);
        part_Ax[1] -= - d_ETA * der2_scalxz<2>(smem_idx, s_Az);
        part_Ax[2] -= - d_ETA * der2_scalxz<1>(smem_idx, s_Az);
        part_Ax[3] += induction<X_AXIS>(smem_idx, stncl, uux, uuy, uuz);
        part_Ax[4] += - d_ETA * der2_scalxz<1>(smem_idx, s_Az);
        part_Ax[5] += - d_ETA * der2_scalxz<2>(smem_idx, s_Az);
        part_Ax[6]  = - d_ETA * der2_scalxz<3>(smem_idx, s_Az);

        part_Ay[0] -= - d_ETA * der2_scalyz<3>(smem_idx, s_Az);
        part_Ay[1] -= - d_ETA * der2_scalyz<2>(smem_idx, s_Az);
        part_Ay[2] -= - d_ETA * der2_scalyz<1>(smem_idx, s_Az);
        part_Ay[3] += induction<Y_AXIS>(smem_idx, stncl, uux, uuy, uuz);
        part_Ay[4] += - d_ETA * der2_scalyz<1>(smem_idx, s_Az);
        part_Ay[5] += - d_ETA * der2_scalyz<2>(smem_idx, s_Az);
        part_Ay[6]  = - d_ETA * der2_scalyz<3>(smem_idx, s_Az);

        part_Az[0] -= - d_ETA * (der2_scalxz<3>(smem_idx, s_Ax) + der2_scalyz<3>(smem_idx, s_Ay));
        part_Az[1] -= - d_ETA * (der2_scalxz<2>(smem_idx, s_Ax) + der2_scalyz<2>(smem_idx, s_Ay));
        part_Az[2] -= - d_ETA * (der2_scalxz<1>(smem_idx, s_Ax) + der2_scalyz<1>(smem_idx, s_Ay));
        part_Az[3] += induction<Z_AXIS>(smem_idx, stncl, uux, uuy, uuz);
        part_Az[4] += - d_ETA * (der2_scalxz<1>(smem_idx, s_Ax) + der2_scalyz<1>(smem_idx, s_Ay));
        part_Az[5] += - d_ETA * (der2_scalxz<2>(smem_idx, s_Ax) + der2_scalyz<2>(smem_idx, s_Ay));
        part_Az[6]  = - d_ETA*(der2_scalxz<3>(smem_idx, s_Ax) + der2_scalyz<3>(smem_idx, s_Ay));


        if (k >= 2*ZBOUND_SIZE && tz + k < zmax) {
            const int write_idx = curr_idx - BOUND_SIZE*d_mxy;
            if (!step_number) {
                d_Ax_dst[write_idx] = r_Ax[0] + BETA*dt*part_Ax[0];
                d_Ay_dst[write_idx] = r_Ay[0] + BETA*dt*part_Ay[0];
                d_Az_dst[write_idx] = r_Az[0] + BETA*dt*part_Az[0];
            } else {
                d_Ax_dst[write_idx] = r_Ax[0] + BETA*(dt*part_Ax[0] 
                                                + ALPHA*(r_Ax[0] - d_Ax_dst[write_idx]) * INVBETAPREV);
                d_Ay_dst[write_idx] = r_Ay[0] + BETA*(dt*part_Ay[0] 
                                                + ALPHA*(r_Ay[0] - d_Ay_dst[write_idx]) * INVBETAPREV);
                d_Az_dst[write_idx] = r_Az[0] + BETA*(dt*part_Az[0] 
                                                + ALPHA*(r_Az[0] - d_Az_dst[write_idx]) * INVBETAPREV);
            }
        }

        #pragma unroll
        for (int i=0; i < Z_PENCIL_LENGTH-1; ++i) {
            r_Ax[i] = r_Ax [i+1];    
            r_Ay[i] = r_Ay [i+1];    
            r_Az[i] = r_Az [i+1];   
        }

        #pragma unroll
        for (int i=0; i < PART_A_SIZE-1; ++i) {
            part_Ax[i] = part_Ax[i+1];
            part_Ay[i] = part_Ay[i+1];
            part_Az[i] = part_Az[i+1];
        }
    }
}
#endif

template<int step_number>
static void rk3_inner_step_cuda_generic(Grid* d_grid, Grid* d_grid_dst, const real dt, CParamConfig* cparams, cudaStream_t hydro_stream, cudaStream_t induct_stream)
//printf("addresses %d %d %d %d\n",d_grid->arr[d_grid->LNRHO], d_grid->arr[d_grid->UUX], d_grid->arr[d_grid->UUY], d_grid->arr[d_grid->UUZ]);
{
//printf("indices %d %d %d %d\n",d_grid->LNRHO, d_grid->UUX, d_grid->UUY, d_grid->UUZ);
    const dim3 tpb((unsigned int) min(RK_THREADS_X, cparams->nx), RK_THREADS_Y, RK_THREADS_Z);
    const dim3 bpg((unsigned int) ceil(cparams->nx / (real)(tpb.x)),
                   (unsigned int) ceil(cparams->ny / (real)(tpb.y)),
                   (unsigned int) ceil((cparams->nz - 2*BOUND_SIZE) / (real)(tpb.z*RK_ELEMS_PER_THREAD)));
    //cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);
    //cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);

//printf("tpb= %d %d %d, bpg= %d %d %d %d \n", tpb.x,tpb.y,tpb.z, bpg.x,bpg.y,bpg.z,cparams->nz - 2*BOUND_SIZE);
    //INTEGRATE
    #ifdef HYDRO
    hydro_step<step_number><<<bpg, tpb, 0, hydro_stream>>>(d_grid->arr[d_grid->LNRHO], 
                                                           d_grid->arr[d_grid->UUX], 
                                                           d_grid->arr[d_grid->UUY], 
                                                           d_grid->arr[d_grid->UUZ],
                                                           d_grid_dst->arr[d_grid_dst->LNRHO], 
                                                           d_grid_dst->arr[d_grid_dst->UUX], 
                                                           d_grid_dst->arr[d_grid_dst->UUY], 
                                                           d_grid_dst->arr[d_grid_dst->UUZ],  
                                                           dt, SEGMENT_MID); 
    CUDA_ERRCHK_KERNEL();
    #endif

    #ifdef MAGNETIC
        induction_step<step_number><<<bpg, tpb, 0, induct_stream>>>(d_grid->arr[d_grid->AAX], 
                                                                    d_grid->arr[d_grid->AAY], 
                                                                    d_grid->arr[d_grid->AAZ], 
                                                                    d_grid->arr[d_grid->UUX], 
                                                                    d_grid->arr[d_grid->UUY], 
                                                                    d_grid->arr[d_grid->UUZ], 
                                                                    d_grid_dst->arr[d_grid_dst->AAX], 
                                                                    d_grid_dst->arr[d_grid_dst->AAY], 
                                                                    d_grid_dst->arr[d_grid_dst->AAZ],
                                                                    dt, SEGMENT_MID);
        CUDA_ERRCHK_KERNEL();
    #else
        (void) induct_stream;      //Suppress warning about unused parameter
    #endif
}

template<int step_number>
static void rk3_outer_step_cuda_generic(Grid* d_grid, Grid* d_grid_dst, const real dt, CParamConfig* cparams, cudaStream_t stream)
{
    const dim3 tpb_fb((unsigned int)min(RK_THREADS_X, cparams->nx), RK_THREADS_Y, RK_THREADS_Z);
    const dim3 bpg_fb((unsigned int) ceil(cparams->nx / (real)(tpb_fb.x)),
                      (unsigned int) ceil(cparams->ny / (real)(tpb_fb.y)),
                      (unsigned int) ceil(BOUND_SIZE / (real)(tpb_fb.z*RK_ELEMS_PER_THREAD)));

    //cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);
    //cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);

    //INTEGRATE

    #ifdef HYDRO
    hydro_step<step_number><<<bpg_fb, tpb_fb, 0, stream>>>(d_grid->arr[d_grid->LNRHO], 
                                                           d_grid->arr[d_grid->UUX], 
                                                           d_grid->arr[d_grid->UUY], 
                                                           d_grid->arr[d_grid->UUZ],
                                                           d_grid_dst->arr[d_grid_dst->LNRHO], 
                                                           d_grid_dst->arr[d_grid_dst->UUX], 
                                                           d_grid_dst->arr[d_grid_dst->UUY], 
                                                           d_grid_dst->arr[d_grid_dst->UUZ],  
                                                           dt, SEGMENT_FRONT); 
    hydro_step<step_number><<<bpg_fb, tpb_fb, 0, stream>>>(d_grid->arr[d_grid->LNRHO], 
                                                           d_grid->arr[d_grid->UUX], 
                                                           d_grid->arr[d_grid->UUY], 
                                                           d_grid->arr[d_grid->UUZ],
                                                           d_grid_dst->arr[d_grid_dst->LNRHO], 
                                                           d_grid_dst->arr[d_grid_dst->UUX], 
                                                           d_grid_dst->arr[d_grid_dst->UUY], 
                                                           d_grid_dst->arr[d_grid_dst->UUZ],  
                                                           dt, SEGMENT_BACK);
    CUDA_ERRCHK_KERNEL();
    #endif

    #ifdef MAGNETIC
        induction_step<step_number><<<bpg_fb, tpb_fb, 0, stream>>>(d_grid->arr[d_grid->AAX], 
                                                                   d_grid->arr[d_grid->AAY], 
                                                                   d_grid->arr[d_grid->AAZ], 
                                                                   d_grid->arr[d_grid->UUX], 
                                                                   d_grid->arr[d_grid->UUY], 
                                                                   d_grid->arr[d_grid->UUZ], 
                                                                   d_grid_dst->arr[AAX], 
                                                                   d_grid_dst->arr[AAY], 
                                                                   d_grid_dst->arr[AAZ],
                                                                   dt, SEGMENT_FRONT);
        induction_step<step_number><<<bpg_fb, tpb_fb, 0, stream>>>(d_grid->arr[d_grid->AAX], 
                                                                   d_grid->arr[d_grid->AAY], 
                                                                   d_grid->arr[d_grid->AAZ], 
                                                                   d_grid->arr[d_grid->UUX], 
                                                                   d_grid->arr[d_grid->UUY], 
                                                                   d_grid->arr[d_grid->UUZ], 
                                                                   d_grid_dst->arr[d_grid_dst->AAX], 
                                                                   d_grid_dst->arr[d_grid_dst->AAY], 
                                                                   d_grid_dst->arr[d_grid_dst->AAZ],
                                                                   dt, SEGMENT_BACK);
        CUDA_ERRCHK_KERNEL();
    #endif
}

/* MR: not needed
//This is just here s.t. we can pass step_number as parameter 
//(easier to interface with Pencil Code without templates)
void rk3_cuda_generic(Grid* d_grid, Grid* d_grid_dst, 
                      const int step_number, const real dt, CParamConfig* cparams, 
                      cudaStream_t stream)
{
    switch(step_number) {
        case 0:
            rk3_step_cuda_generic<0>(d_grid, d_grid_dst, dt, cparams, stream);
            break;
        case 1:
            rk3_step_cuda_generic<1>(d_grid, d_grid_dst, dt, cparams, stream);
            break;
        case 2:
            rk3_step_cuda_generic<2>(d_grid, d_grid_dst, dt, cparams, stream);
            break;
        default :
            CRASH("Invalid step number in rk3_cuda_generic");
    }
}
*/

void rk3_inner_cuda_generic(Grid* d_grid, Grid* d_grid_dst, 
                      const int step_number, const real dt, CParamConfig* cparams, 
                      cudaStream_t hydro_stream, cudaStream_t induct_stream)
{
    switch(step_number) {
        case 0:
            rk3_inner_step_cuda_generic<0>(d_grid, d_grid_dst, dt, cparams, hydro_stream, induct_stream);
            break;
        case 1:
            rk3_inner_step_cuda_generic<1>(d_grid, d_grid_dst, dt, cparams, hydro_stream, induct_stream);
            break;
        case 2:
            rk3_inner_step_cuda_generic<2>(d_grid, d_grid_dst, dt, cparams, hydro_stream, induct_stream);
            break;
        default :
            CRASH("Invalid step number in rk3_cuda_generic");
    }
}

void rk3_outer_cuda_generic(Grid* d_grid, Grid* d_grid_dst, 
                      const int step_number, const real dt, CParamConfig* cparams, 
                      cudaStream_t stream)
{
    switch(step_number) {
        case 0:
            rk3_outer_step_cuda_generic<0>(d_grid, d_grid_dst, dt, cparams, stream);
            break;
        case 1:
            rk3_outer_step_cuda_generic<1>(d_grid, d_grid_dst, dt, cparams, stream);
            break;
        case 2:
            rk3_outer_step_cuda_generic<2>(d_grid, d_grid_dst, dt, cparams, stream);
            break;
        default :
            CRASH("Invalid step number in rk3_cuda_generic");
    }
}
