#include "rk3_entropy.cuh"

#include "gpu/cuda/core/dconsts_core.cuh"
#include "gpu/cuda/core/errorhandler_cuda.cuh"

// Kernel configuration
static const dim3 tpb(4, 4, 4);

// Grid indices
#define IDX_1D(i, j, k) ((i) + (j)*d_mx + (k)*d_mxy)

#ifdef ENTROPY 
static __device__ __forceinline__ real
dot(const real a0, const real a1, const real a2,
    const real b0, const real b1, const real b2)
{
    return a0*b0 + a1*b1 + a2*b2;
}

/** NOTE DANGER!!! Make sure this function is called with lnrho and entropy in
    correct order! */
static __device__ __forceinline__ real
get_lnT(const real lnrho, const real entropy_s)
{

    const real lnT = d_LNTT0 +
                     (entropy_s / d_CV) +
                     (d_GAMMA - real(1.0))*(lnrho - d_LNRHO0);

    // TODO NOTE: This is not the actual lnT, this part accumulates huge
    // arithmetic error so replaced temporarily with a toy solution that's free
    // of any significant loss in precision (see model_rk3.cc also)
    //const real lnT = .5f * (entropy_s + lnrho);
    // NOTE UPDATE: The error comes from d_LNTT0: if it is bigger than the
    // range of the values in the grid, then with addition we inherently lose
    // precision. Therefore special care must be taken when choosing d_LNTT0
    return lnT;
}

static __device__ __forceinline__ real
get_inv_pT(const real lnrho, const real entropy)
{
    const real lnT = get_lnT(lnrho, entropy);
    return real(1.) / exp(lnrho + lnT);
}

/**
*   Returns the derivative at point f3. The input parameters f0..f2 are the
*   points to the "left" from f3, and f4..f6 are on the "right"
*/
static __device__ __forceinline__ real
der_scal(const real f0, const real f1, const real f2,
         const real f4, const real f5, const real f6,
         const real ds)// Grid spacing, f.ex. cparams->dsx
{
    const real res = (          (f6 - f0)
                             + real(9.)  * (f1 - f5)
                             + real(45.) * (f4 - f2)
                            )
                            / (real(60.) * ds);
    return res;
}

static __device__ __forceinline__ real
der2_scal(const real f0, const real f1, const real f2,
          const real f3,
          const real f4, const real f5, const real f6,
          const real ds)// Grid spacing, f.ex. cparams->dsx
{
    const real res = (
	                          real(2.)   * (f0 + f6)
	                        - real(27.)  * (f1 + f5)
	                        + real(270.) * (f2 + f4)
	                        - real(490.) * f3
                            )
	                        / (real(180.)*ds*ds);
    return res;
}

static __device__ __forceinline__ real
der_scalx(const int i, const int j, const int k, const real* __restrict__ scal)
{
    return der_scal(scal[IDX_1D(i-3, j, k)],
                     scal[IDX_1D(i-2, j, k)],
                     scal[IDX_1D(i-1, j, k)],
                     scal[IDX_1D(i+1, j, k)],
                     scal[IDX_1D(i+2, j, k)],
                     scal[IDX_1D(i+3, j, k)],
                     d_DSX);
}

static __device__ real
der_scaly(const int i, const int j, const int k, const real* __restrict__ scal)
{
    return der_scal(scal[IDX_1D(i, j-3, k)],
                     scal[IDX_1D(i, j-2, k)],
                     scal[IDX_1D(i, j-1, k)],
                     scal[IDX_1D(i, j+1, k)],
                     scal[IDX_1D(i, j+2, k)],
                     scal[IDX_1D(i, j+3, k)],
                     d_DSY);
}

static __device__ real
der_scalz(const int i, const int j, const int k, const real* __restrict__ scal)
{
    return der_scal(scal[IDX_1D(i, j, k-3)],
                     scal[IDX_1D(i, j, k-2)],
                     scal[IDX_1D(i, j, k-1)],
                     scal[IDX_1D(i, j, k+1)],
                     scal[IDX_1D(i, j, k+2)],
                     scal[IDX_1D(i, j, k+3)],
                     d_DSZ);
}

static __device__ __forceinline__ real
der2_scalx(const int i, const int j, const int k, const real* __restrict__ scal)
{
    return der2_scal(scal[IDX_1D(i-3, j, k)],
                     scal[IDX_1D(i-2, j, k)],
                     scal[IDX_1D(i-1, j, k)],
                     scal[IDX_1D(i, j, k)],
                     scal[IDX_1D(i+1, j, k)],
                     scal[IDX_1D(i+2, j, k)],
                     scal[IDX_1D(i+3, j, k)],
                     d_DSX);
}

static __device__ real
der2_scaly(const int i, const int j, const int k, const real* __restrict__ scal)
{
    return der2_scal(scal[IDX_1D(i, j-3, k)],
                     scal[IDX_1D(i, j-2, k)],
                     scal[IDX_1D(i, j-1, k)],
                     scal[IDX_1D(i, j, k)],
                     scal[IDX_1D(i, j+1, k)],
                     scal[IDX_1D(i, j+2, k)],
                     scal[IDX_1D(i, j+3, k)],
                     d_DSY);
}

static __device__ real
der2_scalz(const int i, const int j, const int k, const real* __restrict__ scal)
{
    return der2_scal(scal[IDX_1D(i, j, k-3)],
                     scal[IDX_1D(i, j, k-2)],
                     scal[IDX_1D(i, j, k-1)],
                     scal[IDX_1D(i, j, k)],
                     scal[IDX_1D(i, j, k+1)],
                     scal[IDX_1D(i, j, k+2)],
                     scal[IDX_1D(i, j, k+3)],
                     d_DSZ);
}

static __device__ real
vec_dot_nabla_scal(const int i, const int j, const int k,
                   const real* __restrict__ vecx,
                   const real* __restrict__ vecy,
                   const real* __restrict__ vecz,
                   const real* __restrict__ scal)
{
    const int idx = IDX_1D(i, j, k);

    return vecx[idx] * der_scalx(i, j, k, scal) +
           vecy[idx] * der_scaly(i, j, k, scal) +
           vecz[idx] * der_scalz(i, j, k, scal);
}

static __device__ __forceinline__ real
der_scalx_lnT(const int i, const int j, const int k,
               const real* __restrict__ lnrho,
               const real* __restrict__ entropy_s)
{
    const dim3 offset(1, 0, 0);
    const real ds = d_DSX;

    real f[7];
    for (int counter = 0; counter < 7; ++counter) {
        if (counter  == 3)
            continue;
        const size_t idx = IDX_1D(i - (counter-3)*offset.x,
                                  j - (counter-3)*offset.y,
                                  k - (counter-3)*offset.z);
        f[counter] = get_lnT(lnrho[idx], entropy_s[idx]);
    }
    return der_scal(f[0], f[1], f[2], f[4], f[5], f[6], ds);
}

static __device__ real
der_scaly_lnT(const int i, const int j, const int k,
               const real* __restrict__ lnrho,
               const real* __restrict__ entropy_s)
{
    const dim3 offset(0, 1, 0);
    const real ds = d_DSY;

    real f[7];
    for (int counter = 0; counter < 7; ++counter) {
        if (counter  == 3)
            continue;
        const size_t idx = IDX_1D(i - (counter-3)*offset.x,
                                  j - (counter-3)*offset.y,
                                  k - (counter-3)*offset.z);
        f[counter] = get_lnT(lnrho[idx], entropy_s[idx]);
    }
    return der_scal(f[0], f[1], f[2], f[4], f[5], f[6], ds);
}

static __device__ real
der_scalz_lnT(const int i, const int j, const int k,
               const real* __restrict__ lnrho,
               const real* __restrict__ entropy_s)
{
    const dim3 offset(0, 0, 1);
    const real ds = d_DSZ;

    real f[7];
    for (int counter = 0; counter < 7; ++counter) {
        if (counter  == 3)
            continue;
        const size_t idx = IDX_1D(i - (counter-3)*offset.x,
                                  j - (counter-3)*offset.y,
                                  k - (counter-3)*offset.z);
        f[counter] = get_lnT(lnrho[idx], entropy_s[idx]);
    }
    return der_scal(f[0], f[1], f[2], f[4], f[5], f[6], ds);
}

static __device__ __forceinline__ real
der2_scalx_lnT(const int i, const int j, const int k,
               const real* __restrict__ lnrho,
               const real* __restrict__ entropy_s)
{
    const dim3 offset(1, 0, 0);
    const real ds = d_DSX;

    real f[7];
    for (int counter = 0; counter < 7; ++counter) {
        const size_t idx = IDX_1D(i - (counter-3)*offset.x,
                                  j - (counter-3)*offset.y,
                                  k - (counter-3)*offset.z);
        f[counter] = get_lnT(lnrho[idx], entropy_s[idx]);
    }
    return der2_scal(f[0], f[1], f[2], f[3], f[4], f[5], f[6], ds);
}

static __device__ real
der2_scaly_lnT(const int i, const int j, const int k,
               const real* __restrict__ lnrho,
               const real* __restrict__ entropy_s)
{
    const dim3 offset(0, 1, 0);
    const real ds = d_DSY;

    real f[7];
    for (int counter = 0; counter < 7; ++counter) {
        const size_t idx = IDX_1D(i - (counter-3)*offset.x,
                                  j - (counter-3)*offset.y,
                                  k - (counter-3)*offset.z);
        f[counter] = get_lnT(lnrho[idx], entropy_s[idx]);
    }
    return der2_scal(f[0], f[1], f[2], f[3], f[4], f[5], f[6], ds);
}

static __device__ real
der2_scalz_lnT(const int i, const int j, const int k,
               const real* __restrict__ lnrho,
               const real* __restrict__ entropy_s)
{
    const dim3 offset(0, 0, 1);
    const real ds = d_DSZ;

    real f[7];
    for (int counter = 0; counter < 7; ++counter) {
        const size_t idx = IDX_1D(i - (counter-3)*offset.x,
                                  j - (counter-3)*offset.y,
                                  k - (counter-3)*offset.z);
        f[counter] = get_lnT(lnrho[idx], entropy_s[idx]);
    }
    return der2_scal(f[0], f[1], f[2], f[3], f[4], f[5], f[6], ds);
}

static __device__ real
laplace_scal_lnT(const int i, const int j, const int k,
                 const real* __restrict__ lnrho,
                 const real* __restrict__ entropy_s)
{
    const real d2dx2_lnT = der2_scalx_lnT(i, j, k, lnrho, entropy_s);
    const real d2dy2_lnT = der2_scaly_lnT(i, j, k, lnrho, entropy_s);
    const real d2dz2_lnT = der2_scalz_lnT(i, j, k, lnrho, entropy_s);

    return d2dx2_lnT + d2dy2_lnT + d2dz2_lnT;
}

static __device__ real
dot_grad_lnT2(const int i, const int j, const int k,
                 const real* __restrict__ lnrho,
                 const real* __restrict__ entropy_s)
{
    const real ddx_lnT = der_scalx_lnT(i, j, k, lnrho, entropy_s);
    const real ddy_lnT = der_scaly_lnT(i, j, k, lnrho, entropy_s);
    const real ddz_lnT = der_scalz_lnT(i, j, k, lnrho, entropy_s);

    return dot(ddx_lnT, ddy_lnT, ddz_lnT,
               ddx_lnT, ddy_lnT, ddz_lnT);
}

static __device__ real
laplace_scal(const int i, const int j, const int k,
             const real* __restrict__ scal)
{
    return der2_scalx(i, j, k, scal) +
           der2_scaly(i, j, k, scal) +
           der2_scalz(i, j, k, scal);
}

static __device__ void
laplace_vec(const int i, const int j, const int k,
            const real* __restrict__ vecx,
            const real* __restrict__ vecy,
            const real* __restrict__ vecz,
            real* lux, real* luy, real* luz)
{
    *lux = laplace_scal(i, j, k, vecx);
    *luy = laplace_scal(i, j, k, vecy);
    *luz = laplace_scal(i, j, k, vecz);
}

typedef enum { DERNM_XY, DERNM_YZ, DERNM_XZ} DERNM_TYPE;
static __device__ real
dernm_scal(int i, int j, int k, DERNM_TYPE type, const real* __restrict__ scal)
{ // TODO NOTE casts to real
    const real dx = d_DSX;
    const real dy = d_DSY;
    const real dz = d_DSZ;

   if (type == DERNM_XY) {
      const real fac = (1.0/720.0)*(1.0/dx)*(1.0/dy);
      return fac*(
        270.0*( scal[i+1 + (j+1)*d_mx + k*d_mxy]
              - scal[i-1 + (j+1)*d_mx + k*d_mxy]
              + scal[i-1 + (j-1)*d_mx + k*d_mxy]
              - scal[i+1 + (j-1)*d_mx + k*d_mxy] )
        -27.0*( scal[i+2 + (j+2)*d_mx + k*d_mxy]
              - scal[i-2 + (j+2)*d_mx + k*d_mxy]
              + scal[i-2 + (j-2)*d_mx + k*d_mxy]
              - scal[i+2 + (j-2)*d_mx + k*d_mxy] )
        + 2.0*( scal[i+3 + (j+3)*d_mx + k*d_mxy]
              - scal[i-3 + (j+3)*d_mx + k*d_mxy]
              + scal[i-3 + (j-3)*d_mx + k*d_mxy]
              - scal[i+3 + (j-3)*d_mx + k*d_mxy] )
                );
   } else if (type == DERNM_YZ) {
      const real fac = (1.0/720.0)*(1.0/dy)*(1.0/dz);
      return fac*(
        270.0*( scal[i + (j+1)*d_mx + (k+1)*d_mxy]
              - scal[i + (j+1)*d_mx + (k-1)*d_mxy]
              + scal[i + (j-1)*d_mx + (k-1)*d_mxy]
              - scal[i + (j-1)*d_mx + (k+1)*d_mxy] )
        -27.0*( scal[i + (j+2)*d_mx + (k+2)*d_mxy]
              - scal[i + (j+2)*d_mx + (k-2)*d_mxy]
              + scal[i + (j-2)*d_mx + (k-2)*d_mxy]
              - scal[i + (j-2)*d_mx + (k+2)*d_mxy] )
        + 2.0*( scal[i + (j+3)*d_mx + (k+3)*d_mxy]
              - scal[i + (j+3)*d_mx + (k-3)*d_mxy]
              + scal[i + (j-3)*d_mx + (k-3)*d_mxy]
              - scal[i + (j-3)*d_mx + (k+3)*d_mxy] )
                );
   } else if (type == DERNM_XZ) {
      const real fac = (1.0/720.0)*(1.0/dz)*(1.0/dx);
      return fac*(
        270.0*( scal[i+1 + j*d_mx + (k+1)*d_mxy]
              - scal[i-1 + j*d_mx + (k+1)*d_mxy]
              + scal[i-1 + j*d_mx + (k-1)*d_mxy]
              - scal[i+1 + j*d_mx + (k-1)*d_mxy] )
        -27.0*( scal[i+2 + j*d_mx + (k+2)*d_mxy]
              - scal[i-2 + j*d_mx + (k+2)*d_mxy]
              + scal[i-2 + j*d_mx + (k-2)*d_mxy]
              - scal[i+2 + j*d_mx + (k-2)*d_mxy] )
        + 2.0*( scal[i+3 + j*d_mx + (k+3)*d_mxy]
              - scal[i-3 + j*d_mx + (k+3)*d_mxy]
              + scal[i-3 + j*d_mx + (k-3)*d_mxy]
              - scal[i+3 + j*d_mx + (k-3)*d_mxy] )
                );
   } else {
        return NAN;
    }
}

static __device__ void
grad_div_vec(const int i, const int j, const int k,
             const real* __restrict__ vecx,
             const real* __restrict__ vecy,
             const real* __restrict__ vecz,
             real* gdux, real* gduy, real* gduz)
{
   *gdux = der2_scalx(i, j, k, vecx)
         + dernm_scal(i, j, k, DERNM_XY, vecy)
         + dernm_scal(i, j, k, DERNM_XZ, vecz);

   *gduy = dernm_scal(i, j, k, DERNM_XY, vecx)
         + der2_scaly(i, j, k, vecy)
         + dernm_scal(i, j, k, DERNM_YZ, vecz);

   *gduz = dernm_scal(i, j, k, DERNM_XZ, vecx)
         + dernm_scal(i, j, k, DERNM_YZ, vecy)
         + der2_scalz(i, j, k, vecz);
}

static __device__ real
eta_mu0_j_dot_j(const int i, const int j, const int k,
                const real* __restrict__ ax,
                const real* __restrict__ ay,
                const real* __restrict__ az)
{
    real lux, luy, luz;
    laplace_vec(i, j, k, ax, ay, az, &lux, &luy, &luz);

    real gdux, gduy, gduz;
    grad_div_vec(i, j, k, ax, ay, az, &gdux, &gduy, &gduz);

    const real inv_mu0 = real(1.) / d_MU0;

    const real jx = inv_mu0 * (-lux + gdux);
    const real jy = inv_mu0 * (-luy + gduy);
    const real jz = inv_mu0 * (-luz + gduz);

    return d_ETA * d_MU0 * dot(jx, jy, jz, jx, jy, jz);
}

static __device__ real
contract_S_to_scal(const int i, const int j, const int k,
                    const real* __restrict__ uux,
                    const real* __restrict__ uuy,
                    const real* __restrict__ uuz)
{
    real S[3][3];
    real c23, c13;

   c23 = real(2. / 3.); c13 = real(1. / 3.);

   S[0][0] = c23*der_scalx(i, j, k, uux)
         - c13*(der_scaly(i, j, k, uuy)
         + der_scalz(i, j, k, uuz) );
   S[0][1] = real(.5)*( der_scaly(i, j, k, uux)
         + der_scalx(i, j, k, uuy) );
   S[0][2] = real(.5)*( der_scalz(i, j, k, uux)
         + der_scalx(i, j, k, uuz) );

   S[1][0] = S[0][1];
   S[1][1] = c23*der_scaly(i, j, k, uuy)
         - c13*( der_scalx(i, j, k, uux)
         + der_scalz(i, j, k, uuz) );
   S[1][2] = real(.5)*( der_scalz(i, j, k, uuy)
         + der_scaly(i, j, k, uuz) );

   S[2][0] = S[0][2];
   S[2][1] = S[1][2];
   S[2][2] = c23*der_scalz(i, j, k, uuz) //replaced from "c23*der_scalz(i, j, k, uux)"! TODO recheck that ddz_uu_z is the correct one
         - c13*( der_scalx(i, j, k, uux)
         + der_scaly(i, j, k, uuy) );

    real res = 0.;
    for (int j=0; j < 3; ++j) {
        for (int i=0; i < 3; ++i) {
            res += S[i][j]*S[i][j];
        }
    }
    return res;
}

static __device__ real
entropy(const int i, const int j, const int k,
        const real* __restrict__ lnrho,
        const real* __restrict__ uux,
        const real* __restrict__ uuy,
        const real* __restrict__ uuz,
        const real* __restrict__ ax,
        const real* __restrict__ ay,
        const real* __restrict__ az,
        const real* __restrict__ entropy_s)
{
    const int idx = IDX_1D(i, j, k);

    // Convective derivative d/dt s = - (u dot nabla) s
    const real conv_deriv = - vec_dot_nabla_scal(i, j, k, uux, uuy, uuz, entropy_s);

    // Inverse pT
    const real inv_pT = get_inv_pT(lnrho[idx], entropy_s[idx]);

    // Scalar Laplacian (nabla^2 lnT)
    const real nabla2_lnT = laplace_scal_lnT(i, j, k, lnrho, entropy_s);

    // (Grad ln T)^2
    const real _dot_grad_lnT2 = dot_grad_lnT2(i, j, k, lnrho, entropy_s);

    //  eta * mu0 * j^2
    const real _eta_mu0_j_dot_j = eta_mu0_j_dot_j(i, j, k, ax, ay, az);

    // 2 * rho * visc * S contract_op S
    const real S_contract_S = contract_S_to_scal(i, j, k, uux, uuy, uuz);
    const real strain_tensor_term = real(2.) * exp(lnrho[idx]) * d_NU * S_contract_S;

    return conv_deriv + inv_pT*(nabla2_lnT + _dot_grad_lnT2 + _eta_mu0_j_dot_j + strain_tensor_term);
}


template <int step_number>
static __global__ void
entropy_kernel( const real* __restrict__ d_lnrho,  //SOURCE
                const real* __restrict__ d_uux,
                const real* __restrict__ d_uuy,
                const real* __restrict__ d_uuz,
                const real* __restrict__ d_ax,
                const real* __restrict__ d_ay,
                const real* __restrict__ d_az,
                const real* __restrict__ d_entropy_s,
                const real dt,
                real* __restrict__ d_entropy_s_dst)
{
    const real alpha[] = {0.0, -0.53125, -1.1851851851851851};
    const real beta[]  = {0.25, 0.88888888888888884, 0.75};

    const int i = threadIdx.x + blockIdx.x*blockDim.x + BOUND_SIZE;
    const int j = threadIdx.y + blockIdx.y*blockDim.y + BOUND_SIZE;
    const int k = threadIdx.z + blockIdx.z*blockDim.z + BOUND_SIZE;

    // Return if out of bounds
    if (i >= d_nx_max || j >= d_ny_max || k >= d_nz_max)
        return;

    const int ijk = IDX_1D(i, j, k);

    const real entropy_result = entropy(i, j, k,
                                        d_lnrho,
                                        d_uux, d_uuy, d_uuz,
                                        d_ax, d_ay, d_az,
                                        d_entropy_s);

    if (step_number == 0) {
       d_entropy_s_dst[ijk] = d_entropy_s[ijk] +
                              beta[step_number] * entropy_result * dt;
    } else if (step_number >= 1){
        d_entropy_s_dst[ijk] = d_entropy_s[ijk] + beta[step_number] *
                               (
                                    alpha[step_number] *
                                    (
                                        d_entropy_s[ijk] - d_entropy_s_dst[ijk]
                                    )
                                    / beta[step_number-1] + //
                                    entropy_result * dt
                                );
    } else {
        d_entropy_s_dst[ijk] = NAN;
    }
}


void
rk3_entropy_step(const int step_number, const Grid* d_grid, Grid* d_grid_dst, const real dt,
                                   const CParamConfig* cparams,
                                   const cudaStream_t stream)
{
    const size_t nx = cparams->nx;
    const size_t ny = cparams->ny;
    const size_t nz = cparams->nz;

    const dim3 bpg((unsigned int) ceil(nx / tpb.x),
                   (unsigned int) ceil(ny / tpb.y),
                   (unsigned int) ceil(nz / tpb.z));

    #define KERNEL_PARAMETER_LIST d_grid->arr[d_grid->LNRHO],\
                                          d_grid->arr[d_grid->UUX],\
                                          d_grid->arr[d_grid->UUY],\
                                          d_grid->arr[d_grid->UUZ],\
                                          d_grid->arr[d_grid->AAX],\
                                          d_grid->arr[d_grid->AAY],\
                                          d_grid->arr[d_grid->AAZ],\
                                          d_grid->arr[d_grid->SS],\
                                          dt,\
                                          d_grid_dst->arr[d_grid->SS]

    if (step_number == 0)
        entropy_kernel<0><<<bpg, tpb, 0, stream>>>(KERNEL_PARAMETER_LIST);
    else if (step_number == 1)
        entropy_kernel<1><<<bpg, tpb, 0, stream>>>(KERNEL_PARAMETER_LIST);
    else if (step_number == 2)
        entropy_kernel<2><<<bpg, tpb, 0, stream>>>(KERNEL_PARAMETER_LIST);
    else
        CRASH("Invalid step number");
    CUDA_ERRCHK_KERNEL();

    #undef KERNEL_PARAMETER_LIST
}
#endif
