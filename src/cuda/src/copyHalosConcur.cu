//                             copyHaloAsync.cu
//                             --------------------

/* Date:   6-Jun-2017
   Author: M. Rheinhardt
   Description: Copier functions for the different "plates" of the halo and th efull inner data cube with host-device concurrency.
                Load balance yet to be established.
*/

//C libraries
#include <stdio.h>
#include <stdlib.h>

#include "../cparam_c.h"
#include "defines_dims_PC.h"
#define EXTERN extern
#include "dconsts.cuh"

static cudaStream_t strFront=NULL, strBack=NULL, strBot=NULL, strTop=NULL, strLeftRight=NULL;
static int mxy;
static int halo_yz_size;
static float *halo_yz, *d_halo_yz; 

/****************************************************************************************************************/
__global__ void unpackOyzPlates(float* d_grid,float* d_halo_yz)
{
//  unpacks buffer for yz halos in global memory
        
        int halo_ind=threadIdx.x + (threadIdx.y + blockIdx.x*d_COMP_DOMAIN_SIZE_Y)*(2*d_BOUND_SIZE), grid_ind;
        const int start_offset=(d_GRID_Z_OFFSET+d_NX)*d_BOUND_SIZE;

        grid_ind=start_offset + blockIdx.x*d_GRID_Z_OFFSET + threadIdx.y*d_NX + threadIdx.x;
        if (threadIdx.x>=d_BOUND_SIZE) grid_ind+=d_COMP_DOMAIN_SIZE_X;
/*if (threadIdx.y>120 && blockIdx.x==0) printf("threadIdx.y,halo_ind,grid_ind= %d, %d, %d \n",threadIdx.y,halo_ind,grid_ind);        
if (blockIdx.x==0){
if (threadIdx.x==0&& threadIdx.y==0) printf("halo_ind,grid_ind= %d %d, %d \n",blockIdx.x,halo_ind,grid_ind);        
if (threadIdx.x==5&& threadIdx.y==127) printf("halo_ind,grid_ind= %d %d, %d \n",blockIdx.x,halo_ind,grid_ind);        
}
if (blockIdx.x==127){
if (threadIdx.x==0&& threadIdx.y==0) printf("halo_ind,grid_ind= %d %d, %d \n",blockIdx.x,halo_ind,grid_ind);        
if (threadIdx.x==5&& threadIdx.y==127) printf("halo_ind,grid_ind= %d %d, %d \n",blockIdx.x,halo_ind,grid_ind);        
}*/
        d_grid[grid_ind]=d_halo_yz[halo_ind];
}
/****************************************************************************************/
__global__ void packIyzPlates(float* d_grid,float* d_halo_yz)
{
//  packs inner yz halos in buffer d_halo_yz on device

        const int halo_ind=threadIdx.x + (threadIdx.y + blockIdx.x*(d_COMP_DOMAIN_SIZE_Y-2*d_BOUND_SIZE))*(2*d_BOUND_SIZE);
        const int start_offset=((d_GRID_Z_OFFSET+d_NX)*2+1)*d_BOUND_SIZE;

        int grid_ind=start_offset + blockIdx.x*d_GRID_Z_OFFSET + threadIdx.y*d_NX + threadIdx.x;
        if (threadIdx.x>=d_BOUND_SIZE) grid_ind+=d_COMP_DOMAIN_SIZE_X-2*d_BOUND_SIZE;

/*if (blockIdx.x==0){
if (threadIdx.x==0&& threadIdx.y==0) printf("halo_ind,grid_ind= %d %d, %d \n",blockIdx.x,halo_ind,grid_ind);        
if (threadIdx.x==5&& threadIdx.y==127-6) printf("halo_ind,grid_ind= %d %d, %d \n",blockIdx.x,halo_ind,grid_ind);        
}
if (blockIdx.x==127-6){
if (threadIdx.x==0&& threadIdx.y==0) printf("halo_ind,grid_ind= %d %d, %d \n",blockIdx.x,halo_ind,grid_ind);        
if (threadIdx.x==5&& threadIdx.y==127-6) printf("halo_ind,grid_ind= %d %d, %d \n",blockIdx.x,halo_ind,grid_ind);        
}*/
        d_halo_yz[halo_ind] = d_grid[grid_ind];
}
/****************************************************************************************/
//Headers
#include "../cdata_c.h"
//using namespace PC;
/****************************************************************************************************************/
__host__ void initializeCopying()
{ 
        mxy=mx*my;
        halo_yz_size=2*nghost*ny*nz*sizeof(float);      // size of buffer for yz halos

        cudaMalloc(&d_halo_yz,halo_yz_size);            // buffer for yz halos in device
        halo_yz=(float*) malloc(halo_yz_size);          // buffer for yz halos in host
 
        cudaStreamCreate(&strFront);
        cudaStreamCreate(&strBack);
        cudaStreamCreate(&strBot);
        cudaStreamCreate(&strTop);
        cudaStreamCreate(&strLeftRight);
}
/****************************************************************************************************************/
__host__ void finalizeCopying()
{
        cudaFree(&d_halo_yz);
        free(halo_yz);

        cudaStreamDestroy(strFront);
        cudaStreamDestroy(strBack);
        cudaStreamDestroy(strBot);
        cudaStreamDestroy(strTop);
        cudaStreamDestroy(strLeftRight);
}
/****************************************************************************************************************/
__host__ void copyOxyPlates(float* grid, float* d_grid)
{
//  copies outer xy halos from host to device

        const int size=mxy*nghost*sizeof(float);
        const int offset=mxy*(mz-nghost);

        // front plate
        cudaHostRegister(grid, size, cudaHostRegisterDefault);
        cudaMemcpyAsync(d_grid, grid, size, cudaMemcpyHostToDevice, strFront);

        // back plate
        cudaHostRegister(grid+offset, size, cudaHostRegisterDefault);
        cudaMemcpyAsync(d_grid+offset, grid+offset, size, cudaMemcpyHostToDevice, strBack);
}
/****************************************************************************************************************/
__host__ void copyOxzPlates(float* grid, float* d_grid)
{
//  copies outer xz halos from host to device

        const int size=mx*nghost*sizeof(float);

        int offset=mxy*nghost;
        int i;

        // bottom plate
        for (i=0;i<nz;i++)
        {
          cudaHostRegister(grid+offset, size, cudaHostRegisterDefault);
          cudaMemcpyAsync(d_grid+offset, grid+offset, size, cudaMemcpyHostToDevice, strBot);
          offset+=mxy;
        }

        // top plate
        offset=mxy*nghost+mx*(my-nghost);
        for (i=0;i<nz;i++)
        {
          cudaHostRegister(grid+offset, size, cudaHostRegisterDefault);
          cudaMemcpyAsync(d_grid+offset, grid+offset, size, cudaMemcpyHostToDevice, strTop);
          offset+=mxy;
        }
}
/****************************************************************************************************************/
__host__ void copyOyzPlates(float* grid, float* d_grid)
{
//  copies outer yz halos from host to device: they are first packed into the buffer halo_yz, which is then copied 
//  into device buffer d_halo_yz, finally unpacked on device.

        const int size=nghost*sizeof(float);
        const int x_inc=mx-nghost;

        int i,j;
        int halo_ind=0;
        int offset=mx*(my+1)*nghost;

        for (i=0;i<nz;i++)
        {
                for (j=0;j<ny;j++)
                {
                        // left plate
                        cudaMemcpy(halo_yz+halo_ind,grid+offset,size,cudaMemcpyHostToHost);  // also async?
                        halo_ind+=nghost;
                        offset+=x_inc;
                        // right plate
                        cudaMemcpy(halo_yz+halo_ind,grid+offset,size,cudaMemcpyHostToHost);  // also async?
                        halo_ind+=nghost;
                        offset+=nghost;
                }
                offset+=2*mx*nghost;
        }
        cudaHostRegister(halo_yz, halo_yz_size, cudaHostRegisterDefault);
        cudaMemcpyAsync(d_halo_yz, halo_yz, halo_yz_size, cudaMemcpyHostToDevice, strLeftRight);

//  unpacking in global memory; done by GPU kernel in stream strLeftRight

        int numBlocks=nz;
        dim3 threads(2*nghost,ny,1);    // 2*nghost*ny  need to be <=1024 !!!
//printf("halo_yz(0:2)= %f, %f, %f, \n",*(halo_yz),*(halo_yz+1),*(halo_yz+2));
        unpackOyzPlates<<<numBlocks,threads,0,strLeftRight>>>(d_grid,d_halo_yz);
        cudaDeviceSynchronize();
/*float buf[3];
offset=mxy*nghost+mx*nghost;
cudaMemcpy(&buf,d_grid+offset,3*sizeof(float),cudaMemcpyDeviceToHost);
printf("buf(0:2)= %f, %f, %f, \n",buf[0],buf[1],buf[2]); */
}
/****************************************************************************************************************/
__host__ void unlockHostMemOuter(float* grid,float* d_grid)
{
//  after copy of outer halos: synchronizes streams and releases pinned memory

     	// front and back plates
        cudaStreamSynchronize(strFront);
	cudaHostUnregister(grid);	

        cudaStreamSynchronize(strBack);
	cudaHostUnregister(grid+mxy*(mz-nghost));

        int offset=mxy*nghost;
        int i;

        // outer bottom plate
	cudaStreamSynchronize(strBot);
        for (i=0;i<nz;i++)
        {
        	cudaHostUnregister(grid+offset);
        	offset+=mxy;
        }
        // outer top plate
	cudaStreamSynchronize(strTop);
        offset=mxy*nghost+mx*(my-nghost);
        for (i=0;i<nz;i++)
        {
        	cudaHostUnregister(grid+offset);
        	offset+=mxy;
        }
	// left & right plates
        cudaStreamSynchronize(strLeftRight);
        cudaHostUnregister(halo_yz);
}
/****************************************************************************************************************/
__host__ void unlockHostMemInner(float* grid,float* d_grid)
{
//  after copy of inner halos: synchronizes streams and releases pinned memory

        int offset=(mxy+mx+1)*nghost;
        int i;

        cudaStreamSynchronize(strFront);

        // inner front plate
        for (i=0;i<nghost;i++)
        {
          cudaHostUnregister(grid+offset);
          offset+=mxy;
        }

        cudaStreamSynchronize(strBack);

        // inner back plate
        offset=mxy*nz+(mx+1)*nghost;
        for (i=0;i<nghost;i++)
        {
          cudaHostUnregister(grid+offset);
          offset+=mxy;
        }

        cudaStreamSynchronize(strBot);

        // inner bottom plate
        offset=(2*mxy+mx+1)*nghost;
        for (i=0;i<nz-2*nghost;i++)
        {
          cudaHostUnregister(grid+offset);
          offset+=mxy;
        }

        cudaStreamSynchronize(strTop);

        // inner top plate
        offset=2*mxy*nghost+mx*ny+nghost;
        for (i=0;i<nz-2*nghost;i++)
        {
          cudaHostUnregister(grid+offset);
          offset+=mxy;
        }

        cudaStreamSynchronize(strLeftRight);
        cudaHostUnregister(halo_yz);
}
/****************************************************************************************************************/
__host__ void copyOuterHalos(float* grid, float* d_grid)
{
//  copies complete outer halo

        copyOxyPlates(grid, d_grid);
        copyOxzPlates(grid, d_grid);
        copyOyzPlates(grid, d_grid);
 	unlockHostMemOuter(grid, d_grid);
}
/****************************************************************************************/
__host__ void copyIxyPlates(float* grid, float* d_grid)    // or kernel?
{
//  copies inner xy halos from device to host

        const size_t px=mx*sizeof(float);
        const size_t sx=nx*sizeof(float);

        int offset=(mxy+mx+1)*nghost;
        int i;

        // inner front plate
        for (i=0;i<nghost;i++)
        {
          cudaHostRegister(grid+offset, px*ny, cudaHostRegisterDefault);
          cudaMemcpy2DAsync(grid+offset, px, d_grid+offset, px, sx, ny, cudaMemcpyDeviceToHost, strFront);
          offset+=mxy;
        }
        // inner back plate
        offset=mxy*nz+(mx+1)*nghost;
        for (i=0;i<nghost;i++)
        {
          cudaHostRegister(grid+offset, px*ny, cudaHostRegisterDefault);
          cudaMemcpy2DAsync(grid+offset, px, d_grid+offset, px, sx, ny, cudaMemcpyDeviceToHost, strBack);
          offset+=mxy;
        }
}
/****************************************************************************************/
__host__ void copyIxzPlates(float* grid, float* d_grid)    // or __global__?
{
//  copies inner xz halos from device to host

        const int px=mx*sizeof(float);
        const int sx=nx*sizeof(float);

        int offset=(2*mxy+mx+1)*nghost;
        int i;

        // inner bottom plate
        for (i=0;i<nz-2*nghost;i++)
        {
          cudaHostRegister(grid+offset, px*nghost, cudaHostRegisterDefault);
          cudaMemcpy2DAsync( grid+offset, px, d_grid+offset, px, sx, nghost, cudaMemcpyDeviceToHost, strBot);
          offset+=mxy;
        }
        // inner top plate
        offset=2*mxy*nghost+mx*ny+nghost;
        for (i=0;i<nz-2*nghost;i++)
        {
          cudaHostRegister(grid+offset, px*nghost, cudaHostRegisterDefault);
          cudaMemcpy2DAsync( grid+offset, px, d_grid+offset, px, sx, nghost, cudaMemcpyDeviceToHost, strTop);
          offset+=mxy;
        }
}
/****************************************************************************************/
__host__ void copyIyzPlates(float* grid, float* d_grid)
{
//  copies inner yz halos from device to host: they are first packed into the buffer d_halo_yz, which is then copied 
//  into host buffer halo_yz, finally unpacked on host.


        //d_halo_yz has to have at least size (2*nghost)*(ny-2*nghost)*(nz-2*nghost).
        const int size=nghost*sizeof(float);
        const int halo_size=2*nghost*(ny-2*nghost)*(nz-2*nghost)*sizeof(float);
        const int x_inc=nx-nghost;

        int i,j;
        int halo_ind=0;
        int offset=((mxy+mx)*2+1)*nghost;
        dim3 threads(2*nghost,ny-2*nghost,1);

        packIyzPlates<<<nz-2*nghost,threads,0,strLeftRight>>>(d_grid,d_halo_yz);
        cudaHostRegister(halo_yz, halo_size, cudaHostRegisterDefault);
        cudaMemcpyAsync(halo_yz, d_halo_yz, halo_size, cudaMemcpyDeviceToHost,strLeftRight);

// unpack on host side

        for (i=0;i<nz-2*nghost;i++)
        {
                for (j=0;j<ny-2*nghost;j++)
                {
                        // inner left plate
                        cudaMemcpyAsync(grid+offset,halo_yz+halo_ind,size,cudaMemcpyHostToHost,strLeftRight);
                        halo_ind+=nghost;
                        offset+=x_inc;
                        // inner right plate
                        cudaMemcpyAsync(grid+offset,halo_yz+halo_ind,size,cudaMemcpyHostToHost,strLeftRight);
                        halo_ind+=nghost;
                        offset+=3*nghost;
                }
                offset+=4*mx*nghost;
        }
}
/****************************************************************************************************************/
/*__global__ void setIxyPlates(float* d_grid, int mx, int mxy, int nz, int nghost)
{
// sets d_grid[linear_index] = -(linear_index+1) in global memory in inner xy halos

        int start_offset=(mxy+mx+1)*nghost;
        int grid_ind=start_offset + threadIdx.x + threadIdx.y*mx + threadIdx.z*mxy;

        // inner front plate
        d_grid[grid_ind] = (float) (-grid_ind-1);

        // inner back plate
        grid_ind += (nz-nghost)*mxy;
        d_grid[grid_ind] = (float) (-grid_ind-1);
}*/
/****************************************************************************************************************/
/*__global__ void setIxzPlates(float* d_grid, int mx, int mxy, int ny, int nghost)
{
// sets d_grid[linear_index] = -(linear_index+1) in global memory in inner xz halos


        int start_offset=(2*mxy+mx+1)*nghost;
        int grid_ind=start_offset + threadIdx.x + threadIdx.y*mx + threadIdx.z*mxy;

        // inner bottom plate
        d_grid[grid_ind] = (float) (-grid_ind-1);

        // inner top plate
        grid_ind += (ny-nghost)*mx;
        d_grid[grid_ind] = (float) (-grid_ind-1);
}*/
/****************************************************************************************/
/*__global__ void setIyzPlates(float* d_grid,int mx,int nx,int mxy,int nghost)
{
// sets d_grid[linear_index] = -(linear_index+1) in global memory in inner yz halos

        const int start_offset=((mxy+mx)*2+1)*nghost;

        int grid_ind=start_offset + threadIdx.z*mxy + threadIdx.y*mx + threadIdx.x;
        d_grid[grid_ind] = (float)(-grid_ind-1);
        
        grid_ind+=nx-nghost;
        d_grid[grid_ind] = (float)(-grid_ind-1);
}*/
/****************************************************************************************************************/
__host__ void copyInnerHalos(float* grid, float* d_grid)
{
//  copies all inner halos from device to host

/* for testing: sets elements of inner halo to their negative linear index -1.

        dim3 threadsxy(nx,ny,nghost);
        setIxyPlates<<<1,threadsxy>>>(d_grid, mx, mxy, nz, nghost);

        dim3 threadsxz(nx,nghost,nz-2*nghost);
        setIxzPlates<<<1,threadsxz>>>(d_grid, mx, mxy, ny, nghost);

        dim3 threadsyz(nghost,ny-2*nghost,nz-2*nghost);
        setIyzPlates<<<1,threadsyz>>>(d_grid, mx, nx, mxy, nghost);
*/
        copyIxyPlates(grid, d_grid);
        copyIxzPlates(grid, d_grid);
        copyIyzPlates(grid, d_grid);
        unlockHostMemInner(grid, d_grid);
}
/****************************************************************************************************************/
__host__ void copyAll(float* grid, float* d_grid)
{
// copies the full data cube from host to device.

 	int size=mxy*mz*sizeof(float);
	cudaHostRegister(grid,size,cudaHostRegisterDefault);
	cudaMemcpy(d_grid, grid, size, cudaMemcpyHostToDevice);
	cudaHostUnregister(grid);
}
/****************************************************************************************************************/
__host__ void copyInnerAll(float* grid, float* d_grid)
{
// copies the full inner data cube from device to host

        size_t px=mx*sizeof(float);
        size_t sx=nx*sizeof(float);
	const long offset=mxy*nghost;
        long offset_data=offset+(mx+1)*nghost;

        cudaHostRegister(grid+offset,mxy*nz*sizeof(float),cudaHostRegisterDefault);
        for (int nn=0;nn<nz;nn++) {
        	cudaMemcpy2DAsync( grid+offset_data, px, d_grid+offset_data, px, sx, ny, cudaMemcpyDeviceToHost, strFront);
    		offset_data+=mxy;
	}
        cudaStreamSynchronize(strFront);
        cudaHostUnregister(grid+offset);
}
/****************************************************************************************************************/
