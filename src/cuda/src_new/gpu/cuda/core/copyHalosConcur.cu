//                             copyHaloAsync.cu
//                             --------------------

/* Date:   6-Jun-2017
   Author: M. Rheinhardt
   Description: Copier functions for the different "plates" of the halo and the full inner data cube with host-device concurrency.
                Load balance yet to be established.
*/

//C libraries
#include <stdio.h>
#include <stdlib.h>

#include "common/errorhandler.h"

#include "dconsts_core.cuh"
#include "copyHalosConcur.cuh"
#include "errorhandler_cuda.cuh"

cudaStream_t strFront;
const int mxy=mx*my;
static int halo_yz_size=0;
static real *halo_yz=NULL; 

const int bot=0, top=1, tot=2;
int halo_widths_x[3]={nghost,nghost,2*nghost};
int halo_widths_y[3]={nghost,nghost,2*nghost};
int halo_widths_z[3]={nghost,nghost,2*nghost};

/****************************************************************************************************************/
__global__ void setCube(real* d_grid,int offset)
{         
// sets d_grid[linear_index] = -(linear_index+1) in global memory in inner xy halos

        int grid_ind=threadIdx.x + threadIdx.y*mx + blockIdx.x*mxy + offset;
//if (offset==1263 && blockIdx.x==0) printf("grid_ind= %d %d \n",grid_ind);        
      
       	d_grid[grid_ind] = (float) (-grid_ind-1);  
}   
/****************************************************************************************************************/
__global__ void unpackOyzPlates(real* d_grid,real* d_halo_yz)
{
//  unpacks buffer for yz halos in global memory
        
        int halo_ind=threadIdx.x + (threadIdx.y + blockIdx.x*blockDim.y)*d_halo_widths_x[tot];

        int grid_ind=blockIdx.x*mxy + threadIdx.y*mx + threadIdx.x;
        if (threadIdx.x>=d_halo_widths_x[bot]) grid_ind+=mx-d_halo_widths_x[tot];

        d_grid[grid_ind]=d_halo_yz[halo_ind];
}
/****************************************************************************************/
__global__ void packIyzPlates(real* d_grid,real* d_halo_yz)
{
//  packs inner yz halos in buffer d_halo_yz on device
        
        const int halo_ind=threadIdx.x + (threadIdx.y + blockIdx.x*(ny-d_halo_widths_y[tot]))*d_halo_widths_x[tot];

        int grid_ind= blockIdx.x*mxy + threadIdx.y*mx + threadIdx.x;
        if (threadIdx.x>=d_halo_widths_x[bot]) grid_ind+=nx-d_halo_widths_x[tot];
/*
if (blockIdx.x==0){
//if (threadIdx.x==0 && threadIdx.y==0) printf("start_offset,halo_ind,grid_ind= %d %d, %d \n",start_offset,halo_ind,grid_ind);
if (threadIdx.x==0 && threadIdx.y==0) printf("d_halo_widths_x= %d %d, %d \n",d_halo_widths_x[0],d_halo_widths_x[1],d_halo_widths_x[2]);
}
if (blockIdx.x==127-6){
if (threadIdx.x==0&& threadIdx.y==0) printf("halo_ind,grid_ind= %d %d, %d \n",blockIdx.x,halo_ind,grid_ind);        
if (threadIdx.x==5&& threadIdx.y==127-6) printf("halo_ind,grid_ind= %d %d, %d \n",blockIdx.x,halo_ind,grid_ind);        
}*/
        d_halo_yz[halo_ind] = d_grid[grid_ind];
//printf("blockIdx,d_halo_yz[halo_ind],grid_ind,halo_ind= %f %d %d \n",d_halo_yz[halo_ind],grid_ind,halo_ind);        
}
/****************************************************************************************/
__global__ void setInnerAll(real* d_grid, const int start_offset)
{
        int grid_ind= blockIdx.x*mxy + threadIdx.y*mx + threadIdx.x + start_offset;
	//d_grid[grid_ind] = -d_grid[grid_ind];
	d_grid[grid_ind] = -(grid_ind+1);
}
/****************************************************************************************/
void setInnerHalos(const GPUContext & ctx,const int w,bool lfirstGPU,bool llastGPU)
{
// sets d_grid[linear_index] = -(linear_index+1) in global memory in inner xy halos
	dim3 threads(nx,ny,1);
	int offset=(mxy+mx+1)*nghost;
        if (lfirstGPU)
	{
printf("first, d_halo_widths_z, offset= %d %d\n",ctx.d_halo_widths_z[bot],offset);
        // inner front plate
		setCube<<<ctx.d_halo_widths_z[bot],threads>>>(ctx.d_grid.arr[w],offset);
	}
        if (llastGPU)
	{
        // inner back plate 
		offset=mxy*(ctx.d_cparams.mz-ctx.d_halo_widths_z[top]-nghost) + (mx+1)*nghost;
printf("last, d_halo_widths_z, offset= %d %d\n",ctx.d_halo_widths_z[top],offset);
		setCube<<<ctx.d_halo_widths_z[top],threads>>>(ctx.d_grid.arr[w],offset);
	}

        int iza=(lfirstGPU ? ctx.d_halo_widths_z[bot] : 0)+nghost,
            ize=ctx.d_cparams.mz - nghost - (llastGPU ? ctx.d_halo_widths_z[top] : 0);
        offset = iza*mxy + (mx+1)*nghost;

        // inner bottom plate
	threads=dim3(nx,3,1);
   	setCube<<<ize-iza,threads>>>(ctx.d_grid.arr[w],offset);
 
        offset += mx*(ny-halo_widths_y[top]);
        // inner top plate
   	setCube<<<ize-iza,threads>>>(ctx.d_grid.arr[w],offset);

        offset = iza*mxy + (nghost+halo_widths_y[bot])*mx + nghost;

        // inner left plate
//printf("left offset %d\n",offset);
        threads=dim3(halo_widths_x[bot],ny-halo_widths_y[tot],1);
   	setCube<<<ize-iza,threads>>>(ctx.d_grid.arr[w],offset);

        offset += nx-halo_widths_x[top];

//printf("right offset %d\n",offset);
        // inner right plate
        threads=dim3(halo_widths_x[top],ny-halo_widths_y[tot],1);
   	setCube<<<ize-iza,threads>>>(ctx.d_grid.arr[w],offset);
}
/****************************************************************************************/
void setOuterHalos(const GPUContext & ctx,real* h_grid,bool lfirstGPU,bool llastGPU)
{
// sets d_grid[linear_index] = -(linear_index+1) in global memory in outer halos
        int offset=0,i;

        if (lfirstGPU)
        {
//printf("first, offset= %d %d\n",ctx.d_halo_widths_z[bot],offset);
        // outer front plate
                for (i=0; i<mxy*halo_widths_z[bot]; i++) h_grid[i]=-(i+1);
        // outer back plate 
                offset=mxy*(mz-halo_widths_z[top]);
//printf("last, offset= %d %d %d\n",ctx.d_cparams.nz,ctx.d_halo_widths_z[tot],offset);
                for (i=0; i<mxy*ctx.d_halo_widths_z[top]; i++) h_grid[i+offset]=-(i+offset+1);

        // outer bottom plate
        int iz,offset = halo_widths_z[bot]*mxy;
 	for (iz=0; iz<nz; iz++){
                for (i=0; i<mx*halo_widths_y[bot]; i++) h_grid[i+offset]=-(i+offset+1);
		offset+=mxy;
	}
        // outer top plate
        offset = halo_widths_z[bot]*mxy+mx*(my-halo_widths_y[top]);
 	for (iz=0; iz<nz; iz++){
                for (i=0; i<mx*halo_widths_y[bot]; i++) h_grid[i+offset]=-(i+offset+1);
		offset+=mxy;
	}

        // inner left and right plate
//printf("left offset %d\n",offset);
        int iy,offset0 = halo_widths_z[bot]*mxy+mx*halo_widths_y[bot];
 	for (iz=0; iz<nz; iz++){
	        offset=offset0;
 		for (iy=0; iy<ny; iy++){
        	        for (i=0; i<halo_widths_x[bot]; i++){
				 h_grid[i+offset]=-(i+offset+1);
			}
			offset+=mx-halo_widths_x[top];
                	for (i=0; i<halo_widths_x[top]; i++){
				 h_grid[i+offset]=-(i+offset+1);
			}
			offset+=halo_widths_x[top];
		}
		offset0+=mxy;
	}
//printf("right offset %d\n",offset);
	}
}
/****************************************************************************************/
__host__ void initHaloConcur(GPUContext & ctx, bool lfirstGPU, bool llastGPU)
{ 
        for (int i=0; i<3; i++ ){
        	ctx.d_halo_widths_x[i] = halo_widths_x[i];
        	ctx.d_halo_widths_y[i] = halo_widths_y[i];
        	ctx.d_halo_widths_z[i] = halo_widths_z[i];
	}

        if (!lfirstGPU) ctx.d_halo_widths_z[bot]=nghost;
   	if (!llastGPU) ctx.d_halo_widths_z[top]=nghost;

        ctx.d_halo_widths_z[tot]=ctx.d_halo_widths_z[bot]+ctx.d_halo_widths_z[top];

	// size of buffer for yz halos 
        ctx.d_halobuffer_size=halo_widths_x[tot]*(my-halo_widths_y[tot])*(ctx.d_cparams.mz-ctx.d_halo_widths_z[tot])*sizeof(real);
//??
        CUDA_ERRCHK( cudaStreamCreate(&ctx.d_copy_streams[FRONT]) );
        CUDA_ERRCHK( cudaStreamCreate(&ctx.d_copy_streams[BACK]) );
        CUDA_ERRCHK( cudaStreamCreate(&ctx.d_copy_streams[BOT]) );
        CUDA_ERRCHK( cudaStreamCreate(&ctx.d_copy_streams[TOP]) );
        CUDA_ERRCHK( cudaStreamCreate(&ctx.d_copy_streams[LEFTRIGHT]) );

        CUDA_ERRCHK( cudaMemcpyToSymbol(d_halo_widths_x, ctx.d_halo_widths_x, 3*sizeof(int)) );
        CUDA_ERRCHK( cudaMemcpyToSymbol(d_halo_widths_y, ctx.d_halo_widths_y, 3*sizeof(int)) );
        CUDA_ERRCHK( cudaMemcpyToSymbol(d_halo_widths_z, ctx.d_halo_widths_z, 3*sizeof(int)) );
        CUDA_ERRCHK( cudaMalloc(&ctx.d_halobuffer,ctx.d_halobuffer_size));            // buffer for yz halos in device

        halo_yz_size=max(halo_yz_size,ctx.d_halobuffer_size);
}
/****************************************************************************************************************/
__host__ void finalizeCopying()
{
	free(halo_yz);
}
/****************************************************************************************************************/
__host__ void copyOxyPlates(const GPUContext & ctx,GridType w,real *h_grid,bool lfirstGPU,bool llastGPU)
{
//  copies outer xy halos from host to device

//setOuterHalos(ctx,h_grid,lfirstGPU,llastGPU);

        int size;
//printf("in copyOxyPlates");

        // front plate
        if (lfirstGPU){
        	size=mxy*halo_widths_z[bot]*sizeof(real);
//printf("front:size %d \n",size);
        	//!!!cudaHostRegister(h_grid, size, cudaHostRegisterDefault);    // time-critical!
        	cudaMemcpyAsync(ctx.d_grid.arr[w], h_grid, size, cudaMemcpyHostToDevice, ctx.d_copy_streams[FRONT]);
	}
        // back plate
        if (llastGPU){
		int h_offset, d_offset;
	        h_offset=mxy*(mz-halo_widths_z[top]);
                d_offset=mxy*(ctx.d_cparams.mz-halo_widths_z[top]);
        	size=mxy*halo_widths_z[top]*sizeof(real);
//printf("back:size, h_offset, d_offset: %d %d %d\n",size, h_offset, d_offset);
	        //!!!cudaHostRegister(h_grid+h_offset, size, cudaHostRegisterDefault);     // time-critical!
        	cudaMemcpyAsync(ctx.d_grid.arr[w]+d_offset, h_grid+h_offset, size, cudaMemcpyHostToDevice, ctx.d_copy_streams[BACK]);
	}
}
/****************************************************************************************************************/
__host__ void copyOxzPlates(const GPUContext & ctx,int w,real * h_grid,bool lfirstGPU,bool llastGPU)
{
//  copies outer xz halos from host to device

        int size, h_offset, d_offset, h_offset0, d_offset0, i, ia, ie;
//printf("in copyOxzPlates");
        // bottom plate
        size=mx*halo_widths_y[bot]*sizeof(real);
        h_offset0=(lfirstGPU ? halo_widths_z[bot] : ctx.start_idx.z)*mxy;
        d_offset0=(lfirstGPU ? halo_widths_z[bot] : 0)*mxy;
       
        h_offset=h_offset0; d_offset=d_offset0; 
        ia=lfirstGPU ? halo_widths_z[bot] : 0; ie=llastGPU ? ctx.d_cparams.mz-halo_widths_z[top] : ctx.d_cparams.mz;
        for (i=ia;i<ie;i++)
        {
          //!!!cudaHostRegister(h_grid+h_offset, size, cudaHostRegisterDefault);     // time-critical!
          cudaMemcpyAsync(ctx.d_grid.arr[w]+d_offset, h_grid+h_offset, size, cudaMemcpyHostToDevice, ctx.d_copy_streams[BOT]);
          d_offset+=mxy; h_offset+=mxy;
        }

        // top plate
        size=mx*halo_widths_y[top]*sizeof(real);
        h_offset=h_offset0 + mx*(my-halo_widths_y[top]);
        d_offset=d_offset0 + mx*(my-halo_widths_y[top]);
        for (i=ia;i<ie;i++)
        {
          //!!!cudaHostRegister(h_grid+h_offset, size, cudaHostRegisterDefault);      // time-critical!
          cudaMemcpyAsync(ctx.d_grid.arr[w]+d_offset, h_grid+h_offset, size, cudaMemcpyHostToDevice, ctx.d_copy_streams[TOP]);
          d_offset+=mxy; h_offset+=mxy;
        }
}
/****************************************************************************************************************/
__host__ void copyOyzPlates(const GPUContext & ctx,GridType w,real *h_grid,bool lfirstGPU,bool llastGPU)
{
//  copies outer yz halos from host to device: they are first packed into the buffer halo_yz, which is then copied 
//  into device buffer d_halo_yz, finally unpacked on device.
        const int x_inc=mx-halo_widths_x[top];

        int i,j,k,iza,ize;
        int h_offset, h_offset0;
        int halo_ind=0;
        const int lsize=halo_widths_x[bot]*sizeof(real), rsize=halo_widths_x[top]*sizeof(real);

        iza = lfirstGPU ? halo_widths_z[bot] : 0;
        ize = ctx.d_cparams.mz - (llastGPU ? halo_widths_z[top] : 0);

        h_offset0=((lfirstGPU ? halo_widths_z[bot] : 0) + ctx.start_idx.z)*mxy + halo_widths_y[bot]*mx;
        h_offset=h_offset0;
//goto cont; //!!!
// Packing on host side: time-critical!

        for (i=iza;i<ize;i++)
        {
                for (j=0;j<my-halo_widths_y[tot];j++)
                {
                        // left plate
                        //cudaMemcpyAsync(halo_yz+halo_ind,h_grid+h_offset,lsize,cudaMemcpyHostToHost, ctx.d_copy_streams[LEFTRIGHT]); 
                        //for (k=0;k<halo_widths_x[bot];k++) *(halo_yz+halo_ind+k)=*(h_grid+h_offset+k);
                        memcpy(halo_yz+halo_ind,h_grid+h_offset,lsize);
                        halo_ind+=halo_widths_x[bot];
                        h_offset+=x_inc;
                        // right plate
                        //cudaMemcpyAsync(halo_yz+halo_ind,h_grid+h_offset,rsize,cudaMemcpyHostToHost, ctx.d_copy_streams[LEFTRIGHT]);
                        //for (k=0;k<halo_widths_x[top];k++) *(halo_yz+halo_ind+k)=*(h_grid+h_offset+k);
                        memcpy(halo_yz+halo_ind,h_grid+h_offset,rsize);
                        halo_ind+=halo_widths_x[top];
                        h_offset+=halo_widths_x[top];
                }
                h_offset=h_offset0+mxy;
        }
cont:
//printf("halo_yz_size,ctx.d_halobuffer_size= %d %d \n",halo_yz_size,ctx.d_halobuffer_size);
        cudaMemcpyAsync(ctx.d_halobuffer, halo_yz, ctx.d_halobuffer_size, cudaMemcpyHostToDevice, ctx.d_copy_streams[LEFTRIGHT]);

//  unpacking in global memory; done by GPU kernel in stream LEFTRIGHT
        int numBlocks=ize-iza;
        dim3 threads(halo_widths_x[tot],my-halo_widths_y[tot],1);    // 2*nghost*ny  needs to be <=1024 !!!
        int d_offset0=(lfirstGPU ? halo_widths_z[bot] : 0)*mxy + halo_widths_y[bot]*mx;

//printf("halo_yz(0:2)= %f, %f, %f, \n",*(halo_yz),*(halo_yz+1),*(halo_yz+2));
        unpackOyzPlates<<<numBlocks,threads,0,ctx.d_copy_streams[LEFTRIGHT]>>>(ctx.d_grid.arr[w]+d_offset0,ctx.d_halobuffer);
        cudaDeviceSynchronize();

/*float buf[3];
offset=mxy*nghost+mx*nghost;
cudaMemcpy(&buf,d_grid+offset,3*sizeof(real),cudaMemcpyDeviceToHost);
printf("buf(0:2)= %f, %f, %f, \n",buf[0],buf[1],buf[2]); */
}
/****************************************************************************************************************/
__host__ void synchronizeStreams(const GPUContext & ctx,bool lfirstGPU,bool llastGPU)
{
//  after copy of halos: synchronizes streams.

        if (lfirstGPU) cudaStreamSynchronize(ctx.d_copy_streams[FRONT]);
        if (llastGPU) cudaStreamSynchronize(ctx.d_copy_streams[BACK]);
	cudaStreamSynchronize(ctx.d_copy_streams[BOT]);
	cudaStreamSynchronize(ctx.d_copy_streams[TOP]);
        cudaStreamSynchronize(ctx.d_copy_streams[LEFTRIGHT]);
}
/****************************************************************************************************************/
__host__ void unlockHostMemOuter(const GPUContext & ctx,real* h_grid,bool lfirstGPU,bool llastGPU)
{
//  after copy of outer halos: releases pinned memory.
        int h_offset0, h_offset, i, ia, ie;

     	// front and back plates
        if (lfirstGPU) cudaHostUnregister(h_grid);

        if (llastGPU){
                h_offset=mxy*(mz-halo_widths_z[top]);
		cudaHostUnregister(h_grid+h_offset);
	}

        h_offset0=(lfirstGPU ? halo_widths_z[bot] : ctx.start_idx.z)*mxy;
        ia=lfirstGPU ? halo_widths_z[bot] : 0; ie=llastGPU ? ctx.d_cparams.mz-halo_widths_z[top] : ctx.d_cparams.mz;

        // outer bottom plate
        h_offset=h_offset0;
        for (i=ia;i<ie;i++)
        {
          cudaHostUnregister(h_grid+h_offset);
          h_offset+=mxy;
        }

        // outer top plate
        h_offset=h_offset0 + mx*(my-halo_widths_y[top]);
        for (i=ia;i<ie;i++)
        {
          cudaHostUnregister(h_grid+h_offset);
          h_offset+=mxy;
        }
}
/****************************************************************************************************************/
__host__ void lockHostMemyz(){
        if (halo_yz==NULL) halo_yz=(real*) malloc(halo_yz_size*sizeof(real));          // buffer for yz halos in host
        cudaHostRegister(halo_yz, halo_yz_size, cudaHostRegisterDefault);
}
/****************************************************************************************************************/
__host__ void unlockHostMemyz(){
        cudaHostUnregister(halo_yz);
}
/****************************************************************************************************************/
__host__ void unlockHostMemInner(const GPUContext & ctx,real* h_grid,bool lfirstGPU,bool llastGPU)
{
//  after copy of inner halos: and releases pinned memory

        int i, h_offset;

        if (lfirstGPU) {
        // inner front plate
                h_offset=(mxy+mx+1)*nghost;
                for (i=0;i<halo_widths_z[bot];i++)
                {
                  cudaHostUnregister(h_grid+h_offset);
          	  h_offset+=mxy;
        	}
        }

        if (llastGPU) {
        // inner back plate
                h_offset=mxy*(mz-nghost-halo_widths_z[top])+(mx+1)*nghost;
                for (i=0;i<halo_widths_z[top];i++)
                {
                  cudaHostUnregister(h_grid+h_offset);
          	  h_offset+=mxy;
        	}
        }

        int iza=(lfirstGPU ? ctx.d_halo_widths_z[bot] : 0)+nghost,
            ize=ctx.d_cparams.mz - nghost - (llastGPU ? ctx.d_halo_widths_z[top] : 0);
        int h_offset0 = (iza+ctx.start_idx.z)*mxy + (mx+1)*nghost;

        // inner bottom plate
        h_offset=h_offset0;
        for (i=iza;i<ize;i++)
        {
          cudaHostUnregister(h_grid+h_offset);
          h_offset+=mxy;
        }
        // inner top plate
        h_offset = h_offset0+mx*(ny-halo_widths_y[top]);
        for (i=iza;i<ize;i++)
        {
          cudaHostUnregister(h_grid+h_offset);
          h_offset+=mxy;
        }
}
/****************************************************************************************************************/
__host__ void copyIxyPlates(const GPUContext & ctx,GridType w,real *h_grid,bool lfirstGPU,bool llastGPU)    // or kernel?
{
//  copies inner xy halos from device to host

        const size_t px=mx*sizeof(real);
        const size_t sx=nx*sizeof(real);

        int i, h_offset;
//h_offset=(mxy+mx)*nghost;
//dim3 threads(mx,my,1);
//setInnerAll<<<nz,threads>>>(ctx.d_grid.arr[w],h_offset); 
//setInnerHalos(ctx,w,lfirstGPU,llastGPU);
	if (lfirstGPU) {
        // inner front plate
	        h_offset=(mxy+mx+1)*nghost;   
        	for (i=0;i<halo_widths_z[bot];i++)
	        {
        	  //!!!cudaHostRegister(h_grid+h_offset, px*ny, cudaHostRegisterDefault);
	          cudaMemcpy2DAsync(h_grid+h_offset, px, ctx.d_grid.arr[w]+h_offset, px, sx, ny, 
                                    cudaMemcpyDeviceToHost, ctx.d_copy_streams[FRONT]);
	          //cudaMemcpy(h_grid+h_offset,ctx.d_grid.arr[w]+h_offset, sizeof(real)*mxy*nghost,cudaMemcpyDeviceToHost); 
//,printf("firstGPU: i, host_offset, dev_offset= %d %d %d \n", i, h_grid+h_offset, ctx.d_grid.arr[w]+h_offset);
        	  h_offset+=mxy;
 	        }
        }

	if (llastGPU) {
        // inner back plate
	        h_offset=mxy*(mz-nghost-halo_widths_z[top])+(mx+1)*nghost;
        	int d_offset=mxy*(ctx.d_cparams.mz-nghost-halo_widths_z[top])+(mx+1)*nghost;;
        	for (i=0;i<halo_widths_z[top];i++)
	        {
        	  //!!!cudaHostRegister(h_grid+h_offset, px*ny, cudaHostRegisterDefault);
	          cudaMemcpy2DAsync(h_grid+h_offset, px, ctx.d_grid.arr[w]+d_offset, px, sx, ny, 
                                    cudaMemcpyDeviceToHost, ctx.d_copy_streams[BACK]);
//printf("lastGPU: i, host_offset, dev_offset= %d %d %d \n", i, h_grid+h_offset, ctx.d_grid.arr[w]+d_offset);
        	  h_offset+=mxy; d_offset+=mxy;
	        }
        }
}
/****************************************************************************************/
__host__ void copyIxzPlates(const GPUContext & ctx,GridType w,real *h_grid,bool lfirstGPU,bool llastGPU)    // or __global__?
{
//  copies inner xz halos from device to host

        const int px=mx*sizeof(real);
        const int sx=nx*sizeof(real);

        int iza=(lfirstGPU ? ctx.d_halo_widths_z[bot] : 0)+nghost,
            ize=ctx.d_cparams.mz - nghost - (llastGPU ? ctx.d_halo_widths_z[top] : 0);
        int h_offset0 = (iza+ctx.start_idx.z)*mxy + (mx+1)*nghost;
        int d_offset0 = iza*mxy + (mx+1)*nghost;

        // inner bottom plate
        int i;
        int h_offset=h_offset0, d_offset=d_offset0;

        for (i=iza;i<ize;i++)
        {
          //!!!cudaHostRegister(h_grid+h_offset, px*halo_widths_y[bot], cudaHostRegisterDefault);
          cudaMemcpy2DAsync(h_grid+h_offset, px, ctx.d_grid.arr[w]+d_offset, px, sx, halo_widths_y[bot], 
                            cudaMemcpyDeviceToHost,ctx.d_copy_streams[BOT]);
       	  h_offset+=mxy; d_offset+=mxy;
        }
        // inner top plate
        h_offset = h_offset0+mx*(ny-halo_widths_y[top]);
        d_offset = d_offset0+mx*(ny-halo_widths_y[top]);
        for (i=iza;i<ize;i++)
        {
          //!!!cudaHostRegister(h_grid+h_offset, px*halo_widths_y[top], cudaHostRegisterDefault);
          cudaMemcpy2DAsync(h_grid+h_offset, px, ctx.d_grid.arr[w]+d_offset, px, sx, halo_widths_y[top],
                            cudaMemcpyDeviceToHost, ctx.d_copy_streams[TOP]);
       	  h_offset+=mxy; d_offset+=mxy;
        }
}
/****************************************************************************************/
__host__ void copyIyzPlates(const GPUContext & ctx,GridType w,real *h_grid,bool lfirstGPU,bool llastGPU)
{
//  copies inner yz halos from device to host: they are first packed into the buffer d_halo_yz, which is then copied 
//  into host buffer halo_yz, finally unpacked on host.

        //d_halo_yz has to have at least size (2*nghost)*(ny-2*nghost)*(nz-2*nghost).
        const int lsize=halo_widths_x[bot]*sizeof(real), rsize=halo_widths_x[top]*sizeof(real);

        const int iza=(lfirstGPU ? ctx.d_halo_widths_z[bot] : 0)+nghost,
                  ize=ctx.d_cparams.mz - nghost - (llastGPU ? ctx.d_halo_widths_z[top] : 0);

        const int nzu=ize-iza,
                  nyu=ny-halo_widths_y[tot],
                  nxu=halo_widths_x[tot];

        const int halo_size=nxu*nyu*nzu*sizeof(real);
        const int x_inc=nx-halo_widths_x[top];

        int d_offset = iza*mxy + (nghost+halo_widths_y[bot])*mx + nghost;
//printf("d_offset= %d\n", d_offset);
        dim3 threads(nxu,nyu,1);
        packIyzPlates<<<nzu,threads,0,ctx.d_copy_streams[LEFTRIGHT]>>>(ctx.d_grid.arr[w]+d_offset,ctx.d_halobuffer);
        cudaMemcpyAsync(halo_yz, ctx.d_halobuffer, halo_size, cudaMemcpyDeviceToHost,ctx.d_copy_streams[LEFTRIGHT]);
        //cudaMemcpy(halo_yz, ctx.d_halobuffer, halo_size, cudaMemcpyDeviceToHost);
        cudaDeviceSynchronize();
//printf("halo_yz \n");
//for (int i=0; i<nxu*nyu*nzu; i++) printf("%f ",halo_yz[i]);   //*(halo_yz+i));
//printf("\n");
     	// unpack on host side
        int i,j,k,halo_ind=0,h_offset0=d_offset + ctx.start_idx.z*mxy,h_offset;

        for (i=iza;i<ize;i++)
        {
                h_offset=h_offset0;
                for (j=0;j<nyu;j++)
                {
                        // inner left plate
                        //cudaMemcpyAsync(h_grid+h_offset,halo_yz+halo_ind,lsize,cudaMemcpyHostToHost,ctx.d_copy_streams[LEFTRIGHT]);
                        //for (k=0;k<halo_widths_x[bot];k++) *(h_grid+h_offset+k)=*(halo_yz+halo_ind+k);
                        memcpy(h_grid+h_offset,halo_yz+halo_ind,lsize);
                        halo_ind+=halo_widths_x[bot];
                        h_offset+=x_inc;

                        // inner right plate
                        //cudaMemcpyAsync(h_grid+h_offset,halo_yz+halo_ind,rsize,cudaMemcpyHostToHost,ctx.d_copy_streams[LEFTRIGHT]);
                        //for (k=0;k<halo_widths_x[top];k++) *(h_grid+h_offset+k)=*(halo_yz+halo_ind+k);
                        memcpy(h_grid+h_offset,halo_yz+halo_ind,rsize);
                        halo_ind+=halo_widths_x[top];
                        h_offset+=2*nghost+halo_widths_x[top];
                }
                h_offset0=h_offset0+mxy;
        }
}
/****************************************************************************************************************/
/*__global__ void setIxzPlates(real* d_grid, int mx, int mxy, int ny, int nghost)
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
/*__global__ void setIyzPlates(real* d_grid,int mx,int nx,int mxy,int nghost)
{
// sets d_grid[linear_index] = -(linear_index+1) in global memory in inner yz halos

        const int start_offset=((mxy+mx)*2+1)*nghost;

        int grid_ind=start_offset + threadIdx.z*mxy + threadIdx.y*mx + threadIdx.x;
        d_grid[grid_ind] = (float)(-grid_ind-1);
        
        grid_ind+=nx-nghost;
        d_grid[grid_ind] = (float)(-grid_ind-1);
}*/
//Headers
#include "../cdata_c.h"
/****************************************************************************************************************/
__host__ void initializeCopying()
{ 
//printf("lperi= %d %d %d \n", lperi[0],lperi[1],lperi[2]);
	// halo widths for undivided data cube
        if (!lperi[0]){
                if (lfirst_proc_x) halo_widths_x[bot]=nghost+1;
                if (llast_proc_x) halo_widths_x[top]=nghost+1;
        }
        if (!lyinyang) {
                if (!lperi[1]){
                        if (lfirst_proc_y) halo_widths_y[bot]=nghost+1;
                        if (llast_proc_y) halo_widths_y[top]=nghost+1;
                }
                if (!lperi[2]){
                        if (lfirst_proc_z) halo_widths_z[bot]=nghost+1;
                        if (llast_proc_z) halo_widths_z[top]=nghost+1;
                }
        }
//printf("halo_widths_x= %d %d\n",halo_widths_x[bot],halo_widths_x[top]);
//printf("halo_widths_y= %d %d\n",halo_widths_y[bot],halo_widths_y[top]);
//printf("halo_widths_z= %d %d\n",halo_widths_z[bot],halo_widths_z[top]);
        halo_widths_x[tot]=halo_widths_x[bot]+halo_widths_x[top];
        halo_widths_y[tot]=halo_widths_y[bot]+halo_widths_y[top];
        halo_widths_z[tot]=halo_widths_z[bot]+halo_widths_z[top];
}
/****************************************************************************************************************/
//  copies all inner halos from device to host

/* for testing: sets elements of inner halo to their negative linear index -1.

        dim3 threadsxy(nx,ny,nghost);
        setIxyPlates<<<1,threadsxy>>>(d_grid, mx, mxy, nz, nghost);

        dim3 threadsxz(nx,nghost,nz-2*nghost);
        setIxzPlates<<<1,threadsxz>>>(d_grid, mx, mxy, ny, nghost);

        dim3 threadsyz(nghost,ny-2*nghost,nz-2*nghost);
        setIyzPlates<<<1,threadsyz>>>(d_grid, mx, nx, mxy, nghost);
*/
/****************************************************************************************************************/
__host__ void copyInnerAll(real* grid, real* d_grid)
{
// copies the full inner data cube from device to host

        size_t px=mx*sizeof(real);
        size_t sx=nx*sizeof(real);
	const long offset=mxy*nghost;
        long offset_data=offset+(mx+1)*nghost;

        cudaHostRegister(grid+offset,mxy*nz*sizeof(real),cudaHostRegisterDefault);
        for (int nn=0;nn<nz;nn++) {
        	cudaMemcpy2DAsync( grid+offset_data, px, d_grid+offset_data, px, sx, ny, cudaMemcpyDeviceToHost, strFront);
    		offset_data+=mxy;
	}
        cudaStreamSynchronize(strFront);
        cudaHostUnregister(grid+offset);
}
/****************************************************************************************************************/
