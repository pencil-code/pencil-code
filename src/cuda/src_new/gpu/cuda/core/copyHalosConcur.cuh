#include "gpu/cuda/cuda_generic.cuh"

void initializeCopying();
void finalizeCopying();
void initHaloConcur(GPUContext & ctx, bool lfirstGPU, bool llastGPU);
void copyOxyPlates(const GPUContext & ctx,GridType w,real *h_grid,bool lfirstGPU,bool llastGPU);
void copyOxzPlates(const GPUContext & ctx,GridType w,real *h_grid,bool lfirstGPU,bool llastGPU);
void copyOyzPlates(const GPUContext & ctx,GridType w,real *h_grid,bool lfirstGPU,bool llastGPU);
void synchronizeStreams(const GPUContext & ctx,bool lfirstGPU,bool llastGPU);
void unlockHostMemOuter(const GPUContext & ctx,real* h_grid,bool lfirstGPU,bool llastGPU);
void unlockHostMemInner(const GPUContext & ctx,real* h_grid,bool lfirstGPU,bool llastGPU);
void unlockHostMemyz();
void lockHostMemyz();

void copyIxyPlates(const GPUContext & ctx,GridType w,real *h_grid,bool lfirstGPU,bool llastGPU);
void copyIxzPlates(const GPUContext & ctx,GridType w,real *h_grid,bool lfirstGPU,bool llastGPU);
void copyIyzPlates(const GPUContext & ctx,GridType w,real *h_grid,bool lfirstGPU,bool llastGPU);

void copyInnerAll(real* grid, real* d_grid);

