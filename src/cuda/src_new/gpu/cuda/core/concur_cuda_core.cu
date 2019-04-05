#include "concur_cuda_core.cuh"

void init_concur_ctx(ConcurContext* ctx)
{
    int low_prio, high_prio;
    cudaDeviceGetStreamPriorityRange(&low_prio, &high_prio);
    for (int i=0; i < NUM_STREAMS; ++i)
        cudaStreamCreateWithPriority(&ctx->streams[(StreamName)i], cudaStreamDefault, high_prio + i);

    for (int i=0; i < NUM_EVENTS; ++i)
        cudaEventCreate(&ctx->events[(EventName)i]);
}

void destroy_concur_ctx(ConcurContext* ctx)
{
    for (int i=0; i < NUM_STREAMS; ++i)
        cudaStreamDestroy(ctx->streams[(StreamName)i]);

    for (int i=0; i < NUM_EVENTS; ++i)
        cudaEventDestroy(ctx->events[(EventName)i]);
}
