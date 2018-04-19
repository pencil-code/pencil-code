#pragma once

typedef enum {
    STREAM_GLOBAL = 0,
    STREAM_LOCAL_INDUCT,
    STREAM_LOCAL_HYDRO,
    NUM_STREAMS
} StreamName;

typedef enum {
    EVENT_LOCAL_BC_DONE = 0,
    NUM_EVENTS
} EventName;

typedef struct {
    cudaStream_t streams[NUM_STREAMS];
    cudaEvent_t events[NUM_EVENTS];
} ConcurContext;

void init_concur_ctx(ConcurContext* ctx);
void destroy_concur_ctx(ConcurContext* ctx);

