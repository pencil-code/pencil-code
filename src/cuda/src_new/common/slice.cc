#include "errorhandler.h"

#define INCLUDED_FROM_SLICE_NAME_DEFINER
#include "slice.h"

void slice_malloc(Slice* s, CParamConfig *cparams, RunConfig *run_params)
{
    if (run_params->slice_axis != 'z') CRASH("Slice axis other that z not yet supported!");
    const size_t slice_size_bytes = sizeof(real) * cparams->mx * cparams->my;

    for (int i=0; i < NUM_SLICES; ++i) {
        //s->arr[i] = (real*) malloc(slice_size_bytes);
        s->arr[i] = (real*) calloc(sizeof(real), cparams->mx*cparams->my);
        if (s->arr[i] == NULL) CRASH("Malloc fail");
    }
}


void slice_free(Slice* s)
{
    for (int i=0; i < NUM_SLICES; ++i) {
        free(s->arr[i]); s->arr[i] = NULL;
    } 
}
