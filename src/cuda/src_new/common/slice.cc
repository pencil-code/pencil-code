#include "errorhandler.h"
#include "slice.h"

#define SLICE_TYPE(type, name) name
const char* slice_names[] = { SLICE_TYPES };
#undef SLICE_TYPE

void slice_malloc(Slice* s, CParamConfig *cparams, RunConfig *run_params)
{
    if (run_params->slice_axis != 'z') { 
        printf("Got slice axis \"%c\", but slice axis other than z not yet " 
               "supported.\nMake sure that slice_axis is defined correctly in "
               "the local run.conf in the project directory\n"
               "(see <astaroth root dir>/src/configs/run.conf for defaults).\n", 
                run_params->slice_axis); 
                CRASH("Slice axis other that z not yet supported"); 
    }

    for (int i=0; i < NUM_SLICES; ++i) {
        s->arr[i] = (real*) calloc(sizeof(real), cparams->mx*cparams->my);
        if (s->arr[i] == NULL) { CRASH("Malloc fail"); }
    }
}

void slice_free(Slice* s)
{
    for (int i=0; i < NUM_SLICES; ++i) {
        free(s->arr[i]); s->arr[i] = NULL;
    } 
}
