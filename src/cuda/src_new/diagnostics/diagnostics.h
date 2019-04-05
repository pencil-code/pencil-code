#pragma once 
#include "common/grid.h"

// only defined & used if GPU_ASTAROTH is defined

void timeseries_diagnostics(const Grid & h_grid);
void init_diagnostics();

