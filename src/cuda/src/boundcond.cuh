

#pragma once

void boundcond_cuda(float* d_lnrho, float* d_uu_x, float* d_uu_y, float* d_uu_z);

void periodic_boundcond_scal_cuda(float* d_scal);
