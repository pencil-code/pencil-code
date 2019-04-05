

#pragma once

void set_grids_zero(float* lnrho, float* uu_x, float* uu_y, float* uu_z);

void lnrho_init(float* lnrho);

void hydro_init(float* uu_x, float* uu_y, float* uu_z);

void debug_grid_init(float* lnrho, float* uu_x, float* uu_y, float* uu_z, int type);

void debug_clear_grid(float* lnrho, float* uu_x, float* uu_y, float* uu_z);
