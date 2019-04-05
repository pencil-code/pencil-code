

#pragma once

void save_grid_information(float time);

int load_grid_data(float* grid, const char* path);

int save_grid_data(float* grid, const char* path);

//template <char slice_axis>//looks like template function call from .cu doesnt like g++ compiled .cpps
void save_anim_slice(char slice_axis, int slice_step, float* slice_lnrho, float* slice_uu, float* slice_uu_x, float* slice_uu_y, float* slice_uu_z);

void init_ts_file();

void save_ts(float t, float dt, int step, float urms, float uxrms, float uyrms, float uzrms, float uxmax, float uymax, float uzmax, float rhorms, float umax, float rhomax,
             float uxmin, float uymin, float uzmin, float rhomin, float umin);

void init_used_k_file(float kaver);

void save_used_k(float kk_vec_x, float kk_vec_y, float kk_vec_z, float phi, float forcing_kk_part_x, float forcing_kk_part_y, float forcing_kk_part_z);

