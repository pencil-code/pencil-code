

#pragma once


void get_slice_cuda( char slice_axis, float* d_slice_lnrho, float* d_slice_uu, 
		     float* d_slice_uu_x, float* d_slice_uu_y, float* d_slice_uu_z, 
		     float* d_lnrho, float* d_uu_x, float* d_uu_y, float* d_uu_z);
