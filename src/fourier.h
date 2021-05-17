!
!  $Id$  
!
  private
!
  public :: initialize_fourier
  public :: fourier_transform, fourier_transform_xy, fourier_transform_xz
  public :: fourier_transform_x,fourier_transform_y
  public :: fourier_transform_shear, fourier_transform_shear_xy
  public :: fourier_transform_other, fourier_transform_xy_xy
  public :: fourier_transform_y_y
  public :: fourier_shift_yz_y, fourier_shift_y
  public :: fourier_transform_xy_xy_other
  public :: fft_x_parallel, fft_y_parallel, fft_z_parallel
  public :: fft_xy_parallel, fft_xyz_parallel
  public :: setup_extrapol_fact, vect_pot_extrapol_z_parallel, field_extrapol_z_parallel
  public :: fourier_transform_real_1
  public :: fft_xy_parallel_2D_other
!  
  public :: kx_fft, kx_fft2
  public :: ky_fft, ky_fft2
  public :: kz_fft, kz_fft2
!
  real, dimension (nxgrid) :: kx_fft, kx_fft2
  real, dimension (nygrid) :: ky_fft, ky_fft2
  real, dimension (nzgrid) :: kz_fft, kz_fft2
