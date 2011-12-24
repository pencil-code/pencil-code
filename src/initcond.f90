! $Id$
!
!  This module contains code used by the corresponding physics
!  modules to set up various initial conditions (stratitication,
!  perturbations, and other structures). This module is not used
!  during run time (although it is used by the physics modules that
!  are used both during run time and for the initial condition).
!
module Initcond
!
  use Cdata
  use General
  use Messages
!
  implicit none
!
  private
!
  public :: arcade_x, vecpatternxy, bipolar, bipolar_restzero
  public :: soundwave,sinwave,sinwave_phase,coswave,coswave_phase,cos_cos_sin
  public :: hatwave
  public :: acosy
  public :: sph_constb
  public :: gaunoise, posnoise, posnoise_rel
  public :: gaunoise_rprof
  public :: gaussian, gaussian3d, gaussianpos, beltrami, bessel_x, bessel_az_x
  public :: beltrami_complex,bhyperz
  public :: rolls, tor_pert
  public :: jump, bjump, bjumpz, stratification, stratification_x
  public :: modes, modev, modeb, crazy
  public :: trilinear, baroclinic
  public :: triquad, isotdisk
  public :: diffrot, olddiffrot
  public :: const_omega
  public :: powern, power_randomphase
  public :: planet, planet_hc
  public :: random_isotropic_KS
  public :: htanh, htube, htube2, htube_x, hat, hat3d
  public :: htube_erf, xpoint, xpoint2
  public :: wave_uu, wave, parabola, linprof
  public :: sinxsinz, cosx_cosy_cosz, cosx_coscosy_cosz
  public :: x_siny_cosz, x1_siny_cosz, x1_cosy_cosz, lnx_cosy_cosz
  public :: sinx_siny_sinz, cosx_siny_cosz, sinx_siny_cosz
  public :: sin2x_sin2y_cosz, cos2x_cos2y_cos2z, x3_cosy_cosz
  public :: cosx_cosz, cosy_cosz, cosy_sinz
  public :: cosxz_cosz, cosyz_sinz
  public :: halfcos_x, magsupport, vfield
  public :: uniform_x, uniform_y, uniform_z, uniform_phi, phi_comp_over_r
  public :: vfluxlayer, hfluxlayer, hfluxlayer_y
  public :: vortex_2d
  public :: vfield2
  public :: hawley_etal99a
  public :: robertsflow
  public :: const_lou
  public :: corona_init,mdi_init,mag_init,temp_hydrostatic
  public :: innerbox
  public :: couette, couette_rings
  public :: strange,phi_siny_over_r2
  public :: ferriere_uniform_x, ferriere_uniform_y
  public :: rotblob, pre_stellar_cloud
!
  interface posnoise            ! Overload the `posnoise' function
    module procedure posnoise_vect
    module procedure posnoise_scal
  endinterface
!
  interface posnoise_rel        ! Overload the `posnoise' function
    module procedure posnoise_rel_vect
    module procedure posnoise_rel_scal
  endinterface
!
  interface gaunoise            ! Overload the `gaunoise' function
    module procedure gaunoise_vect
    module procedure gaunoise_scal
    module procedure gaunoise_prof_vect
    module procedure gaunoise_prof_scal
  endinterface
!
  interface gaunoise_rprof      ! Overload the `gaunoise_rprof' function
    module procedure gaunoise_rprof_vect
    module procedure gaunoise_rprof_scal
  endinterface
!
  character(LEN=labellen) :: wave_fmt1='(1x,a,4f8.2)'
!
  contains
!***********************************************************************
    subroutine phi_siny_over_r2(f,i)
!
!  A_phi ~ sin(y)/R^2 if R>=0.7 field (in terms of vector potential)
!
!  14-oct-09/irina: coded
!
      integer :: i,j,k
      real, dimension (mx,my,mz,mfarray) :: f
      real :: dmx,radius, dtheta,theta
!
     dmx=(l2-l1)/mx
     dtheta=pi/my
     f(:,:,:,i)=0.0
!
     do j=1,mx
       radius=l1+(j-1)*dmx
     if (radius>=0.7) then
       do k=1,my
        theta=pi*(k-1)*dtheta
        f(j,k,:,i)=sin(theta)/radius**2
       enddo
     endif
    enddo
!
    endsubroutine phi_siny_over_r2
!***********************************************************************
    subroutine sinxsinz(ampl,f,i,kx,ky,kz,KKx,KKy,KKz)
!
!  sinusoidal wave. Note: f(:,:,:,j) with j=i+1 is set.
!
!  26-jul-02/axel: coded
!
      integer :: i,j
      real, dimension (mx,my,mz,mfarray) :: f
      real,optional :: kx,ky,kz,KKx,KKy,KKz
      real :: ampl,kx1=pi/2.,ky1=0.,kz1=pi/2.,KKx1=0.,KKy1=0.,KKz1=0.
!
!  wavenumber k, helicity H=ampl (can be either sign)
!
!  sinx(kx*x)*sin(kz*z)
!
      if (present(kx)) kx1=kx
      if (present(ky)) ky1=ky
      if (present(kz)) kz1=kz
!
!  Gaussian scake heights
!
      if (present(KKx)) KKx1=KKx
      if (present(KKy)) KKy1=KKy
      if (present(KKz)) KKz1=KKz
!
      if (ampl==0) then
        if (lroot) print*,'sinxsinz: ampl=0 in sinx*sinz wave; kx,kz=',kx1,kz1
      else
        if (lroot) print*,'sinxsinz: sinx*sinz wave; ampl,kx,kz=',ampl,kx1,kz1
        j=i+1
        f(:,:,:,j)=f(:,:,:,j)+ampl*(&
          spread(spread(cos(kx1*x)*exp(-.5*(KKx1*x)**2),2,my),3,mz)*&
          spread(spread(cos(ky1*y)*exp(-.5*(KKy1*y)**2),1,mx),3,mz)*&
          spread(spread(cos(kz1*z)*exp(-.5*(KKz1*z)**2),1,mx),2,my))
      endif
!
    endsubroutine sinxsinz
!***********************************************************************
    subroutine sinx_siny_sinz(ampl,f,i,kx,ky,kz)
!
!  sinusoidal wave, adapted from sinxsinz (that routine was already doing
!  this, but under a different name)
!
!  20-jan-07/axel: adapted
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real,optional :: kx,ky,kz
      real :: ampl,kx1=pi/2.,ky1=0.,kz1=pi/2.
!
!  wavenumber k, helicity H=ampl (can be either sign)
!
!  sinx(kx*x)*sin(kz*z)
!
      if (present(kx)) kx1=kx
      if (present(ky)) ky1=ky
      if (present(kz)) kz1=kz
      if (ampl==0) then
        if (lroot) print*,'sinx_siny_sinz: ampl=0'
      else
        if (lroot) write(*,wave_fmt1) 'sinx_siny_sinz: ampl,kx,ky,kz=', &
                                      ampl,kx1,ky1,kz1
        f(:,:,:,i)=f(:,:,:,i)+ampl*(spread(spread(sin(kx1*x),2,my),3,mz)&
                                   *spread(spread(sin(ky1*y),1,mx),3,mz)&
                                   *spread(spread(sin(kz1*z),1,mx),2,my))
      endif
!
    endsubroutine sinx_siny_sinz
!***********************************************************************
    subroutine sinx_siny_cosz(ampl,f,i,kx,ky,kz)
!
!  sinusoidal wave, adapted from sinxsinz (that routine was already doing
!  this, but under a different name)
!
!   2-dec-03/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real,optional :: kx,ky,kz
      real :: ampl,kx1=pi/2.,ky1=0.,kz1=pi/2.
!
!  wavenumber k, helicity H=ampl (can be either sign)
!
!  sinx(kx*x)*sin(kz*z)
!
      if (present(kx)) kx1=kx
      if (present(ky)) ky1=ky
      if (present(kz)) kz1=kz
      if (ampl==0) then
        if (lroot) print*,'sinx_siny_cosz: ampl=0'
      else
        if (lroot) write(*,wave_fmt1) 'sinx_siny_cosz: ampl,kx,ky,kz=', &
                                      ampl,kx1,ky1,kz1
        f(:,:,:,i)=f(:,:,:,i)+ampl*(spread(spread(sin(kx1*x),2,my),3,mz)&
                                   *spread(spread(sin(ky1*y),1,mx),3,mz)&
                                   *spread(spread(cos(kz1*z),1,mx),2,my))
      endif
!
    endsubroutine sinx_siny_cosz
!***********************************************************************
    subroutine x_siny_cosz(ampl,f,i,kx,ky,kz)
!
!  sinusoidal wave, adapted from sinxsinz (that routine was already doing
!  this, but under a different name)
!
!   2-dec-03/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real,optional :: kx,ky,kz
      real :: ampl,kx1=pi/2.,ky1=0.,kz1=pi/2.
!
!  wavenumber k, helicity H=ampl (can be either sign)
!
!  sinx(kx*x)*sin(kz*z)
!
      if (present(kx)) kx1=kx
      if (present(ky)) ky1=ky
      if (present(kz)) kz1=kz
      if (ampl==0) then
        if (lroot) print*,'x_siny_cosz: ampl=0'
      else
        if (lroot) write(*,wave_fmt1) 'x_siny_cosz: ampl,kx,ky,kz=', &
                                      ampl,kx1,ky1,kz1
        f(:,:,:,i)=f(:,:,:,i)+ampl*(spread(spread(   (    x),2,my),3,mz)&
                                   *spread(spread(sin(ky1*y),1,mx),3,mz)&
                                   *spread(spread(cos(kz1*z),1,mx),2,my))
      endif
!
    endsubroutine x_siny_cosz
!***********************************************************************
    subroutine x1_siny_cosz(ampl,f,i,kx,ky,kz,phasey)
!
!  sinusoidal wave, adapted from sinxsinz (that routine was already doing
!  this, but under a different name)
!
!   2-dec-03/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real,optional :: kx,ky,kz,phasey
      real :: ampl,kx1=pi/2.,ky1=0.,kz1=pi/2., phasey1=0.
!
!  wavenumber k, helicity H=ampl (can be either sign)
!
!  sinx(kx*x)*sin(kz*z)
!
      if (present(kx)) kx1=kx
      if (present(ky)) ky1=ky
      if (present(kz)) kz1=kz
      if (present(phasey)) phasey1=phasey
      if (ampl==0) then
        if (lroot) print*,'x1_siny_cosz: ampl=0'
      else
        if (lroot) write(*,wave_fmt1) 'x1_siny_cosz: ampl,kx,ky,kz=', &
                                      ampl,kx1,ky1,kz1
        f(:,:,:,i)=f(:,:,:,i)+ampl*(spread(spread(   ( 1./x),2,my),3,mz)&
                                   *spread(spread(sin(ky1*y+phasey1),1,mx),3,mz)&
                                   *spread(spread(cos(kz1*z),1,mx),2,my))
      endif
!
    endsubroutine x1_siny_cosz
!***********************************************************************
    subroutine x1_cosy_cosz(ampl,f,i,kx,ky,kz)
!
!  sinusoidal wave, adapted from sinxsinz (that routine was already doing
!  this, but under a different name)
!
!   1-jul-07/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real,optional :: kx,ky,kz
      real :: ampl,kx1=pi/2.,ky1=0.,kz1=pi/2.
!
!  wavenumber k, 1/x radial profile
!
      if (present(kx)) kx1=kx
      if (present(ky)) ky1=ky
      if (present(kz)) kz1=kz
      if (ampl==0) then
        if (lroot) print*,'x1_siny_cosz: ampl=0'
      else
        if (lroot) write(*,wave_fmt1) 'x1_siny_cosz: ampl,kx,ky,kz=', &
                                      ampl,kx1,ky1,kz1
        f(:,:,:,i)=f(:,:,:,i)+ampl*(spread(spread(   ( 1./x),2,my),3,mz)&
                                   *spread(spread(cos(ky1*y),1,mx),3,mz)&
                                   *spread(spread(cos(kz1*z),1,mx),2,my))
      endif
!
    endsubroutine x1_cosy_cosz
!***********************************************************************
    subroutine lnx_cosy_cosz(ampl,f,i,kx,ky,kz)
!
!  sinusoidal wave, adapted from sinxsinz (that routine was already doing
!  this, but under a different name)
!
!   2-jul-07/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real,optional :: kx,ky,kz
      real :: ampl,kx1=pi/2.,ky1=0.,kz1=pi/2.
!
!  wavenumber k, 1/x radial profile
!
      if (present(kx)) kx1=kx
      if (present(ky)) ky1=ky
      if (present(kz)) kz1=kz
      if (ampl==0) then
        if (lroot) print*,'lnx_cosy_cosz: ampl=0'
      else
        if (lroot) write(*,wave_fmt1) 'lnx_cosy_cosz: ampl,kx,ky,kz=', &
                                      ampl,kx1,ky1,kz1
        f(:,:,:,i)=f(:,:,:,i)+ampl*(spread(spread(alog(   x),2,my),3,mz)&
                                   *spread(spread(cos(ky1*y),1,mx),3,mz)&
                                   *spread(spread(cos(kz1*z),1,mx),2,my))
      endif
!
    endsubroutine lnx_cosy_cosz
!***********************************************************************
    subroutine cosx_siny_cosz(ampl,f,i,kx,ky,kz)
!
!  sinusoidal wave, adapted from sinxsinz (that routine was already doing
!  this, but under a different name)
!
!  15-mar-07/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real,optional :: kx,ky,kz
      real :: ampl,kx1=pi/2.,ky1=0.,kz1=pi/2.
!
!  wavenumber k, helicity H=ampl (can be either sign)
!
!  sinx(kx*x)*sin(kz*z)
!
      if (present(kx)) kx1=kx
      if (present(ky)) ky1=ky
      if (present(kz)) kz1=kz
      if (ampl==0) then
        if (lroot) print*,'cosx_siny_cosz: ampl=0'
      else
        if (lroot) write(*,wave_fmt1) 'cosx_siny_cosz: ampl,kx,ky,kz=', &
                                      ampl,kx1,ky1,kz1
        f(:,:,:,i)=f(:,:,:,i)+ampl*(spread(spread(cos(kx1*x),2,my),3,mz)&
                                   *spread(spread(sin(ky1*y),1,mx),3,mz)&
                                   *spread(spread(cos(kz1*z),1,mx),2,my))
      endif
!
    endsubroutine cosx_siny_cosz
!***********************************************************************
    subroutine sin2x_sin2y_cosz(ampl,f,i,kx,ky,kz)
!
!  sinusoidal wave, adapted from sinxsinz (that routine was already doing
!  this, but under a different name)
!
!   2-dec-03/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real,optional :: kx,ky,kz
      real :: ampl,kx1=pi/2.,ky1=0.,kz1=pi/2.
!
!  wavenumber k, helicity H=ampl (can be either sign)
!
!  sin2(kx*x)*sin2(kz*z)
!
      if (present(kx)) kx1=kx
      if (present(ky)) ky1=ky
      if (present(kz)) kz1=kz
      if (ampl==0) then
        if (lroot) print*,'sin2x_sin2y_cosz: ampl=0'
      else
        if (lroot) write(*,wave_fmt1) 'sin2x_sin2y_cosz: ampl,kx,ky,kz=',&
                                      ampl,kx1,ky1,kz1
        f(:,:,:,i)=f(:,:,:,i)+ampl*(spread(spread(sin(kx1*x)**2,2,my),3,mz)&
                                   +spread(spread(sin(ky1*y)**2,1,mx),3,mz))&
                                   *spread(spread(cos(kz1*z),1,mx),2,my)
      endif
!
    endsubroutine sin2x_sin2y_cosz
!***********************************************************************
    subroutine cosxz_cosz(ampl,f,i,kx,kz)
!
!  sinusoidal wave, adapted from sinxsinz (that routine was already doing
!  this, but under a different name)
!
!  13-aug-09/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real,optional :: kx,kz
      real :: ampl,kx1=pi/2.,kz1=pi/2.
!
!  wavenumber k, helicity H=ampl (can be either sign)
!
!  sinx(kx*x)*sin(kz*z)
!
      if (present(kx)) kx1=kx
      if (present(kz)) kz1=kz
      if (ampl==0) then
        if (lroot) print*,'cosxz_cosz: ampl=0'
      else
        if (lroot) write(*,wave_fmt1) 'cos(x-cosz): ampl,kx,kz=', &
                                      ampl,kx1,kz1
        f(:,:,:,i)=f(:,:,:,i)+ampl*cos( &
          spread(spread(kx1*x,2,my),3,mz)- &
          spread(spread(cos(kz1*z),1,mx),2,my))
      endif
!
    endsubroutine cosxz_cosz
!***********************************************************************
    subroutine cosyz_sinz(ampl,f,i,ky,kz)
!
!  sinusoidal wave, adapted from sinxsinz (that routine was already doing
!  this, but under a different name)
!
!  13-aug-09/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real,optional :: ky,kz
      real :: ampl,ky1=pi/2.,kz1=pi/2.
!
!  wavenumber k, helicity H=ampl (can be either sign)
!
!  sinx(ky*y)*sin(kz*z)
!
      if (present(ky)) ky1=ky
      if (present(kz)) kz1=kz
      if (ampl==0) then
        if (lroot) print*,'cosyz_sinz: ampl=0'
      else
        if (lroot) write(*,wave_fmt1) 'cos(y-sinz): ampl,ky,kz=', &
                                      ampl,ky1,kz1
        f(:,:,:,i)=f(:,:,:,i)+ampl*cos( &
          spread(spread(ky1*y,1,mx),3,mz)- &
          spread(spread(cos(kz1*z),1,mx),2,my))
      endif
!
    endsubroutine cosyz_sinz
!***********************************************************************
    subroutine cosx_cosy_cosz(ampl,f,i,kx,ky,kz)
!
!  sinusoidal wave, adapted from sinxsinz (that routine was already doing
!  this, but under a different name)
!
!   2-dec-03/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real,optional :: kx,ky,kz
      real :: ampl,kx1=pi/2.,ky1=0.,kz1=pi/2.
!
!  wavenumber k, helicity H=ampl (can be either sign)
!
!  sinx(kx*x)*sin(kz*z)
!
      if (present(kx)) kx1=kx
      if (present(ky)) ky1=ky
      if (present(kz)) kz1=kz
      if (ampl==0) then
        if (lroot) print*,'cosx_cosy_cosz: ampl=0'
      else
        if (lroot) write(*,wave_fmt1) 'cosx_cosy_cosz: ampl,kx,ky,kz=', &
                                      ampl,kx1,ky1,kz1
        f(:,:,:,i)=f(:,:,:,i)+ampl*(spread(spread(cos(kx1*x),2,my),3,mz)&
                                   *spread(spread(cos(ky1*y),1,mx),3,mz)&
                                   *spread(spread(cos(kz1*z),1,mx),2,my))
      endif
!
    endsubroutine cosx_cosy_cosz
!***********************************************************************
    subroutine cosx_cosz(ampl,f,i,kx,kz)
!
!  initial condition for potential field test
!
!  15-aug-09/axel: adapted from cosy_sinz
!
      integer :: i,j
      real, dimension (mx,my,mz,mfarray) :: f
      real,optional :: kx,kz
      real :: ampl,kx1=1.,kz1=pi
!
!  wavenumber k
!
      if (present(kx)) kx1=kx
      if (present(kz)) kz1=kz
!
      j=i; f(:,:,:,j)=f(:,:,:,j)+ampl*spread(spread(cos(kx1*x),2,my),3,mz)&
                                     *spread(spread(cos(kz1*z),1,mx),2,my)
!
    endsubroutine cosx_cosz
!***********************************************************************
    subroutine cosy_cosz(ampl,f,i,ky,kz)
!
!  initial condition for potential field test
!
!   06-oct-06/axel: coded
!   11-oct-06/wolf: modified to only set one component of aa
!
      integer :: i,j
      real, dimension (mx,my,mz,mfarray) :: f
      real,optional :: ky,kz
      real :: ampl,ky1=1.,kz1=pi
!
!  wavenumber k
!
      if (present(ky)) ky1=ky
      if (present(kz)) kz1=kz
!
      j=i; f(:,:,:,j)=f(:,:,:,j)+ampl*spread(spread(cos(ky1*y),1,mx),3,mz)&
                                     *spread(spread(cos(kz1*z),1,mx),2,my)
!
    endsubroutine cosy_cosz
!***********************************************************************
    subroutine cosy_sinz(ampl,f,i,ky,kz)
!
!  initial condition for potential field test
!
!   06-oct-06/axel: coded
!   11-oct-06/wolf: modified to only set one component of aa
!
      integer :: i,j
      real, dimension (mx,my,mz,mfarray) :: f
      real,optional :: ky,kz
      real :: ampl,ky1=1.,kz1=pi
!
!  wavenumber k
!
      if (present(ky)) ky1=ky
      if (present(kz)) kz1=kz
!
      j=i; f(:,:,:,j)=f(:,:,:,j)+ampl*spread(spread(cos(ky1*y),1,mx),3,mz)&
                                     *spread(spread(sin(kz1*z),1,mx),2,my)
!
!  Axel, are you OK with this? If so, I'll eliminate j above.
!
!  Don't do this: we now call this twice from magnetic.f90, and setting
!  just Ax /= 0 makes perfect sense for potential bc tests.
!
!      j=i+1; f(:,:,:,j)=f(:,:,:,j)+ampl*spread(spread(sin(ky1*y),1,mx),3,mz)&
!                                       *spread(spread(cos(kz1*z),1,mx),2,my)
!
    endsubroutine cosy_sinz
!***********************************************************************
    subroutine x3_cosy_cosz(ampl,f,i,ky,kz)
!
!  special initial condition for producing toroidal field (ndynd decay test)
!
!   06-oct-06/axel: coded
!   11-oct-06/wolf: modified to only set one component of aa
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real,optional :: ky,kz
      real :: ampl,ky1=1.,kz1=pi
!
!  wavenumber k
!
      if (present(ky)) ky1=ky
      if (present(kz)) kz1=kz
!
      f(:,:,:,i)=f(:,:,:,i)+ampl*spread(spread(x*(1.-x)*(x-.2),2,my),3,mz)&
                                *spread(spread(cos(ky1*y),1,mx),3,mz)&
                                *spread(spread(cos(kz1*z),1,mx),2,my)
!
    endsubroutine x3_cosy_cosz
!***********************************************************************
    subroutine cosx_coscosy_cosz(ampl,f,i,kx,ky,kz)
!
!  sinusoidal wave, adapted from sinxsinz (that routine was already doing
!  this, but under a different name)
!
!   2-dec-03/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real,optional :: kx,ky,kz
      real :: ampl,kx1=pi/2.,ky1=0.,kz1=pi/2.
!
!  wavenumber k, helicity H=ampl (can be either sign)
!
!  sinx(kx*x)*sin(kz*z)
!
      if (present(kx)) kx1=kx
      if (present(ky)) ky1=ky
      if (present(kz)) kz1=kz
      if (ampl==0) then
        if (lroot) print*,'cosx_cosy_cosz: ampl=0'
      else
        if (lroot) write(*,wave_fmt1) 'cosx_cosy_cosz: ampl,kx,ky,kz=', &
                                      ampl,kx1,ky1,kz1
        f(:,:,:,i)=f(:,:,:,i)+ampl*(spread(spread(cos(kx1*x),2,my),3,mz) &
                                   *spread(spread(-cos(ky1*y)*(          &
                                     ((1./9.)*sin(ky*y)**8)+((8./63.)*   &
                                     sin(ky*y)**6)+((16./105.)*          &
                                     sin(ky*y)**4)+((64./315.)*          &
                                     sin(ky*y)**2)                       &
                                     + (128./315.)                       &
                                    ),1,mx),3,mz)                        &
                                   *spread(spread(cos(kz1*z),1,mx),2,my))
      endif
!
    endsubroutine cosx_coscosy_cosz
!***********************************************************************
    subroutine cos2x_cos2y_cos2z(ampl,f,i,kx,ky,kz)
!
!  sinusoidal wave, adapted from sinxsinz (that routine was already doing
!  this, but under a different name)
!
!  21-feb-08/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real,optional :: kx,ky,kz
      real :: ampl,kx1=1.,ky1=1.,kz1=1.
!
!  wavenumber k, helicity H=ampl (can be either sign)
!
!  sinx(kx*x)*sin(kz*z)
!
      if (present(kx)) kx1=kx
      if (present(ky)) ky1=ky
      if (present(kz)) kz1=kz
      if (ampl==0) then
        if (lroot) print*,'cos2x_cos2y_cos2z: ampl=0'
      else
        if (lroot) write(*,wave_fmt1) 'cos2x_cos2y_cos2z: ampl,kx,ky,kz=', &
                                      ampl,kx1,ky1,kz1
        if (ampl>0.) then
          f(:,:,:,i)=f(:,:,:,i)+ampl*tanh(10.* &
              spread(spread(cos(.5*kx1*x)**2,2,my),3,mz) &
             *spread(spread(cos(.5*ky1*y)**2,1,mx),3,mz) &
             *spread(spread(cos(.5*kz1*z)**2,1,mx),2,my))
        else
          f(:,:,:,i)=f(:,:,:,i)-ampl*(1.-tanh(10.* &
              spread(spread(cos(.5*kx1*x)**2,2,my),3,mz) &
             *spread(spread(cos(.5*ky1*y)**2,1,mx),3,mz) &
             *spread(spread(cos(.5*kz1*z)**2,1,mx),2,my)))
        endif
      endif
!
    endsubroutine cos2x_cos2y_cos2z
!***********************************************************************
    subroutine xpoint(ampl,f,i,x0,y0)
!
! Creates a magnetic X point for a 2D run.
! The vector potential is given by A_z = ampl*(x-x0)*(y -y0)
!
!   9-dec-10/bing: coded
!
      integer, intent(in) :: i
      real, dimension (mx,my,mz,mfarray) ,intent(inout) :: f
      real, intent(in) :: ampl,x0,y0
!
      if (ampl==0) then
        if (lroot) print*,'xpoint: ampl=0'
      else
        f(:,:,n1,i)=ampl*spread(x-x0,2,my)*spread(y-y0,1,mx)
      endif
!
    endsubroutine xpoint
!***********************************************************************
    subroutine xpoint2(ampl,f,i,x0,y0)
!
! Creates a magnetic X point (e.g., for a 2D run);
! The vector potential is given by A_z = (1/2)*ampl*[(x-x0)^2-(y-y0)^2]
!
!  21-may-11/axel: coded
!
      integer, intent(in) :: i
      real, dimension (mx,my,mz,mfarray) ,intent(inout) :: f
      real, intent(in) :: ampl,x0,y0
!
      if (ampl==0) then
        if (lroot) print*,'xpoint: ampl=0'
      else
        f(:,:,n1,i)=.5*ampl*(spread(x-x0,2,my)**2-spread(y-y0,1,mx)**2)
      endif
!
    endsubroutine xpoint2
!***********************************************************************
    subroutine innerbox(ampl,ampl2,f,i,width)
!
!  set background value plus a different value inside a box
!
!   9-mar-08/axel: coded
!
      integer :: i,ll1,ll2,mm1,mm2,nn1,nn2
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,ampl2,width
!
      intent(in)  :: ampl,ampl2,i,width
      intent(out) :: f
!
!  inner box
!
      if (ampl==0.and.ampl2==0) then
        if (lroot) print*,'innerbox: ampl=0'
      else
        if (lroot) write(*,wave_fmt1) 'innerbox: ampl,ampl2,width=', &
                                                 ampl,ampl2,width
        f(:,:,:,i)=ampl
        call find_index_range(x,mx,-width,width,ll1,ll2)
        call find_index_range(y,my,-width,width,mm1,mm2)
        call find_index_range(z,mz,-width,width,nn1,nn2)
        if (lroot) print*,'innerbox: ll1,ll2,mm1,mm2,nn1,nn2=', &
                                     ll1,ll2,mm1,mm2,nn1,nn2
        f(ll1:ll2,mm1:mm2,nn1:nn2,i)=ampl2
      endif
!
    endsubroutine innerbox
!***********************************************************************
    subroutine couette(ampl,mu,f,i)
!
!   couette flow,  Omega = a + b/r^2
!
!   12-jul-07/mgellert: coded
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl, mu, omegao, omegai, rinner, router, a, b
!
      rinner=xyz0(1)
      router=xyz1(1)
      omegai=ampl
      omegao=ampl*mu
!
      if (ampl==0) then
        if (lroot) print*,'couette flow: omegai=0'
      else
        if (lroot) write(*,wave_fmt1) 'couette flow: omegai,omegao=', omegai,omegao
        a = ( (omegao/omegai) - (rinner/router)**2 ) / ( 1. - (rinner/router)**2 ) * omegai
        b = ( ( 1. - (omegao/omegai) ) / ( 1. - (rinner/router)**2 ) ) * rinner**2 * omegai
        f(:,:,:,i)=f(:,:,:,i)+spread(spread((a*x + b/x),2,my),3,mz)
      endif
!
    endsubroutine couette
!***********************************************************************
    subroutine couette_rings(ampl,mu,nr,om_nr,gap,f,i)
!
!   * (finite) couette flow  with top/bottom split in several rings
!   * can be used with the 'freeze' condition for top/bottom BCs
!   * gap is the width of the gap between two rings filled with a tanh smoothing profile
!     (in experiments you usually have something around gap=0.05*D)
!
!   18-jul-07/mgellert: coded
!   21-aug-07/mgellert: made compatible with nprocx>1
!   09-oct-07/mgellert: smoothing profile between neighboring rings added
!
      integer                           :: i, k, l, nr
      real, dimension (nr)               :: om_nr, xsteps
      real, dimension (nr+2)             :: om_all, om_diff
      real, dimension (:), allocatable   :: omx
      real, dimension (mx,my,mz,mfarray) :: f
      real                              :: ampl, mu, gap, omegao, omegai, rinner, router, step
      real                              :: x0, y0
      character(len=20)                 :: unfmt
!
      rinner=xyz0(1)
      router=xyz1(1)
      omegai=ampl
      omegao=ampl*mu
!
      allocate(omx(mx))
      omx=0.
!
      step=(router-rinner)/nr
      do k=1,nr
        xsteps(k)=rinner+k*step-gap/2. ! ring boundaries
      enddo
!
      if (gap>tini) then
        om_all(1)=omegai
        om_all(2:nr+1)=om_nr
        om_all(nr+2)=omegao
        om_diff=0.
        om_diff=om_all(1:nr+1)-om_all(2:nr+2) ! difference in omega from ring to ring
        omx(l1:l2)=omegai
        do k=1,nr+1
          if (k==1) then
            x0=rinner+gap/2.
          else
            x0=rinner+(k-1)*step-gap/2.
          endif
          y0=0.5*om_diff(k)
          omx(l1:l2) = omx(l1:l2) - ( y0*(1.+tanh((x(l1:l2)-x0)/(0.15*gap+tini))) );
        enddo
      else
        do l=l1,l2
          k=nr
          do while (k>0)
            if (x(l)<=xsteps(k)) omx(l)=om_nr(k)*omegai
            k=k-1
          enddo
        enddo
      endif
!
      if (ampl==0) then
        if (lroot) print*,'couette flow with rings: omegai=0'
      else
        write(unfmt,FMT='(A12,1I2.2,A6)') '(A19,I2,A24,',nr,'2F7.3)'
        if (lroot) write(*,FMT=unfmt) 'couette flow with ',nr,' rings: omegai,omegao = ', omegai,omegao
        write(unfmt,FMT='(A12,1I2.2,A5)') '(A19,I2,A24,',nr,'F7.3)'
        if (lroot) write(*,FMT=unfmt) 'couette flow with ',nr,' rings: omega_rings   = ',om_nr*omegai
        if (lroot) write(*,FMT=unfmt) 'couette flow with ',nr,' rings: radius_rings  = ',xsteps
        f(:,:,:,i)=f(:,:,:,i)+spread(spread((omx*x),2,my),3,mz) ! velocity up=omx*x
      endif
!
      if (allocated(omx)) deallocate(omx)
!
    endsubroutine couette_rings
!***********************************************************************
    subroutine hat(ampl,f,i,width,kx,ky,kz)
!
!  hat bump
!
!   2-may-03/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real,optional :: kx,ky,kz
      real :: ampl,width,k=1.,width2,k2
!
!  prepare
!
      width2=width**2
!
!  set x-hat
!
      if (present(kx)) then
        k=kx
        k2=k**2
        if (ampl==0) then
          if (lroot) print*,'hat: ampl=0; kx=',k
        else
          if (lroot) print*,'hat: kx,i,ampl=',k,i,ampl
          f(:,:,:,i)=f(:,:,:,i)+ampl*spread(spread(.5+.5*tanh(k2*(width2-x**2)),2,my),3,mz)
        endif
      endif
!
!  set y-hat
!
      if (present(ky)) then
        k=ky
        k2=k**2
        if (ampl==0) then
          if (lroot) print*,'hat: ampl=0; ky=',k
        else
          if (lroot) print*,'hat: ky,i,ampl=',k,i,ampl
          f(:,:,:,i)=f(:,:,:,i)+ampl*spread(spread(.5+.5*tanh(k2*(width2-y**2)),1,mx),3,mz)
        endif
      endif
!
!  set z-hat
!
      if (present(kz)) then
        k=kz
        k2=k**2
        if (ampl==0) then
          if (lroot) print*,'hat: ampl=0; kz=',k
        else
          if (lroot) print*,'hat: kz,i,ampl=',k,i,ampl
          f(:,:,:,i)=f(:,:,:,i)+ampl*spread(spread(.5+.5*tanh(k2*(width2-z**2)),1,mx),2,my)
        endif
      endif
!
    endsubroutine hat
!***********************************************************************
    subroutine hat3d(ampl,f,i,width,kx,ky,kz)
!
!  Three-dimensional hat bump
!
!   9-nov-04/anders: coded
!
      integer :: i,l
      real, dimension (mx,my,mz,mfarray) :: f
      real :: kx,ky,kz,kx2,ky2,kz2
      real :: ampl,width,width2
!
!  set hat
!
      if (lroot) print*,'hat3d: kx,ky,kz=',kx,ky,kz
      if (ampl==0) then
        if (lroot) print*,'hat3d: ampl=0'
      else
        if (lroot) print*,'hat3d: ampl=',ampl
        width2=width**2
        kx2=kx**2
        ky2=ky**2
        kz2=kz**2
        do l=l1,l2; do m=m1,m2; do n=n1,n2
            f(l,m,n,i) = f(l,m,n,i) + ampl*( &
                (0.5+0.5*tanh(kx2*(width2-x(l)**2))) * &
                (0.5+0.5*tanh(ky2*(width2-y(m)**2))) * &
                (0.5+0.5*tanh(kz2*(width2-z(n)**2))) )
        enddo; enddo; enddo
      endif
!
    endsubroutine hat3d
!***********************************************************************
    subroutine gaussian(ampl,f,i,kx,ky,kz)
!
!  gaussian bump
!
!   2-may-03/axel: coded
!  20-sep-03/axel: added 1/2 factor in defn; hopefully ok with everyone?
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real,optional :: kx,ky,kz
      real :: ampl,k=1.
!
!  wavenumber k
!
!  set x-wave
!
      if (present(kx)) then
        k=kx
        if (ampl==0) then
          if (lroot) print*,'gaussian: ampl=0; kx=',k
        else
          if (lroot) print*,'gaussian: kx,i=',k,i
          f(:,:,:,i)=f(:,:,:,i)+ampl*spread(spread(exp(-.5*(k*x)**2),2,my),3,mz)
        endif
      endif
!
!  set y-wave
!
      if (present(ky)) then
        k=ky
        if (ampl==0) then
          if (lroot) print*,'gaussian: ampl=0; ky=',k
        else
          if (lroot) print*,'gaussian: ky,i=',k,i
          f(:,:,:,i)=f(:,:,:,i)+ampl*spread(spread(exp(-.5*(k*y)**2),1,mx),3,mz)
        endif
      endif
!
!  set z-wave
!
      if (present(kz)) then
        k=kz
        if (ampl==0) then
          if (lroot) print*,'gaussian: ampl=0; kz=',k
        else
          if (lroot) print*,'gaussian: kz,i=',k,i
          f(:,:,:,i)=f(:,:,:,i)+ampl*spread(spread(exp(-.5*(k*z)**2),1,mx),2,my)
        endif
      endif
!
    endsubroutine gaussian
!***********************************************************************
    subroutine gaussian3d(ampl,f,i,radius)
!
!  gaussian 3-D bump
!
!  28-may-03/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,radius,radius21
!
      radius21=1./radius**2
      do n=n1,n2; do m=m1,m2
        f(l1:l2,m,n,i)=f(l1:l2,m,n,i)+ampl*exp(-(x(l1:l2)**2+y(m)**2+z(n)**2)*radius21)
      enddo; enddo
!
    endsubroutine gaussian3d
!***********************************************************************
    subroutine gaussianpos(ampl,f,i,radius,posx,posy,posz)
!
!  gaussian 3-D bump centered in specific position
!
!  21-sep-09/rplasson: coded from gaussian3d
!  Maybe could have been done by extending gaussian3d, but didn't want to interfere
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,radius,posx, posy, posz, radius21
!
      radius21=1./radius**2
      do n=n1,n2; do m=m1,m2
        f(l1:l2,m,n,i)=f(l1:l2,m,n,i)+ampl*exp(-((x(l1:l2)-posx)**2+(y(m)-posy)**2+(z(n)-posz)**2)*radius21)
!        do l=l1,l2
!          print*,"c(",x(l),",",y(m),",",z(n),")=", f(l,m,n,i)
!        enddo
      enddo; enddo
!
    endsubroutine gaussianpos
!***********************************************************************
    subroutine parabola(ampl,f,i,kx,ky,kz)
!
!  gaussian bump
!
!   2-may-03/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real,optional :: kx,ky,kz
      real :: ampl,k=1.
!
!  wavenumber k
!
!  set x-wave
!
      if (present(kx)) then
        k=kx
        if (ampl==0) then
          if (lroot) print*,'parabola: ampl=0; kx=',k
        else
          if (lroot) print*,'parabola: kx,i=',k,i
          f(:,:,:,i)=f(:,:,:,i)+ampl*spread(spread((-(k*x)**2),2,my),3,mz)
        endif
      endif
!
!  set y-wave
!
      if (present(ky)) then
        k=ky
        if (ampl==0) then
          if (lroot) print*,'parabola: ampl=0; ky=',k
        else
          if (lroot) print*,'parabola: ky,i=',k,i
          f(:,:,:,i)=f(:,:,:,i)+ampl*spread(spread((-(k*y)**2),1,mx),3,mz)
        endif
      endif
!
!  set z-wave
!
      if (present(kz)) then
        k=kz
        if (ampl==0) then
          if (lroot) print*,'parabola: ampl=0; kz=',k
        else
          if (lroot) print*,'parabola: kz,i=',k,i
          f(:,:,:,i)=f(:,:,:,i)+ampl*spread(spread((-(k*z)**2),1,mx),2,my)
        endif
      endif
!
    endsubroutine parabola
!***********************************************************************
    subroutine wave(ampl,f,i,kx,ky,kz)
!
!  sinusoidal wave
!
!   6-jul-02/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real,optional :: kx,ky,kz
      real :: ampl,k=1.
!
!  wavenumber k
!
!  set x-wave
!
      if (present(kx)) then
        k=kx
        if (ampl==0) then
          if (lroot) print*,'wave: ampl=0; kx=',k
        else
          if (lroot) print*,'wave: kx,i=',k,i
          f(:,:,:,i)=f(:,:,:,i)+ampl*spread(spread(sin(k*x),2,my),3,mz)
        endif
      endif
!
!  set y-wave
!
      if (present(ky)) then
        k=ky
        if (ampl==0) then
          if (lroot) print*,'wave: ampl=0; ky=',k
        else
          if (lroot) print*,'wave: ky,i=',k,i
          f(:,:,:,i)=f(:,:,:,i)+ampl*spread(spread(sin(k*y),1,mx),3,mz)
        endif
      endif
!
!  set z-wave
!
      if (present(kz)) then
        k=kz
        if (ampl==0) then
          if (lroot) print*,'wave: ampl=0; kz=',k
        else
          if (lroot) print*,'wave: kz,i=',k,i
          f(:,:,:,i)=f(:,:,:,i)+ampl*spread(spread(sin(k*z),1,mx),2,my)
        endif
      endif
!
    endsubroutine wave
!***********************************************************************
    subroutine linprof(ampl,f,i,kx,ky,kz)
!
!  periodic linear profile
!
!  21-sep-09/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real,optional :: kx,ky,kz
      real :: ampl,k=1.
!
!  wavenumber k
!
!  set x-linprof
!
      if (present(kx)) then
        k=kx
        if (ampl==0) then
          if (lroot) print*,'linprof: ampl=0; kx=',k
        else
          if (lroot) print*,'linprof: kx,i=',k,i
          f(:,:,:,i)=f(:,:,:,i)+ampl*spread(spread(asin(sin(k*x)),2,my),3,mz)
        endif
      endif
!
!  set y-linprof
!
      if (present(ky)) then
        k=ky
        if (ampl==0) then
          if (lroot) print*,'linprof: ampl=0; ky=',k
        else
          if (lroot) print*,'linprof: ky,i=',k,i
          f(:,:,:,i)=f(:,:,:,i)+ampl*spread(spread(asin(sin(k*y)),1,mx),3,mz)
        endif
      endif
!
!  set z-linprof
!
      if (present(kz)) then
        k=kz
        if (ampl==0) then
          if (lroot) print*,'linprof: ampl=0; kz=',k
        else
          if (lroot) print*,'linprof: kz,i=',k,i
          f(:,:,:,i)=f(:,:,:,i)+ampl*spread(spread(asin(sin(k*z)),1,mx),2,my)
        endif
      endif
!
    endsubroutine linprof
!***********************************************************************
    subroutine wave_uu(ampl,f,i,kx,ky,kz)
!
!  sinusoidal wave
!
!  14-apr-03/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real,optional :: kx,ky,kz
      real :: ampl,k=1.
!
!  wavenumber k
!
!  set x-wave
!
      if (present(kx)) then
        k=kx
        if (ampl==0) then
          if (lroot) print*,'wave_uu: ampl=0; kx=',k
        else
          if (lroot) print*,'wave_uu: kx,i=',k,i
          f(:,:,:,i)=log(1.+ampl*spread(spread(sin(k*x),2,my),3,mz)*f(:,:,:,iux))
        endif
      endif
!
!  set y-wave
!
      if (present(ky)) then
        k=ky
        if (ampl==0) then
          if (lroot) print*,'wave_uu: ampl=0; ky=',k
        else
          if (lroot) print*,'wave_uu: ky,i=',k,i
          f(:,:,:,i)=log(1.+ampl*spread(spread(sin(k*y),1,mx),3,mz)*f(:,:,:,iuy))
        endif
      endif
!
!  set z-wave
!
      if (present(kz)) then
        k=kz
        if (ampl==0) then
          if (lroot) print*,'wave_uu: ampl=0; kz=',k
        else
          if (lroot) print*,'wave_uu: kz,i=',k,i,iuz
          f(:,:,:,i)=log(1.+ampl*spread(spread(sin(k*z),1,mx),2,my)*f(:,:,:,iuz))
        endif
      endif
!
    endsubroutine wave_uu
!***********************************************************************
    subroutine acosy(ampl,f,i,ky)
!
!  mode
!
!  30-oct-03/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,ky
!
      do n=n1,n2; do m=m1,m2
        f(l1:l2,m,n,i)=ampl*cos(ky*y(m))
        f(l1:l2,m,n,i+1)=0.
        f(l1:l2,m,n,i+2)=0
      enddo; enddo
!
    endsubroutine acosy
!***********************************************************************
    subroutine modes(ampl,coef,f,i,kx,ky,kz)
!
!  mode
!
!  30-oct-03/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      complex :: coef
      complex :: ii=(0.,1.)
      real :: ampl,kx,ky,kz
!
      do n=n1,n2; do m=m1,m2
        f(l1:l2,m,n,i)=ampl*real(coef*exp(ii*(kx*x(l1:l2)+ky*y(m)+kz*z(n))))
      enddo; enddo
!
    endsubroutine modes
!***********************************************************************
    subroutine modev(ampl,coef,f,i,kx,ky,kz)
!
!  mode
!
!  30-oct-03/axel: coded
!
      integer :: i,ivv
      real, dimension (mx,my,mz,mfarray) :: f
      complex, dimension (3) :: coef
      complex :: ii=(0.,1.)
      real :: ampl,kx,ky,kz
!
      do n=n1,n2; do m=m1,m2
        do ivv=0,2
          f(l1:l2,m,n,ivv+i)=ampl*real(coef(ivv+1)*exp(ii*(kx*x(l1:l2)+ky*y(m)+kz*z(n))))
        enddo
      enddo;enddo
!
    endsubroutine modev
!***********************************************************************
    subroutine modeb(ampl,coefb,f,i,kx,ky,kz)
!
!  mode
!
!  30-oct-03/axel: coded
!
      integer :: i,ivv
      real, dimension (mx,my,mz,mfarray) :: f
      complex, dimension (3) :: coef,coefb
      complex :: ii=(0.,1.)
      real :: ampl,kx,ky,kz,k2
!
      print*,'print Ak coefficients from Bk coefficients'
      k2=kx**2+ky**2+kz**2
      coef(1)=(ky*coefb(3)-kz*coefb(2))/k2
      coef(2)=(kz*coefb(1)-kx*coefb(3))/k2
      coef(3)=(kx*coefb(2)-ky*coefb(1))/k2
      print*,'coef=',coef
!
      do n=n1,n2; do m=m1,m2
        do ivv=0,2
          f(l1:l2,m,n,ivv+i)=ampl*real(coef(ivv+1)*exp(ii*(kx*x(l1:l2)+ky*y(m)+kz*z(n))))
        enddo
      enddo; enddo
!
    endsubroutine modeb
!***********************************************************************
    subroutine jump(f,i,fleft,fright,width,dir)
!
!  jump
!
!  19-sep-02/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx) :: profx
      real, dimension (my) :: profy
      real, dimension (mz) :: profz
      real :: fleft,fright,width
      character(len=*) :: dir
!
!  jump; check direction
!
      select case (dir)
!
      case ('x')
        profx=fleft+(fright-fleft)*.5*(1.+tanh(x/width))
        f(:,:,:,i)=f(:,:,:,i)+spread(spread(profx,2,my),3,mz)
!
      case ('y')
        profy=fleft+(fright-fleft)*.5*(1.+tanh(y/width))
        f(:,:,:,i)=f(:,:,:,i)+spread(spread(profy,1,mx),3,mz)
!
      case ('z')
        profz=fleft+(fright-fleft)*.5*(1.+tanh(z/width))
        f(:,:,:,i)=f(:,:,:,i)+spread(spread(profz,1,mx),2,my)
!
      case default
        print*,'jump: no default value'
!
      endselect
!
    endsubroutine jump
!***********************************************************************
    subroutine bjump(f,i,fyleft,fyright,fzleft,fzright,width,dir)
!
!  jump in B-field (in terms of magnetic vector potential)
!
!   9-oct-02/wolf+axel: coded
!  21-apr-05/axel: added possibility of Bz component
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx) :: profy,profz,alog_cosh_xwidth
      real :: fyleft,fyright,fzleft,fzright,width
      character(len=*) :: dir
!
!  jump; check direction
!
      select case (dir)
!
!  use correct signs when calling this routine
!  Ay=+int Bz dx
!  Az=-int By dx
!
!  alog(cosh(x/width)) =
!
      case ('x')
        alog_cosh_xwidth=abs(x/width)+alog(.5*(1.+exp(-2*abs(x/width))))
        profz=.5*(fyright+fyleft)*x &
             +.5*(fyright-fyleft)*width*alog_cosh_xwidth
        profy=.5*(fzright+fzleft)*x &
             +.5*(fzright-fzleft)*width*alog_cosh_xwidth
        f(:,:,:,i+1)=f(:,:,:,i+1)+spread(spread(profy,2,my),3,mz)
        f(:,:,:,i+2)=f(:,:,:,i+2)-spread(spread(profz,2,my),3,mz)
!
      case default
        print*,'bjump: no default value'
!
      endselect
!
    endsubroutine bjump
!***********************************************************************
    subroutine bjumpz(f,i,fxleft,fxright,fyleft,fyright,width,dir)
!
!  jump in B-field (in terms of magnetic vector potential)
!
!   9-oct-02/wolf+axel: coded
!  21-apr-05/axel: added possibility of Bz component
!
      integer :: i,il,im
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mz) :: profx,profy,alog_cosh_zwidth
      real :: fyleft,fyright,fxleft,fxright,width
      character(len=*) :: dir
!
      select case (dir)
      case ('z')
        alog_cosh_zwidth=abs(z/width)+alog(.5*(1.+exp(-2*abs(z/width))))
        profx=.5*(fyright+fyleft)*z &
             +.5*(fyright-fyleft)*width*alog_cosh_zwidth
        profy=.5*(fxright+fxleft)*z &
             +.5*(fxright-fxleft)*width*alog_cosh_zwidth
        do il=1,mx
        do im=1,my
          f(il,im,:,i ) =f(il,im,:,i  )+profx
          f(il,im,:,i+1)=f(il,im,:,i+1)-profy
        enddo
        enddo
      case default
        print*,'bjump: no default value'
!
      endselect
!
    endsubroutine bjumpz
!***********************************************************************
    subroutine beltrami(ampl,f,i,kx,ky,kz,kx2,ky2,kz2,phase)
!
!  Beltrami field (as initial condition)
!
!  19-jun-02/axel: coded
!   5-jul-02/axel: made additive (if called twice), kx,ky,kz are optional
!
      integer :: i,j
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx) :: sfuncx,cfuncx
      real, dimension (my) :: sfuncy,cfuncy
      real, dimension (mz) :: sfuncz,cfuncz
      real,optional :: kx,ky,kz,kx2,ky2,kz2,phase
      real :: ampl,k=1.,ph
!
!  possibility of shifting the Beltrami wave by phase ph
!
      if (present(phase)) then
        if (lroot) print*,'Beltrami: phase=',phase
        ph=phase
      else
        ph=0.
      endif
!
!  wavenumber k, helicity H=ampl (can be either sign)
!
!  set x-dependent Beltrami field
!
      if (present(kx)) then
        k=abs(kx)
        if (k==0) print*,'beltrami: k must not be zero!'
        cfuncx=sign(sqrt(abs(ampl/k)),kx)*cos(k*x+ph)
        sfuncx=sign(sqrt(abs(ampl/k)),kx)*sin(k*x+ph)
        if (present(kx2)) sfuncx=sfuncx*sin(kx2*x+ph)
        if (ampl==0) then
          if (lroot) print*,'beltrami: ampl=0; kx=',k
        elseif (ampl>0) then
          if (lroot) print*,'beltrami: Beltrami field (pos-hel): kx,i=',k,i
          j=i+1; f(:,:,:,j)=f(:,:,:,j)+spread(spread(sfuncx,2,my),3,mz)
          j=i+2; f(:,:,:,j)=f(:,:,:,j)+spread(spread(cfuncx,2,my),3,mz)
        elseif (ampl<0) then
          if (lroot) print*,'beltrami: Beltrami field (neg-hel): kx,i=',k,i
          j=i+1; f(:,:,:,j)=f(:,:,:,j)+spread(spread(cfuncx,2,my),3,mz)
          j=i+2; f(:,:,:,j)=f(:,:,:,j)+spread(spread(sfuncx,2,my),3,mz)
        endif
      endif
!
!  set y-dependent Beltrami field
!
      if (present(ky)) then
        k=abs(ky)
        if (k==0) print*,'beltrami: k must not be zero!'
        cfuncy=sign(sqrt(abs(ampl/k)),ky)*cos(k*y+ph)
        sfuncy=sign(sqrt(abs(ampl/k)),ky)*sin(k*y+ph)
        if (present(ky2)) sfuncy=sfuncy*sin(ky2*y+ph)
        if (ampl==0) then
          if (lroot) print*,'beltrami: ampl=0; ky=',k
        elseif (ampl>0) then
          if (lroot) print*,'beltrami: Beltrami field (pos-hel): ky,i=',k,i
          j=i;   f(:,:,:,j)=f(:,:,:,j)+spread(spread(cfuncy,1,mx),3,mz)
          j=i+2; f(:,:,:,j)=f(:,:,:,j)+spread(spread(sfuncy,1,mx),3,mz)
        elseif (ampl<0) then
          if (lroot) print*,'beltrami: Beltrami field (neg-hel): ky,i=',k,i
          j=i;   f(:,:,:,j)=f(:,:,:,j)+spread(spread(sfuncy,1,mx),3,mz)
          j=i+2; f(:,:,:,j)=f(:,:,:,j)+spread(spread(cfuncy,1,mx),3,mz)
        endif
      endif
!
!  set z-dependent Beltrami field
!
      if (present(kz)) then
        k=abs(kz)
        if (k==0) print*,'beltrami: k must not be zero!'
        cfuncz=sign(sqrt(abs(ampl/k)),kz)*cos(k*z+ph)
        sfuncz=sign(sqrt(abs(ampl/k)),kz)*sin(k*z+ph)
        if (present(kz2)) sfuncz=sfuncz*sin(kz2*z+ph)
        if (ampl==0) then
          if (lroot) print*,'beltrami: ampl=0; kz=',k
        elseif (ampl>0) then
          if (lroot) print*,'beltrami: Beltrami field (pos-hel): kz,i=',k,i
          j=i;   f(:,:,:,j)=f(:,:,:,j)+spread(spread(sfuncz,1,mx),2,my)
          j=i+1; f(:,:,:,j)=f(:,:,:,j)+spread(spread(cfuncz,1,mx),2,my)
        elseif (ampl<0) then
          if (lroot) print*,'beltrami: Beltrami field (neg-hel): kz,i=',k,i
          j=i;   f(:,:,:,j)=f(:,:,:,j)+spread(spread(cfuncz,1,mx),2,my)
          j=i+1; f(:,:,:,j)=f(:,:,:,j)+spread(spread(sfuncz,1,mx),2,my)
        endif
      endif
!
    endsubroutine beltrami
!***********************************************************************
    subroutine bhyperz(ampl,f,i,kz,nfactor)
!
!  Beltrami field with wavelengths that are cube roots of unity
!
!
      integer :: i
      integer :: ix,iy,iz
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,kz,nfactor
!
!
!  set z-dependent Beltrami field
!
      do ix=1,mx; do iy=1,my;do iz=1,mz
        f(ix,iy,iz,i) = ampl*((1-nfactor)*sin(kz*z(iz)) &
            +nfactor*(sin(kz*z(iz)/2.)*cosh(kz*z(iz)*sqrt(3.)/2.) &
                + sqrt(3.)*cos(kz*z(iz)/2.)*sinh(kz*z(iz)*sqrt(3.)/2.) ))
        f(ix,iy,iz,i+1) = ampl*((1-nfactor)*cos(kz*z(iz)) &
            - nfactor*(cos(kz*z(iz)/2.)*cosh(kz*z(iz)*sqrt(3.)/2.) &
              +sqrt(3.)*sin(kz*z(iz)/2.)*sinh(kz*z(iz)*sqrt(3.)/2.) ))
      enddo;enddo;enddo
!
    endsubroutine bhyperz
!***********************************************************************
    subroutine beltrami_complex(ampl,f,i,kx,ky,kz,kx2,ky2,kz2,phase)
!
!  Beltrami field (as initial condition)
!
!  23-sep-10/dhruba: aped from beltrami
!
      integer :: i,j
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx) :: sfuncx,cfuncx
      real, dimension (my) :: sfuncy,cfuncy
      real, dimension (mz) :: sfuncz,cfuncz
      real,optional :: kx,ky,kz,kx2,ky2,kz2,phase
      real :: ampl,k=1.,ph
      complex :: omg,omgsqr
!
! complex cube roots of unity
!
      omg=cmplx(-sqrt(2.),sqrt(3.)/2.)
      omgsqr=conjg(omg)
!
!  possibility of shifting the Beltrami wave by phase ph
!
      if (present(phase)) then
        if (lroot) print*,'Beltrami: phase=',phase
        ph=phase
      else
        ph=0.
      endif
!
!  wavenumber k, helicity H=ampl (can be either sign)
!
!  set x-dependent Beltrami field
!
      if (present(kx)) then
        k=abs(kx)
        if (k==0) print*,'beltrami: k must not be zero!'
        cfuncx=sign(sqrt(abs(ampl/k)),kx)*&
            real(cos(k*x+ph)+omgsqr*cos(omg*k*x+ph)+omg*cos(omgsqr*k*x+ph))
        sfuncx=sign(sqrt(abs(ampl/k)),kx)*&
            real(sin(k*x+ph)+omgsqr*sin(omg*k*x+ph)+omg*sin(omgsqr*k*x+ph))
        if (present(kx2)) sfuncx=sfuncx*sin(kx2*x+ph)
        if (ampl==0) then
          if (lroot) print*,'beltrami: ampl=0; kx=',k
        elseif (ampl>0) then
          if (lroot) print*,'beltrami: Beltrami field (pos-hel): kx,i=',k,i
          j=i+1; f(:,:,:,j)=f(:,:,:,j)+spread(spread(sfuncx,2,my),3,mz)
          j=i+2; f(:,:,:,j)=f(:,:,:,j)+spread(spread(cfuncx,2,my),3,mz)
        elseif (ampl<0) then
          if (lroot) print*,'beltrami: Beltrami field (neg-hel): kx,i=',k,i
          j=i+1; f(:,:,:,j)=f(:,:,:,j)+spread(spread(cfuncx,2,my),3,mz)
          j=i+2; f(:,:,:,j)=f(:,:,:,j)+spread(spread(sfuncx,2,my),3,mz)
        endif
      endif
!
!  set y-dependent Beltrami field
!
      if (present(ky)) then
        k=abs(ky)
        if (k==0) print*,'beltrami: k must not be zero!'
        cfuncy=sign(sqrt(abs(ampl/k)),ky)*&
            real(cos(k*y+ph)+omgsqr*cos(omg*k*y+ph)+omg*cos(omgsqr*k*y+ph))
        sfuncy=sign(sqrt(abs(ampl/k)),ky)*&
            real(sin(k*y+ph)+omgsqr*sin(omg*k*y+ph)+omg*sin(omgsqr*k*y+ph))
        if (present(ky2)) sfuncy=sfuncy*sin(ky2*y+ph)
        if (ampl==0) then
          if (lroot) print*,'beltrami: ampl=0; ky=',k
        elseif (ampl>0) then
          if (lroot) print*,'beltrami: Beltrami field (pos-hel): ky,i=',k,i
          j=i;   f(:,:,:,j)=f(:,:,:,j)+spread(spread(cfuncy,1,mx),3,mz)
          j=i+2; f(:,:,:,j)=f(:,:,:,j)+spread(spread(sfuncy,1,mx),3,mz)
        elseif (ampl<0) then
          if (lroot) print*,'beltrami: Beltrami field (neg-hel): ky,i=',k,i
          j=i;   f(:,:,:,j)=f(:,:,:,j)+spread(spread(sfuncy,1,mx),3,mz)
          j=i+2; f(:,:,:,j)=f(:,:,:,j)+spread(spread(cfuncy,1,mx),3,mz)
        endif
      endif
!
!  set z-dependent Beltrami field
!
      if (present(kz)) then
        k=abs(kz)
        if (k==0) print*,'beltrami: k must not be zero!'
        cfuncz=sign(sqrt(abs(ampl/k)),kz)*&
            real(cos(k*z+ph)-omgsqr*cos(omg*k*z+ph)-omg*cos(omgsqr*k*z+ph))
        sfuncz=sign(sqrt(abs(ampl/k)),kz)*&
            real(sin(k*z+ph)-omgsqr*sin(omg*k*z+ph)-omg*sin(omgsqr*k*z+ph))
        if (present(kz2)) sfuncz=sfuncz*sin(kz2*z+ph)
        if (ampl==0) then
          if (lroot) print*,'beltrami: ampl=0; kz=',k
        elseif (ampl>0) then
          if (lroot) print*,'beltrami: Beltrami field (pos-hel): kz,i=',k,i
          j=i;   f(:,:,:,j)=f(:,:,:,j)+spread(spread(sfuncz,1,mx),2,my)
          j=i+1; f(:,:,:,j)=f(:,:,:,j)+spread(spread(cfuncz,1,mx),2,my)
        elseif (ampl<0) then
          if (lroot) print*,'beltrami: Beltrami field (neg-hel): kz,i=',k,i
          j=i;   f(:,:,:,j)=f(:,:,:,j)+spread(spread(cfuncz,1,mx),2,my)
          j=i+1; f(:,:,:,j)=f(:,:,:,j)+spread(spread(sfuncz,1,mx),2,my)
        endif
      endif
!
    endsubroutine beltrami_complex
!***********************************************************************
    subroutine bessel_x(ampl,f,i,kx)
!
!  Bessel function field (as initial condition)
!
!  12-nov-09/axel: coded
!
      integer :: i,j,l
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx) :: J0,J1
      real :: ampl,kx,kx1
!
!  set x-dependent Bessel function field.
!  For Az, subtract value on boundary so that it matches
!  the perfect conductor boundary condition.
!
      if (ampl==0) then
        if (lroot) print*,'bessel_x: ampl=0; kx=',kx
      else
        if (lroot) print*,'bessel_x: Bessel function field: kx,i=',kx,i
        kx1=1./kx
        do l=1,mx
          J0(l)=kx1*(bessj(0,kx*x(l))-bessj(0,kx*xyz1(1)))
          J1(l)=kx1*bessj(1,kx*x(l))
        enddo
        j=i+1; f(:,:,:,j)=f(:,:,:,j)+ampl*spread(spread(J1,2,my),3,mz)
        j=i+2; f(:,:,:,j)=f(:,:,:,j)+ampl*spread(spread(J0,2,my),3,mz)
      endif
!
    endsubroutine bessel_x
!***********************************************************************
    subroutine bessel_az_x(ampl,f,i,kx)
!
!  Bessel function field (as initial condition)
!
!  12-nov-09/axel: coded
!
      integer :: i,j,l
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx) :: J0
      real :: ampl,kx,kx1
!
!  set x-dependent Bessel function field
!
      if (ampl==0) then
        if (lroot) print*,'bessel_az_x: ampl=0; kx=',kx
      else
        if (lroot) print*,'bessel_az_x: Bessel function field: kx,i=',kx,i
        kx1=1./kx
        do l=1,mx
          J0(l)=kx1*bessj(0,kx*x(l))
        enddo
        j=i+2; f(:,:,:,j)=f(:,:,:,j)+ampl*spread(spread(J0,2,my),3,mz)
      endif
!
    endsubroutine bessel_az_x
!***********************************************************************
    subroutine rolls(ampl,f,i,kx,kz)
!
!  convection rolls (as initial condition)
!
!  21-aug-07/axel: coded
!
      integer :: i,j
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,kx,kz,zbot
!
!  check input parameters
!
      zbot=xyz0(3)
      if (lroot) print*,'rolls: i,kx,kz,zbot=',i,kx,kz,zbot
!
!  set stream function psi=sin(kx*x)*sin(kz*(z-zbot))
!
      j=i
      f(:,:,:,j)=f(:,:,:,j)-ampl*kz*spread(spread(sin(kx*x       ),2,my),3,mz)&
                                   *spread(spread(cos(kz*(z-zbot)),1,mx),2,my)
      j=i+2
      f(:,:,:,j)=f(:,:,:,j)+ampl*kx*spread(spread(cos(kx*x       ),2,my),3,mz)&
                                   *spread(spread(sin(kz*(z-zbot)),1,mx),2,my)
!
    endsubroutine rolls
!***********************************************************************
    subroutine robertsflow(ampl,f,i,relhel)
!
!  Roberts Flow (as initial condition)
!
!   9-jun-05/axel: coded
!
      integer :: i,j
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,k=1.,kf,fac1,fac2,relhel
!
!  prepare coefficients
!
      kf=k*sqrt(2.)
      fac1=sqrt(2.)*ampl*k/kf
      fac2=sqrt(2.)*ampl*relhel
!
      j=i+0; f(:,:,:,j)=f(:,:,:,j)-fac1*spread(spread(cos(k*x),2,my),3,mz)&
                                       *spread(spread(sin(k*y),1,mx),3,mz)
!
      j=i+1; f(:,:,:,j)=f(:,:,:,j)+fac1*spread(spread(sin(k*x),2,my),3,mz)&
                                       *spread(spread(cos(k*y),1,mx),3,mz)
!
      j=i+2; f(:,:,:,j)=f(:,:,:,j)+fac2*spread(spread(cos(k*x),2,my),3,mz)&
                                       *spread(spread(cos(k*y),1,mx),3,mz)
!
    endsubroutine robertsflow
!***********************************************************************
    subroutine vecpatternxy(ampl,f,i,kx,ky,kz)
!
!  horizontal pattern with exponential decay (as initial condition)
!
!   9-jun-05/axel: coded
!
      integer :: i,j
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,kx,ky,kz
!
!  prepare coefficients
!
      j=i+0; f(:,:,:,j)=f(:,:,:,j)-ampl*spread(spread(sin(ky*y),1,mx),3,mz) &
        *spread(spread(exp(-abs(kz*z)),1,mx),2,my)
      j=i+1; f(:,:,:,j)=f(:,:,:,j)+ampl*spread(spread(sin(kx*x),2,my),3,mz) &
        *spread(spread(exp(-abs(kz*z)),1,mx),2,my)
!
    endsubroutine vecpatternxy
!***********************************************************************
    subroutine bipolar(ampl,f,i,kx,ky,kz)
!
!  horizontal pattern with exponential decay (as initial condition)
!
!  24-may-09/axel: coded
!
      integer :: i,j
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,kx,ky,kz
!
!  sets up a nearly force-free bipolar region
!
      j=i+1; f(:,:,:,j)=f(:,:,:,j)+ampl &
        *spread(spread(exp(-(kx*x)**2),2,my),3,mz) &
        *spread(spread(exp(-(ky*y)**2),1,mx),3,mz) &
        *spread(spread(exp(-abs(kz*z)),1,mx),2,my)
      j=i+2; f(:,:,:,j)=f(:,:,:,j)+ampl &
        *spread(spread(exp(-(kx*x)**2)*(-2*x),2,my),3,mz)*kx**2 &
        *spread(spread(exp(-(ky*y)**2),1,mx),3,mz) &
        *spread(spread(exp(-abs(kz*z)),1,mx),2,my)
!
    endsubroutine bipolar
!***********************************************************************
    subroutine bipolar_restzero(ampl,f,i,kx,ky)
!
!  horizontal pattern with exponential decay (as initial condition)
!
!  24-may-09/axel: coded
!
      integer :: i,j
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,kx,ky
!
!  sets up a nearly force-free bipolar region
!
      if (lfirst_proc_z) then
        j=i+1; f(:,:,n1,j)=f(:,:,n1,j)+ampl &
          *spread(exp(-(kx*x)**2),2,my) &
          *spread(exp(-(ky*y)**2),1,mx)
        j=i+2; f(:,:,n1,j)=f(:,:,n1,j)+ampl &
          *spread(exp(-(kx*x)**2)*(-2*x),2,my)*kx**2 &
          *spread(exp(-(ky*y)**2),1,mx)
      endif
!
    endsubroutine bipolar_restzero
!***********************************************************************
    subroutine soundwave(ampl,f,i,kx,ky,kz)
!
!  sound wave (as initial condition)
!
!   2-aug-02/axel: adapted from Beltrami
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real,optional :: kx,ky,kz
      real :: ampl,k=1.,fac
!
!  wavenumber k
!
!  set x-dependent sin wave
!
      if (present(kx)) then
        k=kx; if (k==0) print*,'soundwave: k must not be zero!'; fac=sqrt(abs(ampl/k))
        if (ampl==0) then
          if (lroot) print*,'soundwave: ampl=0; kx=',k
        else
          if (lroot) print*,'soundwave: kx,i=',k,i
          f(:,:,:,i)=f(:,:,:,i)+fac*spread(spread(sin(k*x),2,my),3,mz)
        endif
      endif
!
!  set y-dependent sin wave field
!
      if (present(ky)) then
        k=ky; if (k==0) print*,'soundwave: k must not be zero!'; fac=sqrt(abs(ampl/k))
        if (ampl==0) then
          if (lroot) print*,'soundwave: ampl=0; ky=',k
        else
          if (lroot) print*,'soundwave: ky,i=',k,i
          f(:,:,:,i)=f(:,:,:,i)+fac*spread(spread(sin(k*y),1,mx),3,mz)
        endif
      endif
!
!  set z-dependent sin wave field
!
      if (present(kz)) then
        k=kz; if (k==0) print*,'soundwave: k must not be zero!'; fac=sqrt(abs(ampl/k))
        if (ampl==0) then
          if (lroot) print*,'soundwave: ampl=0; kz=',k
        else
          if (lroot) print*,'soundwave: kz,i=',k,i
          f(:,:,:,i)=f(:,:,:,i)+fac*spread(spread(sin(k*z),1,mx),2,my)
        endif
      endif
!
    endsubroutine soundwave
!***********************************************************************
    subroutine coswave(ampl,f,i,kx,ky,kz)
!
!  cosine wave (as initial condition)
!
!  14-nov-03/axel: adapted from sinwave
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real,optional :: kx,ky,kz
      real :: ampl,k=1.,fac
!
!  wavenumber k
!
!  set x-dependent cos wave
!
      if (present(kx)) then
        k=kx; if (k==0) print*,'coswave: k must not be zero!'; fac=ampl
        if (ampl==0) then
          if (lroot) print*,'coswave: ampl=0; kx=',k
        else
          if (lroot) print*,'coswave: kx,i=',k,i
          f(:,:,:,i)=f(:,:,:,i)+fac*spread(spread(cos(k*x),2,my),3,mz)
        endif
      endif
!
!  set y-dependent cos wave field
!
      if (present(ky)) then
        k=ky; if (k==0) print*,'coswave: k must not be zero!'; fac=ampl
        if (ampl==0) then
          if (lroot) print*,'coswave: ampl=0; ky=',k
        else
          if (lroot) print*,'coswave: ky,i=',k,i
          f(:,:,:,i)=f(:,:,:,i)+fac*spread(spread(cos(k*y),1,mx),3,mz)
        endif
      endif
!
!  set z-dependent cos wave field
!
      if (present(kz)) then
        k=kz; if (k==0) print*,'coswave: k must not be zero!'; fac=ampl
        if (ampl==0) then
          if (lroot) print*,'coswave: ampl=0; kz=',k
        else
          if (lroot) print*,'coswave: kz,i=',k,i
          f(:,:,:,i)=f(:,:,:,i)+fac*spread(spread(cos(k*z),1,mx),2,my)
        endif
      endif
!
    endsubroutine coswave
!***********************************************************************
    subroutine sph_constb(ampl,f,izero)
!
!  Antisymmetric (about equator) toroidal B field and zero poloidal
!  B field in spherical polar coordinate.
!
!  27-aug-09/dhruba: adapted from sinwave
!
      integer :: izero,l,m
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl
!
! set minimum and maximum values for r and theta
!
! set A_{phi} = 0 (izero should be iaa)
!
      f(:,:,:,izero+2)=0.
! set A_r and A_{\theta}
      do l=1,mx
        do m=1,my
          f(l,m,:,izero)=ampl*x(l)*y(m)
          f(l,m,:,izero+1)=ampl*x(l)
        enddo
      enddo
!
    endsubroutine sph_constb
!***********************************************************************
    subroutine hatwave(ampl,f,i,width,kx,ky,kz)
!
!  cosine wave (as initial condition)
!
!   9-jan-08/axel: adapted from coswave
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real,optional :: kx,ky,kz
      real :: ampl,k=1.,fac,width
!
!  wavenumber k
!
!  set x-dependent hat wave
!
      if (present(kx)) then
        k=kx; if (k==0) print*,'hatwave: k must not be zero!'; fac=.5*ampl
        if (ampl==0) then
          if (lroot) print*,'hatwave: ampl=0; kx=',k
        else
          if (lroot) print*,'hatwave: kx,i=',k,i
          f(:,:,:,i)=f(:,:,:,i)+fac*spread(spread(1+tanh(cos(k*x)/width),2,my),3,mz)
        endif
      endif
!
!  set y-dependent hat wave field
!
      if (present(ky)) then
        k=ky; if (k==0) print*,'hatwave: k must not be zero!'; fac=.5*ampl
        if (ampl==0) then
          if (lroot) print*,'hatwave: ampl=0; ky=',k
        else
          if (lroot) print*,'hatwave: ky,i=',k,i
          f(:,:,:,i)=f(:,:,:,i)+fac*spread(spread(1+tanh(5*cos(k*y)),1,mx),3,mz)
        endif
      endif
!
!  set z-dependent hat wave field
!
      if (present(kz)) then
        k=kz; if (k==0) print*,'hatwave: k must not be zero!'; fac=.5*ampl
        if (ampl==0) then
          if (lroot) print*,'hatwave: ampl=0; kz=',k
        else
          if (lroot) print*,'hatwave: kz,i=',k,i
          f(:,:,:,i)=f(:,:,:,i)+fac*spread(spread(1+tanh(5*cos(k*z)),1,mx),2,my)
        endif
      endif
!
    endsubroutine hatwave
!***********************************************************************
    subroutine sinwave(ampl,f,i,kx,ky,kz)
!
!  sine wave (as initial condition)
!
!  14-nov-03/axel: adapted from sound wave
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real,optional :: kx,ky,kz
      real :: ampl,k=1.,fac
!
!  wavenumber k
!
!  set x-dependent sin wave
!
      if (present(kx)) then
        k=kx; if (k==0) print*,'sinwave: k must not be zero!'; fac=ampl
        if (ampl==0) then
          if (lroot) print*,'sinwave: ampl=0; kx=',k
        else
          if (lroot) print*,'sinwave: kx,i=',k,i
          f(:,:,:,i)=f(:,:,:,i)+fac*spread(spread(sin(k*x),2,my),3,mz)
        endif
      endif
!
!  set y-dependent sin wave field
!
      if (present(ky)) then
        k=ky; if (k==0) print*,'sinwave: k must not be zero!'; fac=ampl
        if (ampl==0) then
          if (lroot) print*,'sinwave: ampl=0; ky=',k
        else
          if (lroot) print*,'sinwave: ky,i=',k,i
          f(:,:,:,i)=f(:,:,:,i)+fac*spread(spread(sin(k*y),1,mx),3,mz)
        endif
      endif
!
!  set z-dependent sin wave field
!
      if (present(kz)) then
        k=kz; if (k==0) print*,'sinwave: k must not be zero!'; fac=ampl
        if (ampl==0) then
          if (lroot) print*,'sinwave: ampl=0; kz=',k
        else
          if (lroot) print*,'sinwave: kz,i=',k,i
          f(:,:,:,i)=f(:,:,:,i)+fac*spread(spread(sin(k*z),1,mx),2,my)
        endif
      endif
!
    endsubroutine sinwave
!***********************************************************************
    subroutine sinwave_phase(f,i,ampl,kx,ky,kz,phase)
!
!  Sine wave (as initial condition)
!
!  23-jan-06/anders: adapted from sinwave.
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl, kx, ky, kz, phase
      integer :: i
!
!  Set sin wave
!
      if (lroot) print*, 'sinwave_phase: i, ampl, kx, ky, kz, phase=', &
          i, ampl, kx, ky, kz, phase
!
      do m=m1,m2; do n=n1,n2
        f(l1:l2,m,n,i) = f(l1:l2,m,n,i) + &
            ampl*sin(kx*x(l1:l2)+ky*y(m)+kz*z(n)+phase)
      enddo; enddo
!
    endsubroutine sinwave_phase
!***********************************************************************
    subroutine coswave_phase(f,i,ampl,kx,ky,kz,phase)
!
!  Cosine wave (as initial condition)
!
!  13-jun-06/anders: adapted from sinwave-phase.
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl, kx, ky, kz, phase
      integer :: i
!
!  Set cos wave
!
      if (lroot) print*, 'coswave_phase: i, ampl, kx, ky, kz, phase=', &
          i, ampl, kx, ky, kz, phase
!
      do m=m1,m2; do n=n1,n2
        f(l1:l2,m,n,i) = f(l1:l2,m,n,i) + &
            ampl*cos(kx*x(l1:l2)+ky*y(m)+kz*z(n)+phase)
      enddo; enddo
!
    endsubroutine coswave_phase
!***********************************************************************
    subroutine hawley_etal99a(ampl,f,i,Lxyz)
!
!  velocity perturbations as used by Hawley et al (1999, ApJ,518,394)
!
!  13-jun-05/maurice reyes: sent to axel via email
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx) :: funx
      real, dimension (my) :: funy
      real, dimension (mz) :: funz
      real, dimension (3) :: Lxyz
      real :: k1,k2,k3,k4,phi1,phi2,phi3,phi4,ampl
      integer :: i,l,m,n
!
!  velocity perturbations as used by Hawley et al (1999, ApJ,518,394)
!
      if (lroot) print*,'init_uu: hawley-et-al'
      k1=2.0*pi/(Lxyz(1))
      k2=4.0*pi/(Lxyz(1))
      k3=6.0*pi/(Lxyz(1))
      k4=8.0*pi/(Lxyz(1))
      phi1=k1*0.226818
      phi2=k2*0.597073
      phi3=k3*0.962855
      phi4=k4*0.762091
!
!  use l,m,n as loop variables; inner loop should be on first index
!
      funx=sin(k1*x+phi1)+sin(k2*x+phi2)+sin(k3*x+phi3)+sin(k4*x+phi4)
      funy=sin(k1*y+phi1)+sin(k2*y+phi2)+sin(k3*y+phi3)+sin(k4*y+phi4)
      funz=sin(k1*z+phi1)+sin(k2*z+phi2)+sin(k3*z+phi3)+sin(k4*z+phi4)
      do n=1,mz; do m=1,my; do l=1,mx
        f(l,m,n,i)=ampl*funx(l)*funy(m)*funz(n)
      enddo; enddo; enddo
!
    endsubroutine hawley_etal99a
!***********************************************************************
    subroutine stratification(f,strati_type)
!
!  Read mean stratification from "stratification.dat".
!
!   8-apr-03/axel: coded
!  23-may-04/anders: made structure for other input variables
!
      use EquationOfState, only: eoscalc,ilnrho_lnTT
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, parameter :: ntotal=nz*nprocz,mtotal=nz*nprocz+2*nghost
      real, dimension (mtotal) :: lnrho0,ss0,lnTT0
      real :: tmp,var1,var2
      logical :: exist
      integer :: stat
      character (len=labellen) :: strati_type
!
!  Read mean stratification and write into array.
!  If file is not found in run directory, search under trim(directory).
!
      inquire(file='stratification.dat',exist=exist)
      if (exist) then
        open(19,file='stratification.dat')
      else
        inquire(file=trim(directory)//'/stratification.ascii',exist=exist)
        if (exist) then
          open(19,file=trim(directory)//'/stratification.ascii')
        else
          call fatal_error('stratification','no input file')
        endif
      endif
!
!  Read data - first the entire stratification file.
!
      select case (strati_type)
      case ('lnrho_ss')
        do n=1,mtotal
          read(19,*,iostat=stat) tmp,var1,var2
          if (stat>=0) then
            if (ip<5) print*, 'stratification: z, var1, var2=', tmp, var1, var2
            if (ldensity) lnrho0(n)=var1
            if (lentropy) ss0(n)=var2
          else
            exit
          endif
        enddo
!
      case ('lnrho_lnTT')
        do n=1,mtotal
          read(19,*,iostat=stat) tmp,var1,var2
          if (stat>=0) then
            if (ip<5) print*, 'stratification: z, var1, var2=', tmp, var1, var2
            if (ldensity) lnrho0(n)=var1
            if (ltemperature) lnTT0(n)=var2
            if (lentropy) then
              call eoscalc(ilnrho_lnTT,var1,var2,ss=tmp)
              ss0(n)=tmp
            endif
          else
            exit
          endif
        enddo
      endselect
!
!  Select the right region for the processor afterwards.
!
      select case (n)
!
!  Without ghost zones.
!
      case (ntotal+1)
        if (lentropy) then
          do n=n1,n2
            f(:,:,n,ilnrho)=lnrho0(ipz*nz+n-nghost)
            f(:,:,n,iss)=ss0(ipz*nz+n-nghost)
          enddo
        endif
        if (ltemperature) then
          do n=n1,n2
            f(:,:,n,ilnrho)=lnrho0(ipz*nz+n-nghost)
            f(:,:,n,ilnTT)=lnTT0(ipz*nz+n-nghost)
          enddo
        endif
!
!  With ghost zones.
!
      case (mtotal+1)
        if (lentropy) then
          do n=1,mz
            f(:,:,n,ilnrho)=lnrho0(ipz*nz+n)
            f(:,:,n,iss)=ss0(ipz*nz+n)
          enddo
        endif
        if (ltemperature) then
          do n=1,mz
            f(:,:,n,ilnrho)=lnrho0(ipz*nz+n)
            f(:,:,n,ilnTT)=lnTT0(ipz*nz+n)
          enddo
        endif
!
      case default
        if (lroot) then
          print '(A,I4,A,I4,A,I4,A)','ERROR: The stratification file '// &
                'for this run is allowed to contain either',ntotal, &
                ' lines (without ghost zones) or more than',mtotal, &
                ' lines (with ghost zones). It does contain',n-1, &
                ' lines though.'
        endif
        call fatal_error('','')
!
      endselect
!
      close(19)
!
    endsubroutine stratification
!***********************************************************************
    subroutine stratification_x(f,strati_type)
!
!  read mean stratification from "stratification.dat"
!
!   02-mar-09/petri: adapted from stratification
!
      use EquationOfState, only: eoscalc,ilnrho_lnTT
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, parameter :: ntotal=nx*nprocx,mtotal=nx*nprocx+2*nghost
      real, dimension (mtotal) :: lnrho0,ss0,lnTT0
      real :: tmp,var1,var2
      logical :: exist
      integer :: stat
      character (len=labellen) :: strati_type
!
!  read mean stratification and write into array
!  if file is not found in run directory, search under trim(directory)
!
      inquire(file='stratification.dat',exist=exist)
      if (exist) then
        open(19,file='stratification.dat')
      else
        inquire(file=trim(directory)//'/stratification.ascii',exist=exist)
        if (exist) then
          open(19,file=trim(directory)//'/stratification.ascii')
        else
          call fatal_error('stratification','no input file')
        endif
      endif
!
!  read data
!  first the entire stratification file
!
      select case (strati_type)
      case ('lnrho_ss')
        do n=1,mtotal
          read(19,*,iostat=stat) tmp,var1,var2
          if (stat>=0) then
            if (ip<5) print*,"stratification: ",tmp,var1,var2
            if (ldensity) lnrho0(n)=var1
            if (lentropy) ss0(n)=var2
          else
            exit
          endif
        enddo
!
      case ('lnrho_lnTT')
        do n=1,mtotal
          read(19,*,iostat=stat) tmp,var1,var2
          if (stat>=0) then
            if (ip<5) print*,"stratification: ",tmp,var1,var2
            if (ldensity) lnrho0(n)=var1
            if (ltemperature) lnTT0(n)=var2
            if (lentropy) then
              call eoscalc(ilnrho_lnTT,var1,var2,ss=tmp)
              ss0(n)=tmp
            endif
          else
            exit
          endif
        enddo
      endselect
!
!  select the right region for the processor afterwards
!
      select case (n)
  !
  !  without ghost zones
  !
      case (ntotal+1)
        if (lentropy) then
          do n=l1,l2
            f(n,:,:,ilnrho)=lnrho0(ipx*nx+n-nghost)
            f(n,:,:,iss)=ss0(ipx*nx+n-nghost)
          enddo
        endif
        if (ltemperature) then
          do n=l1,l2
            f(n,:,:,ilnrho)=lnrho0(ipx*nx+n-nghost)
            f(n,:,:,ilnTT)=lnTT0(ipx*nx+n-nghost)
          enddo
        endif
  !
  !  with ghost zones
  !
      case (mtotal+1)
        if (lentropy) then
          do n=1,mx
            f(n,:,:,ilnrho)=lnrho0(ipx*nx+n)
            f(n,:,:,iss)=ss0(ipx*nx+n)
          enddo
        endif
        if (ltemperature) then
          do n=1,mx
            f(n,:,:,ilnrho)=lnrho0(ipx*nx+n)
            f(n,:,:,ilnTT)=lnTT0(ipx*nx+n)
          enddo
        endif
!
      case default
        if (lroot) then
          print '(A,I4,A,I4,A,I4,A)','ERROR: The stratification file '// &
                'for this run is allowed to contain either',ntotal, &
                ' lines (without ghost zones) or more than',mtotal, &
                ' lines (with ghost zones). It does contain',n-1, &
                ' lines though.'
        endif
        call fatal_error('','')
!
      endselect
!
      close(19)
!
    endsubroutine stratification_x
!***********************************************************************
    subroutine planet_hc(ampl,f,eps,radius,gamma,cs20,rho0,width)
!
!  Ellipsoidal planet solution (Goldreich, Narayan, Goodman 1987)
!
!   6-jul-02/axel: coded
!  22-feb-03/axel: fixed 3-D background solution for enthalpy
!  26-Jul-03/anders: Revived from June 1 version
!
      use Mpicomm, only: mpireduce_sum, mpibcast_real
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: hh, xi
      real, dimension (mz) :: hz
      real :: delS,ampl,sigma2,sigma,delta2,delta,eps,radius,a_ell,b_ell,c_ell
      real :: gamma,cs20,gamma_m1,eps2,radius2,width
      real :: lnrhosum_thisbox,rho0
      real, dimension (1) :: lnrhosum_thisbox_tmp,lnrhosum_wholebox
      integer :: l
!
!  calculate sigma
!
      if (lroot) print*,'planet_hc: qshear,eps=',qshear,eps
      eps2=eps**2
      radius2=radius**2
      sigma2=2*qshear/(1.-eps2)
      if (sigma2<0.) then
        if (lroot) print*, &
          'planet_hc: sigma2<0 not allowed; choose another value of eps_planet'
      else
        sigma=sqrt(sigma2)
      endif
!
!  calculate delta
!
      delta2=(2.-sigma)*sigma
      if (lroot) print*,'planet_hc: sigma,delta2,radius=',sigma,delta2,radius
      if (delta2<=0.) then
        if (lroot) print*,'planet_hc: delta2<=0 not allowed'
      else
        delta=sqrt(delta2)
      endif
!
!  calculate gamma_m1
!
      gamma_m1=gamma-1.
      if (lroot) print*,'planet_hc: gamma=',gamma
!
!  ellipse parameters
!
      b_ell = radius
      a_ell = radius/eps
      c_ell = radius*delta
      if (lroot) print*,'planet_hc: Ellipse axes (b_ell,a_ell,c_ell)=', &
          b_ell,a_ell,c_ell
      if (lroot) print*,"planet_hc: integrate hot corona"
!
!  xi=1 inside vortex, and 0 outside
!
      do n=n1,n2; do m=m1,m2
        hh=0.5*delta2*Omega**2*(radius2-x(l1:l2)**2-eps2*y(m)**2)-.5*Omega**2*z(n)**2
        xi=0.5+0.5*tanh(hh/width)
!
!  Calculate velocities (Kepler speed subtracted)
!
        f(l1:l2,m,n,iux)=   eps2*sigma *Omega*y(m)    *xi
        f(l1:l2,m,n,iuy)=(qshear-sigma)*Omega*x(l1:l2)*xi
        if (lentropy) f(l1:l2,m,n,iss)=-log(ampl)*xi
      enddo; enddo
!
      do m=m1,m2; do l=l1,l2
!
!  add continuous vertical stratification to horizontal planet solution
!  NOTE: if width is too small, the vertical integration below may fail.
!
        hz(n2)=1.0  !(initial condition)
        do n=n2-1,n1,-1
          delS=f(l,m,n+1,iss)-f(l,m,n,iss)
          hz(n)=(hz(n+1)*(1.0-0.5*delS)+ &
               Omega**2*0.5*(z(n)+z(n+1))*dz)/(1.0+0.5*delS)
        enddo
!
!  calculate density, depending on what gamma is
!
        if (lentropy) then
          f(l,m,n1:n2,ilnrho)= &
               (log(gamma_m1*hz(n1:n2)/cs20)-gamma*f(l,m,n1:n2,iss))/gamma_m1
          if (lroot) &
            print*,'planet_hc: planet solution with entropy for gamma=',gamma
        else
          if (gamma==1.) then
            f(l,m,n1:n2,ilnrho)=hz(n1:n2)/cs20
            if (lroot) print*,'planet_hc: planet solution for gamma=1'
          else
            f(l,m,n1:n2,ilnrho)=log(gamma_m1*hz(n1:n2)/cs20)/gamma_m1
            if (lroot) print*,'planet_hc: planet solution for gamma=',gamma
          endif
        endif
!
      enddo; enddo
!
      if (gamma_m1<0. .and. lroot) &
          print*,'planet_hc: must have gamma>1 for planet solution'
!
!  Use average density of box as unit density
!
      lnrhosum_thisbox = sum(f(l1:l2,m1:m2,n1:n2,ilnrho))
      if (ip<14) &
        print*,'planet_hc: iproc,lnrhosum_thisbox=',iproc,lnrhosum_thisbox
!
!  Must put sum_thisbox in 1-dimensional array
!
      lnrhosum_thisbox_tmp = (/ lnrhosum_thisbox /)
!
!  Add sum_thisbox up for all processors, deliver to root
!
      call mpireduce_sum(lnrhosum_thisbox_tmp,lnrhosum_wholebox,1)
      if (lroot .and. ip<14) &
        print*,'planet_hc: lnrhosum_wholebox=',lnrhosum_wholebox
!
!  Calculate <rho> and send to all processors
!
      if (lroot) rho0 = exp(-lnrhosum_wholebox(1)/nwgrid)
      call mpibcast_real(rho0,1)
      if (ip<14) print*,'planet_hc: iproc,rho0=',iproc,rho0
!
!  Multiply density by rho0 (divide by <rho>)
!
      f(l1:l2,m1:m2,n1:n2,ilnrho) = f(l1:l2,m1:m2,n1:n2,ilnrho) + log(rho0)
!
    endsubroutine planet_hc
!***********************************************************************
    subroutine planet(rbound,f,eps,radius,gamma,cs20,rho0,width,hh0)
!
!  Cylindrical planet solution (Goldreich, Narayan, Goodman 1987)
!
!   jun-03/anders: coded (adapted from old 'planet', now 'planet_hc')
!
      use Mpicomm, only: mpireduce_sum, mpibcast_real
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: hh, xi, r_ell
      real :: rbound,sigma2,sigma,delta2,delta,eps,radius
      real :: gamma,eps2,radius2,width,a_ell,b_ell,c_ell
      real :: gamma_m1,ztop,cs20,hh0
      real :: lnrhosum_thisbox,rho0
      real, dimension (1) :: lnrhosum_thisbox_tmp,lnrhosum_wholebox
!
!  calculate sigma
!
      if (lroot) print*,'planet: qshear,eps=',qshear,eps
      eps2=eps**2
      radius2=radius**2
      sigma2=2*qshear/(1.-eps2)
      if (sigma2<0. .and. lroot) then
        print*, &
            'planet: sigma2<0 not allowed; choose another value of eps_planet'
      else
        sigma=sqrt(sigma2)
      endif
!
      gamma_m1=gamma-1.
!
!  calculate delta
!
      delta2=(2.-sigma)*sigma
      if (lroot) print*,'planet: sigma,delta2,radius=',sigma,delta2,radius
      if (delta2<0. .and. lroot) then
        print*,'planet: delta2<0 not allowed'
      else
        delta=sqrt(delta2)
      endif
!
      ztop=z0+lz
      if (lroot) print*,'planet: ztop=', ztop
!
!  Cylinder vortex 3-D solution (b_ell along x, a_ell along y)
!
      b_ell = radius
      a_ell = radius/eps
      c_ell = radius*delta
      if (lroot) print*,'planet: Ellipse axes (b_ell,a_ell,c_ell)=', &
          b_ell,a_ell,c_ell
      if (lroot) print*,'planet: width,rbound',width,rbound
!
      do n=n1,n2; do m=m1,m2
        r_ell = sqrt(x(l1:l2)**2/b_ell**2+y(m)**2/a_ell**2)
!
!  xi is 1 inside vortex and 0 outside
!
        xi = 1/(exp((1/width)*(r_ell-rbound))+1.0)
!
!  Calculate enthalpy inside vortex
!
        hh = 0.5*delta2*Omega**2*(radius2-x(l1:l2)**2-eps2*y(m)**2) &
             -0.5*Omega**2*z(n)**2 + 0.5*Omega**2*ztop**2 + hh0
!
!  Calculate enthalpy outside vortex
!
        where (r_ell>1.0) hh=-0.5*Omega**2*z(n)**2 + 0.5*Omega**2*ztop**2 + hh0
!
!  Calculate velocities (Kepler speed subtracted)
!
        f(l1:l2,m,n,iux)=   eps2*sigma *Omega*y(m)*xi
        f(l1:l2,m,n,iuy)=(qshear-sigma)*Omega*x(l1:l2)*xi
!
!  calculate density, depending on what gamma is
!
        if (lentropy) then
          f(l1:l2,m,n,ilnrho)=(log(gamma_m1*hh/cs20)-gamma*f(l1:l2,m,n,iss))/gamma_m1
        else
          if (gamma==1.) then
            f(l1:l2,m,n,ilnrho) = hh/cs20
          else
            f(l1:l2,m,n,ilnrho) = log(gamma_m1*hh/cs20)/gamma_m1
          endif
        endif
      enddo; enddo
!
!  Use average density of box as unit density
!
      lnrhosum_thisbox = sum(f(l1:l2,m1:m2,n1:n2,ilnrho))
      if (ip<14) &
        print*,'planet_hc: iproc,lnrhosum_thisbox=',iproc,lnrhosum_thisbox
!
!  Must put sum_thisbox in 1-dimensional array
!
      lnrhosum_thisbox_tmp = (/ lnrhosum_thisbox /)
!
!  Add sum_thisbox up for all processors, deliver to root
!
      call mpireduce_sum(lnrhosum_thisbox_tmp,lnrhosum_wholebox,1)
      if (lroot .and. ip<14) &
          print*,'planet_hc: lnrhosum_wholebox=',lnrhosum_wholebox
!
!  Calculate <rho> and send to all processors
!
      if (lroot) rho0 = exp(-lnrhosum_wholebox(1)/nwgrid)
      call mpibcast_real(rho0,1)
      if (ip<14) print*,'planet_hc: iproc,rho0=',iproc,rho0
!
!  Multiply density by rho0 (divide by <rho>)
!
      f(l1:l2,m1:m2,n1:n2,ilnrho) = f(l1:l2,m1:m2,n1:n2,ilnrho) + log(rho0)
!
    endsubroutine planet
!***********************************************************************
    subroutine vortex_2d(f,b_ell,width,rbound)
!
!  Ellipsoidal planet solution (Goldreich, Narayan, Goodman 1987)
!
!   8-jun-04/anders: adapted from planet
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: r_ell, xi
      real :: sigma,eps_ell,a_ell,b_ell,width,rbound
!
!  calculate sigma
!
      eps_ell=0.5
      if (lroot) print*,'vortex_2d: qshear,eps_ell=',qshear,eps_ell
      sigma=sqrt(2*qshear/(1.-eps_ell**2))
!
!  ellipse parameters
!
      a_ell = b_ell/eps_ell
      if (lroot) print*,'vortex_2d: Ellipse axes (b_ell,a_ell)=', b_ell,a_ell
!
!  Limit vortex to within r_ell
!
      do n=n1,n2; do m=m1,m2
        r_ell = sqrt(x(l1:l2)**2/b_ell**2+y(m)**2/a_ell**2)
        xi = 1./(exp((1/width)*(r_ell-rbound))+1.)
!
!  Calculate velocities (Kepler speed subtracted)
!
        f(l1:l2,m,n,iux)=eps_ell**2*sigma*Omega*y(m)*xi
        f(l1:l2,m,n,iuy)=(qshear-sigma)  *Omega*x(l1:l2)*xi
      enddo; enddo
!
    endsubroutine vortex_2d
!***********************************************************************
    subroutine baroclinic(f,gamma,dlnrhobdx,co1_ss,co2_ss,cs20)
!
!  Remark 4-dez-09 bing: this routine is nowhere called! Is it still in use?
!
!  Baroclinic shearing sheet initial condition
!  11-nov-03/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: sz,I_int
      real :: gamma,dlnrhobdx,co1_ss,co2_ss,cs20
!
!  Specify vertical entropy and integral of exp(-sz/cp)*z
!
      do n=n1,n2; do m=m1,m2
        if (co1_ss/=0.0 .and. co2_ss==0.0) then
          if (lroot) print*,'baroclinic: sz =', co1_ss,'*abs(z(n))'
          sz = co1_ss*abs(z(n))
          I_int = 1/(co1_ss**2)*( 1 - exp(-sz) * (1+co1_ss*abs(z(n))) )
        elseif (co1_ss==0.0 .and. co2_ss/=0.0) then
          if (lroot) print*,'baroclinic: sz =', co2_ss,'*zz**2'
          sz = co2_ss*z(n)**2
          I_int = -1/(2*co2_ss)*( exp(-co2_ss*z(n)**2)-1 )
        elseif (lroot) then
          print*, 'baroclinic: no valid sz specified'
        endif
!
        f(l1:l2,m,n,iss) = sz
!
!  Solution to hydrostatic equlibrium in the z-direction
!
        f(l1:l2,m,n,ilnrho) = 1/(gamma-1)*log( (1-gamma)/cs20 * I_int + 1 ) - sz
!
!  Toroidal velocity comes from hyd. stat. eq. equ. in the x-direction
!
        f(l1:l2,m,n,iuy) = cs20/(2*Omega)*exp( gamma*f(l1:l2,m,n,iss) + &
            (gamma-1)*f(l1:l2,m,n,ilnrho) ) * dlnrhobdx/gamma
!
      enddo; enddo
!
    endsubroutine baroclinic
!***********************************************************************
    subroutine crazy(ampl,f,i)
!
!  A crazy initial condition
!  (was useful to initialize all points with finite values)
!
!  19-may-02/axel: coded
!
      integer :: i,j
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl
!
      if (lroot) print*, 'crazy: sinusoidal magnetic field: for debugging purposes'
      j=i; f(:,:,:,j)=f(:,:,:,j)+ampl*&
        spread(spread(sin(2*x),2,my),3,mz)*&
        spread(spread(sin(3*y),1,mx),3,mz)*&
        spread(spread(cos(1*z),1,mx),2,my)
      j=i+1; f(:,:,:,j)=f(:,:,:,j)+ampl*&
        spread(spread(sin(5*x),2,my),3,mz)*&
        spread(spread(sin(1*y),1,mx),3,mz)*&
        spread(spread(cos(2*z),1,mx),2,my)
      j=i+2; f(:,:,:,j)=f(:,:,:,j)+ampl*&
        spread(spread(sin(3*x),2,my),3,mz)*&
        spread(spread(sin(4*y),1,mx),3,mz)*&
        spread(spread(cos(2*z),1,mx),2,my)
!
    endsubroutine crazy
!***********************************************************************
    subroutine strange(ampl,f,i)
!
!  A strange initial condition !
!
!  24-april-09/dhruba: coded
!
      integer :: i,ix,iy
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl
!      real :: arg1,arg2
!
      if (lroot) print*, 'strange: magnetic field that satisfies open', &
          'vf condn: for debugging purposes'
      do ix=l1,l2
        do iy=m1,m2
!          arg1= pi*(x(ix)-xyz0(1))/(xyz1(1)-xyz0(1))
!          arg2= pi*(y(iy)-xyz0(2))/(xyz1(2)-xyz0(2))
          f(ix,iy,:,i)=ampl*cos(y(iy))
          f(ix,iy,:,i+1)=0.
          f(ix,iy,:,i+2)=0.
       enddo
     enddo
!
    endsubroutine strange
!***********************************************************************
    subroutine htanh(ampl,f,i1,eps)
!
!  Horizontal flux tanh (for vector potential, or passive scalar)
!
!  30-may-11/axel: tanh layer
!
      integer :: i1
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,eps
!
      if (ampl==0) then
        f(:,:,:,i1)=0
        if (lroot) print*,'htanh: set variable to zero; i1=',i1
      else
        if (lroot) then
          print*,'htanh: implement y-dependent flux tanh in xz-plane; i1=',i1
          print*,'htanh: eps=',eps
        endif
!
! completely quenched "gaussian"
!
        do n=n1,n2; do m=m1,m2
          f(l1:l2,m,n,i1)=ampl*alog(cosh(eps*y(m)))
        enddo; enddo
      endif
!
    endsubroutine htanh
!***********************************************************************
    subroutine htube(ampl,f,i1,i2,radius,eps,center1_x,center1_z)
!
!  Horizontal flux tube (for vector potential, or passive scalar)
!
!   7-jun-02/axel+vladimir: coded
!  11-sep-02/axel: allowed for scalar field (if i1=i2)
!
      integer :: i1,i2
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: tmp,tube_radius_sqr !,modulate
      real :: ampl,radius,eps
      real :: center1_x,center1_z
!
      if (ampl==0) then
        f(:,:,:,i1:i2)=0
        if (lroot) print*,'htube: set variable to zero; i1,i2=',i1,i2
      else
        if (lroot) then
          print*,'htube: implement y-dependent flux tube in xz-plane; i1,i2=',i1,i2
          print*,'htube: radius,eps=',radius,eps
          if (i1==i2) then
            print*,'htube: set scalar'
          else
            print*,'htube: set vector'
          endif
        endif
!
! completely quenched "gaussian"
!
        do n=n1,n2; do m=m1,m2
          tube_radius_sqr=(x(l1:l2)-center1_x)**2+(z(n)-center1_z)**2
!         tmp=.5*ampl/modulate*exp(-tube_radius_sqr)/&
!                   (max((radius*modulate)**2-tube_radius_sqr,1e-6))
          tmp=ampl/(1.+tube_radius_sqr/radius**2)
!
!  check whether vector or scalar
!
          if (i1==i2) then
            f(l1:l2,m,n,i1)=tmp
          elseif (i1+2==i2) then
            f(l1:l2,m,n,i1 )=+(z(n)-center1_z)*tmp
            f(l1:l2,m,n,i1+1)=tmp*eps
            f(l1:l2,m,n,i1+2)=-(x(l1:l2)-center1_x)*tmp
          else
            if (lroot) print*,'htube: bad value of i2=',i2
          endif
!
        enddo; enddo
      endif
!
    endsubroutine htube
!***********************************************************************
    subroutine htube_x(ampl,f,i1,i2,radius,eps,center1_y,center1_z)
!
!  Horizontal flux tube pointing in the x-direction
!  (for vector potential, or passive scalar)
!
!  14-apr-09/axel: adapted from htube, used in paper with Violaine Auger
!
      integer :: i1,i2
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: tmp,modulate,tube_radius_sqr
      real :: ampl,radius,eps,kx
      real :: center1_y,center1_z
!
      if (ampl==0) then
        f(:,:,:,i1:i2)=0
        if (lroot) print*,'htube_x: set variable to zero; i1,i2=',i1,i2
      else
        kx=2*pi/Lx
        if (lroot) then
          print*,'htube_x: implement y-dependent flux tube in xz-plane; i1,i2=',i1,i2
          print*,'htube_x: radius,eps=',radius,eps
        endif
!
!  modulation pattern
!
        if (eps==0.) then
          modulate=1.
        else
          modulate=1.+eps*cos(kx*x(l1:l2))
        endif
!
! completely quenched "gaussian"
!
        do n=n1,n2; do m=m1,m2
          tube_radius_sqr=(y(m)-center1_y)**2+(z(n)-center1_z)**2
          tmp=modulate/(1.+tube_radius_sqr/radius**2)
!
!  check whether vector or scalar
!
          if (i1==i2) then
            if (lroot.and.ip<10) print*,'htube_x: set scalar'
            f(l1:l2,m,n,i1)=tmp
          elseif (i1+2==i2) then
            if (lroot.and.ip<10) print*,'htube_x: set vector'
            f(l1:l2,m,n,i1  )=+0.
            f(l1:l2,m,n,i1+1)=-(z(n)-center1_z)*tmp
            f(l1:l2,m,n,i1+2)=+(y(m)-center1_y)*tmp
         else
            if (lroot) print*,'htube_x: bad value of i2=',i2
          endif
!
        enddo; enddo
      endif
!
    endsubroutine htube_x
!***********************************************************************
    subroutine htube_erf(ampl,f,i1,i2,a,eps,center1_x,center1_z,width)
!
!  Horizontal flux tube (for vector potential) which gives
!  error-function border profile for the magnetic field., or passive scalar
!
!   18-mar-09/dhruba: adapted from htube
!
      use Sub, only: erfunc
!
      integer :: i1,i2,l
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,a,eps,width,tmp,radius,r_minus_a
      real :: center1_x,center1_z
!
      if (ampl==0) then
        f(:,:,:,i1:i2)=0
        if (lroot) print*,'htube_erf: set variable to zero; i1,i2=',i1,i2
      else
        if (lroot) then
          print*,'htube_erf: implement y-dependent flux tube in xz-plane; i1,i2=',i1,i2
          print*,'htube_erf: radius,eps=',a,eps
          if (i1==i2) then
            print*,'htube_erf: set scalar'
          else
            print*,'htube_erf: set vector'
          endif
        endif
!
! An integral of error function.
!
        do n=n1,n2; do l=l1,l2
          radius= sqrt((x(l)-center1_x)**2+(z(n)-center1_z)**2)
          r_minus_a= radius -a
          if (radius > tini) then
!            tmp = (-(exp(-width*a_minus_r**2))/(4.*sqrt(pi)*width) +  &
!                radius*(1+erfunc(width*a_minus_r))/4. + &
!                2*a*(exp(-(a**2)*(width**2)) - exp(-(a_minus_r**2)*(width**2)))/(8.*radius*width) + &
!                (1+2*(a**2)*(width**2))*(erfunc(a*width) - erfunc(width*a_minus_r))/(8.*radius*width**2))/radius
            tmp = (radius**2/2-0.5*(radius+a)*(r_minus_a*erfunc(r_minus_a/width)+width*exp(-r_minus_a**2/width**2)/sqrt(pi))&
                  +0.25*width**2*erfunc(r_minus_a/width)+0.5*a*(-a*erfunc(-a/width)+width*exp(-a**2/width**2)/sqrt(pi))&
                  -0.25*width**2*erfunc(-a/width))/radius**2
          else
            tmp = 0
            write(*,*) 'radius <  tini',radius,tini
          endif
!
!  check whether vector or scalar
!
          if (i1==i2) then
            f(l,:,n,i1)=tmp
          elseif (i1+2==i2) then
            f(l,:,n,i1 )=f(l,:,n,i1)-(z(n)-center1_z)*tmp*ampl
            f(l,:,n,i1+1)=f(l,:,n,i1+1)+tmp*eps
            f(l,:,n,i1+2)=f(l,:,n,i1+2)+(x(l)-center1_x)*tmp*ampl
         else
            if (lroot) print*,'htube_erf: bad value of i2=',i2
          endif
!
        enddo; enddo
      endif
!
    endsubroutine htube_erf
!***********************************************************************
    subroutine htube2(ampl,f,i1,i2,radius,epsilon_nonaxi)
!
!  Horizontal flux tube (for vector potential, or passive scalar)
!
!   7-jun-02/axel+vladimir: coded
!  11-sep-02/axel: allowed for scalar field (if i1=i2)
!
      integer :: i1,i2
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: tmp,modulate
      real :: ampl,radius,epsilon_nonaxi,ky
!
      if (ampl==0) then
        f(:,:,:,i1:i2)=0
        if (lroot) print*,'htube2: set variable to zero; i1,i2=',i1,i2
      else
        ky=2*pi/Ly
        if (lroot) then
          print*,'htube2: implement y-dependent flux tube in xz-plane; i1,i2=',i1,i2
          print*,'htube2: radius,epsilon_nonaxi=',radius,epsilon_nonaxi
        endif
!
!  constant, when epsilon_nonaxi; otherwise modulation about zero
!
        do n=n1,n2; do m=m1,m2
          if (epsilon_nonaxi==0) then
            modulate(:)=1.0
          else
            modulate=epsilon_nonaxi*sin(ky*y(m))
          endif
!
! completely quenched "gaussian"
!
          tmp=.5*ampl*modulate*exp(-(x(l1:l2)**2+z(n)**2)/radius**2)
!
!  check whether vector or scalar
!
          if (i1==i2) then
            if (lroot .and. n==n1 .and. m==m1) print*,'htube2: set scalar'
            f(l1:l2,m,n,i1)=tmp
          elseif (i1+2==i2) then
            if (lroot .and. n==n1 .and. m==m1) print*,'htube2: set vector'
            f(l1:l2,m,n,i1 )=+z(n)*tmp
            f(l1:l2,m,n,i1+1)=0.
            f(l1:l2,m,n,i1+2)=-x(l1:l2)*tmp
          else
            if (lroot) print*,'htube2: bad value of i2=',i2
          endif
        enddo; enddo
      endif
!
    endsubroutine htube2
!***********************************************************************
    subroutine magsupport(ampl,f,gravz,cs0,rho0)
!
!  magnetically supported horizontal flux layer
!  (for aa):  By^2 = By0^2 * exp(-z/H),
!  where H=2*cs20/abs(gravz) and ampl=cs0*sqrt(2*rho0)
!  should be used when invoking this routine.
!  Here, ampl=pmag/pgas.
!
!   7-dec-02/axel: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,H,A0,gravz,cs0,rho0,lnrho0
!
      if (ampl==0) then
        if (lroot) print*,'magsupport: do nothing'
      else
        lnrho0=log(rho0)
        H=(1+ampl)*cs0**2/abs(gravz)
        A0=-2*H*ampl*cs0*sqrt(2*rho0)
        if (lroot) print*,'magsupport: H,A0=',H,A0
        do n=n1,n2; do m=m1,m2
          f(l1:l2,m,n,iaa)=A0*exp(-.5*z(n)/H)
          f(l1:l2,m,n,ilnrho)=lnrho0-z(n)/H
        enddo; enddo
      endif
!
    endsubroutine magsupport
!***********************************************************************
    subroutine hfluxlayer(ampl,f,i,zflayer,width)
!
!  Horizontal flux layer (for vector potential)
!
!  19-jun-02/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,zflayer,width
!
      if (ampl==0) then
        f(:,:,:,i:i+2)=0
        if (lroot) print*,'hfluxlayer: set variable to zero; i=',i
      else
        if (lroot) print*,'hfluxlayer: horizontal flux layer; i=',i
        if ((ip<=16).and.lroot) print*,'hfluxlayer: ampl,width=',ampl,width
        do n=n1,n2; do m=m1,m2
          f(l1:l2,m,n,i  )=0.0
          f(l1:l2,m,n,i+1)=ampl*tanh((z(n)-zflayer)/width)
          f(l1:l2,m,n,i+2)=0.0
        enddo; enddo
      endif
!
    endsubroutine hfluxlayer
!***********************************************************************
    subroutine hfluxlayer_y(ampl,f,i,zflayer,width)
!
!  Horizontal flux layer (for vector potential)
!
!  09-apr-10/piyali: copied from hfluxlayer
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,zflayer,width
!
      if (ampl==0) then
        f(:,:,:,i:i+2)=0
        if (lroot) print*,'hfluxlayer-y: set variable to zero; i=',i
      else
        if (lroot) print*,'hfluxlayer-y: horizontal flux layer; i=',i
        if ((ip<=16).and.lroot) print*,'hfluxlayer-y: ampl,width=',ampl,width
        do n=n1,n2; do m=m1,m2
          f(l1:l2,m,n,i  )=ampl*tanh((z(n)-zflayer)/width)
          f(l1:l2,m,n,i+1)=0.0
          f(l1:l2,m,n,i+2)=0.0
        enddo; enddo
      endif
!
    endsubroutine hfluxlayer_y
!***********************************************************************
    subroutine vfluxlayer(ampl,f,i,xflayer,width)
!
!  Vertical flux layer (for vector potential)
!
!  22-mar-04/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,xflayer,width
!
      if (ampl==0) then
        f(:,:,:,i:i+2)=0
        if (lroot) print*,'hfluxlayer: set variable to zero; i=',i
      else
        if (lroot) print*,'hfluxlayer: horizontal flux layer; i=',i
        if ((ip<=16).and.lroot) print*,'hfluxlayer: ampl,width=',ampl,width
        do n=n1,n2; do m=m1,m2
          f(l1:l2,m,n,i  )=0.0
          f(l1:l2,m,n,i+1)=0.0
          f(l1:l2,m,n,i+2)=-ampl*tanh((x(l1:l2)-xflayer)/width)
        enddo; enddo
      endif
!
    endsubroutine vfluxlayer
!***********************************************************************
    subroutine arcade_x(ampl,f,i,kx,kz)
!
!  Arcade-like structures around x=0
!
!  17-jun-04/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,kx,kz,zmid
!
      if (ampl==0) then
        f(:,:,:,i:i+2)=0
        if (lroot) print*,'expcos_x: set variable to zero; i=',i
      else
        zmid=.5*(xyz0(3)+xyz1(3))
        if ((ip<=16).and.lroot) then
          print*,'arcade_x: i,zmid=',i,zmid
          print*,'arcade_x: ampl,kx,kz=',ampl,kx,kz
        endif
!
        do n=n1,n2; do m=m1,m2
!         f(l1:l2,m,n,i+1)=f(l1:l2,m,n,i+1)+ampl*exp(-.5*(kx*x(l1:l2))**2)* &
!           cos(min(abs(kz*(z(n)-zmid)),.5*pi))
          f(l1:l2,m,n,i+1)=f(l1:l2,m,n,i+1) &
            +ampl*cos(kx*x(l1:l2))*exp(-abs(kz*z(n)))
        enddo; enddo
!
      endif
!
    endsubroutine arcade_x
!***********************************************************************
    subroutine halfcos_x(ampl,f,i)
!
!  Uniform B_x field (for vector potential)
!
!  19-jun-02/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,kz,zbot
!
      if (ampl==0) then
        f(:,:,:,i:i+2)=0
        if (lroot) print*,'halfcos_x: set variable to zero; i=',i
      else
        print*,'halscos_x: half cosine x-field ; i=',i
        kz=0.5*pi/Lz
        zbot=xyz0(3)
        ! ztop=xyz0(3)+Lxyz(3)
        if ((ip<=16).and.lroot) print*,'halfcos_x: ampl,kz=',ampl,kz
        do n=n1,n2; do m=m1,m2
          f(l1:l2,m,n,i  )=0.0
          f(l1:l2,m,n,i+1)=-ampl*sin(kz*(z(n)-zbot))
          f(l1:l2,m,n,i+2)=0.0
        enddo; enddo
      endif
!
    endsubroutine halfcos_x
!
!***********************************************************************
    subroutine uniform_x(ampl,f,i)
!
!  Uniform B_x field (for vector potential)
!
!  19-jun-02/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl
!
      if (ampl==0) then
        f(:,:,:,i:i+2)=0
        if (lroot) print*,'uniform_x: set variable to zero; i=',i
      else
        print*,'uniform_x: uniform x-field ; i=',i
        if ((ip<=16).and.lroot) print*,'uniform_x: ampl=',ampl
        do n=n1,n2; do m=m1,m2
          f(l1:l2,m,n,i  )=0.0
          f(l1:l2,m,n,i+1)=-ampl*z(n)
          f(l1:l2,m,n,i+2)=0.0
        enddo; enddo
      endif
!
    endsubroutine uniform_x
!***********************************************************************
    subroutine uniform_y(ampl,f,i)
!
!  Uniform B_y field (for vector potential)
!
!  27-jul-02/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl
!
      if (ampl==0) then
        f(:,:,:,i:i+2)=0
        if (lroot) print*,'uniform_y: set variable to zero; i=',i
      else
        print*,'uniform_y: uniform y-field ; i=',i
        if ((ip<=16).and.lroot) print*,'uniform_y: ampl=',ampl
        do n=n1,n2; do m=m1,m2
          f(l1:l2,m,n,i  )=ampl*z(n)
          f(l1:l2,m,n,i+1)=0.0
          f(l1:l2,m,n,i+2)=0.0
        enddo; enddo
      endif
!
    endsubroutine uniform_y
!***********************************************************************
    subroutine uniform_z(ampl,f,i)
!
!  Uniform B_z field (for vector potential)
!
!  24-jul-03/axel: adapted from uniform_x
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl
!
      if (ampl==0) then
        f(:,:,:,i:i+2)=0
        if (lroot) print*,'uniform_z: set variable to zero; i=',i
      else
        print*,'uniform_z: uniform z-field ; i=',i
        if ((ip<=16).and.lroot) print*,'uniform_z: ampl=',ampl
        do n=n1,n2; do m=m1,m2
          if (lcartesian_coords) then
            f(l1:l2,m,n,i  )=0.0
            f(l1:l2,m,n,i+1)=+ampl*x(l1:l2)
            f(l1:l2,m,n,i+2)=0.0
          elseif (lcylindrical_coords) then
            f(l1:l2,m,n,i  )=0.0
            f(l1:l2,m,n,i+1)=-ampl*x(l1:l2)*y(m)
            f(l1:l2,m,n,i+2)=0.0
          endif
        enddo; enddo
      endif
!
    endsubroutine uniform_z
!***********************************************************************
    subroutine uniform_phi(ampl,f,i)
!
!  Uniform B_phi field (for vector potential)
!
!  27-jul-02/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: rr
      real :: ampl
!
      if (ampl==0) then
        f(:,:,:,i:i+2)=0
        if (lroot) print*,'uniform_phi: set variable to zero; i=',i
      else
        print*,'uniform_phi: uniform phi-field ; i=',i
        if ((ip<=16).and.lroot) print*,'uniform_phi: ampl=',ampl
        do n=n1,n2; do m=m1,m2
          rr=sqrt(x(l1:l2)**2+y(m)**2)
          f(l1:l2,m,n,i  )=0.0
          f(l1:l2,m,n,i+1)=0.0
          f(l1:l2,m,n,i+2)=-ampl*rr
        enddo; enddo
      endif
!
    endsubroutine uniform_phi
!***********************************************************************
    subroutine phi_comp_over_r(ampl,f,i)
!
!  B_phi ~ 1/R field (in terms of vector potential)
!  meaningful mainly in cylindrical coordinates, otherwise it will be By~1/x
!
!  05-jul-07/mgellert: coded
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl
!
      if (ampl==0) then
        f(:,:,:,i:i+2)=0
        if (lroot) print*,'phi_comp_over_r: set variable to zero; i=',i
      else
        if (coord_system=='cylindric') then
          print*,'phi_comp_over_r: set Bphi ~ 1/r ; i=',i
        else
          print*,'phi_comp_over_r: set By ~ 1/x ; i=',i
        endif
        if ((ip<=16).and.lroot) print*,'phi_comp_over_r: ampl=',ampl
        do n=n1,n2; do m=m1,m2
          f(l1:l2,m,n,i  )=0.0!ampl*z(n)/x(l1:l2)
          f(l1:l2,m,n,i+1)=0.0
          f(l1:l2,m,n,i+2)=-ampl*log(x(l1:l2))
        enddo; enddo
      endif
!
    endsubroutine phi_comp_over_r
!***********************************************************************
    subroutine vfield(ampl,f,i,kx)
!
!  Vertical field, for potential field test
!
!  14-jun-02/axel: coded
!  02-aug-2005/joishi: allowed for arbitrary kx
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl
      real,optional :: kx
      real :: k
!
      if (present(kx)) then
         k = kx
      else
         k = 2*pi/Lx
      endif
!
      if (ampl==0) then
        f(:,:,:,i:i+2)=0
        if (lroot) print*,'vfield: set variable to zero; i=',i
      else
        if (lroot) print*,'vfield: implement x-dependent vertical field'
        if ((ip<=8).and.lroot) print*,'vfield: x-dependent vertical field'
        do n=n1,n2; do m=m1,m2
          f(l1:l2,m,n,i  )=0.0
          f(l1:l2,m,n,i+1)=ampl*sin(k*x(l1:l2))
          f(l1:l2,m,n,i+2)=0.0
        enddo; enddo
      endif
!
    endsubroutine vfield
!***********************************************************************
    subroutine vfield2(ampl,f,i)
!
!  Vertical field, zero on boundaries
!
!  22-jun-04/anders: adapted from vfield
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,kx
!
      if (ampl==0) then
        f(:,:,:,i:i+2)=0
        if (lroot) print*,'vfield2: set variable to zero; i=',i
      else
        kx=2*pi/Lx
        if (lroot) print*,'vfield2: implement x-dependent vertical field'
        do n=n1,n2; do m=m1,m2
          f(l1:l2,m,n,i  )=0.0
          f(l1:l2,m,n,i+1)=ampl*cos(kx*x(l1:l2))
          f(l1:l2,m,n,i+2)=0.0
        enddo; enddo
      endif
!
    endsubroutine vfield2
!***********************************************************************
    subroutine posnoise_vect(ampl,f,i1,i2)
!
!  Add Gaussian noise (= normally distributed) white noise for variables i1:i2
!
!  28-may-04/axel: adapted from gaunoise
!
      integer :: i,i1,i2
      real, dimension (mx) :: tmp
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl
!
!  set gaussian random noise vector
!
      if (ampl==0) then
        if (lroot) print*,'posnoise_vect: ampl=0 for i1,i2=',i1,i2
      else
        if ((ip<=8).and.lroot) print*,'posnoise_vect: i1,i2=',i1,i2
        do n=1,mz; do m=1,my
          do i=i1,i2
            call random_number_wrapper(tmp)
            f(:,m,n,i)=f(:,m,n,i)+ampl*tmp
          enddo
        enddo; enddo
      endif
!
    endsubroutine posnoise_vect
!***********************************************************************
    subroutine posnoise_scal(ampl,f,i)
!
!  Add Gaussian (= normally distributed) white noise for variable i
!
!  28-may-04/axel: adapted from gaunoise
!
      integer :: i
      real, dimension (mx) :: tmp
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl
!
!  set positive random noise vector
!
      if (ampl==0) then
        if (lroot) print*,'posnoise_scal: ampl=0 for i=',i
      else
        if ((ip<=8).and.lroot) print*,'posnoise_scal: i=',i
        if (lroot) print*,'posnoise_scal: variable i=',i
        do n=1,mz; do m=1,my
          call random_number_wrapper(tmp)
          f(:,m,n,i)=f(:,m,n,i)+ampl*tmp
        enddo; enddo
      endif
!
!  Wouldn't the following be equivalent (but clearer)?
!
!  call posnoise_vect(ampl,f,i,i)
!
!
    endsubroutine posnoise_scal
!***********************************************************************
    subroutine posnoise_rel_vect(ampl,ampl_rel,f,i1,i2)
!
!  Add noise (= box distributed) white noise for variables i1:i2
!
!  28-may-04/axel: adapted from gaunoise
!
      integer :: i,i1,i2
      real, dimension (mx) :: tmp
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,ampl_rel
!
!  set random noise vector
!
      if (ampl==0) then
        if (lroot) print*,'posnoise_vect: ampl=0 for i1,i2=',i1,i2
      else
        if ((ip<=8).and.lroot) print*,'posnoise_vect: i1,i2=',i1,i2
          if (lroot) print*,'posnoise_vect: variable i=',i
        do n=1,mz; do m=1,my
          do i=i1,i2
            call random_number_wrapper(tmp)
            f(:,m,n,i)=f(:,m,n,i)+ampl*(1.+ampl_rel*(tmp-0.5))
          enddo
        enddo; enddo
      endif
!
    endsubroutine posnoise_rel_vect
!***********************************************************************
    subroutine posnoise_rel_scal(ampl,ampl_rel,f,i)
!
!  Add noise (= box distributed) white noise for variables i1:i2
!
!  28-may-04/axel: adapted from gaunoise
!
      integer :: i
      real, dimension (mx) :: tmp
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,ampl_rel
!
!  set random noise vector
!
      if (ampl==0) then
        if (lroot) print*,'posnoise_vect: ampl=0 for i=',i
      else
        if ((ip<=8).and.lroot) print*,'posnoise_vect: i=',i
          if (lroot) print*,'posnoise_vect: variable i=',i
        do n=1,mz; do m=1,my
          call random_number_wrapper(tmp)
          f(:,m,n,i)=f(:,m,n,i)+ampl*(1.+ampl_rel*(tmp-0.5))
        enddo; enddo
      endif
!
    endsubroutine posnoise_rel_scal
!***********************************************************************
    subroutine gaunoise_vect(ampl,f,i1,i2)
!
!  Add Gaussian noise (= normally distributed) white noise for variables i1:i2
!
!  23-may-02/axel: coded
!  10-sep-03/axel: result only *added* to whatever f array had before
!
      real :: ampl
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i1,i2
!
      real, dimension (mx) :: r,p,tmp
      integer :: i
!
      intent(in)    :: ampl,i1,i2
      intent(inout) :: f
!
!  set gaussian random noise vector
!
      if (ampl==0) then
        if (lroot) print*,'gaunoise_vect: ampl=0 for i1,i2=',i1,i2
      else
        if ((ip<=8).and.lroot) print*,'gaunoise_vect: i1,i2=',i1,i2
        do n=1,mz; do m=1,my
          do i=i1,i2
            if (lroot.and.m==1.and.n==1) print*,'gaunoise_vect: variable i=',i
            if (modulo(i-i1,2)==0) then
              call random_number_wrapper(r)
              call random_number_wrapper(p)
              tmp=sqrt(-2*log(r))*sin(2*pi*p)
            else
              tmp=sqrt(-2*log(r))*cos(2*pi*p)
            endif
            f(:,m,n,i)=f(:,m,n,i)+ampl*tmp
          enddo
        enddo; enddo
      endif
!
    endsubroutine gaunoise_vect
!***********************************************************************
    subroutine gaunoise_scal(ampl,f,i)
!
!  Add Gaussian (= normally distributed) white noise for variable i
!
!  23-may-02/axel: coded
!  10-sep-03/axel: result only *added* to whatever f array had before
!
      real :: ampl
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i
!
      real, dimension (mx) :: r,p,tmp
!
      intent(in)    :: ampl,i
      intent(inout) :: f
!
!  set gaussian random noise vector
!
      if (ampl==0) then
        if (lroot) print*,'gaunoise_scal: ampl=0 for i=',i
      else
        if ((ip<=8).and.lroot) print*,'gaunoise_scal: i=',i
        if (lroot) print*,'gaunoise_scal: variable i=',i
        do n=1,mz; do m=1,my
          call random_number_wrapper(r)
          call random_number_wrapper(p)
          tmp=sqrt(-2*log(r))*sin(2*pi*p)
          f(:,m,n,i)=f(:,m,n,i)+ampl*tmp
        enddo; enddo
      endif
!
!  Wouldn't the following be equivalent (but clearer)?
!
!  call gaunoise_vect(ampl,f,i,i)
!
!
    endsubroutine gaunoise_scal
!***********************************************************************
    subroutine gaunoise_prof_vect(ampl,f,i1,i2)
!
!  Add Gaussian (= normally distributed) white noise for variables i1:i2
!  with amplitude profile AMPL.
!
! 18-apr-04/wolf: adapted from gaunoise_vect
!
      real, dimension (nz) :: ampl
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i1,i2
!
      real, dimension (mx) :: r,p,tmp
      integer :: i
!
      intent(in)    :: ampl,i1,i2
      intent(inout) :: f
!
      if ((ip<=8).and.lroot) print*,'GAUNOISE_PROF_VECT: i1,i2=',i1,i2
      do n=1,mz; do m=1,my
        do i=i1,i2
          if (lroot.and.m==1.and.n==1) print*,'gaunoise_vect: variable i=',i
          if (modulo(i-i1,2)==0) then
            call random_number_wrapper(r)
            call random_number_wrapper(p)
            tmp=sqrt(-2*log(r))*sin(2*pi*p)
          else
            tmp=sqrt(-2*log(r))*cos(2*pi*p)
          endif
          f(:,m,n,i)=f(:,m,n,i)+ampl(n)*tmp
        enddo
      enddo; enddo
!
    endsubroutine gaunoise_prof_vect
!***********************************************************************
    subroutine gaunoise_prof_scal(ampl,f,i)
!
!  Add Gaussian (= normally distributed) white noise for variable i with
!  amplitude profile AMPL.
!
! 18-apr-04/wolf: coded
!
      real, dimension (nz) :: ampl
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i
!
      intent(in)    :: ampl,i
      intent(inout) :: f
!
      if ((ip<=8).and.lroot) print*,'GAUNOISE_PROF_SCAL: i=',i
      call gaunoise_prof_vect(ampl,f,i,i)
!
    endsubroutine gaunoise_prof_scal
!***********************************************************************
    subroutine gaunoise_rprof_vect(ampl,f,i1,i2,rnoise_int,rnoise_ext)
!
!  Add Gaussian noise within rnoise_int < r < rnoise_ext.
!  Use PROF as buffer variable so we don't need to allocate a large
!  temporary.
!
!  18-apr-04/wolf: coded
!
      use Sub, only: cubic_step, get_radial_distance
!
      real :: ampl,rnoise_int,rnoise_ext
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i1,i2
!
      real, dimension (mx) :: prof, rr, r, p, tmp, rr_cyl, rr_sph
      real :: dr
      integer :: i
!
      intent(in)  :: ampl,i1,i2
      intent(out) :: f
!
!  The default is noise in the range r_int < r < r_ext, but the user
!  is allowed to use a different range by initializing the variables
!  rnoise_int and rnoise_ext
!
      if (rnoise_int == impossible) rnoise_int=r_int
      if (rnoise_ext == impossible) rnoise_ext=r_ext
!
!  set up profile
!
      dr = rnoise_ext-max(0.,rnoise_int)
!
      do n=1,mz; do m=1,my
        call get_radial_distance(rr_sph,rr_cyl)
        if (lcylindrical_coords.or.lcylinder_in_a_box) rr=rr_cyl
        if (lspherical_coords  .or.lsphere_in_a_box)   rr=rr_sph
        prof = 1 - cubic_step(rr,rnoise_ext,0.25*dr,SHIFT=-1.)
!
        if (r_int>0.) then
          prof = prof*cubic_step(rr,rnoise_int,0.25*dr,SHIFT=1.)
        endif
        prof = ampl*prof
!
        do i=i1,i2
          if (lroot.and.m==1.and.n==1) print*,'gaunoise_vect: variable i=',i
          if (modulo(i-i1,2)==0) then
            call random_number_wrapper(r)
            call random_number_wrapper(p)
            tmp=sqrt(-2*log(r))*sin(2*pi*p)
          else
            tmp=sqrt(-2*log(r))*cos(2*pi*p)
          endif
          f(:,m,n,i)=f(:,m,n,i)+prof*tmp
        enddo
!
      enddo; enddo
!
    endsubroutine gaunoise_rprof_vect
!***********************************************************************
    subroutine gaunoise_rprof_scal(ampl,f,i,rnoise_int,rnoise_ext)
!
!  Add Gaussian noise within r_int < r < r_ext.
!  Use PROF as buffer variable so we don't need to allocate a large
!  temporary.
!
!  18-apr-04/wolf: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,rnoise_int,rnoise_ext
      integer :: i
!
      intent(in) :: ampl,i
      intent(out) :: f
!
      call gaunoise_rprof_vect(ampl,f,i,i,rnoise_int,rnoise_ext)
!
    endsubroutine gaunoise_rprof_scal
!***********************************************************************
    subroutine trilinear(f,ivar,amplx,amply,amplz)
!
!  Produce a profile that is linear in any non-periodic direction, but
!  periodic in periodic ones (for testing purposes).
!
!  5-nov-02/wolf: coded
! 23-nov-02/axel: included scaling factor ampl, corrected lperi argument
!
      real :: amplx,amply,amplz
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: ivar
!
      real, dimension (nx) :: tmp
!
      if (lroot) print*, 'trilinear: ivar = ', ivar
!
!  x direction
!
      do n=n1,n2; do m=m1,m2
        if (lperi(1)) then
          tmp = amplx*sin(2*pi/Lx*(x(l1:l2)-xyz0(1)-0.25*Lxyz(1)))
        else
          tmp = amplx*x(l1:l2)
        endif
!
!  y direction
!
        if (lperi(2)) then
          tmp = tmp + amply*sin(2*pi/Ly*(y(m)-xyz0(2)-0.25*Lxyz(2)))
        else
          tmp = tmp + amply*y(m)
        endif
!
!  z direction
!
        if (lperi(3)) then
          tmp = tmp + amplz*sin(2*pi/Lz*(z(n)-xyz0(3)-0.25*Lxyz(3)))
        else
          tmp = tmp + amplz*z(n)
        endif
!
        f(l1:l2,m,n,ivar) = tmp
!
      enddo; enddo
!
    endsubroutine trilinear
!***********************************************************************
    subroutine triquad(ampl,f,ivar,kx,ky,kz,kxx,kyy,kzz)
!
!  Produce a profile that is quadratic in any non-periodic direction, but
!  constant in periodic ones, with optional coefficients.
!
!  20-jul-09/hubbard: coded
!
      real :: ampl
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: ivar
!
      real, dimension (nx) :: tmp
!
      real,optional :: kx,ky,kz,kxx,kyy,kzz
!
      if (lroot) print*, 'triquad: ivar = ', ivar
!
!  linear portion
!
!  x direction
!
      do n=n1,n2; do m=m1,m2
        tmp = 0.
        if (.not. lperi(1)) then
          if (present(kx)) then
            tmp = kx*x(l1:l2)
          endif
        endif
!
!  y direction
!
        if (.not. lperi(2)) then
          if (present(ky)) then
            tmp = tmp + ky*y(m)
          endif
        endif
!
!  z direction
!
        if (.not. lperi(3)) then
          if (present(kz)) then
            tmp = tmp + kz*z(n)
          endif
        endif
!
        f(l1:l2,m,n,ivar) = ampl*tmp
!
      enddo; enddo
!
!  quadratic portion
!  negative ampl is because parabola is based on -x**2
!
  if (present(kxx)) call parabola(-ampl*abs(kxx)/kzz,f,icc,kx=kxx)
  if (present(kyy)) call parabola(-ampl*abs(kyy)/kzz,f,icc,ky=kyy)
  if (present(kzz)) call parabola(-ampl*abs(kzz)/kzz,f,icc,kz=kzz)
!
    endsubroutine triquad
!***********************************************************************
    subroutine isotdisk(powerlr,f,ivar,zoverh,hoverr)
!
!  eek
!
! 21-jul-09/hubbard: coded
!
      real :: powerlr
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: ivar
!
      real :: zoverh, hoverr
!
      if (lroot) print*, 'isotdisk: ivar = ', ivar
!
!
      do n=1,mz; do m=1,my
        f(1:mx,m,n,ivar) = f(1:mx,m,n,ivar) + ((1.5*zoverh**2-powerlr)*x(1:mx) &
            -zoverh/(hoverr)*z(n) &
            +(powerlr/2-3*zoverh**2)*x(1:mx)**2 &
            -0.5/(hoverr**2)*z(n)**2 &
            +3.*zoverh/hoverr*z(n)*x(1:mx))
!
      enddo; enddo
!
    endsubroutine isotdisk
!***********************************************************************
    subroutine cos_cos_sin(ampl,f,ivar)
!
!  Produce a profile that is linear in any non-periodic direction, but
!  periodic in periodic ones (for testing purposes).
!
!  7-dec-02/axel: coded
!
      integer :: ivar
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,kx,ky,kz
!
      if (lroot) print*, 'cos_cos_sin: ivar = ', ivar
!
      kx=2*pi/Lx*3
      ky=2*pi/Ly*3
      kz=pi/Lz
      do n=n1,n2; do m=m1,m2
        f(l1:l2,m,n,ivar) = ampl*cos(kx*x(l1:l2))*cos(ky*y(m))*sin(kz*z(n))
      enddo; enddo
!
    endsubroutine cos_cos_sin
!***********************************************************************
    subroutine tor_pert(ampl,f,ivar)
!
!  Produce a profile that is periodic in the y- and z-directions.
!  For testing the Balbus-Hawley instability of a toroidal magnetic field
!
!  12-feb-03/ulf: coded
!
      integer :: ivar
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,ky,kz
!
      if (lroot) print*, 'tor_pert: sinusoidal modulation of ivar = ', ivar
!
      ky=2*pi/Ly
      kz=2.*pi/Lz
      do n=n1,n2; do m=m1,m2
        f(l1:l2,m,n,ivar) = ampl*cos(ky*y(m))*cos(kz*z(n))
      enddo; enddo
!
    endsubroutine tor_pert
!***********************************************************************
    subroutine const_omega(ampl,f,ivar)
!
!  Set up profile for differential rotation
!
!  16-jul-03/axel: coded
!
      real :: ampl
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: ivar
!
      if (lroot) print*, 'const_omega: constant angular velcoity  = ', ivar
!
      do n=n1,n2; do m=m1,m2
        f(l1:l2,m,n,ivar) = ampl*x(l1:l2)*sinth(m)
      enddo; enddo
!
    endsubroutine const_omega
!***********************************************************************
    subroutine diffrot(ampl,f,ivar)
!
!  Set up profile for differential rotation
!
!  16-jul-03/axel: coded
!
      real :: ampl
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: ivar
!
      if (lroot) print*, 'diffrot: sinusoidal modulation of ivar = ', ivar
!
      do n=n1,n2; do m=m1,m2
        f(l1:l2,m,n,ivar) = ampl*cos(x(l1:l2))*cos(z(n))
      enddo; enddo
!
    endsubroutine diffrot
!***********************************************************************
    subroutine olddiffrot(ampl,f,ivar)
!
!  Set up profile for differential rotation
!
!  16-jul-03/axel: coded
!
      real :: ampl
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: ivar
!
      real :: kx,kz
!
      if (lroot) print*, 'olddiffrot: sinusoidal modulation of ivar = ', ivar
!
      kx=.5*pi/Lx
      kz=.5*pi/Lz
      do n=n1,n2; do m=m1,m2
        f(l1:l2,m,n,ivar) = ampl*sin(kx*x(l1:l2))*cos(kz*z(n))
      enddo; enddo
!
    endsubroutine olddiffrot
!***********************************************************************
    subroutine powern(ampl,initpower,cutoff,f,i1,i2)
!
!   Produces k^initpower*exp(-k**2/cutoff**2)  spectrum.
!   Still just one processor (but can be remeshed afterwards).
!
!   07-may-03/tarek: coded
!
      use Fourier, only: fourier_transform
!
      real :: ampl,initpower,cutoff
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i1,i2
!
      real, dimension (:,:,:), allocatable :: k2, u_re, u_im
      real, dimension (nx) :: k2x
      real, dimension (ny) :: k2y
      real, dimension (nz) :: k2z
!
      integer :: i, stat
!
!  Allocate memory for arrays.
!
      allocate(k2(nx,ny,nz),stat=stat)
      if (stat>0) call fatal_error('powern','Could not allocate memory for k2')
      allocate(u_re(nx,ny,nz),stat=stat)
      if (stat>0) call fatal_error('powern','Could not allocate memory for u_re')
      allocate(u_im(nx,ny,nz),stat=stat)
      if (stat>0) call fatal_error('powern','Could not allocate memory for u_im')
!
      if (ampl==0) then
        f(:,:,:,i1:i2)=0
        if (lroot) print*,'powern: set variable to zero; i1,i2=',i1,i2
      else
        call gaunoise_vect(ampl,f,i1,i2) ! which has a k^2. spectrum
!
        if ((initpower/=2.).or.(cutoff/=0.)) then
!
          k2x = cshift((/(i-(nx+1)/2,i=0,nx-1)/),+(nx+1)/2)*2*pi/Lx
          k2 =      (spread(spread(k2x,2,ny),3,nz))**2
!
          k2y = cshift((/(i-(ny+1)/2,i=0,ny-1)/),+(ny+1)/2)*2*pi/Ly
          k2 = k2 + (spread(spread(k2y,1,nx),3,nz))**2
!
          k2z = cshift((/(i-(nz+1)/2,i=0,nz-1)/),+(nz+1)/2)*2*pi/Lz
          k2 = k2 + (spread(spread(k2z,1,nx),2,ny))**2
!
          k2(1,1,1) = 1.  ! Avoid division by zero
!
          do i=i1,i2
            u_re=f(l1:l2,m1:m2,n1:n2,i)
            u_im=0.
            !  fft of gausian noise w/ k^2 spectrum
            call fourier_transform(u_re,u_im)
            ! change to k^n spectrum
            u_re =(k2)**(.25*initpower-.5)*u_re
            u_im =(k2)**(.25*initpower-.5)*u_im
            ! cutoff (changed to hyperviscous cutoff filter)
            if (cutoff /= 0.) then
              u_re = u_re*exp(-(k2/cutoff**2.)**2)
              u_im = u_im*exp(-(k2/cutoff**2.)**2)
            endif
            ! back to real space
            call fourier_transform(u_re,u_im,linv=.true.)
            f(l1:l2,m1:m2,n1:n2,i)=u_re
!
            if (lroot) then
              if (cutoff==0) then
                print*,'powern: k^',initpower,' spectrum : var  i=',i
              else
                print*,'powern: with cutoff : k^n*exp(-k^4/k0^4) w/ n=', &
                    initpower,', k0 =',cutoff,' : var  i=',i
              endif
            endif
          enddo !i
        endif !(initpower/=2.).or.(cutoff/=0.)
!
      endif !(ampl==0)
!
!  Deallocate arrays.
!
      if (allocated(k2)) deallocate(k2)
      if (allocated(u_re)) deallocate(u_im)
!
    endsubroutine powern
!***********************************************************************
    subroutine power_randomphase(ampl,initpower,cutoff,f,i1,i2,lscale_tobox)
!
!   Produces k^initpower*exp(-k**2/cutoff**2)  spectrum.
!   Still just one processor (but can be remeshed afterwards).
!
!   07-may-03/tarek: coded
!   08-may-08/nils: adapted to work on multiple processors
!   06-jul-08/nils+andre: Fixed problem when running on
!      mult. procs (thanks to Andre Kapelrud for finding the bug)
!
      use Fourier, only: fourier_transform
!
      logical, intent(in), optional :: lscale_tobox
      logical :: lscale_tobox1
      integer :: i,i1,i2,ikx,iky,ikz,stat
      real, dimension (:,:,:), allocatable :: k2, u_re, u_im, r
      real, dimension (:), allocatable :: kx, ky, kz
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,initpower,mhalf,cutoff,scale_factor
!
      if (present(lscale_tobox)) then
        lscale_tobox1 = lscale_tobox
      else
        lscale_tobox1 = .false.
      endif
!
!  Allocate memory for arrays.
!
      allocate(k2(nx,ny,nz),stat=stat)
      if (stat>0) call fatal_error('powern','Could not allocate memory for k2')
      allocate(u_re(nx,ny,nz),stat=stat)
      if (stat>0) call fatal_error('powern','Could not allocate memory for u_re')
      allocate(u_im(nx,ny,nz),stat=stat)
      if (stat>0) call fatal_error('powern','Could not allocate memory for u_im')
      allocate(r(nx,ny,nz),stat=stat)
      if (stat>0) call fatal_error('powern','Could not allocate memory for r')
      allocate(kx(nxgrid),stat=stat)
      if (stat>0) call fatal_error('powern', &
          'Could not allocate memory for kx')
      allocate(ky(nygrid),stat=stat)
      if (stat>0) call fatal_error('powern', &
          'Could not allocate memory for ky')
      allocate(kz(nzgrid),staT=stat)
      if (stat>0) call fatal_error('powern', &
          'Could not allocate memory for kz')
!
      if (ampl==0) then
        f(:,:,:,i1:i2)=0
        if (lroot) print*,'power_randomphase: set variable to zero; i1,i2=',i1,i2
      else
!
!  calculate k^2
!
        scale_factor=1
        if (lscale_tobox1) scale_factor=2*pi/Lx
        kx=cshift((/(i-(nxgrid+1)/2,i=0,nxgrid-1)/),+(nxgrid+1)/2)*scale_factor
!
        scale_factor=1
        if (lscale_tobox1) scale_factor=2*pi/Ly
        ky=cshift((/(i-(nygrid+1)/2,i=0,nygrid-1)/),+(nygrid+1)/2)*scale_factor
!
        scale_factor=1
        if (lscale_tobox1) scale_factor=2*pi/Lz
        kz=cshift((/(i-(nzgrid+1)/2,i=0,nzgrid-1)/),+(nzgrid+1)/2)*scale_factor
!
!  integration over shells
!
        if (lroot .AND. ip<10) &
             print*,'power_randomphase:fft done; now integrate over shells...'
        do ikz=1,nz
          do iky=1,ny
            do ikx=1,nx
              k2(ikx,iky,ikz)=kx(ikx)**2+ky(iky+ipy*ny)**2+kz(ikz+ipz*nz)**2
            enddo
          enddo
        enddo
        if (lroot) k2(1,1,1) = 1.  ! Avoid division by zero
!
!  To get shell integrated power spectrum E ~ k^n, we need u ~ k^m
!  and since E(k) ~ u^2 k^2 we have n=2m+2, so m=n/2-1.
!  Further, since we operate on k^2, we need m/2 (called mhalf below)
!
        mhalf=.5*(.5*initpower-1)
!
!  generate all 3 velocity components separately
!
        do i=i1,i2
          ! generate k^n spectrum with random phase (between -pi and pi)
          call random_number_wrapper(r); u_re=ampl*k2**mhalf*cos(pi*(2*r-1))
          call random_number_wrapper(r); u_im=ampl*k2**mhalf*sin(pi*(2*r-1))
          ! cutoff (changed to hyperviscous cutoff filter)
          if (cutoff /= 0.) then
            u_re = u_re*exp(-(k2/cutoff**2.)**2)
            u_im = u_im*exp(-(k2/cutoff**2.)**2)
          endif
          ! back to real space
          call fourier_transform(u_re,u_im,linv=.true.)
          f(l1:l2,m1:m2,n1:n2,i)=u_re
          if (lroot) then
            if (cutoff==0) then
              print*,'powern: k^',initpower,' spectrum : var  i=',i
            else
              print*,'powern: with cutoff : k^n*exp(-k^4/k0^4) w/ n=', &
                  initpower,', k0 =',cutoff,' : var  i=',i
            endif
          endif
        enddo !i
!
      endif !(ampl==0)
!
!  Deallocate arrays.
!
      if (allocated(k2))   deallocate(k2)
      if (allocated(u_re)) deallocate(u_re)
      if (allocated(u_im)) deallocate(u_im)
      if (allocated(r))  deallocate(r)
      if (allocated(kx)) deallocate(kx)
      if (allocated(ky)) deallocate(ky)
      if (allocated(kz)) deallocate(kz)
!
    endsubroutine power_randomphase
!***********************************************************************
    subroutine random_isotropic_KS(initpower,f,i1,N_modes)
!
!   produces random, isotropic field from energy spectrum following the
!   KS method (Malik and Vassilicos, 1999.)
!
!   more to do; unsatisfactory so far - at least for a steep power-law energy spectrum
!
!   24-sept-04/snod: coded first attempt
!
      use Sub, only: cross, dot, dot2
!
      integer :: modeN,N_modes,l,n,m,i1
      real, dimension (mx,my,mz,mfarray) :: f
!
! how many wavenumbers?
      real, dimension (3,1024) :: kk,RA,RB !or through whole field for each wavenumber?
      real, dimension (3) :: k_unit,vec,ee,e1,e2,field
      real :: initpower,kmin,ps,k,phi,theta,alpha,beta,dk
      real :: ex,ey,ez,norm,kdotx,r
!
!    minlen=Lxyz(1)/(nx-1)
!    kmax=2.*pi/minlen
!    N_modes=int(0.5*(nx-1))
!    hh=Lxyz(1)/(nx-1)
!    pta=(nx)**(1.0/(nx-1))
!    do modeN=1,N_modes
!       ggt=(kkmax-kkmin)/(N_modes-1)
!       ggt=(kkmax/kkmin)**(1./(N_modes-1))
!        k(modeN)=kmin+(ggt*(modeN-1))
!        k(modeN)=(modeN+3)*2*pi/Lxyz(1)
!       k(modeN)=kkmin*(ggt**(modeN-1)
!    enddo
!
!    do modeN=1,N_modes
!       if (modeN==1)delk(modeN)=(k(modeN+1)-K(modeN))
!       if (modeN==N_modes)delk(modeN)=(k(modeN)-k(modeN-1))
!       if (modeN>1.and.modeN<N_modes)delk(modeN)=(k(modeN+1)-k(modeN-2))/2.0
!    enddo
!          mk=(k2*k2)*((1.0 + (k2/(bk_min*bk_min)))**(0.5*initpower-2.0))
!
!  set kmin
!
      kmin=2.*pi/Lxyz(1)
!
      do modeN=1,N_modes
!
!  pick wavenumber
!
        k=modeN*kmin
!
!  calculate dk
!
        dk=1.0*kmin
!
!   pick 4 random angles for each mode
!
!
        call random_number_wrapper(r); theta=pi*(2*r - 0.)
        call random_number_wrapper(r); phi=pi*(2*r - 0.)
        call random_number_wrapper(r); alpha=pi*(2*r - 0.)
        call random_number_wrapper(r); beta=pi*(2*r - 0.)
!       call random_number_wrapper(r); gamma(modeN)=pi*(2*r - 0.)  ! random phase?
!
!
!   make a random unit vector by rotating fixed vector to random position
!   (alternatively make a random transformation matrix for each k)
!
        k_unit(1)=sin(theta)*cos(phi)
        k_unit(2)=sin(theta)*sin(phi)
        k_unit(3)=cos(theta)
!
!   make a vector kk of length k from the unit vector for each mode
!
        kk(:,modeN)=k*k_unit(:)
!
!   construct basis for plane having rr normal to it
!   (bit of code from forcing to construct x', y')
!
        if ((k_unit(2)==0).and.(k_unit(3)==0)) then
          ex=0.; ey=1.; ez=0.
        else
          ex=1.; ey=0.; ez=0.
        endif
        ee = (/ex, ey, ez/)
        call cross(k_unit(:),ee,e1)
        call dot2(e1,norm); e1=e1/sqrt(norm) ! e1: unit vector perp. to kk
        call cross(k_unit(:),e1,e2)
        call dot2(e2,norm); e2=e2/sqrt(norm) ! e2: unit vector perp. to kk, e1
!
!   make two random unit vectors RB and RA in the constructed plane
!
        RA(:,modeN) = cos(alpha)*e1 + sin(alpha)*e2
        RB(:,modeN) = cos(beta)*e1  + sin(beta)*e2
!
!   define the power spectrum (ps=sqrt(2.*power_spectrum(k)*delta_k/3.))
!
        ps=(k**(initpower/2.))*sqrt(dk*2./3.)
!
!   give RA and RB length ps
!
        RA(:,modeN)=ps*RA(:,modeN)
        RB(:,modeN)=ps*RB(:,modeN)
!
!   form RA = RA x k_unit and RB = RB x k_unit
!
        call cross(RA(:,modeN),k_unit(:),RA(:,modeN))
        call cross(RB(:,modeN),k_unit(:),RB(:,modeN))
!
      enddo
!
!   make the field
!
      do n=n1,n2
        do m=m1,m2
          do l=l1,l2
            field=0.  ! initialize field
            vec(1)=x(l)    ! local coordinates?
            vec(2)=y(m)
            vec(3)=z(n)
            do modeN=1,N_modes  ! sum over N_modes modes
              call dot(kk(:,modeN),vec,kdotx)
              field = field + cos(kdotx)*RA(:,modeN) + sin(kdotx)*RB(:,modeN)
            enddo
            f(l,m,n,i1)   = f(l,m,n,i1)   + field(1)
            f(l,m,n,i1+1) = f(l,m,n,i1+1) + field(2)
            f(l,m,n,i1+2) = f(l,m,n,i1+2) + field(3)
          enddo
        enddo
      enddo
!
    endsubroutine random_isotropic_KS
!***********************************************************************
    subroutine corona_init(f)
!
!  Initialize the density for a given temperature profile
!  in the vertical (z) direction by solving for hydrostatic
!  equilibrium.
!  The temperature is hard coded as three polynomials.
!
!  07-dec-05/bing : coded.
!
      use EquationOfState, only: lnrho0,gamma,gamma_m1,cs20,cs2top,cs2bot
      use Mpicomm, only: mpibcast_real
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: tmp,ztop,zbot
      integer, parameter :: prof_nz=150
      real, dimension (prof_nz) :: prof_lnT,prof_lnrho,prof_z
      integer :: i,lend,j
!
      ! file location settings
      character (len=*), parameter :: lnrho_dat = 'driver/b_lnrho.dat'
      character (len=*), parameter :: lnT_dat = 'driver/b_lnT.dat'
!
      ! temperature given as function lnT(z) in SI units
      ! [T] = K   &   [z] = Mm   & [rho] = kg/m^3
      if (pretend_lnTT) print*,'corona_init: not implemented for pretend_lnTT=T'
!
      if (lroot) then
        inquire(IOLENGTH=lend) tmp
        open (10,file=lnT_dat,form='unformatted',status='unknown',recl=lend*prof_nz)
        read (10) prof_lnT
        read (10) prof_z
        close (10)
!
        open (10,file=lnrho_dat,form='unformatted',status='unknown',recl=lend*prof_nz)
        read (10) prof_lnrho
        close (10)
      endif
!
      call mpibcast_real(prof_lnT,prof_nz)
      call mpibcast_real(prof_z,prof_nz)
      call mpibcast_real(prof_lnrho,prof_nz)
!
      prof_z = prof_z*1.e6/unit_length
      prof_lnT = prof_lnT - alog(real(unit_temperature))
      prof_lnrho = prof_lnrho - alog(real(unit_density))
!
! simple linear interpolation
!
      do j=n1,n2
         do i=1,prof_nz-1
            if (z(j) >= prof_z(i) .and. z(j) < prof_z(i+1) ) then
               f(:,:,j,ilnrho) = (prof_lnrho(i)*(prof_z(i+1) - z(j)) +   &
                    prof_lnrho(i+1)*(z(j)-prof_z(i)) ) / (prof_z(i+1)-prof_z(i))
!
               tmp =  (prof_lnT(i)*(prof_z(i+1) - z(j)) +   &
                    prof_lnT(i+1)*(z(j)-prof_z(i)) ) / (prof_z(i+1)-prof_z(i))
!
               if (ltemperature) then
                  f(:,:,j,ilnTT) = tmp
               elseif (lentropy) then
                  f(:,:,j,iss) = (alog(gamma_m1/cs20)+tmp- &
                       gamma_m1*(f(l1,m1,j,ilnrho)-lnrho0))/gamma
               endif
               exit
            endif
         enddo
         if (z(j) >= prof_z(prof_nz)) then
            f(:,:,j,ilnrho) = prof_lnrho(prof_nz)
!
            tmp =  prof_lnT(prof_nz)
!
            if (ltemperature) then
               f(:,:,j,ilnTT) = tmp
            elseif (lentropy) then
               f(:,:,j,iss) = (alog(gamma_m1/cs20)+tmp- &
                    gamma_m1*(f(l1,m1,j,ilnrho)-lnrho0))/gamma
            endif
         endif
      enddo
!
      ztop=xyz0(3)+Lxyz(3)
      zbot=xyz0(3)
!
      do i=1,prof_nz-1
         if (ztop >= prof_z(i) .and. ztop < prof_z(i+1) ) then
!
            tmp =  (prof_lnT(i)*(prof_z(i+1) - ztop) +   &
                 prof_lnT(i+1)*(ztop-prof_z(i)) ) / (prof_z(i+1)-prof_z(i))
            cs2top = gamma_m1*exp(tmp)
!
         elseif (ztop >= prof_z(prof_nz)) then
            cs2top = gamma_m1*exp(prof_lnT(prof_nz))
         endif
         if (zbot >= prof_z(i) .and. zbot < prof_z(i+1) ) then
!
            tmp =  (prof_lnT(i)*(prof_z(i+1) - zbot) +   &
                 prof_lnT(i+1)*(zbot-prof_z(i)) ) / (prof_z(i+1)-prof_z(i))
            cs2bot = gamma_m1*exp(tmp)
!
         endif
      enddo
!
    endsubroutine corona_init
!***********************************************************************
    subroutine mdi_init(f,periodic)
!
!  Intialize the vector potential
!  by potential field extrapolation
!  of a mdi magnetogram
!
!  13-dec-05/bing : coded.
!
      use Fourier, only: fourier_transform_other
      use Mpicomm, only: mpibcast_real,stop_it_if_any
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
!
      real, dimension (:,:), allocatable :: kx,ky,k2,Bz0_i,Bz0_r,A_r,A_i
      logical, intent (in) :: periodic
      real :: zref,Bzflux
      logical :: exists
      integer :: i,j,idx2,idy2,stat,iostat,lend
      integer :: nxinit,nyinit
!
      ! file location settings
      character (len=*), parameter :: mag_field_dat = 'driver/mag_field.dat'
!
!  Allocate memory for arrays.
!
      if (.not.lequidist(1).or..not.lequidist(2)) call fatal_error('mdi_init', &
          'not yet implemented for non-equidistant grids')
!
      iostat = 0
      if (periodic) then
        nxinit=nxgrid
        nyinit=nygrid
      else
        nxinit=2*nxgrid
        nyinit=2*nygrid
      endif
!
      allocate(kx(nxinit,nyinit),stat=stat);     iostat=max(stat,iostat)
      allocate(ky(nxinit,nyinit),stat=stat);     iostat=max(stat,iostat)
      allocate(k2(nxinit,nyinit),stat=stat);     iostat=max(stat,iostat)
      allocate(Bz0_i(nxinit,nyinit),stat=stat);  iostat=max(stat,iostat)
      allocate(Bz0_r(nxinit,nyinit),stat=stat);  iostat=max(stat,iostat)
      allocate(A_r(nxinit,nyinit),stat=stat);    iostat=max(stat,iostat)
      allocate(A_i(nxinit,nyinit),stat=stat);    iostat=max(stat,iostat)
!
      call stop_it_if_any((iostat>0),'mdi_init: '// &
          'Could not allocate memory for variables, please check')
!
!  Auxiliary quantities:
!
!  idx2 and idy2 are essentially =2, but this makes compilers
!  complain if nyinit=1 (in which case this is highly unlikely to be
!  correct anyway), so we try to do this better:
      idx2 = min(2,nxinit)
      idy2 = min(2,nyinit)
!
      if (periodic) then
        kx = spread(kx_fft,2,nyinit)
        ky = spread(ky_fft,1,nxinit)
      else
        kx=spread( &
            cshift((/(i-(nxinit+1)/2,i=0,nxinit-1)/),+(nxinit+1)/2)*pi/Lx &
            ,2,nyinit)
        ky=spread( &
            cshift((/(i-(nyinit+1)/2,i=0,nyinit-1)/),+(nyinit+1)/2)*pi/Ly &
            ,1,nxinit)
      endif
!
      k2 = kx*kx + ky*ky
!
      if (lroot) then
        inquire(file=mag_field_dat,exist=exists)
        call stop_it_if_any(.not.exists, &
            'mdi_init: Magnetogram file not found: "'//trim(mag_field_dat)//'"')
        inquire(IOLENGTH=lend) unit_magnetic
        open (11,file=mag_field_dat,form='unformatted',status='unknown', &
            recl=lend*nxgrid*nygrid,access='direct')
        read (11,rec=1) Bz0_r(1:nxgrid,1:nygrid)
        close (11)
        if (.not.periodic) then
          do i=1,nxgrid
            do j=1,nygrid
              Bz0_r(nxgrid+i,j)=Bz0_r(nxgrid+1-i,j)
              Bz0_r(i,nygrid+j)=Bz0_r(i,nygrid+1-j)
              Bz0_r(nxgrid+i,nygrid+j)=Bz0_r(nxgrid+1-i,nygrid+1-j)
            enddo
          enddo
        endif
      else
        call stop_it_if_any(.false.,'')
      endif
      call mpibcast_real(Bz0_r,(/nxinit,nyinit/))
!
      if (lroot) then
        if (nxgrid==1.and.nygrid/=1) then
          Bzflux =  sum(abs(Bz0_r * 1e-4))*dy*unit_length
          write (*,'(A,E10.2)') 'Bz flux sum(|B|)*dl [Tm] :',Bzflux
        elseif (nxgrid/=1.and.nygrid==1) then
          Bzflux =  sum(abs(Bz0_r * 1e-4))*dx*unit_length
          write (*,'(A,E10.2)') 'Bz flux sum(|B|)*dl [Tm] :',Bzflux
        elseif (nxgrid/=1.and.nygrid/=1) then
          Bzflux =  sum(abs(Bz0_r * 1e-4))*dx*unit_length*dy*unit_length
          write (*,'(A,E10.2)') 'Bz flux sum(|B|)*dA [Tm^2] :',Bzflux
        endif
      endif
!
      Bz0_i = 0.
      Bz0_r = Bz0_r * 1e-4 / unit_magnetic ! Gauss to Tesla and SI to PENCIL units
!
!  Fourier Transform of Bz0:
!
      call fourier_transform_other(Bz0_r,Bz0_i)
!
      do i=1,mz
!
!  Calculate transformed vector potential for every z layer
!
        zref = z(i) - xyz0(3)
!
        if (nygrid > 1) then
          where (k2 /= 0 )
            A_r = -Bz0_i*ky/k2*exp(-sqrt(k2)*zref )
            A_i =  Bz0_r*ky/k2*exp(-sqrt(k2)*zref )
          elsewhere
            A_r = -Bz0_i*ky/ky(1,idy2)*exp(-sqrt(k2)*zref )
            A_i =  Bz0_r*ky/ky(1,idy2)*exp(-sqrt(k2)*zref )
          endwhere
!
          call fourier_transform_other(A_r,A_i,linv=.true.)
!
          f(l1:l2,m1:m2,i,iax)=A_r(ipx*nx+1:(ipx+1)*nx,ipy*ny+1:(ipy+1)*ny)
        else
          f(l1:l2,m1:m2,i,iax)=0.
        endif
!
        if (nxgrid > 1) then
          where (k2 /= 0 )
            A_r =  Bz0_i*kx/k2*exp(-sqrt(k2)*zref )
            A_i = -Bz0_r*kx/k2*exp(-sqrt(k2)*zref )
          elsewhere
            A_r =  Bz0_i*kx/kx(idx2,1)*exp(-sqrt(k2)*zref )
            A_i = -Bz0_r*kx/kx(idx2,1)*exp(-sqrt(k2)*zref )
          endwhere
!
          call fourier_transform_other(A_r,A_i,linv=.true.)
!
          f(l1:l2,m1:m2,i,iay)=A_r(ipx*nx+1:(ipx+1)*nx,ipy*ny+1:(ipy+1)*ny)
        else
          f(l1:l2,m1:m2,i,iay)=0.
        endif
!
        f(l1:l2,m1:m2,i,iaz)=0.
      enddo
!
!  Deallocate arrays.
!
      if (allocated(kx)) deallocate(kx)
      if (allocated(ky)) deallocate(ky)
      if (allocated(k2)) deallocate(k2)
      if (allocated(Bz0_i)) deallocate(Bz0_i)
      if (allocated(Bz0_r)) deallocate(Bz0_r)
      if (allocated(A_r)) deallocate(A_r)
      if (allocated(A_i)) deallocate(A_i)
!
    endsubroutine mdi_init
!***********************************************************************
    subroutine mag_init(f)
!
!  Intialize the vector potential with a potential field extrapolation.
!
!  29-Mar-2011/Bourdin.KIS : coded, adapted parts from 'mdi_init'.
!
      use Fourier, only: setup_extrapol_fact, field_extrapol_z_parallel
      use Mpicomm, only: stop_it_if_any, mpisend_real, mpirecv_real, sum_xy
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      real, dimension (:,:), allocatable :: Bz
      real, dimension (:,:,:), allocatable :: exp_fact
      integer, parameter :: bnx=nxgrid, bny=ny/nprocx ! data in pencil shape
      integer, parameter :: enx=nygrid, eny=nx/nprocy ! transposed data in pencil shape
      integer, parameter :: unit=11
      integer, parameter :: tag_xy=131, tag_z=132
      integer :: py, pz, partner
      real :: Bz_flux, Bz_flux_local
      logical :: exists
      integer :: alloc_err, rec_len
      real, parameter :: reduce_factor=0.25
!
      ! file location settings
      character (len=*), parameter :: mag_field_dat = 'driver/mag_field.dat'
!
      if (.not. lperi(1) .or. .not. lperi(2)) call fatal_error ('mag_init', &
          'Currently only implemented for xy-periodic setups!')
      if (.not. lequidist(1) .or. .not. lequidist(2)) call fatal_error ('mag_init', &
          'not yet implemented for non-equidistant grids')
      if (mod (nygrid, nprocxy) /= 0) call fatal_error ('mag_init', &
          'nygrid needs to be an integer multiple of nprocx*nprocy')
      if (mod (nxgrid, nprocxy) /= 0) call fatal_error ('mag_init', &
          'nxgrid needs to be an integer multiple of nprocx*nprocy')
!
!  Allocate memory for arrays.
!
      allocate (Bz(bnx,bny), stat=alloc_err)
      if (alloc_err > 0) call fatal_error('mag_init', &
          'Could not allocate memory for Bz', .true.)
      allocate (exp_fact(enx,eny,mz), stat=alloc_err)
      if (alloc_err > 0) call fatal_error('mag_init', &
          'Could not allocate memory for exp_fact', .true.)
!
      if (lroot) then
        inquire (file=mag_field_dat, exist=exists)
        call stop_it_if_any(.not. exists, &
            'mag_init: Magnetogram file not found: "'//trim(mag_field_dat)//'"')
        inquire (iolength=rec_len) unit_magnetic
        open (unit, file=mag_field_dat, form='unformatted', recl=rec_len*bnx*bny, access='direct')
        do py = 1, nprocxy-1
          partner = py + ipz*nprocxy
          ! read Bz data for remote processors
          read (unit, rec=1+py) Bz
          ! send Bz data to remote
          call mpisend_real (Bz, (/ bnx, bny /), partner, tag_xy)
        enddo
        ! read local Bz data
        read (unit, rec=1) Bz
        close (unit)
      else
        call stop_it_if_any(.false.,'')
        if (lfirst_proc_z) then
          ! receive Bz data
          partner = ipy + ipz*nprocxy
          call mpirecv_real (Bz, (/ bnx, bny /), 0, tag_xy)
        endif
      endif
!
      if (nprocz > 1) then
        ! distribute Bz along the z-direction
        if (lfirst_proc_z) then
          do pz = 1, nprocz-1
            partner = ipx + ipy*nprocx + pz*nprocxy
            call mpisend_real (Bz, (/ bnx, bny /), partner, tag_z)
          enddo
        else
          partner = ipx + ipy*nprocx
          call mpirecv_real (Bz, (/ bnx, bny /), partner, tag_z)
        endif
      endif
!
      ! Gauss to Tesla and SI to PENCIL units
      Bz = Bz * 1e-4 / unit_magnetic
!
      if (lfirst_proc_z) then
        if ((nxgrid==1).and.(nygrid/=1)) then
          Bz_flux_local = sum(abs(Bz)) * dy * unit_magnetic*unit_length
          call sum_xy (Bz_flux_local, Bz_flux)
          if (lroot) write (*,*) 'Total vertical flux: sum(|Bz|)*dy [T*m] =', Bz_flux
        elseif ((nxgrid/=1).and.(nygrid==1)) then
          Bz_flux_local = sum(abs(Bz)) * dx * unit_magnetic*unit_length
          call sum_xy (Bz_flux_local, Bz_flux)
          if (lroot) write (*,*) 'Total vertical flux: sum(|Bz|)*dx [T*m] =', Bz_flux
        elseif ((nxgrid/=1).and.(nygrid/=1)) then
          Bz_flux_local = sum(abs(Bz)) * dx*dy * unit_magnetic*unit_length**2
          call sum_xy (Bz_flux_local, Bz_flux)
          if (lroot) write (*,*) 'Total vertical flux: sum(|Bz|)*(dx*dy) [T*m^2] =', Bz_flux
        endif
      endif
!
!  Perform the field extrapolation in parallel.
!
      call setup_extrapol_fact (z, xyz0(3), exp_fact, reduce_factor)
      call field_extrapol_z_parallel (Bz, f(l1:l2,m1:m2,:,iax:iay), exp_fact)
!
      deallocate (Bz, exp_fact)
!
    endsubroutine mag_init
!***********************************************************************
    subroutine temp_hydrostatic(f,rho0)
!
! 07-dec-05/bing : coded.
!      intialize the density for a given temperature profile in vertical
!      z direction by solving the hydrostatic equilibrium:
!      dlnrho = - dlnTT + (cp-cv)/T g dz
!
      use EquationOfState, only: lnrho0,gamma,cs2top,cs2bot
      use Gravity, only: gravz
      use Mpicomm, only: mpibcast_real
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: tmp,ztop,zbot
      integer, parameter :: prof_nz=150
      real, dimension (prof_nz) :: prof_lnT,prof_z
      real :: tmprho,tmpT,tmpdT,tmpz,dz_step,lnrho_0,rho0
      integer :: i,lend,j
!
      ! file location settings
      character (len=*), parameter :: lnT_dat = 'driver/b_lnT.dat'
!
      ! temperature given as function lnT(z) in SI units
      ! [T] = K   &   [z] = Mm   & [rho] = kg/m^3
      if (pretend_lnTT) print*,'temp_hydrostatic: not implemented for pretend_lnTT=T'
!
      if (lnrho0 > 0.99*alog(impossible)) then
        call warning("lnrho0 from eos module not useful","use rho_const from density instead")
        lnrho_0 = alog(rho0)
      else
        lnrho_0 = lnrho0
      endif
!
      if (lroot) then
        inquire(IOLENGTH=lend) tmp
        open (10,file=lnT_dat,form='unformatted',status='unknown',recl=lend*prof_nz)
        read (10) prof_lnT
        read (10) prof_z
        close (10)
      endif
      !
      call mpibcast_real(prof_lnT,prof_nz)
      call mpibcast_real(prof_z,prof_nz)
      !
      prof_z = prof_z*1.e6/unit_length
      prof_lnT = prof_lnT - alog(real(unit_temperature))
      !
      ! get step width
      ! should be smaler than grid width and
      ! data width
      !
      dz_step = min((prof_z(2)-prof_z(1)),minval(1./dz_1))
      dz_step = dz_step/10.
      !
      do j=n1,n2
         tmprho = lnrho_0
         tmpT = prof_lnT(1)
         tmpz = prof_z(1)
         !
         ztop=xyz0(3)+Lxyz(3)
         zbot=xyz0(3)
         !
         do while (tmpz <= ztop)
            if (abs(tmpz-zbot) < dz_step) cs2bot = (gamma-1.)*exp(tmpT)
            if (abs(tmpz-ztop) < dz_step) cs2top = (gamma-1.)*exp(tmpT)
            if (abs(tmpz-z(j)) <= dz_step) then
               f(:,:,j,ilnrho) = tmprho
               f(:,:,j,ilnTT)  = tmpT
            endif
            ! new z coord
            tmpz = tmpz+dz_step
            ! get T at new z
            do i=1,prof_nz-1
               if (tmpz >= prof_z(i)  .and. tmpz < prof_z(i+1) ) then
                  tmpdT = (prof_lnT(i+1)-prof_lnT(i))/(prof_z(i+1)-prof_z(i)) * (tmpz-prof_z(i)) + prof_lnT(i) -tmpT
                  tmpT = tmpT + tmpdT
               elseif (tmpz >= prof_z(prof_nz)) then
                  tmpdT = prof_lnT(prof_nz) - tmpT
                  tmpT = tmpT + tmpdT
               endif
            enddo
            tmprho = tmprho - tmpdT + gamma/(gamma-1.)*gravz*exp(-tmpT) * dz_step
         enddo
      enddo
!
    endsubroutine temp_hydrostatic
!***********************************************************************
    subroutine const_lou(ampl,f,i)
!
!  PLEASE ADD A DESCRIPTION
!
!  5-nov-05/weezy: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl
      integer::i
!
      do n=n1,n2; do m=m1,m2
        f(l1:l2,m,n,i  )=ampl*cos(2.*pi*y(m))/32.*pi
        f(l1:l2,m,n,i+1)=ampl*cos(2.*pi*z(n))/32.*pi
        f(l1:l2,m,n,i+2)=ampl*cos(2.*pi*x(l1:l2))/32.*pi
      enddo; enddo
!
    endsubroutine const_lou
!***********************************************************************
    subroutine ferriere_uniform_x(ampl,f,i)
!
!  Uniform B_x field propto rho (for vector potential)
!
!  This routine sets up an initial magnetic field x-parallel with a
!  magnitude directly proportional to the density. In entropy.f90 we require
!  Galactic-hs or Ferriere-hs to be set for init_ss, in density.f90
!  Galactic-hs should be set for initlnrho and in gravity_simple.f90 use
!  Ferriere for gravz_profile
!
!  09-jan-10/fred: coded
!
      use Mpicomm, only: mpireduce_sum, mpibcast_real
!
      integer :: i,icpu
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,tmp1
      real, dimension(1)::tmp3
      real, dimension(ncpus)::sumtmp,tmp2
!
      tmp2(:)=0.0
      sumtmp(:)=0.0
      if (ldensity_nolog) then
        tmp1=sum(f(l1:l2,m1,n1:n2,irho))
      else
        tmp1=sum(exp(f(l1:l2,m1,n1:n2,ilnrho)))
      endif
!
!  Calculate the total mass on each processor tmp1 and identify it with the
!  appropriate processor in the array tmp2
!
      do icpu=1,ncpus
        tmp3=tmp1
        call mpibcast_real(tmp3,1,icpu-1)
        tmp2(icpu)=tmp3(1)
      enddo
!
!  If nprocz is 1 then start summing mass below from zero (sumtmp above).
!  Otherwise sum the masses on the processors below from which to start
!  summing the mass on this processor.
!
      if (ncpus>nprocy) then
        do icpu=nprocy+1,ncpus
          sumtmp(icpu)=sumtmp(icpu-nprocy)+tmp2(icpu-nprocy)
        enddo
      endif
      if (lroot) print*,'sumtmp =',sumtmp
      print*,'sumtmp on iproc =',sumtmp(iproc+1),iproc
      if (ampl==0) then
        f(:,:,:,i:i+2)=0
        if (lroot) print*,'ferriere_uniform_x: set B field to zero; i=',i
      else
        print*,'ferriere_uniform_x: uniform x-field propto rho; i=',i
        if ((ip<=16).and.lroot) print*,'uniform_x: ampl=',ampl
        do n=n1,n2
        do m=m1,m2
          f(l1:l2,m,n,i  )=0.0
          if (ldensity_nolog) then
            f(l1:l2,m,n,i+1)=ampl*(sumtmp(iproc+1)+&
                sum(f(l1:l2,m,n1:n,irho)))*dx*dz
          else
            f(l1:l2,m,n,i+1)=ampl*(sumtmp(iproc+1)+&
                sum(exp(f(l1:l2,m,n1:n,ilnrho))))*dx*dz
          endif
          f(l1:l2,m,n,i+1)=-ampl*(sumtmp(iproc+1)+&
              sum(exp(f(l1:l2,m,n1:n,ilnrho))))*dx*dz
          f(l1:l2,m,n,i+2)=0.0
        enddo
        enddo
      endif
!
    endsubroutine ferriere_uniform_x
!*****************************************************************************
    subroutine ferriere_uniform_y(ampl,f,i)
!
!  Uniform B_y field (for vector potential)
!
!  22-jan-10/fred
!
!  This routine sets up an initial magnetic field y-parallel(azimuthal) with a
!  magnitude directly proportional to the density. In entropy.f90 we require
!  Galactic-hs or Ferriere-hs to be set for init_ss, in density.f90
!  Galactic-hs should be set for initlnrho and in gravity_simple.f90 use
!  Ferriere for gravz_profile
!
      use Mpicomm, only: mpireduce_sum, mpibcast_real
!
      integer :: i ,icpu
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,tmp1
      real, dimension(1)::tmp3
      real, dimension(ncpus)::sumtmp,tmp2
!
      tmp2(:)=0.0
      sumtmp(:)=0.0
      if (ldensity_nolog) then
        tmp1=sum(f(l1:l2,m1,n1:n2,irho))
      else
        tmp1=sum(exp(f(l1:l2,m1,n1:n2,ilnrho)))
      endif
!
!  Calculate the total mass on each processor tmp1 and identify it with the
!  appropriate processor in the array tmp2
!
      do icpu=1,ncpus
        tmp3=tmp1
        call mpibcast_real(tmp3,1,icpu-1)
        tmp2(icpu)=tmp3(1)
      enddo
!
!  If nprocz is 1 then start summing mass below from zero (sumtmp above).
!  Otherwise sum the masses on the processors below from which to start
!  summing the mass on this processor.
!
      if (ncpus>nprocy) then
        do icpu=nprocy+1,ncpus
          sumtmp(icpu)=sumtmp(icpu-nprocy)+tmp2(icpu-nprocy)
        enddo
      endif
      if (lroot) print*,'sumtmp =',sumtmp
      print*,'sumtmp on iproc =',sumtmp(iproc+1),iproc
      if (ampl==0) then
        f(:,:,:,i:i+2)=0
        if (lroot) print*,'ferriere_uniform_y: set B field to zero; i=',i
      else
        print*,'ferriere_uniform_y: uniform y-field approx rho; i=',i
        if ((ip<=16).and.lroot) print*,'uniform_y: ampl=',ampl
        do n=n1,n2
        do m=m1,m2
          if (ldensity_nolog) then
            f(l1:l2,m,n,i)=ampl*(sumtmp(iproc+1)+&
                sum(f(l1:l2,m,n1:n,irho)))*dx*dz
          else
            f(l1:l2,m,n,i)=ampl*(sumtmp(iproc+1)+&
                sum(exp(f(l1:l2,m,n1:n,ilnrho))))*dx*dz
          endif
          f(l1:l2,m,n,i+1)=0.0
          f(l1:l2,m,n,i+2)=0.0
        enddo
        enddo
      endif
!
    endsubroutine ferriere_uniform_y
!*****************************************************************************
    subroutine rotblob(ampl,incl_alpha,f,i,radius,xsphere,ysphere,zsphere)
!
!  Rigid rotating sphere initial velocity
!  with inclination  angle alpha.
!
!  18-feb-10/mvaisala: coded
!
      integer :: i,j,l
      real, dimension (mx,my,mz,mfarray) :: f
      real, optional :: xsphere,ysphere,zsphere
      real :: omega, radius, phi, rr_rot, ampl, x01=0.
      real :: y01=0., z01=0., x_real, y_real, z_real
      real :: incl_alpha, theta, vel_phi
!
!  Rotating sphere
!
      omega = ampl/radius
      if (present(xsphere)) x01=xsphere
      if (present(ysphere)) y01=ysphere
      if (present(zsphere)) z01=zsphere
      if (omega==0) then
        if (lroot) print*,'The sphere does not rotate!'
      else
        do n = n1,n2
          do m = m1,m2
            do l = l1,l2
              x_real = x(l) - x01
              y_real = y(m) - y01
              z_real = z(n) - z01
              rr_rot = sqrt((x_real)**2+(y_real)**2+(z_real)**2)
              theta  = atan(abs(z_real)/sqrt((x_real)**2+(y_real)**2))
              phi    = atan(y_real/x_real)
              vel_phi = omega*rr_rot*cos(theta)
              if (rr_rot <= radius) then
                j = i
                f(l,m,n,j) = vel_phi*sin(phi)
                j = i+1
                f(l,m,n,j) = -vel_phi*cos(phi)*cos(incl_alpha)
                j = i+2
                f(l,m,n,j) = vel_phi*cos(phi)*sin(incl_alpha)
                if (x_real < 0.0) then
                  j = i
                  f(l,m,n,j) = -f(l,m,n,j)
                  j = i+1
                  f(l,m,n,j) = -f(l,m,n,j)
                endif
              endif
            enddo
          enddo
        enddo
      endif
!
    endsubroutine rotblob
!***********************************************************************
    subroutine pre_stellar_cloud(f, datafile, mass_cloud,  &
        cloud_mode, T_cloud_out_rel, dens_coeff, &
        temp_coeff, temp_trans, temp_coeff_out)
!
!  Creates a isothermal or modified Bonnor-Ebert Sphere to be used as
!  a prestellar cloud.
!
!  datafile: The file that includes the density and temp distribution
!  BE_resolution: The resolution of isothemal BE-solutions
!  mass_cloud: mass of the cloud in solar masses
!
!  13-jul-10/mvaisala: created
!
      use EquationOfState, only: eoscalc,ilnrho_lnTT
!
      integer :: jj, test, n, m, l, len_file
      integer, parameter :: BE_resolution = 2000
      logical :: exist !
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (BE_resolution) :: lnTT_file, lnrho_r, r_rho
      real :: tmp, var1, var2, var3
      real :: bigr, mass_cloud
      real :: x01 = 0.0, y01 = 0.0, z01 = 0.0
      real :: x_real, y_real, z_real, rr_box, counter
      real :: M_sun = 1.98892e30 ! kg (SI)
      real :: T_cloud_out_rel, lnTTpoint
      real :: lnTTpoint0, dens_coeff, temp_coeff, temp_coeff_out
      real :: x_wave, wavelength, QQQ, temp_trans, T_cloud_out_rel0
!
      character (len=labellen) :: datafile, cloud_mode
!
      write (*,*) 'Solar mass:', M_sun
      write (*,*) 'unit_mass:', unit_mass
      write (*,*) 'unit_temperature', unit_temperature
      mass_cloud = (mass_cloud * M_sun) / unit_mass
!
      select case (cloud_mode)
!
         case ('read_modified')
!
! Read a radial density and temperature distribution from a file with a chosen
! name. The file must be in the format:
! 1st line: BIGR (cm)    BE_RESOLUTION (int)
! others:   RADIUS (cm)  DENSITY (g/cm^3)  TEMPERATURE (K)
!
           inquire(file=datafile, exist=exist)
           if (exist) then
             open(19, file=datafile)
           else
             call fatal_error('Bonnor-Ebert Sphere', 'No input file')
           end if
           read(19,*) var1, len_file
           bigr = var1/unit_length
           write (*,*) '(Modified BE-sphere) R = ', var1, 'cm =',&
                        bigr, 'pc_units'
           do jj = 1, len_file
             read(19,*) var1, var2, var3
             r_rho(jj) = var1/unit_length
             lnrho_r(jj) = log(var2/unit_density)
             lnTT_file(jj) = log(var3/unit_temperature)
           end do
!
           counter = 0
           do n = n1,n2
             do m = m1,m2
               do l = l1,l2
                 x_real = x(l) - x01
                 y_real = y(m) - y01
                 z_real = z(n) - z01
                 rr_box = sqrt((x_real)**2+(y_real)**2+(z_real)**2)
                 test = 0
                 if (rr_box <= bigr) then
                   do while (r_rho(test+1) <= rr_box)
                      test = test + 1
                   end do
                   f(l,m,n,ilnrho) = lnrho_r(test)
                   call eoscalc(ilnrho_lnTT, lnrho_r(test), lnTT_file(test),&
                                ss=tmp)
                   f(l,m,n,iss) = tmp
                   counter = counter+1
                 else
                   do while (r_rho(test+1) <= bigr)
                      test = test + 1
                   end do
                   if (rr_box <= temp_trans*bigr) then
                      x_wave = rr_box-bigr
                      wavelength = 2.0*bigr*(temp_trans-1)
                      QQQ = 0.5*(sin(2.0*pi*x_wave/wavelength - pi/2.0) + 1.0)
                      T_cloud_out_rel0 = 1 + QQQ*(T_cloud_out_rel-1.0)
                   else
                      T_cloud_out_rel0 = T_cloud_out_rel
                   end if
                   f(l,m,n,ilnrho) = log(exp(lnrho_r(test))/T_cloud_out_rel0)
                   lnTTpoint = log(exp(lnTT_file(test))* &
                       T_cloud_out_rel0*temp_coeff_out)
                   call eoscalc(ilnrho_lnTT, f(l,m,n,ilnrho), &
                                lnTTpoint, ss=tmp)
                   f(l,m,n,iss) = tmp
                 end if
               end do
             end do
           end do
           write (*,*) 'Covers:', counter, '/', &
                       (abs(n1-n2)*abs(m1-m2)*abs(l1-l2)), 'cells'
!
         case ('read_isothermal')
!
! Read a radial density distribution from a file with a chosen
! name. The file must be in the format:
! 1st line: BIGR (cm)    BE_RESOLUTION (int)   TEMPERATURE (K)
! others:   RADIUS (cm)  DENSITY (g/cm^3)
!
           inquire(file=datafile, exist=exist)
           if (exist) then
             open(19, file=datafile)
           else
             write (*,*) datafile
             call fatal_error('Bonnor-Ebert Sphere', 'No input file')
           end if
           read(19,*) var1, len_file, var2
           bigr = var1/unit_length
           lnTTpoint0 = log(var2*temp_coeff/unit_temperature)
           write (*,*) '(Isothermal BE-sphere) R =', var1, 'cm =',&
                        bigr, 'pc_units, T =', var2, 'K'
           do jj = 1, len_file
             read(19,*) var1, var2
             r_rho(jj) = var1/unit_length
             lnrho_r(jj) = log(var2*dens_coeff/unit_density)
           end do
           write (*,*) 'Temperature, lnTT = ', lnTTpoint0
!
           counter = 0
           do n = n1,n2
             do m = m1,m2
               do l = l1,l2
                 lnTTpoint = lnTTpoint0
                 x_real = x(l) - x01
                 y_real = y(m) - y01
                 z_real = z(n) - z01
                 rr_box = sqrt((x_real)**2+(y_real)**2+(z_real)**2)
                 test = 0
                 if (rr_box <= bigr) then
                   do while (r_rho(test+1) <= rr_box)
                     test = test + 1
                   end do
                   f(l,m,n,ilnrho) = lnrho_r(test)
                   call eoscalc(ilnrho_lnTT, f(l,m,n,ilnrho), lnTTpoint, ss=tmp)
                   f(l,m,n,iss) = tmp
                   counter = counter+1
                 else
                   do while (r_rho(test+1) <= bigr)
                     test = test + 1
                   end do
                   if (rr_box <= temp_trans*bigr) then
                      x_wave = rr_box-bigr
                      wavelength = 2.0*bigr*(temp_trans-1)
                      QQQ = 0.5*(sin(2.0*pi*x_wave/wavelength - pi/2.0) + 1.0)
                      T_cloud_out_rel0 = 1 + QQQ*(T_cloud_out_rel-1.0)
                   else
                      T_cloud_out_rel0 = T_cloud_out_rel
                   end if
                   !WRITE (*,*) 'T_cloud_out_rel:', T_cloud_out_rel0
                   f(l,m,n,ilnrho) = log(exp(lnrho_r(test))/T_cloud_out_rel0)
                   lnTTpoint = log(exp(lnTTpoint0)* &
                       T_cloud_out_rel0*temp_coeff_out)
                   call eoscalc(ilnrho_lnTT, f(l,m,n,ilnrho), lnTTpoint, ss=tmp)
                   f(l,m,n,iss) = tmp
                 end if
               end do
             end do
           end do
           write (*,*) 'Covers:', counter, '/', &
                       (abs(n1-n2)*abs(m1-m2)*abs(l1-l2)), 'cells'
           write (*,*) (temp_trans*bigr-bigr), wavelength/2.0, 4.0*wavelength/2.0
!
         case default
      end select
!
    endsubroutine pre_stellar_cloud
!***********************************************************************
endmodule Initcond
