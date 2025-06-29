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
  use Mpicomm, only : ipx, ipy, ipz
!
  implicit none
!
  private
!
  include 'eos_params.h'

  public :: arcade_x, vecpatternxy, bipolar, bipolar_restzero
  public :: soundwave,sinwave,sinwave_phase,coswave,coswave_phase,cos_cos_sin
  public :: hatwave
  public :: acosy
  public :: sph_constb,tanh_hyperbola
  public :: gaunoise, posnoise, posnoise_rel
  public :: gaunoise_rprof
  public :: gaussian, gaussian3d, gaussianpos
  public :: ABC_field, beltrami, bessel_x, bessel_az_x
  public :: beltramik_general, beltrami_general, beltrami_complex
  public :: beltrami_old, bhyperz, bihelical
  public :: straining, rolls, tor_pert
  public :: jump, bjump, bjumpz, stratification, stratification_x
  public :: stratification_xz
  public :: modes, modev, modeb, crazy, exponential
  public :: trilinear, baroclinic
  public :: triquad, isotdisk
  public :: diffrot, olddiffrot
  public :: const_omega
  public :: powern, power_randomphase, power_randomphase_hel, bunch_davies
  public :: planet, planet_hc
  public :: random_isotropic_KS, random_isotropic_shell
  public :: htanh, vtube, vtube_peri, htube, htube2, htube_x, hat, hat3d
  public :: htube_erf, xpoint, xpoint2
  public :: htube2_x
  public :: Gaussian_By_z
  public :: wave_uu, wave, parabola, linprof
  public :: sinxsinz, cosx_cosy_cosz, cosx_coscosy_cosz
  public :: x_siny_cosz, x1_siny_cosz, x32_siny_cosz, x1_cosy_cosz, lnx_cosy_cosz
  public :: sinx_siny_sinz, sinx_cosy_cosz, cosx_siny_cosz, sinx_siny_cosz
  public :: sin2x_sin2y_cosz, cos2x_cos2y_cos2z, x3_cosy_cosz, x3_siny_cosz
  public :: cosx_cosz, cosy_cosz, cosy_sinz
  public :: cosxz_cosz, cosyz_sinz
  public :: halfcos_x, halfcos_z, magsupport, vfield
  public :: uniform_x, uniform_y, uniform_z, uniform_phi, phi_comp_over_r
  public :: vfluxlayer, hfluxlayer, hfluxlayer_y, hfluxlayer_y_theta
  public :: vortex_2d
  public :: vfield2
  public :: hawley_etal99a
  public :: robertsflow, rotated_robertsflow
  public :: const_lou
  public :: corona_init,mdi_init,mag_init,mag_Az_init,file_init,temp_hydrostatic
  public :: innerbox
  public :: couette, couette_rings
  public :: strange,phi_siny_over_r2
  public :: ferriere_uniform_x, ferriere_uniform_y
  public :: rotblob, rotblob_yz, pre_stellar_cloud, dipole, dipole_tor, switchback
  public :: read_outside_scal_array, read_outside_vec_array
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
      real, optional :: kx,ky,kz,KKx,KKy,KKz
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
!  Gaussian scale heights
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
!
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
    subroutine x3_siny_cosz(ampl,f,i,x1,x2,ky,kz)
!
!  sinusoidal wave, adapted from sinxsinz (that routine was already doing
!  this, but under a different name)
!
!   8-may-12/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real,optional :: ky,kz
      real :: ampl, ky1=0., kz1=pi/2., x1, x2
!
!  wavenumber k, helicity H=ampl (can be either sign)
!
!  (x-x1)*(x2-x)*sin(kz*z)
!
      if (present(ky)) ky1=ky
      if (present(kz)) kz1=kz
      if (ampl==0) then
        if (lroot) print*,'x3_siny_cosz: ampl=0'
      else
        if (lroot) write(*,wave_fmt1) 'x3_siny_cosz: ampl,ky,kz=', &
                                      ampl,ky1,kz1
        f(:,:,:,i)=f(:,:,:,i)+ampl*(spread(spread((x-x1)*(x2-x),2,my),3,mz)&
                                   *spread(spread(sin(ky1*y),1,mx),3,mz)&
                                   *spread(spread(cos(kz1*z),1,mx),2,my))
      endif
!
    endsubroutine x3_siny_cosz
!***********************************************************************
    subroutine x_siny_cosz(ampl,f,i,kx,ky,kz,xbot, nexp)
!
!  sinusoidal wave, adapted from sinxsinz (that routine was already doing
!  this, but under a different name)
!
!   2-dec-03/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real, optional :: kx, ky, kz, xbot, nexp
      real :: ampl, kx1=pi/2., ky1=0., kz1=pi/2., xbot_, nexp_
!
!  wavenumber k, helicity H=ampl (can be either sign)
!
!  sinx(kx*x)^nexp*sin(kz*z)
!
      if (present(kx)) kx1=kx
      if (present(ky)) ky1=ky
      if (present(kz)) kz1=kz
      if (present(xbot)) then
        xbot_=xbot
      else
        xbot_=0.
      endif
!
!  If nexp=0, we have the usual sin(theta) dependence.
!  To concentrate the field more strongly toward the equator,
!  we can ut nexp=1 or larger. (The actual exponent must be even.)
!
      if (present(nexp)) then
        nexp_=nexp
      else
        nexp_=0.
      endif
!
      if (ampl==0) then
        if (lroot) print*,'x_siny_cosz: ampl=0'
      else
        if (lroot) write(*,wave_fmt1) 'x_siny_cosz: ampl,kx,ky,kz=', &
                                      ampl,kx1,ky1,kz1
        f(:,:,:,i)=f(:,:,:,i)+ampl*(spread(spread((x-xbot_),2,my),3,mz)&
                                   *spread(spread(sin(ky1*y)**(1.+2.*nexp_),1,mx),3,mz)&
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
    subroutine x32_siny_cosz(ampl,f,i,kx,ky,kz,phasey)
!
!  sinusoidal wave, adapted from sinxsinz (that routine was already doing
!  this, but under a different name)
!
!  29-may-24/axel: coded
!
      use General, only: roptest
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real,optional :: kx,ky,kz,phasey
!
      real :: ampl,kx1,ky1,kz1,phasey1
!
!  wavenumber k, helicity H=ampl (can be either sign)
!
!  sinx(kx*x)*sin(kz*z)
!
      if (ampl==0) then
        if (lroot) print*,'x32_siny_cosz: ampl=0'
      else
        kx1 = roptest(kx,pi/2.)
        ky1 = roptest(ky,0.)
        kz1 = roptest(kz,pi/2.)
        phasey1 = roptest(phasey,0.)

        if (lroot) write(*,wave_fmt1) 'x32_siny_cosz: ampl,kx,ky,kz=',ampl,kx1,ky1,kz1
        f(:,:,:,i)=f(:,:,:,i)+ampl*(spread(spread(1./x**3-x**2,2,my),3,mz) &
                                   *spread(spread(sin(ky1*y+phasey1),1,mx),3,mz) &
                                   *spread(spread(cos(kz1*z),1,mx),2,my))
      endif
!
    endsubroutine x32_siny_cosz
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
        f(:,:,:,i)=f(:,:,:,i)+ampl*(spread(spread(log(    x),2,my),3,mz)&
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
    subroutine sinx_cosy_cosz(ampl,f,i,kx,ky,kz)
!
!  sinusoidal wave, adapted from cosx_siny_cosz
!
!  26-mar-25/TP: coded
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
        if (lroot) print*,'sinx_cosy_cosz: ampl=0'
      else
        if (lroot) write(*,wave_fmt1) 'sinx_cosy_cosz: ampl,kx,ky,kz=', &
                                      ampl,kx1,ky1,kz1
        f(:,:,:,i)=f(:,:,:,i)+ampl*(spread(spread(sin(kx1*x),2,my),3,mz)&
                                   *spread(spread(cos(ky1*y),1,mx),3,mz)&
                                   *spread(spread(cos(kz1*z),1,mx),2,my))
      endif
!
    endsubroutine sinx_cosy_cosz
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
        om_diff(1:nr+1)=om_all(1:nr+1)-om_all(2:nr+2) ! difference in omega from ring to ring
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
    subroutine gaussianpos(ampl,f,i,radius,posx,posy,posz,prefac,slope)
!
!  gaussian 3-D bump centered in specific position
!
!  21-sep-09/rplasson: coded from gaussian3d
!  Maybe could have been done by extending gaussian3d, but didn't want to interfere
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,radius,posx, posy, posz, radius21, alp
      real, dimension (nx) :: phase_factor
      real, optional :: slope
      character(len=*), optional :: prefac
!
!  compute slope
!
      if (present(slope)) then
        alp=slope
      else
        alp=-1.
      endif
!
!  compute phase_factor
!
      radius21=1./radius**2
      do n=n1,n2; do m=m1,m2
        if (present(prefac)) then
          select case (prefac)
          case ('cosx'); phase_factor=+cos(alp*x(l1:l2))
          case ('sinx'); phase_factor=-sin(alp*x(l1:l2))
          case default; phase_factor=1.
          endselect
        else
          phase_factor=1.
        endif
!
!  add Gaussian blob to f(l1:l2,m,n,i) with phase_factor
!
        f(l1:l2,m,n,i)=f(l1:l2,m,n,i)+ampl*phase_factor &
            *exp(-((x(l1:l2)-posx)**2+(y(m)-posy)**2+(z(n)-posz)**2)*radius21)
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
      real, optional :: kx,ky,kz
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
      real, optional :: kx,ky,kz
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
      real, optional :: kx,ky,kz
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
      real, optional :: kx,ky,kz
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
    subroutine jump(f,i,fleft,fright,width,xmid,ymid,zmid,dir)
!
!  jump
!
!  19-sep-02/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my) :: profxy
      real, dimension (mx) :: profx
      real, dimension (my) :: profy
      real, dimension (mz) :: profz
      real :: fleft,fright,width
      real :: xmid,ymid,zmid
      character(len=*) :: dir
      integer :: l,m
!
!  jump; check direction
!
      select case (dir)
!
      case ('x')
        profx=fleft+(fright-fleft)*.5*(1.+tanh((x-xmid)/width))
        f(:,:,:,i)=f(:,:,:,i)+spread(spread(profx,2,my),3,mz)
!
      case ('y')
        profy=fleft+(fright-fleft)*.5*(1.+tanh((y-ymid)/width))
        f(:,:,:,i)=f(:,:,:,i)+spread(spread(profy,1,mx),3,mz)
!
      case ('z')
        profz=fleft+(fright-fleft)*.5*(1.+tanh((z-zmid)/width))
        f(:,:,:,i)=f(:,:,:,i)+spread(spread(profz,1,mx),2,my)
!
!  2-D shocks
!
      case ('xy')
        do l=1,mx
        do m=1,my
          profxy(l,m)=fleft+(fright-fleft)*.25* &
            (1.+tanh((x(l)-xmid)/width))*(1.+tanh((y(m)-ymid)/width))
        enddo
        enddo
        f(:,:,:,i)=f(:,:,:,i)+spread(profxy,3,mz)
!
!  2-D diagonal shocks
!
      case ('x-y')
        do l=1,mx
        do m=1,my
          profxy(l,m)=fleft+(fright-fleft)*.5*(1.+tanh(((x(l)-xmid)+(y(m)-ymid))/width))
        enddo
        enddo
        f(:,:,:,i)=f(:,:,:,i)+spread(profxy,3,mz)
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
        alog_cosh_xwidth=abs(x/width)+log(.5*(1.+exp(-2*abs(x/width))))
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
        alog_cosh_zwidth=abs(z/width)+log(.5*(1.+exp(-2*abs(z/width))))
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
    subroutine beltrami_old(ampl,f,i,kx,ky,kz,kx2,ky2,kz2,phase)
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
      real, optional :: kx,ky,kz,kx2,ky2,kz2,phase
      real :: ampl,k=1.,ph
!
!  This routine should be removed by 2020
!
      call fatal_error('beltrami_old','email Axel if you need this routine')
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
        if (k==0) print*,'beltrami_old: k must not be zero!'
        cfuncx=sign(sqrt(abs(ampl/k)),kx)*cos(k*x+ph)
        sfuncx=sign(sqrt(abs(ampl/k)),kx)*sin(k*x+ph)
        if (present(kx2)) sfuncx=sfuncx*sin(kx2*x+ph)
        if (ampl==0) then
          if (lroot) print*,'beltrami_old: ampl=0; kx=',k
        elseif (ampl>0) then
          if (lroot) print*,'beltrami_old: Beltrami field (pos-hel): kx,i=',k,i
          j=i+1; f(:,:,:,j)=f(:,:,:,j)+spread(spread(sfuncx,2,my),3,mz)
          j=i+2; f(:,:,:,j)=f(:,:,:,j)+spread(spread(cfuncx,2,my),3,mz)
        elseif (ampl<0) then
          if (lroot) print*,'beltrami_old: Beltrami field (neg-hel): kx,i=',k,i
          j=i+1; f(:,:,:,j)=f(:,:,:,j)+spread(spread(cfuncx,2,my),3,mz)
          j=i+2; f(:,:,:,j)=f(:,:,:,j)+spread(spread(sfuncx,2,my),3,mz)
        endif
      endif
!
!  set y-dependent Beltrami field
!
      if (present(ky)) then
        k=abs(ky)
        if (k==0) print*,'beltrami_old: k must not be zero!'
        cfuncy=sign(sqrt(abs(ampl/k)),ky)*cos(k*y+ph)
        sfuncy=sign(sqrt(abs(ampl/k)),ky)*sin(k*y+ph)
        if (present(ky2)) sfuncy=sfuncy*sin(ky2*y+ph)
        if (ampl==0) then
          if (lroot) print*,'beltrami_old: ampl=0; ky=',k
        elseif (ampl>0) then
          if (lroot) print*,'beltrami_old: Beltrami field (pos-hel): ky,i=',k,i
          j=i;   f(:,:,:,j)=f(:,:,:,j)+spread(spread(cfuncy,1,mx),3,mz)
          j=i+2; f(:,:,:,j)=f(:,:,:,j)+spread(spread(sfuncy,1,mx),3,mz)
        elseif (ampl<0) then
          if (lroot) print*,'beltrami_old: Beltrami field (neg-hel): ky,i=',k,i
          j=i;   f(:,:,:,j)=f(:,:,:,j)+spread(spread(sfuncy,1,mx),3,mz)
          j=i+2; f(:,:,:,j)=f(:,:,:,j)+spread(spread(cfuncy,1,mx),3,mz)
        endif
      endif
!
!  set z-dependent Beltrami field
!
      if (present(kz)) then
        k=abs(kz)
        if (k==0) print*,'beltrami_old: k must not be zero!'
        cfuncz=sign(sqrt(abs(ampl/k)),kz)*cos(k*z+ph)
        sfuncz=sign(sqrt(abs(ampl/k)),kz)*sin(k*z+ph)
        if (present(kz2)) sfuncz=sfuncz*sin(kz2*z+ph)
        if (ampl==0) then
          if (lroot) print*,'beltrami_old: ampl=0; kz=',k
        elseif (ampl>0) then
          if (lroot) print*,'beltrami_old: Beltrami field (pos-hel): kz,i=',k,i
          j=i;   f(:,:,:,j)=f(:,:,:,j)+spread(spread(sfuncz,1,mx),2,my)
          j=i+1; f(:,:,:,j)=f(:,:,:,j)+spread(spread(cfuncz,1,mx),2,my)
        elseif (ampl<0) then
          if (lroot) print*,'beltrami_old: Beltrami field (neg-hel): kz,i=',k,i
          j=i;   f(:,:,:,j)=f(:,:,:,j)+spread(spread(cfuncz,1,mx),2,my)
          j=i+1; f(:,:,:,j)=f(:,:,:,j)+spread(spread(sfuncz,1,mx),2,my)
        endif
      endif
!
    endsubroutine beltrami_old
!***********************************************************************
    subroutine ABC_field(f,i,kx,ky,kz,ABC,x0,y0,z0,width,sigma)
!
!  ABC field (as initial condition)
!
!  24-aug-19/axel: coded
!
      use Sub, only: cubic_step
!
      integer :: i,j,l,m,n
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx) :: sfuncx,cfuncx,xprof
      real, dimension (my) :: sfuncy,cfuncy,yprof
      real, dimension (mz) :: sfuncz,cfuncz,zprof
      real, dimension (3) :: ABC, width
      real :: kx,ky,kz, x0,y0,z0, prof, sigma1
      real, optional :: sigma
!
      if (present(sigma)) then
        sigma1=sigma
      else
        sigma1=1.
      endif
!
!  Envelope profile (can turn it off by setting x0=0.).
!
      if (lroot) print*,'ABC_field: x0,width=',x0,width(1)
      if (x0>0.) then
        xprof=cubic_step(x,-x0,width(1))-cubic_step(x,x0,width(1))
        if (lroot) print*,'xprof=',xprof
      else
        xprof=1.
      endif
!
!  Envelope profile (can turn it off by setting y0=0.).
!
      if (lroot) print*,'ABC_field: y0,width=',y0,width(2)
      if (y0>0.) then
        yprof=cubic_step(y,-y0,width(2))-cubic_step(y,y0,width(2))
        if (lroot) print*,'yprof=',yprof
      else
        yprof=1.
      endif
!
!  Envelope profile (can turn it off by setting z0=0.).
!
      if (lroot) print*,'ABC_field: z0,width=',z0,width(3)
      if (z0>0.) then
        zprof=cubic_step(z,-z0,width(3))-cubic_step(z,z0,width(3))
        if (lroot) print*,'zprof=',zprof
      else
        zprof=1.
      endif
!
!  Set x-dependent part of ABC field.
!
      sfuncx=ABC(2)*sin(kx*x)*xprof
      cfuncx=ABC(2)*cos(kx*x)*xprof
!
!  Set y-dependent part of ABC field.
!
      sfuncy=ABC(3)*sin(ky*y)*yprof
      cfuncy=ABC(3)*cos(ky*y)*yprof
!
!  Set z-dependent part of ABC field.
!
      sfuncz=ABC(1)*sin(kz*z)*zprof
      cfuncz=ABC(1)*cos(kz*z)*zprof
!
      do n=n1,n2
      do m=m1,m2
      do l=l1,l2
        prof=xprof(l)*yprof(m)*zprof(n)
        j=i  ; f(l,m,n,j)=f(l,m,n,j)+prof*(sigma1*sfuncz(n)+cfuncy(m))
        j=i+1; f(l,m,n,j)=f(l,m,n,j)+prof*(sigma1*sfuncx(l)+cfuncz(n))
        j=i+2; f(l,m,n,j)=f(l,m,n,j)+prof*(sigma1*sfuncy(m)+cfuncx(l))
      enddo
      enddo
      enddo
!
    endsubroutine ABC_field
!***********************************************************************
    subroutine beltrami(ampl,f,i,kx,ky,kz,kx2,ky2,kz2,phase,sigma,z0,width)
!
!  Beltrami field (as initial condition)
!
!  19-jun-02/axel: coded
!   5-jul-02/axel: made additive (if called twice), kx,ky,kz are optional
!
      use Sub, only: cubic_step
!
      integer :: i,j
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx) :: sfuncx,cfuncx
      real, dimension (my) :: sfuncy,cfuncy
      real, dimension (mz) :: sfuncz,cfuncz,zprof
      real, optional :: kx,ky,kz,kx2,ky2,kz2,phase,sigma,z0,width
      real :: ampl,k,ph,sig
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
!  possibility of a fractional helicity
!
      if (present(sigma)) then
        if (lroot) print*,'Beltrami: sigma=',sigma
        sig=sigma
      else
        sig=1.
      endif
!
!  possibility of envelope profile
!
      if (present(z0)) then
        if (.not.present(width)) call fatal_error('beltrami','width?')
        if (lroot) print*,'Beltrami: z0,width=',z0,width
        if (z0>0.) then
          zprof=cubic_step(z,-z0,width)-cubic_step(z,z0,width)
        else
          zprof=1.
        endif
      else
        zprof=1.
      endif
!
!  wavenumber k, helicity H=ampl (can be either sign)
!
!  set x-dependent Beltrami field
!
      if (present(kx)) then
        k=kx
        if (k==0) print*,'beltrami: k must not be zero!'
        cfuncx=ampl*cos(k*x+ph)
        sfuncx=ampl*sin(k*x+ph)
        if (present(kx2)) sfuncx=sfuncx*sin(kx2*x+ph)
        if (ampl==0) then
          if (lroot) print*,'beltrami: ampl=0; kx=',k
        elseif (ampl>0) then
          if (lroot) print*,'beltrami: Beltrami field (pos-hel): kx,i=',k,i
          j=i+1; f(:,:,:,j)=f(:,:,:,j)+spread(spread(sfuncx,2,my),3,mz)*sig
          j=i+2; f(:,:,:,j)=f(:,:,:,j)+spread(spread(cfuncx,2,my),3,mz)
        elseif (ampl<0) then
          if (lroot) print*,'beltrami: Beltrami field (neg-hel): kx,i=',k,i
          j=i+1; f(:,:,:,j)=f(:,:,:,j)+spread(spread(cfuncx,2,my),3,mz)
          j=i+2; f(:,:,:,j)=f(:,:,:,j)+spread(spread(sfuncx,2,my),3,mz)*sig
        endif
      endif
!
!  set y-dependent Beltrami field
!
      if (present(ky)) then
        k=ky
        if (k==0) print*,'beltrami: k must not be zero!'
        cfuncy=ampl*cos(k*y+ph)
        sfuncy=ampl*sin(k*y+ph)
        if (present(ky2)) sfuncy=sfuncy*sin(ky2*y+ph)
        if (ampl==0) then
          if (lroot) print*,'beltrami: ampl=0; ky=',k
        elseif (ampl>0) then
          if (lroot) print*,'beltrami: Beltrami field (pos-hel): ky,i=',k,i
          j=i;   f(:,:,:,j)=f(:,:,:,j)+spread(spread(cfuncy,1,mx),3,mz)
          j=i+2; f(:,:,:,j)=f(:,:,:,j)+spread(spread(sfuncy,1,mx),3,mz)*sig
        elseif (ampl<0) then
          if (lroot) print*,'beltrami: Beltrami field (neg-hel): ky,i=',k,i
          j=i;   f(:,:,:,j)=f(:,:,:,j)+spread(spread(sfuncy,1,mx),3,mz)*sig
          j=i+2; f(:,:,:,j)=f(:,:,:,j)+spread(spread(cfuncy,1,mx),3,mz)
        endif
      endif
!
!  set z-dependent Beltrami field
!
      if (present(kz)) then
        k=kz
        if (k==0) print*,'beltrami: k must not be zero!'
        cfuncz=ampl*cos(k*z+ph)*zprof
        sfuncz=ampl*sin(k*z+ph)*zprof
        if (present(kz2)) sfuncz=sfuncz*sin(kz2*z+ph)
        if (ampl==0) then
          if (lroot) print*,'beltrami: ampl=0; kz=',k
        elseif (ampl>0) then
          if (lroot) print*,'beltrami: Beltrami field (pos-hel): kz,i=',k,i
          j=i;   f(:,:,:,j)=f(:,:,:,j)+spread(spread(sfuncz,1,mx),2,my)*sig
          j=i+1; f(:,:,:,j)=f(:,:,:,j)+spread(spread(cfuncz,1,mx),2,my)
        elseif (ampl<0) then
          if (lroot) print*,'beltrami: Beltrami field (neg-hel): kz,i=',k,i
          j=i;   f(:,:,:,j)=f(:,:,:,j)+spread(spread(cfuncz,1,mx),2,my)
          j=i+1; f(:,:,:,j)=f(:,:,:,j)+spread(spread(sfuncz,1,mx),2,my)*sig
        endif
      endif
!
    endsubroutine beltrami
!***********************************************************************
    subroutine beltramik_general(ampl,f,i,kx,ky,kz,phase)
!
!  Beltrami field (as initial condition)
!  Currently not additive
!
!  16-apr-21/axel: coded
!
      integer :: i,j
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: kx, ky, kz
      real :: ampl, phase
      complex :: phase_factor_x
!
!  preparations; j is the first index of the array for the complex part.
!  phase_factor is the phase factor.
!
      j=i+3
      phase_factor_x=exp(cmplx(0.,phase+xyz0(1)))
!
      if (kx>0.and.kx<=nx-1) then
        if (ipx==0) then
          f(l1+kx,m1,n1,i+1)=+ampl*real(phase_factor_x)
          f(l1+kx,m1,n1,j+2)=-ampl*real(phase_factor_x)
        endif
!
        if (ipx==nprocx-1) then
          f(l2+1-kx,m1,n1,i+1)=+ampl*aimag(phase_factor_x)
          f(l2+1-kx,m1,n1,j+2)=+ampl*aimag(phase_factor_x)
        endif
      endif
!
    endsubroutine beltramik_general
!***********************************************************************
    subroutine beltrami_general(ampl,f,i,kx,ky,kz,phase)
!
!  Beltrami field (as initial condition)
!
!  19-jun-02/axel: coded
!   5-jul-02/axel: made additive (if called twice), kx,ky,kz are optional
!
      integer :: i,j,l,m,n
      real, dimension (mx,my,mz,mfarray) :: f
      real :: kx, ky, kz, phase, k, cfunc, sfunc
      !real :: ex=.1, ey=.33, ez=.58
      !real :: ex=1., ey=.0, ez=1.
      !real :: ex=1., ey=.0, ez=0.
      real :: ex=0., ey=.0, ez=1.
      real :: kxe_x, kxe_y, kxe_z, kxkxe_x, kxkxe_y, kxkxe_z
      real :: ampl
!
!  wavenumber k, helicity H=ampl (can be either sign)
!
      kxe_x=ky*ez-kz*ey
      kxe_y=kz*ex-kx*ez
      kxe_z=kx*ey-ky*ex
!
      kxkxe_x=ky*kxe_z-kz*kxe_y
      kxkxe_y=kz*kxe_x-kx*kxe_z
      kxkxe_z=kx*kxe_y-ky*kxe_x
!
      k=sqrt(kx**2+ky**2+kz**2)
!
      do n=1,mz
      do m=1,my
      do l=1,mx
        cfunc=abs(ampl)*cos(kx*x(l)+ky*y(m)+kz*z(n)+phase)
        sfunc=    ampl *sin(kx*x(l)+ky*y(m)+kz*z(n)+phase)
        j=i  ; f(l,m,n,j)=f(l,m,n,j)+kxkxe_x*cfunc+k*kxe_x*sfunc
        j=i+1; f(l,m,n,j)=f(l,m,n,j)+kxkxe_y*cfunc+k*kxe_y*sfunc
        j=i+2; f(l,m,n,j)=f(l,m,n,j)+kxkxe_z*cfunc+k*kxe_z*sfunc
      enddo
      enddo
      enddo
!
    endsubroutine beltrami_general
!***********************************************************************
    subroutine bihelical(ampl,f,i,kx,ky,kz,kx2,ky2,kz2,phase,sym)
!
!  Bihelical field (as initial condition)
!
!  19-jun-02/axel: coded
!   5-jul-02/axel: made additive (if called twice), kx,ky,kz are optional
!
      integer :: i,j
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx) :: sfuncx,cfuncx
      real, dimension (my) :: sfuncy,cfuncy
      real, dimension (mz) :: sfuncz,cfuncz
      logical, optional :: sym
      real, optional :: kx,ky,kz,kx2,ky2,kz2,phase
      real :: ampl,k=1.,kp,km,ph
!
!  possibility of shifting the Bihelical wave by phase ph
!
      if (present(phase)) then
        if (lroot) print*,'Bihelical: phase=',phase
        ph=phase
      else
        ph=0.
      endif
!
!  wavenumber k, helicity H=ampl (can be either sign)
!
!  set x-dependent Bihelical field
!
      if (present(kx)) then
        if (present(sym)) then
          km=kx-.5
          kp=kx+.5
        else
          km=kx
          kp=kx+1.
        endif
        if (k==0) print*,'bihelical: k must not be zero!'
        cfuncx=ampl*cos(km*x+ph)
        sfuncx=ampl*sin(kp*x+ph)
        if (present(kx2)) sfuncx=sfuncx*sin(kx2*x+ph)
        if (ampl==0) then
          if (lroot) print*,'bihelical: ampl=0; kx=',k
        elseif (ampl>0) then
          if (lroot) print*,'bihelical: Bihelical field (pos-hel): kx,i=',k,i
          j=i+1; f(:,:,:,j)=f(:,:,:,j)+spread(spread(sfuncx,2,my),3,mz)
          j=i+2; f(:,:,:,j)=f(:,:,:,j)+spread(spread(cfuncx,2,my),3,mz)
        elseif (ampl<0) then
          if (lroot) print*,'bihelical: Bihelical field (neg-hel): kx,i=',k,i
          j=i+1; f(:,:,:,j)=f(:,:,:,j)+spread(spread(cfuncx,2,my),3,mz)
          j=i+2; f(:,:,:,j)=f(:,:,:,j)+spread(spread(sfuncx,2,my),3,mz)
        endif
      endif
!
!  set y-dependent Bihelical field
!
      if (present(ky)) then
        if (present(sym)) then
          km=ky-.5
          kp=ky+.5
        else
          km=ky
          kp=ky+1.
        endif
        if (k==0) print*,'bihelical: k must not be zero!'
        cfuncy=ampl*cos(km*y+ph)
        sfuncy=ampl*sin(kp*y+ph)
        if (present(ky2)) sfuncy=sfuncy*sin(ky2*y+ph)
        if (ampl==0) then
          if (lroot) print*,'bihelical: ampl=0; ky=',k
        elseif (ampl>0) then
          if (lroot) print*,'bihelical: Bihelical field (pos-hel): ky,i=',k,i
          j=i;   f(:,:,:,j)=f(:,:,:,j)+spread(spread(cfuncy,1,mx),3,mz)
          j=i+2; f(:,:,:,j)=f(:,:,:,j)+spread(spread(sfuncy,1,mx),3,mz)
        elseif (ampl<0) then
          if (lroot) print*,'bihelical: Bihelical field (neg-hel): ky,i=',k,i
          j=i;   f(:,:,:,j)=f(:,:,:,j)+spread(spread(sfuncy,1,mx),3,mz)
          j=i+2; f(:,:,:,j)=f(:,:,:,j)+spread(spread(cfuncy,1,mx),3,mz)
        endif
      endif
!
!  set z-dependent Bihelical field
!
      if (present(kz)) then
        if (present(sym)) then
          km=kz-.5
          kp=kz+.5
        else
          km=kz
          kp=kz+1.
        endif
        if (k==0) print*,'bihelical: k must not be zero!'
        cfuncz=ampl*cos(km*z+ph)
        sfuncz=ampl*sin(kp*z+ph)
        if (present(kz2)) sfuncz=sfuncz*sin(kz2*z+ph)
        if (ampl==0) then
          if (lroot) print*,'bihelical: ampl=0; kz=',k
        elseif (ampl>0) then
          if (lroot) print*,'bihelical: Bihelical field (pos-hel): kz,i=',k,i
          j=i;   f(:,:,:,j)=f(:,:,:,j)+spread(spread(sfuncz,1,mx),2,my)
          j=i+1; f(:,:,:,j)=f(:,:,:,j)+spread(spread(cfuncz,1,mx),2,my)
        elseif (ampl<0) then
          if (lroot) print*,'bihelical: Bihelical field (neg-hel): kz,i=',k,i
          j=i;   f(:,:,:,j)=f(:,:,:,j)+spread(spread(cfuncz,1,mx),2,my)
          j=i+1; f(:,:,:,j)=f(:,:,:,j)+spread(spread(sfuncz,1,mx),2,my)
        endif
      endif
!
    endsubroutine bihelical
!***********************************************************************
    subroutine bhyperz(ampl,f,i,kz,nfactor)
!
!  Beltrami field with wavelengths that are cube roots of unity
!
      integer :: i
      integer :: ix,iy,iz
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,kz,nfactor,maxAx,maxAy,zmax,zmin,zzm !,maxA
      complex :: omega,Ax,Ay
!
!  set z-dependent Beltrami field
!
      omega=cmplx(1./2.,sqrt(3.)/2.)
      zmin=xyz0(3)
      zmax=xyz1(3)
      zzm=max(abs(zmin),abs(zmax))
      maxAx= abs(sin(omega*kz*zzm)+sin(omega*omega*kz*zzm))
      maxAy= abs(cos(omega*kz*zzm)+cos(omega*omega*kz*zzm))
!      maxA = max(maxAx,maxAy)
      do ix=1,mx; do iy=1,my;do iz=1,mz
        Ax = (1.-nfactor)*cos(kz*z(iz))+ nfactor*real(sin(omega*kz*z(iz)))
        Ay = (1.-nfactor)*sin(kz*z(iz))+ nfactor*real(cos(omega*kz*z(iz)))
!        Ax = real(sin(omega*kz*z(iz)))
!        Ay = real(cos(omega*kz*z(iz)))
        f(ix,iy,iz,i) = ampl*real(Ax)/kz
        f(ix,iy,iz,i+1) = ampl*real(Ay)/kz
      enddo;enddo;enddo
!
    endsubroutine bhyperz
!***********************************************************************
    subroutine beltrami_complex(ampl,f,i,kx,ky,kz,kx2,ky2,kz2,phase)
!
!  Beltrami field (as initial condition)
!
!  23-sep-10/dhruba: adapted from beltrami
!
      integer :: i,j
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx) :: sfuncx,cfuncx
      real, dimension (my) :: sfuncy,cfuncy
      real, dimension (mz) :: sfuncz,cfuncz
      real, optional :: kx,ky,kz,kx2,ky2,kz2,phase
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
    subroutine straining(ampl,f,i,kx,ky,kz,dimensionality)
!
!  convection straining (as initial condition)
!
!  23-mar-16/axel: coded
!
      integer :: i,j,dimensionality
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,kx,ky,kz
!
!  check input parameters
!
      if (lroot) print*,'straining: i,kx,ky,kz=',i,kx,ky,kz
!
!  set stream function psi=sin(kx*x)*sin(kz*(z-zbot))
!
      j=i
      f(:,:,:,j)=f(:,:,:,j)+ampl*spread(spread(sin(kx*x),2,my),3,mz)&
                                *spread(spread(cos(ky*y),1,mx),3,mz)&
                                *spread(spread(cos(kz*z),1,mx),2,my)
      j=i+1
      f(:,:,:,j)=f(:,:,:,j)+ampl*spread(spread(cos(kx*x),2,my),3,mz)&
                                *spread(spread(sin(ky*y),1,mx),3,mz)&
                                *spread(spread(cos(kz*z),1,mx),2,my)
      j=i+2
      f(:,:,:,j)=f(:,:,:,j)+ampl*spread(spread(cos(kx*x),2,my),3,mz)&
                                *spread(spread(cos(ky*y),1,mx),3,mz)&
                                *spread(spread(sin(kz*z),1,mx),2,my)&
                                *(-1.)*(dimensionality-1)
!
    endsubroutine straining
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
    subroutine robertsflow(ampl,f,i,relhel,kx,flow)
!
!  Roberts Flow (as initial condition)
!
!   9-jun-05/axel: coded
!
      integer :: i,j
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,k=1.,kf,fac1,fac2,relhel
      real, optional :: kx
      character (len=labellen) :: flowtype='I'
      character (len=labellen), optional :: flow
!
!  Possibility of changing the wavenumber
!
      if (present(kx)) then
        k=kx
      endif
!
!  Possibility of changing the flow
!
      if (present(flow)) then
        flowtype=flow
      endif
!
!  prepare coefficients
!
      kf=k*sqrt(2.)
      fac1=sqrt(2.)*ampl*k/kf
      fac2=sqrt(2.)*ampl*relhel
!
      if (flowtype=='I-shift' .or. flowtype=='II-shift') then
!
!  shifted by 90 degrees in the x and y directions
!
        j=i+0; f(:,:,:,j)=f(:,:,:,j)+fac1*spread(spread(sin(k*x),2,my),3,mz)&
                                         *spread(spread(cos(k*y),1,mx),3,mz)
!
        j=i+1; f(:,:,:,j)=f(:,:,:,j)-fac1*spread(spread(cos(k*x),2,my),3,mz)&
                                         *spread(spread(sin(k*y),1,mx),3,mz)
!
        if (flowtype=='I-shift') then
          j=i+2; f(:,:,:,j)=f(:,:,:,j)+fac2*spread(spread(sin(k*x),2,my),3,mz)&
                                           *spread(spread(sin(k*y),1,mx),3,mz)
        elseif (flowtype=='II-shift') then
          j=i+2; f(:,:,:,j)=f(:,:,:,j)+fac2*spread(spread(cos(k*x),2,my),3,mz)&
                                           *spread(spread(cos(k*y),1,mx),3,mz)
        else
          call fatal_error('robertsflow','no such flowtype')
        endif
      else
!
!  original, where field = curl(phi*zz)+phi*zz and curl(phi*zz)+tilde(phi)*zz
!  with phi=cosk0x*cosk0y for flows I and II, respectively.
!
        j=i+0; f(:,:,:,j)=f(:,:,:,j)-fac1*spread(spread(cos(k*x),2,my),3,mz)&
                                         *spread(spread(sin(k*y),1,mx),3,mz)
!
        j=i+1; f(:,:,:,j)=f(:,:,:,j)+fac1*spread(spread(sin(k*x),2,my),3,mz)&
                                         *spread(spread(cos(k*y),1,mx),3,mz)
!
        if (flowtype=='I') then
          j=i+2; f(:,:,:,j)=f(:,:,:,j)+fac2*spread(spread(cos(k*x),2,my),3,mz)&
                                           *spread(spread(cos(k*y),1,mx),3,mz)
        elseif (flowtype=='II') then
          j=i+2; f(:,:,:,j)=f(:,:,:,j)+fac2*spread(spread(sin(k*x),2,my),3,mz)&
                                           *spread(spread(sin(k*y),1,mx),3,mz)
        else
          call fatal_error('robertsflow','no such flowtype')
        endif
      endif
!
    endsubroutine robertsflow
!***********************************************************************
    subroutine rotated_robertsflow(ampl,f,i,relhel,kx,flow)
!
!  By 45 degrees rotated Roberts Flow (as initial condition)
!
!   6-jan-25/axel: coded
!
      integer :: i,j
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,k=1.,kf,fac1,fac2,relhel
      real, optional :: kx
      character (len=labellen) :: flowtype='I'
      character (len=labellen), optional :: flow
!
!  Possibility of changing the wavenumber
!
      if (present(kx)) then
        k=kx
      endif
!
!  Possibility of changing the flow
!
      if (present(flow)) then
        flowtype=flow
      endif
!
!  prepare coefficients
!
      kf=k*sqrt(2.)
      fac1=sqrt(2.)*ampl*k/kf
      fac2=sqrt(2.)*ampl*relhel
!
      j=i+0; f(:,:,:,j)=f(:,:,:,j)+fac1*spread(spread(sin(k*y),1,mx),3,mz)
!
      j=i+1; f(:,:,:,j)=f(:,:,:,j)+fac1*spread(spread(sin(k*x),2,my),3,mz)
!
      if (flowtype=='I') then
        j=i+2; f(:,:,:,j)=f(:,:,:,j)+fac1*(spread(spread(cos(k*x),2,my),3,mz)&
                                          -spread(spread(cos(k*y),1,mx),3,mz))
      elseif (flowtype=='II') then
        j=i+2; f(:,:,:,j)=f(:,:,:,j)+fac1*(spread(spread(cos(k*x),2,my),3,mz)&
                                          +spread(spread(cos(k*y),1,mx),3,mz))
      else
        call fatal_error('robertsflow','no such flowtype')
      endif
!
    endsubroutine rotated_robertsflow
!***********************************************************************
    subroutine exponential(ampl,f,j,KKz)
!
!  horizontal pattern with exponential decay (as initial condition)
!
!   2-mar-13/axel: coded
!
      integer :: j
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,KKz,fact
!
!  By(z)=B0*exp(-z/2H), and put KKz=1/2H, so
!  Ax(z)=-2*H*B0*exp(-z/2H)
!
      fact=-ampl/KKz
      f(:,:,:,j)=f(:,:,:,j)+fact*spread(spread(exp(-KKz*z),1,mx),2,my)
!
    endsubroutine exponential
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
    subroutine soundwave(ampl,f,i,kx,ky,kz,width)
!
!  sound wave (as initial condition)
!
!   2-aug-02/axel: adapted from Beltrami
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx) :: envelope_x
      real, optional :: kx,ky,kz,width
      real :: ampl, k=1., fac
!
!  wavenumber k
!
!  set x-dependent sin wave
!
      if (present(kx)) then
        envelope_x=1.
        if (present(width)) then
          if (width/=0.) envelope_x=exp(-.5*(x/width)**2)
        endif
        k=kx; if (k==0) print*,'soundwave: k must not be zero!'; fac=sqrt(abs(ampl/k))
        if (ampl==0) then
          if (lroot) print*,'soundwave: ampl=0; kx=',k
        else
          if (lroot) print*,'soundwave: kx,i=',k,i
          f(:,:,:,i)=f(:,:,:,i)+fac*spread(spread(envelope_x*sin(k*x),2,my),3,mz)
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
      real, optional :: kx, ky, kz
      real :: ampl, k=1., fac
      real :: kx1, ky1, kz1
!
!  wavenumber k
!
!  set x-dependent cos wave
!
      if (present(kx)) then
        k=kx; if (k==0) print*,'coswave: k must not be zero!'
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
!* 26-sep-12/ccyang: This routine does not work with my machine when not
!*                   all k's are present.
!    subroutine coswave(ampl,f,i,kx,ky,kz)
!!
!!  cosine wave (as initial condition)
!!
!!  14-nov-03/axel: adapted from sinwave
!!
!      integer :: i
!      real, dimension (mx,my,mz,mfarray) :: f
!      real,optional :: kx,ky,kz
!      real :: ampl
!!
!      f(:,:,:,i)=f(:,:,:,i)+ampl*cos( &
!        spread(spread(kx*x,2,my),3,mz)+ &
!        spread(spread(ky*y,1,mx),3,mz)+ &
!        spread(spread(kz*z,1,mx),2,my))
!!
!    endsubroutine coswave
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
    subroutine hatwave(ampl,f,i,width,kx,ky,kz,power,pos)
!
!  cosine wave (as initial condition)
!
!   9-jan-08/axel: adapted from coswave
!  22-jul-22/axel: added keyword pos for positive values only
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real, optional :: kx, ky, kz, power, pos
      real :: ampl,k=1.,fac,width,pow=1.
!
!  wavenumber k
!
      if (present(power)) pow=power
!
!  set x-dependent hat wave
!
      if (present(kx)) then
        k=kx; if (k==0) print*,'hatwave: k must not be zero!'; fac=.5*ampl
        if (ampl==0) then
          if (lroot) print*,'hatwave: ampl=0; kx=',k
        else
          if (lroot) print*,'hatwave: kx,i=',k,i
          !f(:,:,:,i)=f(:,:,:,i)+fac*spread(spread(1+tanh((cos(k*x)**pow)/width),2,my),3,mz)
          if (present(pos)) then
            f(:,:,:,i)=f(:,:,:,i)+.5*fac*(1.+spread(spread(tanh((cos(k*x)**pow)/width),2,my),3,mz))
          else
            f(:,:,:,i)=f(:,:,:,i)+fac*spread(spread(tanh((cos(k*x)**pow)/width),2,my),3,mz)
          endif
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
      real, optional :: kx,ky,kz
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
    subroutine sinwave_phase(f,i,ampl,kx,ky,kz,phase,lnorm_kk)
!
!  Sine wave (as initial condition)
!
!  23-jan-06/anders: adapted from sinwave.
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl, kx, ky, kz, phase, fact, k2
      integer :: i
      logical, optional :: lnorm_kk
!
!  possibility to normalize by 1/k if the initial amplitude
!  is meant to be an amplitude for the B-field.
!
      if (present(lnorm_kk)) then
        k2=kx**2+ky**2+kz**2
        if (lnorm_kk .and. k2/=0.) then
          fact=ampl/sqrt(k2)
        else
          fact=ampl
        endif
      else
        fact=ampl
      endif
!
!  Set sin wave
!
      if (lroot) print*, 'sinwave_phase: i, fact, kx, ky, kz, phase=', &
          i, fact, kx, ky, kz, phase
!
      do m=m1,m2; do n=n1,n2
        f(l1:l2,m,n,i) = f(l1:l2,m,n,i) + &
            fact*sin(kx*x(l1:l2)+ky*y(m)+kz*z(n)+phase)
      enddo; enddo
!
    endsubroutine sinwave_phase
!***********************************************************************
    subroutine coswave_phase(f,i,ampl,kx,ky,kz,phase,lnorm_kk)
!
!  Cosine wave (as initial condition)
!
!  13-jun-06/anders: adapted from sinwave-phase.
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl, kx, ky, kz, phase, fact, k2
      integer :: i
      logical, optional :: lnorm_kk
!
!  possibility to normalize by 1/k if the initial amplitude
!  is meant to be an amplitude for the B-field.
!
      if (present(lnorm_kk)) then
        k2=kx**2+ky**2+kz**2
        if (lnorm_kk .and. k2/=0.) then
          fact=ampl/sqrt(k2)
        else
          fact=ampl
        endif
      else
        fact=ampl
      endif
!
!  Set cos wave
!
      if (lroot) print*, 'coswave_phase: i, fact, kx, ky, kz, phase=', &
          i, fact, kx, ky, kz, phase
!
      do m=m1,m2; do n=n1,n2
        f(l1:l2,m,n,i) = f(l1:l2,m,n,i) + &
            fact*cos(kx*x(l1:l2)+ky*y(m)+kz*z(n)+phase)
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
!  30-apr-16/axel: adapted for polytropic eos
!  28-jun-19/nishant: added an option to use stratification file without
!                     the ghost cells; use lnoghost_strati=T in init_pars
!
      use EquationOfState, only: eoscalc
      use Sub, only: write_zprof
      use Cdata, only: lnoghost_strati
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mzgrid) :: lnrho0,ss0,lnTT0,acc0
      real, dimension (mz) :: lnrho_mz,ss_mz,lnTT_mz
      real :: tmp,var1,var2,var3
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
!
      case ('lnrho_ss')
!
!NS: added this switch for avoiding ghost cells from stratification file
!
       if (lnoghost_strati) then
          if (lroot) print*,'ghost cells are not needed in stratification.dat'
          do n=1,nzgrid
            read(19,*,iostat=stat) tmp,var1,var2
            if (stat==0) then
             if (ip<5) print*, 'stratification: z, var1, var2=', tmp, var1, var2
             if (ldensity) lnrho0(n)=var1
             if (lentropy) ss0(n)=var2
            else
             call fatal_error('stratification','file invalid or too big - ghost cells may have been included')
            endif
          enddo
       else
          if (lroot) print*,'ghost cells are needed in stratification.dat'
          do n=1,mzgrid
            read(19,*,iostat=stat) tmp,var1,var2
            if (stat==0) then
             if (ip<5) print*, 'stratification: z, var1, var2=', tmp, var1, var2
             if (ldensity) lnrho0(n)=var1
             if (lentropy) ss0(n)=var2
            else
             call fatal_error('stratification','file invalid or too short; ghost cells may be missing')
            endif
          enddo
       endif
!
      case ('lnrho_lnTT')
        do n=1,mzgrid
          read(19,*,iostat=stat) tmp,var1,var2
          if (stat==0) then
            if (ip<5) print*, 'stratification: z, var1, var2=', tmp, var1, var2
            if (ldensity) lnrho0(n)=var1
            if (ltemperature) lnTT0(n)=var2
            if (lentropy) then
              call eoscalc(ilnrho_lnTT,var1,var2,ss=tmp)
              ss0(n)=tmp
            endif
          else
            call fatal_error('stratification','file invalid or too short - ghost cells may be missing')
          endif
        enddo
!
      case ('lnrho_lnTT_acc')
        do n=1,mzgrid
          read(19,*,iostat=stat) tmp,var1,var2,var3
          if (stat==0) then
            if (ip<5) print*, 'stratification: z, var1, var2, var3=', tmp, var1, var2, var3
            if (ldensity) lnrho0(n)=var1
            if (ltemperature) lnTT0(n)=var2
            if (lentropy) then
              call eoscalc(ilnrho_lnTT,var1,var2,ss=tmp)
              ss0(n)=tmp
            endif
            if (lascalar) acc0(n)=var3
          else
            call fatal_error('stratification','file invalid or too short - ghost cells may be missing')
          endif
        enddo
!
      case ('lnrho')
         print*,'NS1:'   !do n=1,nzgrid
        do n=1,mzgrid
          read(19,*,iostat=stat) tmp,var1
          if (stat==0) then
            if (ip<5) print*, 'stratification: z, var1=', tmp, var1
            if (ldensity) lnrho0(n)=var1
          else
            call fatal_error('stratification','file invalid or too short - ghost cells may be missing')
          endif
!
        enddo
      endselect
!
!  Select the right region for the processor afterwards.
!
!--   select case (n)
!
!  Without ghost zones.
!
!--   case (nzgrid+1)
      if (n==(nzgrid+1)) then
        if (lentropy) then
          do n=n1,n2
            f(:,:,n,ilnrho)=lnrho0(ipz*nz+(n-nghost))
            f(:,:,n,iss)=ss0(ipz*nz+(n-nghost))
          enddo
        endif
        if (ltemperature) then
          do n=n1,n2
            f(:,:,n,ilnrho)=lnrho0(ipz*nz+(n-nghost))
            f(:,:,n,ilnTT)=lnTT0(ipz*nz+(n-nghost))
          enddo
        endif
        if (ltemperature.and.lascalar) then
          do n=n1,n2
            f(:,:,n,ilnrho)=lnrho0(ipz*nz+(n-nghost))
            f(:,:,n,ilnTT)=lnTT0(ipz*nz+(n-nghost))
            f(:,:,n,iacc)=acc0(ipz*nz+(n-nghost))
          enddo
        endif
        if (.not.lentropy.and..not.ltemperature) then
          do n=n1,n2
            f(:,:,n,ilnrho)=lnrho0(ipz*nz+(n-nghost))
          enddo
        endif
!
!  With ghost zones.
!
!--   case (mzgrid+1)
      elseif (n==(mzgrid+1)) then
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
        if (.not.lentropy.and..not.ltemperature) then
          do n=1,mz
            f(:,:,n,ilnrho)=lnrho0(ipz*nz+n)
          enddo
        endif
!
!--   case default
      else
        if (lroot) then
          print '(A,I4,A,I4,A,I4,A)','ERROR: The stratification file '// &
                'for this run is allowed to contain either',nzgrid, &
                ' lines (without ghost zones) or more than',mzgrid, &
                ' lines (with ghost zones). It does contain',n-1, &
                ' lines though.'
        endif
        call fatal_error('','')
!
!--   endselect
      endif
!
!  occupy profile arrays
!
      if (lentropy) then
        do n=1,mz
          lnrho_mz(n)=lnrho0(ipz*nz+n)
          ss_mz(n)=ss0(ipz*nz+n)
        enddo
        if (lcooling_ss_mz) call write_zprof('ss_mz',ss_mz)
      endif
      if (ltemperature) then
        do n=1,mz
          lnrho_mz(n)=lnrho0(ipz*nz+n)
          lnTT_mz(n)=lnTT0(ipz*nz+n)
        enddo
        if (lcooling_ss_mz) call write_zprof('lnTT_mz',lnTT_mz)
      endif
      if (.not.lentropy.and..not.ltemperature) then
        do n=1,mz
          lnrho_mz(n)=lnrho0(ipz*nz+n)
        enddo
      endif
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
      use EquationOfState, only: eoscalc
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mxgrid) :: lnrho0,ss0,lnTT0
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
        do n=1,mxgrid
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
        do n=1,mxgrid
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
!--   select case (n)
  !
  !  without ghost zones
  !
!--   case (nxgrid+1)
      if (n==nxgrid+1) then
        if (lentropy) then
          do n=l1,l2
            f(n,:,:,ilnrho)=lnrho0(ipx*nx+(n-nghost))
            f(n,:,:,iss)=ss0(ipx*nx+(n-nghost))
          enddo
        endif
        if (ltemperature) then
          do n=l1,l2
            f(n,:,:,ilnrho)=lnrho0(ipx*nx+(n-nghost))
            f(n,:,:,ilnTT)=lnTT0(ipx*nx+(n-nghost))
          enddo
        endif
  !
  !  with ghost zones
  !
!--   case (mxgrid+1)
      elseif (n==mxgrid+1) then
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
!--   case default
      else
        if (lroot) then
          print '(A,I4,A,I4,A,I4,A)','ERROR: The stratification file '// &
                'for this run is allowed to contain either',nxgrid, &
                ' lines (without ghost zones) or more than',mxgrid, &
                ' lines (with ghost zones). It does contain',n-1, &
                ' lines though.'
        endif
        call fatal_error('','')
!
!--   endselect
      endif
!
      close(19)
!
    endsubroutine stratification_x
!***********************************************************************
    subroutine stratification_xz(f,strati_type)
!
!  read mean stratification from "stratification_xz.dat"
!
!  09-aug-14/axel: adapted from stratification
!
      use EquationOfState, only: eoscalc
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nxgrid,nzgrid,mvar) :: slice
      logical :: exist
      integer :: stat
      character (len=labellen) :: strati_type
!
!  read mean stratification and write into array
!  if file is not found in run directory, search under trim(directory)
!
      inquire(file='stratification_xz.dat',exist=exist)
      if (exist) then
        open(19,file='stratification_xz.dat')
      else
        call fatal_error('stratification_xz','no input file')
      endif
!
!  read data
!  first the entire stratification file
!
      select case (strati_type)
      case ('5variables')
        read(19,"(8e10.3)",iostat=stat) slice
        if (stat>=0) then
          if (lroot) print*,"stratification_xz: ",slice(1,1,:)
        else
          call fatal_error('stratification_xz','error reading input file')
        endif
!
      endselect
      close(19)
!
!  select the right region for the processor
!
      do m=m1,m2
        f(l1:l2,m,n1:n2,1:mvar)=slice(ipx*nx+1:(ipx+1)*nx,ipz*nz+1:(ipz+1)*nz,:)
      enddo
      if (lroot) print*,"var file initialized with slice data"
!
    endsubroutine stratification_xz
!***********************************************************************
    subroutine planet_hc(ampl,f,eps,radius,gamma,cs20,rho0,width)
!
!  Ellipsoidal planet solution (Goldreich, Narayan, Goodman 1987)
!
!   6-jul-02/axel: coded
!  22-feb-03/axel: fixed 3-D background solution for enthalpy
!  26-Jul-03/anders: Revived from June 1 version
!
      use Mpicomm, only: mpireduce_sum, mpibcast_real, MPI_COMM_WORLD
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: hh, xi
      real, dimension (mz) :: hz
      real :: delS,ampl,sigma2,sigma,delta2,delta,eps,radius,a_ell,b_ell,c_ell
      real :: gamma,cs20,gamma_m1,eps2,radius2,width
      real :: lnrhosum_thisbox,rho0
      real :: lnrhosum_thisbox_tmp,lnrhosum_wholebox
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
      lnrhosum_thisbox_tmp = lnrhosum_thisbox
!
!  Add sum_thisbox up for all processors, deliver to root
!
      call mpireduce_sum(lnrhosum_thisbox_tmp,lnrhosum_wholebox)
      if (lroot .and. ip<14) &
        print*,'planet_hc: lnrhosum_wholebox=',lnrhosum_wholebox
!
!  Calculate <rho> and send to all processors
!
      if (lroot) rho0 = exp(-lnrhosum_wholebox/nwgrid)
      call mpibcast_real(rho0,comm=MPI_COMM_WORLD)
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
      use Mpicomm, only: mpireduce_sum, mpibcast_real, MPI_COMM_WORLD
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: hh, xi, r_ell
      real :: rbound,sigma2,sigma,delta2,delta,eps,radius
      real :: gamma,eps2,radius2,width,a_ell,b_ell,c_ell
      real :: gamma_m1,ztop,cs20,hh0
      real :: lnrhosum_thisbox,rho0
      real :: lnrhosum_thisbox_tmp,lnrhosum_wholebox
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
      lnrhosum_thisbox_tmp = lnrhosum_thisbox
!
!  Add sum_thisbox up for all processors, deliver to root
!
      call mpireduce_sum(lnrhosum_thisbox_tmp,lnrhosum_wholebox)
      if (lroot .and. ip<14) &
          print*,'planet_hc: lnrhosum_wholebox=',lnrhosum_wholebox
!
!  Calculate <rho> and send to all processors
!
      if (lroot) rho0 = exp(-lnrhosum_wholebox/nwgrid)
      call mpibcast_real(rho0,comm=MPI_COMM_WORLD)
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
          f(l1:l2,m,n,i1)=ampl*log(cosh(eps*y(m)))
        enddo; enddo
      endif
!
    endsubroutine htanh
!***********************************************************************
    subroutine vtube(ampl,f,i1,i2,radius)
!
!  Vertical flux tube (for vector potential)
!  Note: in Cartesian coords this cannot be used for periodic xy mesh
!
!   1-apr-13/axel+MR: coded for cylindrical coordinates
!  14-jul-13/axel: adapted for Cartesian coordinates
!
      integer :: i1,i2
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: rr,rr2,pp,Ap
      real :: ampl,radius
!
!  Aphi=[R^6-(R^2-r^2)^3]/6r, gives Bz=(R^2-r^2)^2
!
      if (ampl==0) then
        f(:,:,:,i1:i2)=0
        if (lroot) print*,'vtube: set variable to zero; i1,i2=',i1,i2
      else
        if (lroot) then
          print*,'vtube: r-dep tube in z-direction; ampl,radius=',ampl,radius
          if (i1==i2) then
            print*,'vtube: set scalar'
          else
            print*,'vtube: set vector'
          endif
        endif
!
!  Start with n=1 to set ghost zone (for freeze-ghost condition).
!  Br is close to zero.
!
        do n=1,n2; do m=m1,m2
!
!  check whether vector or scalar
!
          if (i1==i2) then
            f(l1:l2,m,n,i1)=0.
          elseif (i1+2==i2) then
            if (lcartesian_coords) then
              rr2=x(l1:l2)**2+y(m)**2; rr=sqrt(rr2); pp=atan2(y(m),x(l1:l2))
              Ap=ampl/(6.*rr)*(radius**6-max(radius**2-rr2,0.)**3)
              f(l1:l2,m,n,i1  )=f(l1:l2,m,n,i1  )-Ap*sin(pp)
              f(l1:l2,m,n,i1+1)=f(l1:l2,m,n,i1+1)+Ap*cos(pp)
            elseif (lcylindrical_coords) then
              Ap=ampl*rcyl_mn1/6.*(radius**6-max(radius**2-x(l1:l2)**2,0.)**3)
              f(l1:l2,m,n,i1+1)=f(l1:l2,m,n,i1+1)+Ap
            endif
          else
            if (lroot) print*,'vtube: bad value of i2=',i2
          endif
!
        enddo; enddo
      endif
!
    endsubroutine vtube
!***********************************************************************
    subroutine vtube_peri(ampl,f,i1,i2,radius)
!
!  Vertical flux tube (for vector potential) for periodic xy mesh
!
!  14-jul-13/axel: coded
!
      integer :: i1,i2
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: rr2,tmp
      real :: ampl,radius
!
!  Aphi=[R^6-(R^2-r^2)^3]/6r, gives Bz=(R^2-r^2)^2
!
      if (ampl==0) then
        f(:,:,:,i1:i2)=0
        if (lroot) print*,'vtube: set variable to zero; i1,i2=',i1,i2
      else
        if (lroot) then
          print*,'vtube: r-dep tube in z-direction; ampl,radius=',ampl,radius
          if (i1==i2) then
            print*,'vtube: set scalar'
          else
            print*,'vtube: set vector'
          endif
        endif
!
!  Start with n=1 to set ghost zone (for freeze-ghost condition).
!  Br is close to zero.
!
        do n=1,n2; do m=m1,m2
!
!  check whether vector or scalar
!
          if (i1==i2) then
            f(l1:l2,m,n,i1)=0.
          elseif (i1+2==i2) then
            rr2=x(l1:l2)**2+y(m)**2
            tmp=ampl*exp(-.5*rr2/radius**2)/radius**2
            f(l1:l2,m,n,i1  )=f(l1:l2,m,n,i1  )-tmp*y(m)
            f(l1:l2,m,n,i1+1)=f(l1:l2,m,n,i1+1)+tmp*x(l1:l2)
          else
            if (lroot) print*,'vtube: bad value of i2=',i2
          endif
!
        enddo; enddo
      endif
!
    endsubroutine vtube_peri
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
!                  in Vermersch & Brandenburg (2009, AN 330, 797–806).
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
    subroutine htube2_x(ampl,f,i1,i2,radius,epsilon_nonaxi,qtube,center1_y,center1_z,scale_y)
!
!
!  Horizontal twisted flux tube in the x-direction (for vector potential, or passive scalar)
!  y axis can be stretched to form a flattened tube
!
!  17-may-15/piyali.chatterjee: coded from  htube2
!
      integer :: i1,i2
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx) :: tmp,modulate
      real :: ampl,radius,epsilon_nonaxi,kx,qtube
      real :: center1_y,center1_z,rhorad
      real, optional :: scale_y
!
        kx=2*pi/Lx
        if (lroot) then
          print*,'htube2_x: implement x-dependent flux tube in yz-plane; i1,i2=',i1,i2
          print*,'htube2_x: radius,epsilon_nonaxi=',radius,epsilon_nonaxi
        endif
!
!  constant, when epsilon_nonaxi; otherwise modulation about zero
!
        do n=1,mz; do m=1,my
          if (epsilon_nonaxi==0) then
            modulate(:)=1.0
          else
            modulate(:)=epsilon_nonaxi*sin(kx*x(:))
          endif
!
! completely quenched "gaussian"
          if (present(scale_y))  then
            rhorad=sqrt((scale_y*(y(m)-center1_y))**2+(z(n)-center1_z)**2)
            tmp=.5*ampl*modulate*radius**2*&
                (1.-exp(-(rhorad/radius)**2))/rhorad
          else
            rhorad=sqrt((y(m)-center1_y)**2+(z(n)-center1_z)**2)
            tmp=.5*ampl*modulate*radius**2*&
                (1.-exp(-(rhorad/radius)**2))/rhorad
          endif
!
!  check whether vector or scalar
!
          if (i1==i2) then
            if (lroot.and.ip<10) print*,'htube2_x: set scalar'
            f(:,m,n,i1)=f(:,m,n,i1)+tmp
          elseif (i1+2==i2) then
            if (lroot.and.ip<10) print*,'htube2_x: set vector'
            f(:,m,n,i1  )=f(:,m,n,i1)+qtube*tmp*rhorad
            f(:,m,n,i1+1)=f(:,m,n,i1+1)-(z(n)-center1_z)*tmp/rhorad
            f(:,m,n,i1+2)=f(:,m,n,i1+2)+(y(m)-center1_y)*tmp/rhorad
         else
            if (lroot) print*,'htube2_x: bad value of i2=',i2
          endif
!
        enddo; enddo
!
    endsubroutine htube2_x
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
    subroutine hfluxlayer_y(ampl,f,i,zflayer,width,ladd_bb)
!
!  Horizontal flux layer (for vector potential)
!
!  09-apr-10/piyali: copied from hfluxlayer
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,zflayer,width
      logical, intent(in), optional :: ladd_bb
!
      if (ampl==0) then
        f(:,:,:,i:i+2)=0
        if (lroot) print*,'hfluxlayer-y: set variable to zero; i=',i
      else
        if (lroot) print*,'hfluxlayer-y: horizontal flux layer; i=',i
        if ((ip<=16).and.lroot) print*,'hfluxlayer-y: ampl,width=',ampl,width
        if (ladd_bb) then
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,i  )=f(l1:l2,m,n,i  )+ampl*tanh((z(n)-zflayer)/width)
            f(l1:l2,m,n,i+1)=f(l1:l2,m,n,i+1)
            f(l1:l2,m,n,i+2)=f(l1:l2,m,n,i+2)
          enddo; enddo
        else
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,i  )=ampl*tanh((z(n)-zflayer)/width)
            f(l1:l2,m,n,i+1)=0.0
            f(l1:l2,m,n,i+2)=0.0
          enddo; enddo
        endif
      endif
!
    endsubroutine hfluxlayer_y
!***********************************************************************
    subroutine hfluxlayer_y_theta(ampl,f,i)
!
!  Horizontal flux layer (for vector potential)
!
!  19-jul-13/axel+illa: copied from hfluxlayer_y
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl
!
      if (ampl==0) then
        f(:,:,:,i:i+2)=0
        if (lroot) print*,'hfluxlayer-y-theta: set variable to zero; i=',i
      else
        if (lroot) print*,'hfluxlayer-y-theta: horizontal flux layer; i=',i
        if ((ip<=16).and.lroot) print*,'hfluxlayer-y-theta: ampl=',ampl
        do n=n1,n2; do m=m1,m2
          f(l1:l2,m,n,i  )=-ampl*max(-z(n),0.)
          f(l1:l2,m,n,i+1)=0.0
          f(l1:l2,m,n,i+2)=0.0
        enddo; enddo
      endif
!
    endsubroutine hfluxlayer_y_theta
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
!  Half-cosine (i.e., like gaussian bump, but periodic) for Bx field
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
        if ((ip<=16).and.lroot) print*,'halfcos_x: ampl,kz=',ampl,kz
        do n=n1,n2; do m=m1,m2
          f(l1:l2,m,n,i  )=0.0
          f(l1:l2,m,n,i+1)=-ampl*sin(kz*(z(n)-zbot))
          f(l1:l2,m,n,i+2)=0.0
        enddo; enddo
      endif
!
    endsubroutine halfcos_x
!***********************************************************************
    subroutine halfcos_z(ampl,f,i)
!
!  Half-cosine (i.e., like gaussian bump, but periodic) for Bz field
!
!  19-jun-02/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,kx
!
      if (ampl==0) then
        f(:,:,:,i:i+2)=0
        if (lroot) print*,'halfcos_z: set variable to zero; i=',i
      else
        print*,'halscos_z: half cosine z-field ; i=',i
        kx=2*pi/Lx
        if ((ip<=16).and.lroot) print*,'halfcos_z: ampl,kz=',ampl,kx
        do n=n1,n2; do m=m1,m2
          f(l1:l2,m,n,i  )=0.0
          f(l1:l2,m,n,i+1)=ampl/kx*sin(kx*x(l1:l2))
          f(l1:l2,m,n,i+2)=0.0
        enddo; enddo
      endif
!
    endsubroutine halfcos_z
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
        do n=1,mz; do m=1,my
          f(:,m,n,i  )=0.0
          f(:,m,n,i+1)=-ampl*z(n)
          f(:,m,n,i+2)=0.0
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
    subroutine Gaussian_By_z(ampl,f,i,z0_gaussian,width_gaussian)
!
!  Gaussian B_y field (for vector potential)
!
!  29-mar-23/nishant+rajesh: coded
!
      use Sub, only: erfunc

      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl, z0_gaussian, width_gaussian
!
      if (ampl==0) then
        f(:,:,:,i:i+2)=0
        if (lroot) print*,'uniform_y: set variable to zero; i=',i
      else
        print*,'uniform_y: uniform y-field ; i=',i
        if ((ip<=16).and.lroot) print*,'uniform_y: ampl=',ampl
        do n=n1,n2; do m=m1,m2
          f(l1:l2,m,n,i  )=ampl*erfunc((z(n)-z0_gaussian)/width_gaussian)
          f(l1:l2,m,n,i+1)=0.0
          f(l1:l2,m,n,i+2)=0.0
        enddo; enddo
      endif
!
    endsubroutine Gaussian_By_z
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
      real, optional :: kx
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
      intent(in)    :: ampl,i1,i2
      intent(inout) :: f
!
      real, dimension (mx) :: r,p,tmp
      integer :: i,n,m
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
      do n=n1,n2; do m=1,my
        do i=i1,i2
          if (lroot.and.m==1.and.n==n1) print*,'gaunoise_vect: variable i=',i
          if (modulo(i-i1,2)==0) then
            call random_number_wrapper(r)
            call random_number_wrapper(p)
            tmp=sqrt(-2.*log(r))*sin(2.*pi*p)
          else
            tmp=sqrt(-2.*log(r))*cos(2.*pi*p)
          endif
          f(:,m,n,i)=f(:,m,n,i)+ampl(n-nghost)*tmp
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
      intent(inout) :: f
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
      real, optional :: kx,ky,kz,kxx,kyy,kzz
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
    subroutine power_randomphase(ampl,initpower,kgaussian,kpeak,cutoff, &
      f,i1,i2,lscale_tobox)
!
!  Produces k^initpower*exp(-k**2/cutoff**2) spectrum.
!  However, initpower=-3 produces a k^{-1} spectrum.
!
!  07-may-03/tarek: coded
!  08-may-08/nils: adapted to work on multiple processors
!  06-jul-08/nils+andre: Fixed problem when running on
!      mult. procs (thanks to Andre Kapelrud for finding the bug)
!  08-feb-21/jennifer: added a Gaussian (adapted from power_randomphase_hel
!      but with a shift in k)
!
      use Fourier, only: fft_xyz_parallel
!
      logical, intent(in), optional :: lscale_tobox
      logical :: lscale_tobox1
      integer :: i,i1,i2,ikx,iky,ikz,stat
      real, dimension (:,:,:), allocatable :: k2, u_re, u_im, r
      real, dimension (:,:,:), allocatable :: k2mkpeak
      real, dimension (:), allocatable :: kx, ky, kz
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,initpower,mhalf,cutoff,scale_factor
      real :: nfact=4.,kpeak,kpeak1,kpeak21,nexp1,nexp2,kgaussian,fact
!
      if (present(lscale_tobox)) then
        lscale_tobox1 = lscale_tobox
      else
        lscale_tobox1 = .true.
      endif
!
!  Allocate memory for arrays.
!
      allocate(k2(nx,ny,nz),stat=stat)
      if (stat>0) call fatal_error('powern','Could not allocate memory for k2')
      if (kgaussian /= 0.) then
         allocate(k2mkpeak(nx,ny,nz),stat=stat)
         if (stat>0) call fatal_error('powern','Could not allocate memory for k2mkpeak')
      endif
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
      allocate(kz(nzgrid),stat=stat)
      if (stat>0) call fatal_error('powern', &
          'Could not allocate memory for kz')
!
      if (ampl==0) then
        if (lroot) print*,'power_randomphase: zero, nothing added to variables i1,i2=',i1,i2
      else
!
!  calculate k^2
!
        scale_factor=1
        if (.not.lscale_tobox1) scale_factor=2*pi/Lx
        kx=cshift((/(i-(nxgrid+1)/2,i=0,nxgrid-1)/),+(nxgrid+1)/2)*scale_factor
!
        scale_factor=1
        if (.not.lscale_tobox1) scale_factor=2*pi/Ly
        ky=cshift((/(i-(nygrid+1)/2,i=0,nygrid-1)/),+(nygrid+1)/2)*scale_factor
!
        scale_factor=1
        if (.not.lscale_tobox1) scale_factor=2*pi/Lz
        kz=cshift((/(i-(nzgrid+1)/2,i=0,nzgrid-1)/),+(nzgrid+1)/2)*scale_factor
!
!  integration over shells
!
        if (lroot .AND. ip<10) &
               print*,'power_randomphase:fft done; now integrate over shells...'
!  In 2-D
        if (nz==1) then
          ikz=1
          do iky=1,ny
            do ikx=1,nx
              k2(ikx,iky,ikz)=kx(ikx+ipx*nx)**2+ky(iky+ipy*ny)**2
              if (kgaussian /= 0.) then
                k2mkpeak(ikx,iky,ikz)=(kx(ikx+ipx*nx)-kpeak)**2+(ky(iky+ipy*ny)-kpeak)**2
              endif
            enddo
          enddo
!  In 3-D
        else
          do ikz=1,nz
            do iky=1,ny
              do ikx=1,nx
                k2(ikx,iky,ikz)=kx(ikx+ipx*nx)**2+ky(iky+ipy*ny)**2+kz(ikz+ipz*nz)**2
                if (kgaussian /= 0.) then
                  k2mkpeak(ikx,iky,ikz)=(kx(ikx+ipx*nx)-kpeak)**2+(ky(iky+ipy*ny)-kpeak)**2 &
                  +(kz(ikz+ipz*nz)-kpeak)**2
                endif
              enddo
            enddo
          enddo
        endif
!
        if (lroot) k2(1,1,1) = 1.  ! Avoid division by zero
!
!  To get the shell integrated power spectrum E ~ k^n, we need u ~ k^m
!  and since (in 3-D): E(k) ~ u^2 k^2 we have n=2m+2, so m=n/2-1.
!  but in 2-D, we have: E(k) ~ u^2 k, so we have n=2m+1, so m=(n-1)/2.
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
!
! Gaussian (adapted from power_randomphase_hel)
!
          if (kgaussian /= 0.) then
            nexp1=.25*nfact*initpower
            nexp2=1./nfact
            kpeak1=1./kpeak
            kpeak21=1./kpeak**2
!
!  Multiply by kpeak1**1.5 to eliminate scaling with kpeak,
!  which comes from a kpeak^3 factor in the k^2 dk integration.
!  The 1/2 factor comes from setting uk (as opposed to uk^2).
!
            fact=(kpeak1*scale_factor)**1.5
            fact=fact*kgaussian**(-.5*(initpower+1.))
            r=fact*((k2*kpeak21)**mhalf)/(1.+(k2*kpeak21)**nexp1)**nexp2
!
!  for kpeak /= 0, we have a gaussian centered around kpeak.
!  There are two options, neither of them are currently well normalized.
!
            if (kpeak /= 0.) then
               !r=r*exp(-.5*(((sqrt(k2)-kpeak)/kgaussian)**2))
               r=r*exp(-.5*(((k2-kpeak**2)/kgaussian**2)**2))
            else
               r=r*exp(-.25*(k2/kgaussian**2.-1.))
            endif
            u_re = r*u_re
            u_im = r*u_im
          endif
!
          ! back to real space
          call fft_xyz_parallel(u_re,u_im,linv=.true.)
          f(l1:l2,m1:m2,n1:n2,i)=f(l1:l2,m1:m2,n1:n2,i)+u_re
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
      if (allocated(k2mkpeak)) deallocate(k2mkpeak)
      if (allocated(u_re)) deallocate(u_re)
      if (allocated(u_im)) deallocate(u_im)
      if (allocated(r))  deallocate(r)
      if (allocated(kx)) deallocate(kx)
      if (allocated(ky)) deallocate(ky)
      if (allocated(kz)) deallocate(kz)
!
    endsubroutine power_randomphase
!***********************************************************************
    subroutine power_randomphase_hel(ampl,initpower,initpower2, &
      cutoff,ncutoff,kpeak,f,i1,i2,relhel,kgaussian, &
      lskip_projection,lvectorpotential,lscale_tobox, lsquash, &
      k1hel, k2hel,lremain_in_fourier,lpower_profile_file,qexp, &
      lno_noise,nfact0,lfactors0,compk0,llogbranch0,initpower_med0, &
      kpeak_log0,kbreak0,ldouble0,nfactd0,qirro,lsqrt_qirro,time, &
      cs,lreinit,ltime_old,ltime_new,lrho_nonuni,ilnr,l2d, &
      lnot_amp, lrandom_ampl, lfixed_phase)
!
!  Produces helical (q**n * (1+q)**(N-n))*exp(-k**l/cutoff**l) spectrum
!  when kgaussian=0, where q=k/kpeak, n=initpower, N=initpower2,
!  and l=2*ncutoff.
!  The relative helicity is relhel.
!
!  initpower=0 -> k^2  (if used for vector potential)
!           =2 -> k^4
!          =-3 -> k^{-1}
!          =-3.67 k^{-5/3}
!          =-5 -> k^{-3}
!
!  Note that the rms value is unchanged if kpeak is changed.
!  Thus, to reproduce an unchanged subinertial range after increasing kpeak,
!  for example, one has to *increase* ampl by a factor 1/kpeak**2.5,
!  e.g., from 3e-4 to 17e-4.
!
!  08-sep-14/axel: adapted from power_randomphase
!  17-sep-18/axel: added optional wavenumber interval for helical field.
!  28-jan-19/axel: special treatment of 2-D case (nz==1)
!  14-aug-20/axel: adapted for fft_xyz_parallel
!  7-aug-22/alberto: added optional logarithmic branch and double
!                    broken power law
!
      use Fourier, only: fft_xyz_parallel
      use General, only: loptest, roptest
!
      logical, intent(in), optional :: lscale_tobox, lsquash, lremain_in_fourier, ltime_old
      logical, intent(in), optional :: ltime_new, lrho_nonuni, lnot_amp, lrandom_ampl, lfixed_phase
      logical, intent(in), optional :: lpower_profile_file, lno_noise, lfactors0
      logical, intent(in), optional :: llogbranch0,ldouble0, lreinit, l2d, lsqrt_qirro
      logical :: lvectorpotential, lscale_tobox1, lsquash1, lremain_in_fourier1, lno_noise1
      logical :: lskip_projection,lfactors,llogbranch,ldouble, ltime, ltime_old1
      logical :: ltime_new1, lrho_nonuni1, l2d1, lsqrt_qirro1, lnot_amp1, lrandom_ampl1
      logical :: lfixed_phase1
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i, i1, i2, ikx, iky, ikz, stat, ik, nk, ilnr1
      integer, intent(in), optional :: ilnr
      real, intent(in), optional :: k1hel, k2hel, qexp, nfact0, compk0
      real, intent(in), optional :: initpower_med0, kpeak_log0, kbreak0
      real, intent(in), optional :: nfactd0, qirro, time, cs
      real, dimension (:,:,:,:), allocatable :: u_re, u_im
      real, dimension(3) :: v_re, v_im
      real, dimension (:,:,:), allocatable :: k2, r
      real, dimension (:), allocatable :: kx, ky, kz
      real, dimension (:), allocatable :: kk, lgkk, power_factor, lgff
      real :: ampl,initpower,initpower2,mhalf,cutoff,kpeak,scale_factor,relhel
      real :: nfact, kpeak1, kpeak21, nexp1,nexp2,ncutoff,kgaussian,fact
      real :: lgk0, dlgk, lgf, lgk, lgf2, lgf1, lgk2, lgk1, D1, D2, D3, compk
      real :: kpeak_log, kbreak, kbreak1, kbreak2, kbreak21, initpower_med, initpower_log
      real :: nfactd,nexp3,nexp4
      real :: qexp1, qexp11, qirro1, p, p2, time1, cs1, om, ctime, stime
!
      if (ampl==0.) then
        if (lroot) print*,'power_randomphase: set variable to zero; i1,i2=',i1,i2
        f(:,:,:,i1:i2) = 0.
        return
      endif
!
!  By default, don't scale wavenumbers to the box size.
!
      lscale_tobox1 = loptest(lscale_tobox,.true.)
!
      lsquash1 = loptest(lsquash,.true.)
!
!  Check whether or not we want to remain in Fourier space
!
      lremain_in_fourier1 = loptest(lremain_in_fourier)
!
!  Check whether or not we want ltime_old
!
      ltime_old1 = loptest(ltime_old)
!
!  Check whether or not we want ltime_new
!
      ltime_new1 = loptest(ltime_new)
!
!  Check whether or not we want lrho_nonuni
!
      lrho_nonuni1 = loptest(lrho_nonuni)
      if (lrho_nonuni1) then
        ilnr1 = ioptest(ilnr)
        if (ilnr1 <= 0) call fatal_error('power_randomphase_hel','must provide ilnr')
      endif
!
!  Check whether we want random amplitudes or not
!
      lrandom_ampl1 = loptest(lrandom_ampl)
!
!  Check whether we want no_noise or not
!
      lno_noise1 = loptest(lno_noise)
!
!  Check whether we want lfixed_phase or not
!
      lfixed_phase1 = loptest(lfixed_phase)
!
!  Check whether we want l2d or not
!
      l2d1 = loptest(l2d)
!
!  qexp for q-exponential, qexp=1 by default
!
     qexp1 = roptest(qexp, rdef=1.)
!
!  qirro, is the vortical contribution, qirro=1 for fully irrotational.
!
     qirro1 = roptest(qirro)
!
! added option to use sqrt(qirro) to have correct normalization for trace
! of two-point correlation (added by antonino, madeline, alberto)
!
     lsqrt_qirro1 = loptest(lsqrt_qirro)
     !print*,'value of qirro1',lsqrt_qirro1
!
!  time
!
     time1 = roptest(time)
     ltime = time1/=0.
!
!  cs
!
     cs1 = roptest(cs,1.)
!
!  alberto: added option to compensate spectral shape by a power of k
!
     compk = roptest(compk0)
!
!  alberto: added option to use different values of nfact
!  Here, nfact is the exponent on k/k0, with an 1/nfact outside [1+(k/k0)^n]^(1/nfact).
!  By default, we use a large nfact=4 to have a sharp transition.
!
     nfact = roptest(nfact0,4.)
     lfactors = loptest(lfactors0)
     lnot_amp1 = loptest(lnot_amp)
!
!  alberto: added option to include additional logarithmic branch
!
      llogbranch = loptest(llogbranch0)
!
!  alberto: added option to use double smoothed broken power laws
!           logbranch is set to False if ldouble is used
!
      ldouble = loptest(ldouble0)
      if (ldouble) llogbranch = .false.
!
!  Check if the parameters of the logarithmic branch or the
!  double broken power law are given
!
      if (llogbranch.or.ldouble) then

        initpower_med = roptest(initpower_med0,1.)
        kbreak = roptest(kbreak0,.5)
        kbreak1=1./kbreak
        kbreak2=kbreak**2
        kbreak21=kbreak1**2

        if (llogbranch) kpeak_log = roptest(kpeak_log0,1.)
        if (ldouble) nfactd = roptest(nfactd0,4.)

      endif
!
!  Debug output
!
      if (lroot.and.ip<10) then
        print*,'i1,i2=',i1,i2
        print*,'ampl,initpower,initpower2=',ampl,initpower,initpower2
        print*,'cutoff,ncutoff,kpeak,i1,i2,relhel,kgaussian=',cutoff,ncutoff,kpeak,i1,i2,relhel,kgaussian
        !print*,'lskip_projection,lvectorpotential,lscale_tobox1=',lskip_projection,lvectorpotential,lscale_tobox1
        print*,'lskip_projection,lscale_tobox1,lsquash1=',lskip_projection,lscale_tobox1,lsquash1
        !print*,'k1hel,k2hel,lremain_in_fourier,lpower_profile_file,qexp=',k1hel,k2hel,lremain_in_fourier,lpower_profile_file,qexp
        print*,'lremain_in_fourier1',lremain_in_fourier1
  !     print*,'lno_noise,nfact,lfactors=',lno_noise,nfact,lfactors
  !     print*,'compk,llogbranch,initpower_med=',compk,llogbranch,initpower_med
  !     print*,'kpeak_log,kbreak,ldouble,nfactd,qirro1=',kpeak_log,kbreak,ldouble,nfactd,qirro1
      endif
!
!  Allocate memory for arrays.
!
      allocate(k2(nx,ny,nz),stat=stat)
      if (stat>0) call fatal_error('power_randomphase_hel','Could not allocate k2')
      allocate(r(nx,ny,nz),stat=stat)
      if (stat>0) call fatal_error('power_randomphase_hel','Could not allocate r')
!
      if (i2==i1) then
        allocate(u_re(nx,ny,nz,1),stat=stat)
        if (stat>0) call fatal_error('power_randomphase_hel','Could not allocate u_re')
        allocate(u_im(nx,ny,nz,1),stat=stat)
        if (stat>0) call fatal_error('power_randomphase_hel','Could not allocate u_im')
      else
        allocate(u_re(nx,ny,nz,3),stat=stat)
        if (stat>0) call fatal_error('power_randomphase_hel','Could not allocate u_re')
        allocate(u_im(nx,ny,nz,3),stat=stat)
        if (stat>0) call fatal_error('power_randomphase_hel','Could not allocate u_im')
      endif
!
      allocate(kx(nxgrid),stat=stat)
      if (stat>0) call fatal_error('power_randomphase_hel','Could not allocate kx')
      allocate(ky(nygrid),stat=stat)
      if (stat>0) call fatal_error('power_randomphase_hel','Could not allocate ky')
      allocate(kz(nzgrid),stat=stat)
      if (stat>0) call fatal_error('power_randomphase_hel','Could not allocate kz')
!
!  calculate k^2. If lscale_tobox1=T, then k is calculated as if the
!  box is 2pi in all three directions.
!  If lsquash_aa=T, then the box size in the y and z directions is the
!  same as in the x direction. It is called that way, because, if the
!  aspect ratio is is different from unity, it makes the field squashed.
!
      scale_factor=1
      if (.not.lscale_tobox1) scale_factor=2*pi/Lx
      kx=cshift((/(i-nxgrid/2,i=0,nxgrid-1)/),nxgrid/2)*scale_factor
      if (lroot.and.ip<10) print*,'AXEL: kx=',kx
!
      if (.not. lsquash1) then
        scale_factor=1
        if (.not.lscale_tobox1) scale_factor=2*pi/Ly
      endif
      ky=cshift((/(i-nygrid/2,i=0,nygrid-1)/),nygrid/2)*scale_factor
      if (lroot.and.ip<10) print*,'AXEL: ky=',ky
!
      if (.not. lsquash1) then
        scale_factor=1
        if (.not.lscale_tobox1) scale_factor=2*pi/Lz
      endif
      kz=cshift((/(i-nzgrid/2,i=0,nzgrid-1)/),nzgrid/2)*scale_factor
      if (lroot.and.ip<10) print*,'AXEL: kz=',kz
!
!  Generate flat spectrum with random phase (between -pi and pi).
!  When lrandom_ampl=T, we multiply by a sqrt(-2*alog(k2)) factor, as in the
!  gaunoise routine; the variable k2 is borrowed before it is used for k^2.
!  For qexp /= 1, we use the q-logarithm, the inverse of the q-exponential,
!  alog_q(x) = [x^(1-q)-1]/(1-q), and define qexp11 = 1-q for clearer notation.
!
      if (lno_noise1) then
        u_re=ampl
        u_im=0.
      else
        do i=1,i2-i1+1
          call random_number_wrapper(r)
          if (lrandom_ampl1) then
            call random_number_wrapper(k2)
            if (qexp1==1.) then
              u_im(:,:,:,i)=ampl*sqrt(-2*alog(k2))
            else
              qexp11=1.-qexp1
              u_im(:,:,:,i)=ampl*sqrt(-2*(k2**qexp11-1.)/qexp11)
            endif
!           if (lfixed_phase1) then
!             u_re(:,:,:,i)=u_im(:,:,:,i)
!           else
            if (.not. lfixed_phase1) then
              u_re(:,:,:,i)=u_im(:,:,:,i)*cos(pi*(2*r-1))
              u_im(:,:,:,i)=u_im(:,:,:,i)*sin(pi*(2*r-1))
            endif
          else
            u_re(:,:,:,i)=ampl*cos(pi*(2*r-1))
            u_im(:,:,:,i)=ampl*sin(pi*(2*r-1))
          endif
        enddo
      endif
!
!  Set k^2 array.
!
!  Note: for multiple processors, there may still be a problem in 1-D
!
!  In 1-D
      if (nzgrid==1.and.nygrid==1) then
        ikz=1
        iky=1
        do ikx=1,nx
          k2(ikx,iky,ikz)=kx(ikx+ipx*nx)**2
        enddo
!  In 2-D
      elseif (nzgrid==1) then
        ikz=1
        do iky=1,ny
          do ikx=1,nx
            k2(ikx,iky,ikz)=kx(ikx+ipx*nx)**2+ky(iky+ipy*ny)**2
          enddo
        enddo
!  In 3-D
      else
        do ikz=1,nz
          do iky=1,ny
              do ikx=1,nx
                k2(ikx,iky,ikz)=kx(ikx+ipx*nx)**2+ky(iky+ipy*ny)**2+kz(ikz+ipz*nz)**2
              enddo
          enddo
        enddo
      endif
      if (lroot) k2(1,1,1) = 1.  ! Avoid division by zero
!
!  Do phase correction
!
      if (lfixed_phase1) then
        do i=1,i2-i1+1
          u_re(:,:,:,i)=u_im(:,:,:,i)*cos(-sqrt(k2)*t)
          u_im(:,:,:,i)=u_im(:,:,:,i)*sin(-sqrt(k2)*t)
        enddo
      endif
!
!  To get the shell integrated power spectrum E ~ k^n, we need u ~ k^m
!  and since E(k) ~ u^2 k^2 we have n=2m+2, so m=n/2-1 in 3-D.
!  In general, E(k) ~ u^2 k^(D-1), so we have n=2m+D-1, so m=(n-D+1)/2
!  Further, since we operate on k^2, we need m/2 (called mhalf below)
!
      mhalf=.25*(initpower-dimensionality+1)
!
!  generate all 3 velocity components separately
!  generate k^n spectrum with random phase (between -pi and pi)
!
      nexp1=.25*nfact*(initpower-initpower2)
      !
      !  alberto: adapt input power law exponents to parameters that appear in the
      !           broken power law when a logarithmic branch is included or when
      !           a double broken power law is included
      !
      if (llogbranch.or.ldouble) then
        nexp1=.25*nfact*(initpower_med-initpower2)
        initpower_log=.5*(initpower-initpower_med)
      endif
      if (ldouble) then
        mhalf=.25*(initpower_med-dimensionality+1)
        nexp3=1./nfactd
        nexp4=.5*nfactd*initpower_log
      endif
      nexp2=1./nfact
      kpeak1=1./kpeak
      kpeak21=1./kpeak**2
!
!  alberto: added option to compensate amplitude and peak using D1 and D2
!           such that the maximum of the spectrum is located at kpeak
!           or at the plateau (if present)
!
      D1=0.
      D2=1.
      D3=1.
      fact=1.
      if (lfactors) then
        if (llogbranch) then
          !  alberto: changing sign of nfact allows to use spectral shapes with
          !           initpower - initpower2 > 0
          if (initpower_med<initpower2) then
            nexp1=-nexp1
            nexp2=-nexp2
            nfact=-nfact
          endif
          if (initpower2/=0) then
            if (initpower_med/=0) then
              if (initpower2*initpower_med<0) then
                D1=-initpower_med/initpower2
              else
                D1=1.
              endif
              D2=D1
            endif
            fact=fact*log(1+kpeak_log*kpeak1)**(-initpower_log)
          else
            if (initpower_med==0) fact=fact*(1+D2)**nexp2
            fact=fact*(kpeak_log*kpeak1)**(-initpower_log)
          endif
          if ((initpower>=0).and.(initpower_med<0)) then
            fact=(1+D1)**(-nexp2)*(kbreak*kpeak1)**(-.5*initpower)
            fact=fact*(1+D2*(kbreak*kpeak1)**(.5*nfact*(initpower_med-initpower2)))**(-nexp2)
            fact=fact*log(1+kpeak_log*kbreak1)**(-initpower_log)
            if (initpower_med==initpower2) fact=fact*(1+D2)**(2*nexp2)
          endif
        elseif (ldouble) then
          !  alberto: changing sign of nfact allows to use spectral shapes with
          !           initpower - initpower2 > 0
          if (initpower_med<initpower2) then
            nexp1=-nexp1
            nexp2=-nexp2
          endif
          if (initpower<initpower_med) then
            nexp3=-nexp3
            nexp4=-nexp4
          endif
          if (initpower2/=0) then
            if (initpower_med/=0) then
              if (initpower_med*initpower2<0) then
                D1=-initpower_med/initpower2
              else
                D1=1.
              endif
              D2=D1
            endif
            fact=fact*(1+D3*(kpeak*kbreak1)**(-2*nexp4))**nexp3
          else
            if (initpower_med==0) then
              fact=fact*(1+D2)**nexp2
              if (initpower==0) then
                fact=fact*(1+D1)**nexp2*(1+D3)**nexp3
              endif
            endif
          endif
          if ((initpower>=0).and.(initpower_med<0)) then
            fact=1.
            if ((initpower_med*initpower<0).and.(initpower2/=0)) then
              D1=-initpower_med/initpower2
              D2=D1
            endif
            if (initpower2*initpower_med>0) D2=1.
            if (initpower/=0) then
              D3=-initpower_med/initpower
              fact=fact*(1 + D3)**nexp3
            endif
            fact=fact*(1+D1)**(-nexp2)*(kbreak*kpeak1)**(-.5*initpower_med)
            fact=fact*(1+D2*(kbreak*kpeak1)**(2*nexp1))**nexp2
          endif
        else  ! (ldouble)
          !  alberto: changing sign of nfact allows to use spectral shapes with
          !           initpower - initpower2 > 0
          if ((initpower-initpower2)<0) then
            nexp1=-nexp1
            nexp2=-nexp2
          endif
          if ((initpower*initpower2)<0) then
            D1=-initpower/initpower2
            D2=D1
          elseif ((initpower*initpower2)==0) then
            if ((initpower==0).and.(initpower2)==0) then
              fact=fact*(1+D2)**nexp2
            elseif (initpower2==0) then
              fact=fact*D2**nexp2
            endif
          else
            fact=fact*(1+D2)**nexp2
          endif
          !D1=abs(initpower/initpower2)
          !D2=-1/initpower2
          !if (initpower/=0) then
          !  D2=D2*abs(initpower)
          !endif
        endif
      endif  ! (lfactors)
!
      fact=fact*(1+D1)**nexp2
!
!  Multiply by kpeak1**1.5 to eliminate scaling with kpeak,
!  which comes from a kpeak^3 factor in the k^2 dk integration.
!  The 1/2 factor comes from setting uk (as opposed to uk^2).
!
!  alberto: if lfactor is chosen, the amplitude of the spectrum is
!           independent of kpeak instead of the integrated energy,
!           and it takes the value ampl at kpeak (tested for hydro
!           fields, for magnetic one can choose compk_aa=-.5 to obtain
!           a magnetic field spectrum with value ampl at kpeak)
!
      if (lfactors) then
        !print*,'TEST Lx, kpeak',Lx,kpeak
        if (lnot_amp1) then
          fact=fact*(2*pi/Lx)**0.5
        else
          fact=fact*kpeak1*scale_factor**1.5/(2*pi*ampl)**0.5
          !
          ! alberto: compensate for contribution due to helicity
          ! taking into account generic qirro (only correct when
          ! lsqrt_qirro is true)
          !
          fact=fact/(1 + relhel**2*(1 - qirro1))**0.5
          ! compensate when compk is non-zero
          fact=fact*kpeak21**compk
          !fact=fact*(kpeak1*scale_factor)
        endif
      else
        fact=fact*(kpeak1*scale_factor)**1.5
      endif
!
!  If lvectorpotential, we given slopes are meant to be for the magnetic
!  vector potential, so we use (initpower+3.) instead of (initpower+1.).
!
      if (lvectorpotential) then
        fact=fact*kpeak1
        if (kgaussian /= 0.) fact=fact*kgaussian**(-.5*(initpower+3.))
      else
        if (kgaussian /= 0.) fact=fact*kgaussian**(-.5*(initpower+1.))
      endif
      r=fact*((k2*kpeak21)**mhalf)/(1.+D2*(k2*kpeak21)**nexp1)**nexp2
      if (lroot.and.ip<10) print*,'kpeak,mhalf,nexp1,nexp2=',kpeak,mhalf,nexp1,nexp2
!
!  Examples: for initpower=4., initpower2=-2., get mhalf,nexp1,nexp2 = 0.75, 6.0, 0.25
!  while for Jani Dahl setup, we put: initpower=3., initpower2=-3., nfact_uu=.6666667
!  and get mhalf, nexp1, nexp2 = 0.5, 1.0, 1.5.
!
!  alberto: added possibility to multiply spectrum by a power of k (useful for GW for example)
!           the final spectrum will be compensated by k^(4*compk)
!
      r=r*k2**compk
!
!  alberto: added model for double smoothed broken power law
!
      if (ldouble) r=r/(1.+D3*(k2*kbreak21)**(-nexp4))**nexp3
!
!  alberto: multiply the broken power law by the logarithmic branch using
!           the condition k < kbreak
!
!  Produces a 'double' broken power law for GW spectrum using a logarithmic branch based
!  on Roper Pol, Caprini, Neronov, Semikoz, 2022 (arXiv:2201.05630) by multiplying
!  TauGW to the spectrum, where TauGW is:
!  ln (1 + kpeak_log/kbreak)**n when k <= kbreak and
!  ln (1 + kpeak_log/k)**n when k > kbreak
!
!  The resulting spectrum asymptotically has slopes k^initpower at k < kbreak,
!  k^initpower_med at kbreak < k < kpeak, k^(-initpower2) at k > kpeak
!
      if (llogbranch) then
        where (k2<kbreak2)
          r = r*log(1+kpeak_log*kbreak1)**initpower_log
        elsewhere
          r = r*log(1+kpeak_log/sqrt(k2))**initpower_log
        endwhere
      endif
!
!  cutoff (changed to hyperviscous cutoff filter).
!  The 1/2 factor is needed to account for the fact that the
!  spectrum (velocity squared) is proportional to exp(-k^2/cutoff^2).
!  For ncutoff/=1, one has exp(-(k^2/cutoff^2)^ncutoff).
!
      if (cutoff /= 0.) r=r*exp(-.5*(k2/cutoff**2.)**ncutoff)
!
!  apply Gaussian on top of everything
!
      if (kgaussian /= 0.) r=r*exp(-.25*(k2/kgaussian**2.-1.))
!
!  apply additional profile read from file;
!  first check whether or not we want to read from file
!
      if (loptest(lpower_profile_file)) then
        open(9,file='power_profile.dat',status='old')
        read(9,*) nk,lgk0,dlgk
        if (lroot) print*,'power_randomphase_hel: nk,lgk0,dlgk=',nk,lgk0,dlgk
        if (allocated(kk)) deallocate(kk,power_factor,lgkk,lgff)
        allocate(kk(nk),power_factor(nk),lgkk(nk),lgff(nk))
        do ik=1,nk
          read(9,*) kk(ik),power_factor(ik)
        enddo
        close(9)
        lgkk=alog10(kk)
        lgff=alog10(power_factor)
!
!  Go through each mesh point and interpolate logarithmically to the
!  correct scaling factor for the power (energy) spectrum.
!
        do ikz=1,nz
          do iky=1,ny
            do ikx=1,nx
              lgk=alog10(sqrt(k2(ikx,iky,ikz)))
              ik=int((lgk-lgk0)/dlgk)+1
              if (ik<1.or.ik>nk) then
                print*,'ikz,iky,ikx,lgk,ik,k2=',ikz,iky,ikx,lgk,ik,k2(ikx,iky,ikz)
                call fatal_error('power_randomphase_hel','ik<1.or.ik>nk')
              endif
              lgk1=lgkk(ik)
              lgk2=lgkk(ik+1)
              lgf1=lgff(ik)
              lgf2=lgff(ik+1)
              lgf=lgf1+(lgk-lgk1)*(lgf2-lgf1)/(lgk2-lgk1)
              r(ikx,iky,ikz)=r(ikx,iky,ikz)*10**lgf
              if (ip<14) print*,'iproc,lgk1,lgk,lgk2=',iproc,lgk1,lgk,lgk2
            enddo
          enddo
        enddo
      endif
!
!  Scale with r: allow for special case with *scalars* here.
!  Note that this adds to anything that was set previously.
!
      if (i2==i1) then
        u_re(:,:,:,1)=r*u_re(:,:,:,1)
        u_im(:,:,:,1)=r*u_im(:,:,:,1)
        if (lremain_in_fourier1) then
          f(l1:l2,m1:m2,n1:n2,i1  )=f(l1:l2,m1:m2,n1:n2,i1  )+u_re(:,:,:,1)
          f(l1:l2,m1:m2,n1:n2,i1+1)=f(l1:l2,m1:m2,n1:n2,i1+1)+u_im(:,:,:,1)
        else
          call fft_xyz_parallel(u_re(:,:,:,1),u_im(:,:,:,1),linv=.true.)
          f(l1:l2,m1:m2,n1:n2,i1)=f(l1:l2,m1:m2,n1:n2,i1)+u_re(:,:,:,1)
        endif
      else
        do i=1,3
          u_re(:,:,:,i)=r*u_re(:,:,:,i)
          u_im(:,:,:,i)=r*u_im(:,:,:,i)
        enddo
!
!  Apply projection operator
!  Use r=1/k^2 for normalization in khat_i * khat_j = ki*kj/k2.
!  Remember that for the return transform, data have to be
!  arranged in the order (kz,kx,ky).
!  To allow also for the possibility of longitudinal initial fields, we write
!  (1-q)*(delij-kikj) + q*kikj = (1-q)*delij - (1-2*q)*kikj.
!
!
!  Here, u_re, u_im is what was v_re, v_im before!
!
        if (.not.lskip_projection) then
!
!  Allow for possibility of irrotational contributions of fraction q,
!  so the vortical fraction is (1-q) == p.
!
          if (lsqrt_qirro1) then
            p=sqrt(1.-qirro1)
            p2=sqrt(1.-qirro1)-sqrt(2.*qirro1)
          else
            p=1.-qirro1
            p2=1.-2.*qirro1
          endif
!
!  In 2-D
!
          if (nz==1) then
            ikz=1
            do iky=1,ny
              do ikx=1,nx
!
!  Real part of (ux, uy, uz) -> vx, vy, vz
!  (kk.uu)/k2, ==> vi = ui - ki kj uj, but now we write:
!  (kk.uu)/k2, ==> vi = (1-q)*ui - (1-2q) ki kj uj
!
                !r(ikx,iky,ikz)=(1.-2.*qirro1) * (kx(ikx+ipx*nx)*u_re(ikx,iky,ikz,1) &
                r(ikx,iky,ikz)=p2 * (kx(ikx+ipx*nx)*u_re(ikx,iky,ikz,1) &
                                                +ky(iky+ipy*ny)*u_re(ikx,iky,ikz,2))/k2(ikx,iky,ikz)
                u_re(ikx,iky,ikz,1)=p*u_re(ikx,iky,ikz,1)-kx(ikx+ipx*nx)*r(ikx,iky,ikz)
                u_re(ikx,iky,ikz,2)=p*u_re(ikx,iky,ikz,2)-ky(iky+ipy*ny)*r(ikx,iky,ikz)
                u_re(ikx,iky,ikz,3)=p*u_re(ikx,iky,ikz,3)
!
!  Imaginary part of (ux, uy, uz) -> vx, vy, vz
!  (kk.uu)/k2, vi = ui - ki kj uj
!
                !r(ikx,iky,ikz)=(1.-2.*qirro1) * (kx(ikx+ipx*nx)*u_im(ikx,iky,ikz,1) &
                r(ikx,iky,ikz)=p2 * (kx(ikx+ipx*nx)*u_im(ikx,iky,ikz,1) &
                                                +ky(iky+ipy*ny)*u_im(ikx,iky,ikz,2))/k2(ikx,iky,ikz)
                u_im(ikx,iky,ikz,1)=p*u_im(ikx,iky,ikz,1)-kx(ikx+ipx*nx)*r(ikx,iky,ikz)
                u_im(ikx,iky,ikz,2)=p*u_im(ikx,iky,ikz,2)-ky(iky+ipy*ny)*r(ikx,iky,ikz)
                u_im(ikx,iky,ikz,3)=p*u_im(ikx,iky,ikz,3)
              enddo
            enddo
!  In 3-D
          else
            do ikz=1,nz
              do iky=1,ny
                do ikx=1,nx
!
!  Real part of (ux, uy, uz) -> vx, vy, vz
!  (kk.uu)/k2, vi = ui - ki kj uj
!
                  !r(ikx,iky,ikz)=(1.-2.*qirro1)* &
                  r(ikx,iky,ikz)=p2 * &
                      (kx(ikx+ipx*nx)*u_re(ikx,iky,ikz,1) &
                      +ky(iky+ipy*ny)*u_re(ikx,iky,ikz,2) &
                      +kz(ikz+ipz*nz)*u_re(ikx,iky,ikz,3))/k2(ikx,iky,ikz)
!
!  Possibility of a kinematic time dependence.
!  Commented out for now. This does not seem correct here and would overwrite r
!  when the keyword is set. There is another ltime implementation below which
!  seems to be the correct one.
!  AB: Out-commented part to be deleted.
!
!                    if (ltime) then
!                      om=cs1*sqrt(k2(ikx,iky,ikz))
!                      r(ikx,iky,ikz)=r(ikx,iky,ikz)*sin(om*time1)
!                    endif
!
                  u_re(ikx,iky,ikz,1)=p*u_re(ikx,iky,ikz,1)-kx(ikx+ipx*nx)*r(ikx,iky,ikz)
                  u_re(ikx,iky,ikz,2)=p*u_re(ikx,iky,ikz,2)-ky(iky+ipy*ny)*r(ikx,iky,ikz)
                  u_re(ikx,iky,ikz,3)=p*u_re(ikx,iky,ikz,3)-kz(ikz+ipz*nz)*r(ikx,iky,ikz)
!
!  Imaginary part of (ux, uy, uz) -> vx, vy, vz
!  (kk.uu)/k2, vi = ui - ki kj uj
!
                  !r(ikx,iky,ikz)=(1.-2.*qirro1)* &
                  r(ikx,iky,ikz)=p2 * &
                      (kx(ikx+ipx*nx)*u_im(ikx,iky,ikz,1) &
                      +ky(iky+ipy*ny)*u_im(ikx,iky,ikz,2) &
                      +kz(ikz+ipz*nz)*u_im(ikx,iky,ikz,3))/k2(ikx,iky,ikz)
                  u_im(ikx,iky,ikz,1)=p*u_im(ikx,iky,ikz,1)-kx(ikx+ipx*nx)*r(ikx,iky,ikz)
                  u_im(ikx,iky,ikz,2)=p*u_im(ikx,iky,ikz,2)-ky(iky+ipy*ny)*r(ikx,iky,ikz)
                  u_im(ikx,iky,ikz,3)=p*u_im(ikx,iky,ikz,3)-kz(ikz+ipz*nz)*r(ikx,iky,ikz)
                enddo
              enddo
            enddo
          endif  ! 3D
        endif  !(.not.lskip_projection)
!
!  Make it helical, i.e., multiply by delta_ij + epsilon_ijk ikhat_k*sigma.
!  Use r=sigma/k for normalization of sigma*khat_i = sigma*ki/sqrt(k2).
!  Put r(k=0)=0, but this is only true for the root processor.
!
!  From now on, k2 := sqrt(k2)!
!
        k2=sqrt(k2)
        r=relhel/k2
        if (lroot) r(1,1,1)=0.
!
!  put sigma=0 outside [r1hel,r2hel]
!
        if (present(k1hel) .and. present(k2hel)) then
          if (k1hel>0. .and. k2hel<max_real) then
            where (k2<k1hel .or. k2>k2hel) r = 0.
          endif
        endif
!
!  In 2-D
        if (nz==1) then
          do iky=1,ny
            do ikx=1,nx
              ikz=1
!
!  (vx, vy, vz) -> ux, but put ux=uy=0 then l2d=T. This option only
!  makes sense for the nagnetic vector potential, but the same effect
!  can there be achieved by setting lset_AxAy_zero=T in magnetic.
!  In hydro, we would instead put lset_uz_zero=T.
!
              v_re = u_re(ikx,iky,ikz,:); v_im = u_im(ikx,iky,ikz,:)
              if (l2d1) then
                u_re(ikx,iky,ikz,1:2)=0.
                u_im(ikx,iky,ikz,1:2)=0.
              else
!
                u_re(ikx,iky,ikz,1)=v_re(1) - ky(iky+ipy*ny)*v_im(3)*r(ikx,iky,ikz)
                u_im(ikx,iky,ikz,1)=v_im(1) + ky(iky+ipy*ny)*v_re(3)*r(ikx,iky,ikz)
!
!  (vx, vy, vz) -> uy
!
                u_re(ikx,iky,ikz,2)=v_re(2) + kx(ikx+ipx*nx)*v_im(3)*r(ikx,iky,ikz)
                u_im(ikx,iky,ikz,2)=v_im(2) - kx(ikx+ipx*nx)*v_re(3)*r(ikx,iky,ikz)
              endif
!
!  (vx, vy, vz) -> uz
!
              u_re(ikx,iky,ikz,3)=v_re(3) + ky(iky+ipy*ny)*v_im(1)*r(ikx,iky,ikz) &
                                          - kx(ikx+ipx*nx)*v_im(2)*r(ikx,iky,ikz)
              u_im(ikx,iky,ikz,3)=v_im(3) - ky(iky+ipy*ny)*v_re(1)*r(ikx,iky,ikz) &
                                          + kx(ikx+ipx*nx)*v_re(2)*r(ikx,iky,ikz)
            enddo
          enddo
!  In 3-D
        else
          do ikz=1,nz
            do iky=1,ny
              do ikx=1,nx
!
                v_re = u_re(ikx,iky,ikz,:); v_im = u_im(ikx,iky,ikz,:)
!
!  (vx, vy, vz) -> ux
!
                u_re(ikx,iky,ikz,1)=v_re(1) + kz(ikz+ipz*nz)*v_im(2)*r(ikx,iky,ikz) &
                                            - ky(iky+ipy*ny)*v_im(3)*r(ikx,iky,ikz)
                u_im(ikx,iky,ikz,1)=v_im(1) - kz(ikz+ipz*nz)*v_re(2)*r(ikx,iky,ikz) &
                                            + ky(iky+ipy*ny)*v_re(3)*r(ikx,iky,ikz)
!
!  (vx, vy, vz) -> uy
!
                u_re(ikx,iky,ikz,2)=v_re(2) + kx(ikx+ipx*nx)*v_im(3)*r(ikx,iky,ikz) &
                                            - kz(ikz+ipz*nz)*v_im(1)*r(ikx,iky,ikz)
                u_im(ikx,iky,ikz,2)=v_im(2) - kx(ikx+ipx*nx)*v_re(3)*r(ikx,iky,ikz) &
                                            + kz(ikz+ipz*nz)*v_re(1)*r(ikx,iky,ikz)
!
!  (vx, vy, vz) -> uz
!
                u_re(ikx,iky,ikz,3)=v_re(3) + ky(iky+ipy*ny)*v_im(1)*r(ikx,iky,ikz) &
                                            - kx(ikx+ipx*nx)*v_im(2)*r(ikx,iky,ikz)
                u_im(ikx,iky,ikz,3)=v_im(3) - ky(iky+ipy*ny)*v_re(1)*r(ikx,iky,ikz) &
                                            + kx(ikx+ipx*nx)*v_re(2)*r(ikx,iky,ikz)
!
              enddo
            enddo
          enddo
        endif
!
!  Possibility of a kinematic time dependence.
!
        if (ltime) then
          if (ltime_new1.and..not.lskip_projection) &
            call fatal_error('power_randomphase_hel','must have lskip_projection=T')
          do ikz=1,nz
            do iky=1,ny
              do ikx=1,nx
                  om=cs1*k2(ikx,iky,ikz)
                  ctime=cos(om*time1); stime=sin(om*time1)
                  if (ltime_new1) then
!
!  Pretend that u_re(ikx,iky,ikz,1) and u_im(ikx,iky,ikz,1) is a complex random number,
!  then, after multiplying by i, we have i*z = i*(z'+i*z") = -z" + i*z'.
!
                  u_re(ikx,iky,ikz,1)=-kx(ikx+ipx*nx)*u_im(ikx,iky,ikz,1)*ctime
                  u_im(ikx,iky,ikz,1)=+kx(ikx+ipx*nx)*u_re(ikx,iky,ikz,1)*ctime
                  u_re(ikx,iky,ikz,2)=-ky(iky+ipy*ny)*u_im(ikx,iky,ikz,1)*ctime
                  u_im(ikx,iky,ikz,2)=+ky(iky+ipy*ny)*u_re(ikx,iky,ikz,1)*ctime
                  u_re(ikx,iky,ikz,3)=-kz(ikz+ipz*nz)*u_im(ikx,iky,ikz,1)*ctime
                  u_im(ikx,iky,ikz,3)=+kz(ikz+ipz*nz)*u_re(ikx,iky,ikz,1)*ctime
                elseif (ltime_old1) then
                  u_re(ikx,iky,ikz,1:3)=+u_re(ikx,iky,ikz,1:3)*stime
                  u_im(ikx,iky,ikz,1:3)=-u_im(ikx,iky,ikz,1:3)*ctime
                else
                  u_re(ikx,iky,ikz,1:3)=+u_re(ikx,iky,ikz,1:3)*ctime
                  u_im(ikx,iky,ikz,1:3)=+u_im(ikx,iky,ikz,1:3)*ctime
                endif
              enddo
            enddo
          enddo
        endif
!
!  back to real space, unless lremain_in_fourier=T, in which case
!  we assume that the imaginary part is just next to the real part.
!
        if (lremain_in_fourier1) then
          f(l1:l2,m1:m2,n1:n2,i1  :i2  )=f(l1:l2,m1:m2,n1:n2,i1  :i2  )+u_re
          f(l1:l2,m1:m2,n1:n2,i1+3:i2+3)=f(l1:l2,m1:m2,n1:n2,i1+3:i2+3)+u_im
        else
          do i=1,3
            call fft_xyz_parallel(u_re(:,:,:,i),u_im(:,:,:,i),linv=.true.)
          enddo
          if (loptest(lreinit)) then
            f(l1:l2,m1:m2,n1:n2,i1:i2)=u_re
          else
            f(l1:l2,m1:m2,n1:n2,i1:i2)=f(l1:l2,m1:m2,n1:n2,i1:i2)+u_re
          endif
!
          if (lrho_nonuni1) then
!
!  u_re, u_im now used for lnrho.
!
            do ikz=1,nz
              do iky=1,ny
                do ikx=1,nx
                  u_re(ikx,iky,ikz,1)=(kx(ikx+ipx*nx)*u_re(ikx,iky,ikz,1) &
                                      +ky(iky+ipy*ny)*u_re(ikx,iky,ikz,2) &
                                      +kz(ikz+ipz*nz)*u_re(ikx,iky,ikz,3) &
                                      )/(cs1*k2(ikx,iky,ikz))
                  u_im(ikx,iky,ikz,1)=(kx(ikx+ipx*nx)*u_im(ikx,iky,ikz,1) &
                                      +ky(iky+ipy*ny)*u_im(ikx,iky,ikz,2) &
                                      +kz(ikz+ipz*nz)*u_im(ikx,iky,ikz,3) &
                                      )/(cs1*k2(ikx,iky,ikz))
                enddo
              enddo
            enddo
            call fft_xyz_parallel(u_re(:,:,:,1),u_im(:,:,:,1),linv=.true.)
            f(l1:l2,m1:m2,n1:n2,ilnr1)=u_re(:,:,:,1)
          endif
        endif
      endif   !(i2==i1)
!
!  notification
!
      if (lroot.and..not.ltime) then
        if (cutoff==0) then
          print*,'power_randomphase_hel: k^',initpower,' spectrum: var  i=',i
        else
          print*,'power_randomphase_hel: with cutoff : k^n*exp(-k^4/k0^4) w/ n=', &
                 initpower,', k0 =',cutoff,' : var  i=',i
        endif
      endif
!
!  Deallocate arrays.
!
      if (allocated(k2)) deallocate(k2)
      if (allocated(r))  deallocate(r)
      if (allocated(u_re)) deallocate(u_re)
      if (allocated(u_im)) deallocate(u_im)
      if (allocated(kx)) deallocate(kx)
      if (allocated(ky)) deallocate(ky)
      if (allocated(kz)) deallocate(kz)
!
    endsubroutine power_randomphase_hel
!***********************************************************************
    subroutine bunch_davies(f,i1a,i1b,i2a,i2b,ampl,kpeak,deriv_prefactor)
!
!  21-mar-25/axel: adapted from power_randomphase_hel
!  21-may-25/axel: when kpeak<0, interpret is as sharp cutoff; powerlaw otherwise.
!
      use Fourier, only: fft_xyz_parallel
      use General, only: loptest, roptest
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i, i1a, i1b, i2a, i2b, ikx, iky, ikz, stat, ik, nk
      real, dimension (:,:,:,:), allocatable :: u_re, u_im, v_re, v_im
      real, dimension (:,:,:), allocatable :: k1, r
      real, dimension (:), allocatable :: kx, ky, kz
      real, dimension (:), allocatable :: kk
      real :: ampl, kpeak, deriv_prefactor, scale_factor=1.,ksteepness=5.
!
      if (ampl==0.) then
        if (lroot) print*,'bunch_davies: set variables to zero; i1a,i1b,i2a,i2b=',i1a,i1b,i2a,i2b
        !f(:,:,:,i1a:i1b) = 0.
        !f(:,:,:,i2a:i2b) = 0.
        return
      endif
!
!  Allocate memory for arrays r and k1.
!
      allocate(k1(nx,ny,nz),stat=stat)
      if (stat>0) call fatal_error('bunch_davies','Could not allocate k1')
      allocate(r(nx,ny,nz),stat=stat)
      if (stat>0) call fatal_error('bunch_davies','Could not allocate r')
!
!  Complex auxiliary arrary (u_re,u_im) and (v_re,v_im) 
!
      allocate(u_re(nx,ny,nz,3),stat=stat)
      if (stat>0) call fatal_error('bunch_davies','Could not allocate u_re')
      allocate(u_im(nx,ny,nz,3),stat=stat)
      if (stat>0) call fatal_error('bunch_davies','Could not allocate u_im')
      allocate(v_re(nx,ny,nz,3),stat=stat)
      if (stat>0) call fatal_error('bunch_davies','Could not allocate v_re')
      allocate(v_im(nx,ny,nz,3),stat=stat)
      if (stat>0) call fatal_error('bunch_davies','Could not allocate v_im')
!
!  One-dimensional wavevector arrays.
!
      allocate(kx(nxgrid),stat=stat)
      if (stat>0) call fatal_error('bunch_davies','Could not allocate kx')
      allocate(ky(nygrid),stat=stat)
      if (stat>0) call fatal_error('bunch_davies','Could not allocate ky')
      allocate(kz(nzgrid),stat=stat)
      if (stat>0) call fatal_error('bunch_davies','Could not allocate kz')
!
!  Scale factors if box size is not 2*pi
!
      scale_factor=2*pi/Lx
      kx=cshift((/(i-nxgrid/2,i=0,nxgrid-1)/),nxgrid/2)*scale_factor
      if (lroot.and.ip<10) print*,'AXEL: kx=',kx
!
      ky=cshift((/(i-nygrid/2,i=0,nygrid-1)/),nygrid/2)*scale_factor
      if (lroot.and.ip<10) print*,'AXEL: ky=',ky
!
      kz=cshift((/(i-nzgrid/2,i=0,nzgrid-1)/),nzgrid/2)*scale_factor
      if (lroot.and.ip<10) print*,'AXEL: kz=',kz
!
!  Set the 3 components of v_im to Gaussian-distributed random values.
!  If it is a scalar, then there should be only one entry.
!
      do i=1,1+i1b-i1a
        call random_number_wrapper(r)
        call random_number_wrapper(k1)
        v_re(:,:,:,i)=sqrt(-2*log(r))*cos(2*pi*k1)
        v_im(:,:,:,i)=sqrt(-2*log(r))*sin(2*pi*k1)
      enddo
!
!  Set k1 = |k| array.
!
!  In 1-D
      if (nzgrid==1.and.nygrid==1) then
        ikz=1
        iky=1
        do ikx=1,nx
          k1(ikx,iky,ikz)=sqrt(kx(ikx+ipx*nx)**2)
        enddo
!  In 2-D
      elseif (nzgrid==1) then
        ikz=1
        do iky=1,ny
          do ikx=1,nx
            k1(ikx,iky,ikz)=sqrt(kx(ikx+ipx*nx)**2+ky(iky+ipy*ny)**2)
          enddo
        enddo
!  In 3-D
      else
        do ikz=1,nz
          do iky=1,ny
              do ikx=1,nx
                k1(ikx,iky,ikz)=sqrt(kx(ikx+ipx*nx)**2+ky(iky+ipy*ny)**2+kz(ikz+ipz*nz)**2)
              enddo
          enddo
        enddo
      endif
      if (lroot) k1(1,1,1) = 1.  ! To avoid division by zero.
!
!  Put cutoff at kpeak in v_im.
!
!     if (kpeak<0.) then
!       where(k1>=abs(kpeak))
!         v_im(:,:,:,1)=0.
!         v_im(:,:,:,2)=0.
!         v_im(:,:,:,3)=0.
!       endwhere
!     else
!       where(k1>=kpeak)
!         v_im(:,:,:,1)=v_im(:,:,:,1)*(kpeak/k1)**3
!         v_im(:,:,:,2)=v_im(:,:,:,2)*(kpeak/k1)**3
!         v_im(:,:,:,3)=v_im(:,:,:,3)*(kpeak/k1)**3
!       endwhere
!     endif
!
!  Compute Bunch-Davies vacuum, A = e^(-i*k*eta)/sqrt(2*k), so
!  E = -dA/deta = +i*k*e^(-i*k*eta)/sqrt(2*k) = i*e^(-i*k*eta)*sqrt(k/2)
!  Here, v_im serves as a temporary array until the last line.
!  The correct prefactor ampl of H = Hscript/a should be applied in the
!  call tp this routine.
!  exp(-i*k1) = cos(-i*k1) + i*sin(-i*k1)
!
      do i=1,1+i1b-i1a
        !u_re(:,:,:,i)=+ampl*v_im(:,:,:,i)*cos(-k1)/sqrt(k1*2.)
        !u_im(:,:,:,i)=+ampl*v_im(:,:,:,i)*sin(-k1)/sqrt(k1*2.)
        u_re(:,:,:,i)=+ampl*v_re(:,:,:,i)/sqrt(2.*k1)*.5*(1.-tanh(ksteepness*(k1/kpeak-1.)))
        u_im(:,:,:,i)=+ampl*v_im(:,:,:,i)/sqrt(2.*k1)*.5*(1.-tanh(ksteepness*(k1/kpeak-1.)))
        v_re(:,:,:,i)=-k1*u_im(:,:,:,i)
        v_im(:,:,:,i)=+k1*u_re(:,:,:,i)
print*,'AXEL2a=',iproc,i,sum(u_re(:,:,:,i)**2+u_im(:,:,:,i)**2)
print*,'AXEL2b=',iproc,i,sum(v_re(:,:,:,i)**2+v_im(:,:,:,i)**2)
      enddo
!
!  Fourier transform to real space.
!
      do i=1,1+i1b-i1a
        call fft_xyz_parallel(u_re(:,:,:,i),u_im(:,:,:,i),linv=.true.)
        call fft_xyz_parallel(v_re(:,:,:,i),v_im(:,:,:,i),linv=.true.)
print*,'AXEL3a=',iproc,i,sum(u_re(:,:,:,i)**2+u_im(:,:,:,i)**2)/nwgrid
print*,'AXEL4a=',iproc,i,sum(u_re(:,:,:,i)**2)/nxgrid
print*,'AXEL3b=',iproc,i,sum(v_re(:,:,:,i)**2+v_im(:,:,:,i)**2)/nwgrid
print*,'AXEL4b=',iproc,i,sum(v_re(:,:,:,i)**2)/nxgrid
print*,'AXEL5 =',iproc,nwgrid
      enddo
!
!  Use real parts of u and v for A and E.
!
      f(l1:l2,m1:m2,n1:n2,i1a:i1b)=f(l1:l2,m1:m2,n1:n2,i1a:i1b)+u_re
      f(l1:l2,m1:m2,n1:n2,i2a:i2b)=f(l1:l2,m1:m2,n1:n2,i2a:i2b)+v_re*deriv_prefactor
!
!  Deallocate arrays.
!
      if (allocated(k1)) deallocate(k1)
      if (allocated(r))  deallocate(r)
      if (allocated(u_re)) deallocate(u_re)
      if (allocated(u_im)) deallocate(u_im)
      if (allocated(v_re)) deallocate(v_re)
      if (allocated(v_im)) deallocate(v_im)
      if (allocated(kx)) deallocate(kx)
      if (allocated(ky)) deallocate(ky)
      if (allocated(kz)) deallocate(kz)
!
    endsubroutine bunch_davies
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
    subroutine random_isotropic_shell(f,jf,ampl0,z1,z2)
!
!   random_isotropic_shell
!
!   25-jun-15/axel: coded
!
      use Sub, only: cross, dot, dot2
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      character (len=1) :: dummy
      complex :: ii=(0.,1.)
      integer :: jf,i,nvect,ivect
      real, dimension (3) :: ee,kk,exk
      real, dimension (nx) :: kdotx
      real :: exk2,phik,ampl,ampl0,prof,zeta,z1,z2
!
!  read header
!
      open(10,file='kvect.dat')
      read(10,*) nvect
      read(10,*) dummy
!
!  loop through all vectors
!
      f(:,:,:,jf:jf+2)=0.
      do ivect=1,nvect
        read(10,*) kk,ee,phik,ampl
        call cross(ee,kk,exk)
        call dot2(exk,exk2)
        exk=exk/sqrt(exk2)
!
!  allow for the possibility of a profile
!
        do n=n1,n2
          if (z1/=z2) then
            if (z(n)>z1.and.z(n)<z2) then
              zeta=(z(n)-z1)/(z2-z1)-.5
              prof=cos(pi*zeta)**2
            else
              prof=0.
            endif
          else
            prof=1.
          endif
          prof=prof*ampl0*ampl
!
          do m=m1,m2
            kdotx=kk(1)*x(l1:l2)+kk(2)*y(m)+kk(3)*z(n)
            do i=1,3
              f(l1:l2,m,n,jf+i-1)=f(l1:l2,m,n,jf+i-1) &
                +prof*exk(i)*real(exp(ii*(kdotx+phik)))
            enddo
          enddo
        enddo
      enddo
      close(10)

    endsubroutine random_isotropic_shell
!***********************************************************************
    subroutine corona_init(f,gamma)
!
!  Initialize the density for a given temperature profile
!  in the vertical (z) direction by solving for hydrostatic
!  equilibrium.
!  The temperature is hard coded as three polynomials.
!
!  07-dec-05/bing : coded.
!
      use EquationOfState, only: lnrho0,cs20,cs2top,cs2bot
      use Mpicomm, only: mpibcast_real
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: tmp,ztop,zbot,gamma,gamma_m1,dummy=1.
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
      if (pretend_lnTT) call not_implemented('corona_init','for pretend_lnTT=T')
!
      if (lroot) then
        inquire(IOLENGTH=lend) dummy
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

      gamma_m1=gamma-1.
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
                  f(:,:,j,iss) = (log(gamma_m1/cs20)+tmp- &
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
               f(:,:,j,iss) = (log(gamma_m1/cs20)+tmp- &
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
    subroutine mdi_init(f,periodic,z0aa)
!
!  Intialize the vector potential
!  by potential field extrapolation
!  of a mdi magnetogram
!
!  13-dec-05/bing : coded.
!
      use Fourier, only: fourier_transform_other, kx_fft, ky_fft
      use Mpicomm, only: mpibcast_real,stop_it_if_any
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
!
      real, dimension (:,:), allocatable :: kx,ky,k2,Bz0_i,Bz0_r,A_r,A_i
      logical, intent (in) :: periodic
      real, intent (in) :: z0aa
      real :: zref,Bzflux
      logical :: exists
      integer :: i,j,idx2,idy2,stat,iostat,lend
      integer :: nxinit,nyinit,itmp
!
      ! file location settings
      character (len=*), parameter :: mag_field_dat = 'driver/mag_field.dat'
!
!  Allocate memory for arrays.
!
      if (.not.lequidist(1).or..not.lequidist(2)) call not_implemented('mdi_init', &
          'for non-equidistant grids')
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
        zref = max(0.,z(i) - z0aa)
        if (abs(z(i)-z0aa) <= dz/2. ) then
          itmp = i
          print*,'Height of reference irefz for run.in: ',itmp
        endif
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
!  Save initial vector potential
      if (ipz==0) then
        open (11,file=trim(directory_snap)//trim('/Ax_init.dat'),form='unformatted')
        write (11) f(l1:l2,m1:m2,itmp,iax)
        close (11)
        open (11,file=trim(directory_snap)//trim('/Ay_init.dat'),form='unformatted')
        write (11) f(l1:l2,m1:m2,itmp,iay)
        close (11)
        open (11,file=trim(directory_snap)//trim('/Az_init.dat'),form='unformatted')
        write (11) f(l1:l2,m1:m2,itmp,iaz)
        close (11)
      endif
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
      integer(KIND=ikind8) :: rlen
      real, parameter :: reduce_factor=0.25
!
      ! file location settings
      character (len=*), parameter :: mag_field_dat = 'driver/mag_field.dat'
!
      if (.not. lperi(1) .or. .not. lperi(2)) call fatal_error ('mag_init', &
          'Currently only implemented for xy-periodic setups!')
      if (.not. lequidist(1) .or. .not. lequidist(2)) call not_implemented('mag_init', &
          'for non-equidistant grids')
      if (mod (nygrid, nprocxy) /= 0) call fatal_error ('mag_init', &
          'nygrid needs to be an integer multiple of nprocx*nprocy')
      if (mod (nxgrid, nprocxy) /= 0) call fatal_error ('mag_init', &
          'nxgrid needs to be an integer multiple of nprocx*nprocy')
!
!  Allocate memory for arrays.
!
      allocate (Bz(bnx,bny), stat=alloc_err)
      if (alloc_err > 0) call fatal_error('mag_init','Could not allocate Bz', .true.)
      allocate (exp_fact(enx,eny,mz), stat=alloc_err)
      if (alloc_err > 0) call fatal_error('mag_init','Could not allocate exp_fact', .true.)
!
      if (lroot) then
        inquire (file=mag_field_dat, exist=exists)
        call stop_it_if_any(.not. exists, &
            'mag_init: Magnetogram file not found: "'//trim(mag_field_dat)//'"')
        inquire (iolength=rec_len) 1.0d0
        rlen=rec_len*bnx*bny
        open (unit, file=mag_field_dat, form='unformatted', recl=rlen, access='direct')
        do py = 1, nprocxy-1
          partner = find_proc(modulo(py,nprocx),py/nprocx,0)
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
          call mpirecv_real (Bz, (/ bnx, bny /), 0, tag_xy)
        endif
      endif
!
      if (nprocz > 1) then
        ! distribute Bz along the z-direction
        if (lfirst_proc_z) then
          do pz = 1, nprocz-1
            partner = find_proc(ipx,ipy,pz)
            call mpisend_real (Bz, (/ bnx, bny /), partner, tag_z)
          enddo
        else
          partner = find_proc(ipx,ipy,0)
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
!  Save initial vector potential
      if (ipz==0) then
        open (11,file=trim(directory_snap)//trim('/Ax_init.dat'),form='unformatted')
        write (11) f(l1:l2,m1:m2,n1,iax)
        close (11)
        open (11,file=trim(directory_snap)//trim('/Ay_init.dat'),form='unformatted')
        write (11) f(l1:l2,m1:m2,n1,iay)
        close (11)
        open (11,file=trim(directory_snap)//trim('/Az_init.dat'),form='unformatted')
        write (11) f(l1:l2,m1:m2,n1,iaz)
        close (11)
      endif
!
    endsubroutine mag_init
!***********************************************************************
    subroutine mag_Az_init(f)
!
!  Intialize the vector potential A with a potential field extrapolation from the bottom boundary read from a file.
!
!  30-Mar-2018/Bourdin.KIS : coded.
!
      use Mpicomm, only: stop_it_if_any, distribute_xy
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      real, dimension (:,:), allocatable :: A_global, A_local
      integer, parameter :: unit=11
      logical :: exists
      integer :: alloc_err, rec_len
!
      ! file location settings
      character (len=*), parameter :: Az_init_dat = 'Az_init.dat'
!
!  Allocate memory for arrays.
!
      if (lroot) then
        allocate (A_global(nxgrid,nygrid), stat=alloc_err)
        if (alloc_err > 0) call fatal_error('mag_Az_init', &
            'Could not allocate memory for A_global', .true.)
      endif
      if (lfirst_proc_z) then
        allocate (A_local(nx,ny), stat=alloc_err)
        if (alloc_err > 0) call fatal_error('mag_Az_init', &
            'Could not allocate memory for A_local', .true.)
      endif
!
      call mag_init(f)
!
      if (lroot) then
        inquire (file=Az_init_dat, exist=exists)
        call stop_it_if_any(.not. exists, &
            'mag_Az_init: file not found: "'//trim(Az_init_dat)//'"')
        inquire (iolength=rec_len) 1.0d0
        open (unit, file=Az_init_dat, form='unformatted', recl=rec_len*nxgrid*nygrid, access='direct')
        ! read A components for each ipz layer
        read (unit, rec=1) A_global
        close (unit)
        ! distribute A component
        call distribute_xy(A_local, A_global)
        f(l1:l2,m1:m2,1,iaz) = A_local
      elseif (lfirst_proc_z) then
        call stop_it_if_any(.false.,'')
        call distribute_xy(A_local)
        f(l1:l2,m1:m2,1,iaz) = A_local
      else
        call stop_it_if_any(.false.,'')
      endif
!
      if (lroot) deallocate (A_global)
      if (lfirst_proc_z) deallocate (A_local)
!
    endsubroutine mag_Az_init
!***********************************************************************
    subroutine file_init(f)
!
!  Intialize the vector potential with a vector potential A from a file.
!
!  30-Mar-2018/Bourdin.KIS : coded.
!
      use Mpicomm, only: stop_it_if_any, distribute_xy, mpisend_real, mpirecv_real
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      real, dimension (:,:,:), allocatable :: A_global, A_local
      integer, parameter :: unit=11
      integer, parameter :: tag_z=133
      logical :: exists
      integer :: partner, pz, comp, alloc_err, rec_len
!
      ! file location settings
      character (len=*), parameter :: A_init_dat = 'Axyz_init.dat'
!
!  Allocate memory for arrays.
!
      if (lfirst_proc_xy) then
        allocate (A_global(nxgrid,nygrid,nz), stat=alloc_err)
        if (alloc_err > 0) call fatal_error('file_init', &
            'Could not allocate memory for A_global', .true.)
      endif
      allocate (A_local(nx,ny,nz), stat=alloc_err)
      if (alloc_err > 0) call fatal_error('file_init', &
          'Could not allocate memory for A_local', .true.)
!
      if (lroot) then
        inquire (file=A_init_dat, exist=exists)
        call stop_it_if_any(.not. exists, &
            'file_init: file not found: "'//trim(A_init_dat)//'"')
        inquire (iolength=rec_len) 1.0d0
        open (unit, file=A_init_dat, form='unformatted', recl=rec_len*nxgrid*nygrid*nz, access='direct')
        ! read A components for each ipz layer
        do comp = 1, 3
          do pz = 0, nprocz-1
            read (unit, rec=1+pz+nprocz*(comp-1)) A_global
            ! distribute A component
            if (pz == ipz) then
              call distribute_xy(A_local, A_global)
              f(l1:l2,m1:m2,n1:n2,iax+(comp-1)) = A_local
              if (iglobal_ax_ext/=0 .and. comp == 1) then
                f(:,:,:, iglobal_ax_ext)=0.
                f(l1:l2,m1:m2,n1:n2,iglobal_ax_ext) = A_local
              endif
              if (iglobal_ay_ext/=0 .and. comp == 2) then
                f(:,:,:, iglobal_ay_ext)=0.
                f(l1:l2,m1:m2,n1:n2,iglobal_ay_ext) = A_local
              endif
              if (iglobal_az_ext/=0 .and. comp == 3) then
                f(:,:,:, iglobal_az_ext)=0.
                f(l1:l2,m1:m2,n1:n2,iglobal_az_ext) = A_local
              endif
            else
              partner = find_proc(ipx,ipy,pz)
              call mpisend_real (A_global, (/ nxgrid, nygrid, nz /), partner, tag_z)
            endif
          enddo
        enddo
        close (unit)
      else
        call stop_it_if_any(.false.,'')
        do comp = 1, 3
          if (lfirst_proc_xy) then
            partner = find_proc(ipx,ipy,0)
            call mpirecv_real (A_global, (/ nxgrid, nygrid, nz /), partner, tag_z)
            call distribute_xy(A_local, A_global)
          else
            call distribute_xy(A_local)
          endif
          f(l1:l2,m1:m2,n1:n2,iax+(comp-1)) = A_local
          if (iglobal_ax_ext/=0 .and. comp == 1) then
            f(:,:,:, iglobal_ax_ext)=0.
            f(l1:l2,m1:m2,n1:n2,iglobal_ax_ext) = A_local
          endif
          if (iglobal_ay_ext/=0 .and. comp == 2) then
            f(:,:,:, iglobal_ay_ext)=0.
            f(l1:l2,m1:m2,n1:n2,iglobal_ay_ext) = A_local
          endif
          if (iglobal_az_ext/=0 .and. comp == 3) then
            f(:,:,:, iglobal_az_ext)=0.
            f(l1:l2,m1:m2,n1:n2,iglobal_az_ext) = A_local
          endif
        enddo
      endif
!
      if (lfirst_proc_xy) deallocate (A_global)
      deallocate (A_local)
!
    endsubroutine file_init
!***********************************************************************
    subroutine temp_hydrostatic(f,rho0,gamma)
!
! 07-dec-05/bing : coded.
!      intialize the density for a given temperature profile in vertical
!      z direction by solving the hydrostatic equilibrium:
!      dlnrho = - dlnTT + (cp-cv)/T g dz
!
      use EquationOfState, only: lnrho0,cs2top,cs2bot
      use Gravity, only: gravz
      use Mpicomm, only: mpibcast_real
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ztop,zbot,dummy=1.
      integer, parameter :: prof_nz=150
      real, dimension (prof_nz) :: prof_lnT,prof_z
      real :: tmprho,tmpT,tmpdT,tmpz,dz_step,lnrho_0,rho0,gamma
      integer :: i,lend,j
!
      ! file location settings
      character (len=*), parameter :: lnT_dat = 'driver/b_lnT.dat'
!
      ! temperature given as function lnT(z) in SI units
      ! [T] = K   &   [z] = Mm   & [rho] = kg/m^3
      if (pretend_lnTT) call not_implemented('temp_hydrostatic','for pretend_lnTT=T')
!
      if (lnrho0 > 0.99*log(impossible)) then
        call warning("temp_hydrostatic", &
                     "lnrho0 from eos module not useful, use rho_const from density instead")
        lnrho_0 = log(rho0)
      else
        lnrho_0 = lnrho0
      endif
!
      if (lroot) then
        inquire(IOLENGTH=lend) dummy
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
!  Galactic-hs or Ferriere-hs to be set for initss, in density.f90
!  Galactic-hs should be set for initlnrho and in gravity_simple.f90 use
!  Ferriere for gravz_profile
!
!  09-jan-10/fred: coded
!
      use Mpicomm, only: mpireduce_sum, mpibcast_real
!
      integer :: i,icpu
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,tmp1,tmp3
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
        call mpibcast_real(tmp3,icpu-1)
        tmp2(icpu)=tmp3
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
!  Galactic-hs or Ferriere-hs to be set for initss, in density.f90
!  Galactic-hs should be set for initlnrho and in gravity_simple.f90 use
!  Ferriere for gravz_profile
!
      use Mpicomm, only: mpireduce_sum, mpibcast_real
!
      integer :: i ,icpu
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,tmp1,tmp3
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
        call mpibcast_real(tmp3,icpu-1)
        tmp2(icpu)=tmp3
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
!  with inclination angle alpha.
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
    subroutine rotblob_yz(ampl,f,i,radius,xsphere,ysphere,zsphere)
!
!  Rigid rotating sphere initial velocity
!  with inclination angle alpha.
!
!  18-feb-10/mvaisala: coded
!
      integer :: i,j,l
      real, dimension (mx,my,mz,mfarray) :: f
      real, optional :: xsphere,ysphere,zsphere
      real :: omega, radius, phi, rr_rot, ampl, x01=0.
      real :: y01=0., z01=0., x_real, y_real, z_real
      real :: theta, vel_phi
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
              x_real = x(l) - x01 !
              y_real = y(m) - y01 !
              z_real = z(n) - z01 !
              rr_rot = sqrt((x_real)**2+(y_real)**2+(z_real)**2) !
              theta  = atan(abs(x_real)/sqrt((x_real)**2+(z_real)**2)) !
              phi    = atan(y_real/z_real) !
              vel_phi = omega*rr_rot*cos(theta) !
              if (rr_rot <= radius) then
                j = i+1
                f(l,m,n,j) = -vel_phi*cos(phi)
                j = i+2
                f(l,m,n,j) = vel_phi*sin(phi)
                if (z_real < 0.0) then
                  j = i+1
                  f(l,m,n,j) = -f(l,m,n,j)
                  j = i+2
                  f(l,m,n,j) = -f(l,m,n,j)
                endif
              endif
            enddo
          enddo
        enddo
      endif
!
    endsubroutine rotblob_yz
!***********************************************************************
    subroutine dipole(f,ix,amp,r_inner_,r_outer_)
!
!  initial vector potential for a
!  purely poloidal axisymmetric field \vec{B} \sim [f_r \cos\theta, f_\theta \sin\theta, 0]
!
!  18-jun-13/MR: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: ix
      real :: amp, rpart, rr2, pom2, r_inner, r_outer
      real, optional :: r_inner_, r_outer_

      integer :: l,m,n
!
      if (present(r_inner_).and.present(r_outer_)) then
        r_inner=r_inner_
        r_outer=r_outer_
      else
        r_inner=xyz0(1)
        r_outer=xyz1(1)
      endif
!
!  Dipolar in the z direction in Cartesian coordinates
!
      if (lcartesian_coords) then
        do n=n1,n2
        do m=m1,m2
        do l=l1,l2
          pom2=x(l)**2+y(m)**2
          rr2=pom2+z(n)**2
          if (rr2<=1.) then
            f(l,m,n,ix+1)=amp*sqrt(pom2)
          else
            f(l,m,n,ix+1)=amp*sqrt(pom2/rr2**3)
          endif
        enddo
        enddo
        enddo
!
!  Dipole in spherical coordinantes
!
      elseif (lspherical_coords) then
        !do m = m1,m2
        !  do l = l1,l2
        do m = 1,my
          do l = 1,mx
            !rpart = amp*(xyz0(1)-x(l))*(xyz1(1)-x(l))
            rpart = amp*(r_inner-x(l))*(r_outer-x(l))
            f(l,m,:,ix:ix+1) = 0.
            !if (y(m1)==0.and.y(m2)==pi) then
              f(l,m,:,ix+2) = rpart*sin(y(m))
            !else
            !  f(l,m,:,ix+2) = rpart*sin(pi*(y(m)-xyz0(2))/Lxyz(2))
            !endif
          enddo
        enddo
      endif
!
    endsubroutine dipole
!***********************************************************************
    subroutine switchback(f,ix,amp,amp2,r_inner_,r_outer_)
!
!  initial vector potential for a
!  purely poloidal axisymmetric field \vec{B} \sim [f_r \cos\theta, f_\theta \sin\theta, 0]
!
!  18-jun-13/MR: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: ix
      real :: amp, amp2, rpart, rr2, pom2, r_inner, r_outer
      real, optional :: r_inner_, r_outer_

      integer :: l,m,n
!
      if (present(r_inner_).and.present(r_outer_)) then
        r_inner=r_inner_
        r_outer=r_outer_
      else
        r_inner=xyz0(1)
        r_outer=xyz1(1)
      endif
!
!  Dipolar in the z direction in Cartesian coordinates
!
      if (lcartesian_coords) then
        do n=n1,n2
        do m=m1,m2
        do l=l1,l2
          pom2=x(l)**2+y(m)**2
          rr2=pom2+z(n)**2
          if (rr2<=1.) then
            f(l,m,n,ix+1)=amp*sqrt(pom2)
          else
            f(l,m,n,ix+1)=amp*sqrt(pom2/rr2**3)
          endif
        enddo
        enddo
        enddo
!
!  Dipole in spherical coordinantes
!
      elseif (lspherical_coords) then
        !do m = m1,m2
        !  do l = l1,l2
        do m = 1,my
          do l = 1,mx
            if (x(l)>=r_inner.and.x(l)<=r_outer) then
              rpart=amp*sin(2.*pi*(x(l)-r_inner)/(r_outer-r_inner))
            else
              rpart=0.
            endif
            f(l,m,:,ix:ix+1)=0.
            f(l,m,:,ix+2)=rpart+amp2/x(l)
          enddo
        enddo
      endif
!
    endsubroutine switchback
!***********************************************************************
    subroutine dipole_tor(f,ix,amp,ladd)
!
!  initial vector potential for a
!  purely toroidal axisymmetric field B_\phi \sim \sin\theta.
!  Radial dependence parabolic and vanishing at boundaries.
!
!  18-jun-13/MR: coded
!
      use General, only: loptest

      real, dimension (mx,my,mz,mfarray) :: f
      integer :: ix
      real :: amp, rpart
      logical, optional :: ladd

      integer :: m,l

      if (lspherical_coords) then
        do m = 1,my
          do l = 1,mx
            rpart = amp*(xyz0(1)-x(l))*(xyz1(1)-x(l))
            if (loptest(ladd)) then
              f(l,m,:,ix)   = f(l,m,:,ix)  + 2.*rpart*cos(y(m))
              f(l,m,:,ix+1) = f(l,m,:,ix+1)+    rpart*sin(y(m))
            else
              f(l,m,:,ix)   = 2.*rpart*cos(y(m))
              f(l,m,:,ix+1) =    rpart*sin(y(m))
              f(l,m,:,ix+2) = 0.
            endif
          enddo
        enddo
      else
        call not_implemented('dipole_tor', &
        'dipolar toroidal field for non-spherical coordinates')
      endif
!
    endsubroutine dipole_tor
!***********************************************************************
    subroutine tanh_hyperbola(amp,f,ix,yzero,delta,kk)
!
!  initial vector potential as used by Moffatt and Hunt, "A model for magnetic
!  reconnection", 2002.

!
!  23 June 2016/dhruba.mitra
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: amp,yzero,delta,kk
      integer, intent(in) :: ix
      integer :: m,l
      do m = m1,m2
          do l = l1,l2
             f(l,m,:,ix)   = 0.
            f(l,m,:,ix+1) = 0.
            f(l,m,:,ix+2) = 0.5*amp*tanh((y(m)*y(m) - yzero*yzero-kk*kk*x(l)*x(l))/(delta*delta))
          enddo
        enddo

!
    endsubroutine
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
      use EquationOfState, only: eoscalc
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
      character (len=fnlen) :: datafile
      character (len=labellen) :: cloud_mode
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
           endif
           read(19,*) var1, len_file
           bigr = var1/unit_length
           write (*,*) '(Modified BE-sphere) R = ', var1, 'cm =',&
                        bigr, 'pc_units'
           do jj = 1, len_file
             read(19,*) var1, var2, var3
             r_rho(jj) = var1/unit_length
             lnrho_r(jj) = log(var2/unit_density)
             lnTT_file(jj) = log(var3/unit_temperature)
           enddo
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
                   enddo
                   f(l,m,n,ilnrho) = lnrho_r(test)
                   call eoscalc(ilnrho_lnTT, lnrho_r(test), lnTT_file(test),&
                                ss=tmp)
                   f(l,m,n,iss) = tmp
                   counter = counter+1
                 else
                   do while (r_rho(test+1) <= bigr)
                      test = test + 1
                   enddo
                   if (rr_box <= temp_trans*bigr) then
                      x_wave = rr_box-bigr
                      wavelength = 2.0*bigr*(temp_trans-1)
                      QQQ = 0.5*(sin(2.0*pi*x_wave/wavelength - pi/2.0) + 1.0)
                      T_cloud_out_rel0 = 1 + QQQ*(T_cloud_out_rel-1.0)
                   else
                      T_cloud_out_rel0 = T_cloud_out_rel
                   endif
                   f(l,m,n,ilnrho) = log(exp(lnrho_r(test))/T_cloud_out_rel0)
                   lnTTpoint = log(exp(lnTT_file(test))* &
                       T_cloud_out_rel0*temp_coeff_out)
                   call eoscalc(ilnrho_lnTT, f(l,m,n,ilnrho), &
                                lnTTpoint, ss=tmp)
                   f(l,m,n,iss) = tmp
                 endif
               enddo
             enddo
           enddo
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
           endif
           read(19,*) var1, len_file, var2
           bigr = var1/unit_length
           lnTTpoint0 = log(var2*temp_coeff/unit_temperature)
           write (*,*) '(Isothermal BE-sphere) R =', var1, 'cm =',&
                        bigr, 'pc_units, T =', var2, 'K'
           do jj = 1, len_file
             read(19,*) var1, var2
             r_rho(jj) = var1/unit_length
             lnrho_r(jj) = log(var2*dens_coeff/unit_density)
           enddo
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
                   enddo
                   f(l,m,n,ilnrho) = lnrho_r(test)
                   call eoscalc(ilnrho_lnTT, f(l,m,n,ilnrho), lnTTpoint, ss=tmp)
                   f(l,m,n,iss) = tmp
                   counter = counter+1
                 else
                   do while (r_rho(test+1) <= bigr)
                     test = test + 1
                   enddo
                   if (rr_box <= temp_trans*bigr) then
                      x_wave = rr_box-bigr
                      wavelength = 2.0*bigr*(temp_trans-1)
                      QQQ = 0.5*(sin(2.0*pi*x_wave/wavelength - pi/2.0) + 1.0)
                      T_cloud_out_rel0 = 1 + QQQ*(T_cloud_out_rel-1.0)
                   else
                      T_cloud_out_rel0 = T_cloud_out_rel
                   endif
                   !WRITE (*,*) 'T_cloud_out_rel:', T_cloud_out_rel0
                   f(l,m,n,ilnrho) = log(exp(lnrho_r(test))/T_cloud_out_rel0)
                   lnTTpoint = log(exp(lnTTpoint0)* &
                       T_cloud_out_rel0*temp_coeff_out)
                   call eoscalc(ilnrho_lnTT, f(l,m,n,ilnrho), lnTTpoint, ss=tmp)
                   f(l,m,n,iss) = tmp
                 endif
               enddo
             enddo
           enddo
           write (*,*) 'Covers:', counter, '/', &
                       (abs(n1-n2)*abs(m1-m2)*abs(l1-l2)), 'cells'
           write (*,*) (temp_trans*bigr-bigr), wavelength/2.0, 4.0*wavelength/2.0
!
         case default
      end select
!
    endsubroutine pre_stellar_cloud
!***********************************************************************
    subroutine read_outside_scal_array(f, datafile, i)
!
!   Reads a pre-generated setup of a scalar variable (i.e. lnrho) from an ascii file,
!   generated by an external script.
!   13-may-13/mvaisala: created
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      logical :: exfile
      integer :: l,m,n, lfile, mfile, nfile, io_status, i
      real :: value
      character(len=*) :: datafile
!
      inquire (file=datafile, exist=exfile)
      if (.NOT. exfile) then
        call fatal_error(datafile,' no input file')
      else
        open(19,file=datafile)
      end if

!
!   Every processor reads from the same datafile. Depending on the indices,
!   a specific data point in written into proper domain.
!   Format of the data:
!   l  m  n  value
!   The datafiles do not include boundary zones!
!
      io_status = 0
      do while (io_status == 0)
        read (19,*, iostat=io_status) lfile, mfile, nfile, value
        l = lfile - ipx*nx + 3
        m = mfile - ipy*ny + 3
        n = nfile - ipz*nz + 3
        if (l >= l1 .AND. m >= m1 .AND. n >= n1 .AND. l <= l2 .AND. m <= m2 .AND. n <= n2) then
          f(l,m,n,i) = value
        end if
      end do
!
      close(19)
!
    endsubroutine read_outside_scal_array
!***********************************************************************
    subroutine read_outside_vec_array(f, datafile, i, lbinary,ampl)
!
!   Reads a pre-generated setup of a vector variable (i.e. uu) from an ascii file,
!   generated by an external script.
!   13-may-13/mvaisala: created
!
      use IO, only: input_snap, input_snap_finalize
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx,my,mz,3) :: apot
      logical :: exfile, lbin=.false.
      logical, optional :: lbinary
      integer :: l,m,n, lfile, mfile, nfile, io_status, i
      real :: value1, value2, value3, scale_aa=1.0
      real, optional :: ampl
      character(len=*) :: datafile
!
      inquire (file=datafile, exist=exfile)
      if (present(lbinary)) lbin=lbinary
      if (present(ampl)) scale_aa=ampl
      if (.NOT. exfile) then
        call fatal_error(datafile,' no input file')
      else
        if (lbin) then
          call input_snap(datafile,apot,3,0)
          call input_snap_finalize()
!          apot(:,:,:,1)=0.0
!          apot(:,:,:,2)=spread(spread(0.1*(alog(exp((x+2)/0.2)+exp(-(x+2)/0.2))-&
!                       alog(exp((x+1)/0.2)+exp(-(x+1)/0.2))),2,my),3,mz)-&
!                       spread(spread(0.1*(alog(exp((x-1)/0.2)+exp(-(x-1)/0.2))-&
!                       alog(exp((x-2)/0.2)+exp(-(x-2)/0.2))),2,my),3,mz)
!          apot(:,:,:,2)=spread(spread(0.1*(alog(exp((x+5)/0.2)+exp(-(x+5)/0.2))-&
!                       alog(exp((x+4)/0.2)+exp(-(x+4)/0.2))),2,my),3,mz)-&
!                       spread(spread(0.1*(alog(exp((x+1)/0.2)+exp(-(x+1)/0.2))-&
!                       alog(exp((x-0)/0.2)+exp(-(x-0)/0.2))),2,my),3,mz) + &
!                       -spread(spread(0.1*(alog(exp((x+0)/0.2)+exp(-(x+0)/0.2))-&
!                       alog(exp((x-1)/0.2)+exp(-(x-1)/0.2))),2,my),3,mz)-&
!                       spread(spread(0.1*(alog(exp((x-5)/0.2)+exp(-(x-5)/0.2))-&
!                       alog(exp((x-4)/0.2)+exp(-(x-4)/0.2))),2,my),3,mz)
!          apot(:,:,:,3)=0.0
          f(:,:,:,i  ) = f(:,:,:,i  )+scale_aa*apot(:,:,:,1)
          f(:,:,:,i+1) = f(:,:,:,i+1)+scale_aa*apot(:,:,:,2)
          f(:,:,:,i+2) = f(:,:,:,i+2)+scale_aa*apot(:,:,:,3)
          return
        else
          open(19,file=datafile)
        endif
      end if
!
!   Every processor reads from the same datafile. Depending on the indices,
!   a specific data point in written into proper domain.
!   Format of the data:
!   l  m  n  value_x  value_y  value_z
!   The datafiles do not include boundary zones!
!
      io_status = 0
      do while (io_status == 0)
        read (19,*, iostat=io_status) lfile, mfile, nfile, value1, value2, value3
        l = lfile - ipx*nx + 3
        m = mfile - ipy*ny + 3
        n = nfile - ipz*nz + 3
        if (l >= l1 .AND. m >= m1 .AND. n >= n1 .AND. l <= l2 .AND. m <= m2 .AND. n <= n2) then
          f(l,m,n,i) = value1
          f(l,m,n,i+1) = value2
          f(l,m,n,i+2) = value3
        end if
      end do
!
      close(19)
!
    endsubroutine read_outside_vec_array
!***********************************************************************
endmodule Initcond
