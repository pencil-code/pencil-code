! $Id: initcond.f90 13608 2010-04-09 11:56:42Z mppiyali $
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
  use Mpicomm
  use Sub
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
  public :: gaussian, gaussian3d, gaussianpos, beltrami
  public :: rolls, tor_pert
  public :: jump, bjump, bjumpz
  public :: modes, modev, modeb, crazy
  public :: trilinear, triquad
  public :: diffrot, olddiffrot
  public :: planet, planet_hc
  public :: random_isotropic_KS
  public :: htube, htube2, htube_x, hat, hat3d
  public :: htube_erf
  public :: wave_uu, wave, parabola, linprof
  public :: sinxsinz, cosx_cosy_cosz, cosx_coscosy_cosz
  public :: x_siny_cosz, x1_siny_cosz, x1_cosy_cosz, lnx_cosy_cosz
  public :: sinx_siny_sinz, cosx_siny_cosz, sinx_siny_cosz
  public :: sin2x_sin2y_cosz, cos2x_cos2y_cos2z, x3_cosy_cosz
  public :: cosx_cosz, cosy_cosz, cosy_sinz
  public :: cosxz_cosz, cosyz_sinz
  public :: halfcos_x, magsupport, vfield
  public :: uniform_x, uniform_y, uniform_z, uniform_phi
  public :: vfluxlayer, hfluxlayer, hfluxlayer_y
  public :: vortex_2d
  public :: vfield2
  public :: hawley_etal99a
  public :: robertsflow
  public :: innerbox
  public :: strange,phi_siny_over_r2
  public :: ferriere_uniform_x, ferriere_uniform_y
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
    subroutine beltrami(ampl,f,i,kx,ky,kz,phase)
!
!  Beltrami field (as initial condition)
!
!  19-jun-02/axel: coded
!   5-jul-02/axel: made additive (if called twice), kx,ky,kz are optional
!
      integer :: i,j
      real, dimension (mx,my,mz,mfarray) :: f
      real,optional :: kx,ky,kz,phase
      real :: ampl,k=1.,fac,ph
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
        fac=sign(sqrt(abs(ampl/k)),kx)
        if (iproc==0) print*,'beltrami: fac=',fac
        if (ampl==0) then
          if (lroot) print*,'beltrami: ampl=0; kx=',k
        elseif (ampl>0) then
          if (lroot) print*,'beltrami: Beltrami field (pos-hel): kx,i=',k,i
          j=i+1; f(:,:,:,j)=f(:,:,:,j)+fac*spread(spread(sin(k*x+ph),2,my),3,mz)
          j=i+2; f(:,:,:,j)=f(:,:,:,j)+fac*spread(spread(cos(k*x+ph),2,my),3,mz)
        elseif (ampl<0) then
          if (lroot) print*,'beltrami: Beltrami field (neg-hel): kx,i=',k,i
          j=i+1; f(:,:,:,j)=f(:,:,:,j)+fac*spread(spread(cos(k*x+ph),2,my),3,mz)
          j=i+2; f(:,:,:,j)=f(:,:,:,j)+fac*spread(spread(sin(k*x+ph),2,my),3,mz)
        endif
      endif
!
!  set y-dependent Beltrami field
!
      if (present(ky)) then
        k=abs(ky)
        if (k==0) print*,'beltrami: k must not be zero!'
        fac=sign(sqrt(abs(ampl/k)),ky)
        if (ampl==0) then
          if (lroot) print*,'beltrami: ampl=0; ky=',k
        elseif (ampl>0) then
          if (lroot) print*,'beltrami: Beltrami field (pos-hel): ky,i=',k,i
          j=i;   f(:,:,:,j)=f(:,:,:,j)+fac*spread(spread(cos(k*y+ph),1,mx),3,mz)
          j=i+2; f(:,:,:,j)=f(:,:,:,j)+fac*spread(spread(sin(k*y+ph),1,mx),3,mz)
        elseif (ampl<0) then
          if (lroot) print*,'beltrami: Beltrami field (neg-hel): ky,i=',k,i
          j=i;   f(:,:,:,j)=f(:,:,:,j)+fac*spread(spread(sin(k*y+ph),1,mx),3,mz)
          j=i+2; f(:,:,:,j)=f(:,:,:,j)+fac*spread(spread(cos(k*y+ph),1,mx),3,mz)
        endif
      endif
!
!  set z-dependent Beltrami field
!
      if (present(kz)) then
        k=abs(kz)
        if (k==0) print*,'beltrami: k must not be zero!'
        fac=sign(sqrt(abs(ampl/k)),kz)
        if (ampl==0) then
          if (lroot) print*,'beltrami: ampl=0; kz=',k
        elseif (ampl>0) then
          if (lroot) print*,'beltrami: Beltrami field (pos-hel): kz,i=',k,i
          j=i;   f(:,:,:,j)=f(:,:,:,j)+fac*spread(spread(sin(k*z+ph),1,mx),2,my)
          j=i+1; f(:,:,:,j)=f(:,:,:,j)+fac*spread(spread(cos(k*z+ph),1,mx),2,my)
        elseif (ampl<0) then
          if (lroot) print*,'beltrami: Beltrami field (neg-hel): kz,i=',k,i
          j=i;   f(:,:,:,j)=f(:,:,:,j)+fac*spread(spread(cos(k*z+ph),1,mx),2,my)
          j=i+1; f(:,:,:,j)=f(:,:,:,j)+fac*spread(spread(sin(k*z+ph),1,mx),2,my)
        endif
      endif
!
    endsubroutine beltrami
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
    subroutine robertsflow(ampl,f,i)
!
!  Roberts Flow (as initial condition)
!
!   9-jun-05/axel: coded
!
      integer :: i,j
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,k=1.,kf,fac1,fac2
!
!  prepare coefficients
!
      kf=k*sqrt(2.)
      fac1=sqrt(2.)*ampl*k/kf
      fac2=sqrt(2.)*ampl
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
      if (ipz==0) then
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
      if (lroot) rho0 = exp(-lnrhosum_wholebox(1)/(nxgrid*nygrid*nzgrid))
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
      if (lroot) rho0 = exp(-lnrhosum_wholebox(1)/(nxgrid*nygrid*nzgrid))
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
          f(ix,iy,:,i)=ampl*cos(y(iy))
          f(ix,iy,:,i+1)=0.
          f(ix,iy,:,i+2)=0.
       enddo
     enddo
!
    endsubroutine strange
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
      real, dimension (nx) :: tmp,modulate,tube_radius_sqr
      real :: ampl,radius,eps
      real :: center1_x,center1_z
!
! please remove the variable if not needed anymore
      call keep_compiler_quiet(modulate)
!
      if (ampl==0) then
        f(:,:,:,i1:i2)=0
        if (lroot) print*,'htube: set variable to zero; i1,i2=',i1,i2
      else
        if (lroot) then
          print*,'htube: implement y-dependent flux tube in xz-plane; i1,i2=',i1,i2
          print*,'htube: radius,eps=',radius,eps
        endif
!
! completely quenched "gaussian"
!
        do n=n1,n2; do m=m1,m2
          tube_radius_sqr=(x(l1:l2)-center1_x)**2+(z(n)-center1_z)**2
!         tmp=.5*ampl/modulate*exp(-tube_radius_sqr)/&
!                   (max((radius*modulate)**2-tube_radius_sqr,1e-6))
          tmp=1./(1.+tube_radius_sqr/radius**2)
!
!  check whether vector or scalar
!
          if (i1==i2) then
            if (lroot) print*,'htube: set scalar'
            f(l1:l2,m,n,i1)=tmp
          elseif (i1+2==i2) then
            if (lroot) print*,'htube: set vector'
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
!  Horizontal flux tube (for vector potential) which gives error-function border profile
! for the magnetic field. , or passive scalar)
!
!   18-mar-09/dhruba: aped from htube
!
      integer :: i1,i2,l
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: modulate
      real :: ampl,a,eps,width,tmp,radius,a_minus_r
      real :: center1_x,center1_z
!
! please remove the variable if not needed anymore
      call keep_compiler_quiet(modulate)
!
      if (ampl==0) then
        f(:,:,:,i1:i2)=0
        if (lroot) print*,'htube: set variable to zero; i1,i2=',i1,i2
      else
        if (lroot) then
          print*,'htube: implement y-dependent flux tube in xz-plane; i1,i2=',i1,i2
          print*,'htube: radius,eps=',radius,eps
        endif
!
! An integral of error function.
!
        do n=n1,n2; do m=m1,m2;do l=l1,l2
          radius= sqrt((x(l)-center1_x)**2+(z(n)-center1_z)**2)
          a_minus_r= a - radius
          if (radius .gt. tini) then
             tmp = (-(exp(-width*a_minus_r**2))/(4.*sqrt(pi)*width) +  &
                  radius*(1+erfunc(width*a_minus_r))/4. + &
                  2*a*(exp(-(a**2)*(width**2)) - exp(-(a_minus_r**2)*(width**2)))/(8.*radius*width) + &
                  (1+2*(a**2)*(width**2))*(erfunc(a*width) - erfunc(width*a_minus_r))/(8.*radius*width**2))/radius
          else
             tmp = 0
             write(*,*) 'wrong place:radius,tini',radius,tini
          endif
          write(*,*) 'Dhruba:radius,tini,a,a_minus_r,width,tmp',radius,tini,a,a_minus_r,width,tmp
!         tmp=.5*ampl/modulate*exp(-tube_radius_sqr)/&
!                   (max((radius*modulate)**2-tube_radius_sqr,1e-6))
!
!  check whether vector or scalar
!
          if (i1==i2) then
            if (lroot) print*,'htube: set scalar'
            f(l,m,n,i1)=tmp
          elseif (i1+2==i2) then
            if (lroot) print*,'htube: set vector'
            f(l,m,n,i1 )=-(z(n)-center1_z)*tmp*ampl
            f(l,m,n,i1+1)=tmp*eps
            f(l,m,n,i1+2)=+(x(l)-center1_x)*tmp*ampl
         else
            if (lroot) print*,'htube: bad value of i2=',i2
          endif
!
        enddo; enddo;enddo
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
            if (lroot) print*,'htube2: set scalar'
            f(l1:l2,m,n,i1)=tmp
          elseif (i1+2==i2) then
            if (lroot) print*,'htube2: set vector'
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
          f(l1:l2,m,n,i  )=0.0
          f(l1:l2,m,n,i+1)=+ampl*x(l1:l2)
          f(l1:l2,m,n,i+2)=0.0
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
          if (lroot) print*,'posnoise_vect: variable i=',i
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
          print*, m, n
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
    subroutine trilinear(ampl,f,ivar)
!
!  Produce a profile that is linear in any non-periodic direction, but
!  periodic in periodic ones (for testing purposes).
!
!  5-nov-02/wolf: coded
! 23-nov-02/axel: included scaling factor ampl, corrected lperi argument
!
      real :: ampl
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
          tmp = sin(2*pi/Lx*(x(l1:l2)-xyz0(1)-0.25*Lxyz(1)))
        else
          tmp = x(l1:l2)
        endif
!
!  y direction
!
        if (lperi(2)) then
          tmp = tmp + 10*sin(2*pi/Ly*(y(m)-xyz0(2)-0.25*Lxyz(2)))
        else
          tmp = tmp + 10*y(m)
        endif
!
!  z direction
!
        if (lperi(3)) then
          tmp = tmp + 100*sin(2*pi/Lz*(z(n)-xyz0(3)-0.25*Lxyz(3)))
        else
          tmp = tmp + 100*z(n)
        endif
!
        f(l1:l2,m,n,ivar) = ampl*tmp
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
!       if (modeN.eq.1)delk(modeN)=(k(modeN+1)-K(modeN))
!       if (modeN.eq.N_modes)delk(modeN)=(k(modeN)-k(modeN-1))
!       if (modeN.gt.1.and.modeN.lt.N_modes)delk(modeN)=(k(modeN+1)-k(modeN-2))/2.0
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
        if ((k_unit(2).eq.0).and.(k_unit(3).eq.0)) then
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
!      double precision :: g_B
!      double precision, parameter :: g_B_cgs=6.172D20
!      g_B=g_b_cgs/unit_length
      tmp2(:)=0.0
      sumtmp(:)=0.0
          tmp1=sum(exp(f(l1:l2,m1,n1:n2,ilnrho)))
          do icpu=1,ncpus
          tmp3=tmp1
          call mpibcast_real(tmp3,1,icpu-1)
          tmp2(icpu+nprocy)=tmp3(1)
        enddo
        if (ncpus.gt.nprocy) then
          do icpu=nprocy+1,ncpus
            sumtmp(icpu)=sumtmp(icpu-nprocy)+tmp2(icpu)
          enddo
        endif
        if (lroot) print*,'sumtmp =',sumtmp
        print*,'sumtmp on iproc =',sumtmp(iproc+1)
      if (ampl==0) then
        f(:,:,:,i:i+2)=0
        if (lroot) print*,'ferriere_uniform_x: set variable to zero; i=',i
      else
        print*,'ferriere_uniform_x: uniform x-field approx propto rho ; i=',i
        if ((ip<=16).and.lroot) print*,'uniform_x: ampl=',ampl
        do n=n1,n2; do m=m1,m2
          f(l1:l2,m,n,i  )=0.0
          f(l1:l2,m,n,i+1)=-ampl*(sumtmp(iproc+1)+&
              sum(exp(f(l1:l2,m,n1:n,ilnrho))))*dx*dz
!          f(l1:l2,m,n,i+1)=-ampl*g_B*tanh(z(n)/g_B)
          f(l1:l2,m,n,i+2)=0.0
        enddo; enddo
      endif
!
    endsubroutine ferriere_uniform_x
!***********************************************************************
    subroutine ferriere_uniform_y(ampl,f,i)
!
!  Uniform B_y field (for vector potential)
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
      integer :: i,icpu
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,tmp1
      real, dimension(1)::tmp3
      real, dimension(ncpus)::sumtmp,tmp2
!      double precision :: g_B
!      double precision, parameter :: g_B_cgs=6.172D20
!      g_B=g_b_cgs/unit_length
      tmp2(:)=0.0
      sumtmp(:)=0.0
          tmp1=sum(exp(f(l1:l2,m1,n1:n2,ilnrho)))
          do icpu=1,ncpus
          tmp3=tmp1
          call mpibcast_real(tmp3,1,icpu-1)
          tmp2(icpu+nprocy)=tmp3(1)
        enddo
        if (ncpus.gt.nprocy) then
          do icpu=nprocy+1,ncpus
            sumtmp(icpu)=sumtmp(icpu-nprocy)+tmp2(icpu)
          enddo
        endif
        if (lroot) print*,'sumtmp =',sumtmp
        print*,'sumtmp on iproc =',sumtmp(iproc+1)
      if (ampl==0) then
        f(:,:,:,i:i+2)=0
        if (lroot) print*,'ferriere_uniform_y: set variable to zero; i=',i
      else
        print*,'ferriere_uniform_y: uniform y-field approx propto rho ; i=',i
        if ((ip<=16).and.lroot) print*,'uniform_y: ampl=',ampl
        do n=n1,n2; do m=m1,m2
          f(l1:l2,m,n,i)=ampl*(sumtmp(iproc+1)+&
              sum(exp(f(l1:l2,m,n1:n,ilnrho))))*dx*dz
!          f(l1:l2,m,n,i)=ampl*g_B*tanh(z(n)/g_B)
          f(l1:l2,m,n,i+1)=0.0
          f(l1:l2,m,n,i+2)=0.0
        enddo; enddo
      endif
!
    endsubroutine ferriere_uniform_y
!***********************************************************************
endmodule Initcond
