! $Id: pscalar_nolog.f90,v 1.20 2004-04-30 09:30:50 ajohan Exp $

!  This modules solves the passive scalar advection equation
!  Solves for c, not lnc. Keep ilncc and other names involving "ln"
!  and pretend they are *generic* names. A better generic name
!  might be "pscalar", so ipscalar instead of ilncc.

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Pscalar

  use Cparam
  use Cdata

  implicit none

  character (len=labellen) :: initlncc='zero', initlncc2='zero'
  character (len=40) :: tensor_pscalar_file
  logical :: nopscalar=.false.,reinitalize_lncc=.false.

  ! input parameters
  real :: ampllncc=.1, widthlncc=.5, cc_min=0., lncc_min
  real :: ampllncc2=0.,kx_lncc=1.,ky_lncc=1.,kz_lncc=1.,radius_lncc=0.,epsilon_lncc=0.
  real :: eps_ctog=0.01
  real :: unit_rhocc=1.
  real, dimension(3) :: gradC0=(/0.,0.,0./)

  namelist /pscalar_init_pars/ &
       initlncc,initlncc2,ampllncc,ampllncc2,kx_lncc,ky_lncc,kz_lncc, &
       radius_lncc,epsilon_lncc,widthlncc,cc_min,eps_ctog

  ! run parameters
  real :: pscalar_diff=0.,tensor_pscalar_diff=0.
  real :: rhoccm=0., cc2m=0., gcc2m=0.

  namelist /pscalar_run_pars/ &
       pscalar_diff,nopscalar,tensor_pscalar_diff,gradC0, &
       reinitalize_lncc

  ! other variables (needs to be consistent with reset list below)
  integer :: i_rhoccm=0,i_ccmax=0,i_lnccm=0,i_lnccmz=0
  integer :: i_ucm=0,i_uudcm=0,i_Cz2m=0,i_Cz4m=0,i_Crmsm=0
  integer :: i_cc1m=0,i_cc2m=0,i_cc3m=0,i_cc4m=0,i_cc5m=0
  integer :: i_cc6m=0,i_cc7m=0,i_cc8m=0,i_cc9m=0,i_cc10m=0
  integer :: i_gcc1m=0,i_gcc2m=0,i_gcc3m=0,i_gcc4m=0,i_gcc5m=0
  integer :: i_gcc6m=0,i_gcc7m=0,i_gcc8m=0,i_gcc9m=0,i_gcc10m=0

  contains

!***********************************************************************
    subroutine register_pscalar()
!
!  Initialise variables which should know that we solve for passive
!  scalar: ilncc; increase nvar accordingly
!
!  6-jul-02/axel: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_lncc called twice')
      first = .false.
!
      lpscalar = .true.
      ilncc = nvar+1            ! index to access lncc
      nvar = nvar+1             ! added 1 variable
!
      if ((ip<=8) .and. lroot) then
        print*, 'Register_lncc:  nvar = ', nvar
        print*, 'ilncc = ', ilncc
      endif
!
!  Put variable names in array
!
      varname(ilncc) = 'cc'
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: pscalar_nolog.f90,v 1.20 2004-04-30 09:30:50 ajohan Exp $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('Register_lncc: nvar > mvar')
      endif
!
    endsubroutine register_pscalar
!***********************************************************************
    subroutine initialize_pscalar(f)
!
!  Perform any necessary post-parameter read initialization
!  Since the passive scalar is often used for diagnostic purposes
!  one may want to reinitialize it to its initial distribution.
!
!  24-nov-02/tony: coded
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!
!  set to zero and then call the same initial condition
!  that was used in start.csh
!
      if(reinitalize_lncc) then
        f(:,:,:,ilncc)=0.
        call init_lncc_simple(f)
      endif
!
    endsubroutine initialize_pscalar
!***********************************************************************
    subroutine init_lncc_simple(f)
!
!  initialise passive scalar field; called from start.f90
!
!   6-jul-2001/axel: coded
!
      use Cdata
      use Mpicomm
      use Density
      use Sub
      use Initcond
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!
!  identify module
!
      if (lroot) print*,'init_lncc_simple; initlncc=',initlncc
!
      select case(initlncc)
        case('zero'); f(:,:,:,ilncc)=0.
        case('constant'); f(:,:,:,ilncc) = eps_ctog
        case('hat-x'); call hat(ampllncc,f,ilncc,widthlncc,kx=kx_lncc)
        case('hat-y'); call hat(ampllncc,f,ilncc,widthlncc,ky=ky_lncc)
        case('hat-z'); call hat(ampllncc,f,ilncc,widthlncc,kz=kz_lncc)
        case('gaussian-x'); call gaussian(ampllncc,f,ilncc,kx=kx_lncc)
        case('gaussian-y'); call gaussian(ampllncc,f,ilncc,ky=ky_lncc)
        case('gaussian-z'); call gaussian(ampllncc,f,ilncc,kz=kz_lncc)
        case('parabola-x'); call parabola(ampllncc,f,ilncc,kx=kx_lncc)
        case('parabola-y'); call parabola(ampllncc,f,ilncc,ky=ky_lncc)
        case('parabola-z'); call parabola(ampllncc,f,ilncc,kz=kz_lncc)
        case('gaussian-noise'); call gaunoise(ampllncc,f,ilncc,ilncc)
        case('wave-x'); call wave(ampllncc,f,ilncc,kx=kx_lncc)
        case('wave-y'); call wave(ampllncc,f,ilncc,ky=ky_lncc)
        case('wave-z'); call wave(ampllncc,f,ilncc,kz=kz_lncc)
        case('propto-ux'); call wave_uu(ampllncc,f,ilncc,kx=kx_lncc)
        case('propto-uy'); call wave_uu(ampllncc,f,ilncc,ky=ky_lncc)
        case('propto-uz'); call wave_uu(ampllncc,f,ilncc,kz=kz_lncc)
        case('cosx_cosy_cosz'); call cosx_cosy_cosz(ampllncc,f,ilncc,kx_lncc,ky_lncc,kz_lncc)
        case default; call stop_it('init_lncc: bad initlncc='//trim(initlncc))
      endselect
!
!  superimpose something else
!
      select case(initlncc2)
        case('wave-x'); call wave(ampllncc2,f,ilncc,ky=5.)
      endselect
!
!  add floor value if cc_min is set
!
      if(cc_min/=0.) then
        lncc_min=alog(cc_min)
        if(lroot) print*,'set floor value for cc; cc_min=',cc_min
        f(:,:,:,ilncc)=amax1(lncc_min,f(:,:,:,ilncc))
      endif
!
    endsubroutine init_lncc_simple
!***********************************************************************
    subroutine init_lncc(f,xx,yy,zz)
!
!  initialise passive scalar field; called from start.f90
!
!   6-jul-2001/axel: coded
!
      use Cdata
      use Mpicomm
      use Density
      use Sub
      use Initcond
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz)      :: xx,yy,zz,prof
!
      select case(initlncc)
        case('zero'); f(:,:,:,ilncc)=0.
        case('constant'); f(:,:,:,ilncc) = eps_ctog
        case('hat-x'); call hat(ampllncc,f,ilncc,widthlncc,kx=kx_lncc)
        case('hat-y'); call hat(ampllncc,f,ilncc,widthlncc,ky=ky_lncc)
        case('hat-z'); call hat(ampllncc,f,ilncc,widthlncc,kz=kz_lncc)
        case('gaussian-x'); call gaussian(ampllncc,f,ilncc,kx=kx_lncc)
        case('gaussian-y'); call gaussian(ampllncc,f,ilncc,ky=ky_lncc)
        case('gaussian-z'); call gaussian(ampllncc,f,ilncc,kz=kz_lncc)
        case('parabola-x'); call parabola(ampllncc,f,ilncc,kx=kx_lncc)
        case('parabola-y'); call parabola(ampllncc,f,ilncc,ky=ky_lncc)
        case('parabola-z'); call parabola(ampllncc,f,ilncc,kz=kz_lncc)
        case('gaussian-noise'); call gaunoise(ampllncc,f,ilncc,ilncc)
        case('wave-x'); call wave(ampllncc,f,ilncc,kx=kx_lncc)
        case('wave-y'); call wave(ampllncc,f,ilncc,ky=ky_lncc)
        case('wave-z'); call wave(ampllncc,f,ilncc,kz=kz_lncc)
        case('propto-ux'); call wave_uu(ampllncc,f,ilncc,kx=kx_lncc)
        case('propto-uy'); call wave_uu(ampllncc,f,ilncc,ky=ky_lncc)
        case('propto-uz'); call wave_uu(ampllncc,f,ilncc,kz=kz_lncc)
        case('cosx_cosy_cosz'); call cosx_cosy_cosz(ampllncc,f,ilncc,kx_lncc,ky_lncc,kz_lncc)
        case('sound-wave'); f(:,:,:,ilncc)=-ampllncc*cos(kx_lncc*xx)
        case('tang-discont-z')
           print*,'init_lncc: widthlncc=',widthlncc
        prof=.5*(1.+tanh(zz/widthlncc))
        f(:,:,:,ilncc)=-1.+2.*prof
        case('hor-tube'); call htube2(ampllncc,f,ilncc,ilncc,xx,yy,zz,radius_lncc,epsilon_lncc)
        case default; call stop_it('init_lncc: bad initlncc='//trim(initlncc))
      endselect

!
!  superimpose something else
!
      select case(initlncc2)
        case('wave-x'); call wave(ampllncc2,f,ilncc,ky=5.)
      endselect
!
!  add floor value if cc_min is set
!
      if(cc_min/=0.) then
        lncc_min=alog(cc_min)
        if(lroot) print*,'set floor value for cc; cc_min=',cc_min
        f(:,:,:,ilncc)=amax1(lncc_min,f(:,:,:,ilncc))
      endif
!
      if(ip==0) print*,xx,yy,zz !(prevent compiler warnings)
    endsubroutine init_lncc
!***********************************************************************
    subroutine dlncc_dt(f,df,uu,glnrho)
!
!  passive scalar evolution
!  calculate dc/dt=-uu.gcc + pscaler_diff*[del2cc + glnrho.gcc]
!
!  20-may-03/axel: coded
!
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: uu,gcc,glnrho
      real, dimension (nx) :: cc,rho,ugcc,diff_op,del2cc
      real, dimension (nx) :: cc1,gcc1,gcc2
      integer :: j
!
      intent(in)  :: f,uu,glnrho
      intent(out) :: df
!
!  identify module and boundary conditions
!
      if (nopscalar) then
        if (headtt.or.ldebug) print*,'not SOLVED: dcc_dt'
      else
        if (headtt.or.ldebug) print*,'SOLVE dcc_dt'
      endif
      if (headtt) call identify_bcs('cc',ilncc)
!
!  gradient of passive scalar
!  allow for possibility to turn off passive scalar
!  without changing file size and recompiling everything.
!
      if (.not. nopscalar) then ! i.e. if (pscalar)
        call grad(f,ilncc,gcc)
        call dot_mn(uu,gcc,ugcc)
!
!  passive scalar equation
!
        if(lhydro) df(l1:l2,m,n,ilncc)=df(l1:l2,m,n,ilncc)-ugcc
!
!  diffusion operator
!
        if (pscalar_diff/=0.) then
          if(headtt) print*,'dlncc_dt: pscalar_diff=',pscalar_diff
          call dot_mn(glnrho,gcc,diff_op)
          call del2(f,ilncc,del2cc)
          diff_op=diff_op+del2cc
          df(l1:l2,m,n,ilncc)=df(l1:l2,m,n,ilncc)+pscalar_diff*diff_op
        endif
!
!  add diffusion of imposed constant gradient of c
!  restrict ourselves (for the time being) to z-gradient only
!  makes sense really only for periodic boundary conditions
!
        do j=1,3
          if (gradC0(j)/=0.) then
            df(l1:l2,m,n,ilncc)=df(l1:l2,m,n,ilncc)-gradC0(j)*uu(:,j)
          endif
        enddo
!
!  tensor diffusion (but keep the isotropic one)
!
        if (tensor_pscalar_diff/=0.) call tensor_diff(f,df,tensor_pscalar_diff,gcc)
!
!  For the timestep calculation, need maximum diffusion
!
        if (lfirst.and.ldt) then
          call max_for_dt(pscalar_diff,maxdiffus)
          call max_for_dt(tensor_pscalar_diff,maxdiffus)
        endif
!
      endif
!
!  diagnostics
!
!  output for double and triple correlators (assume z-gradient of cc)
!  <u_k u_j d_j c> = <u_k c uu.gradlncc>
!
      if (ldiagnos) then
        cc=f(l1:l2,m,n,ilncc)
        rho=exp(f(l1:l2,m,n,ilnrho))
        cc1=rho*abs(cc)
        call dot2_mn(gcc,gcc2); gcc1=sqrt(gcc2)
        if (i_rhoccm/=0) call sum_mn_name(rho*cc/unit_rhocc,i_rhoccm)
        if (i_ccmax/=0) call max_mn_name(cc,i_ccmax)
        if (i_lnccmz/=0) call xysum_mn_name_z(cc,i_lnccmz)
        if (i_ucm/=0) call sum_mn_name(uu(:,3)*cc,i_ucm)
        if (i_uudcm/=0) call sum_mn_name(uu(:,3)*ugcc,i_uudcm)
        if (i_Cz2m/=0) call sum_mn_name(rho*cc*z(n)**2,i_Cz2m)
        if (i_Cz4m/=0) call sum_mn_name(rho*cc*z(n)**4,i_Cz4m)
        if (i_Crmsm/=0) call sum_mn_name((rho*cc)**2,i_Crmsm,lsqrt=.true.)
        if (i_cc1m/=0) call sum_mn_name(cc1   ,i_cc1m)
        if (i_cc2m/=0) call sum_mn_name(cc1**2,i_cc2m)
        if (i_cc3m/=0) call sum_mn_name(cc1**3,i_cc3m)
        if (i_cc4m/=0) call sum_mn_name(cc1**4,i_cc4m)
        if (i_cc5m/=0) call sum_mn_name(cc1**5,i_cc5m)
        if (i_cc6m/=0) call sum_mn_name(cc1**6,i_cc6m)
        if (i_cc7m/=0) call sum_mn_name(cc1**7,i_cc7m)
        if (i_cc8m/=0) call sum_mn_name(cc1**8,i_cc8m)
        if (i_cc9m/=0) call sum_mn_name(cc1**9,i_cc9m)
        if (i_cc10m/=0)call sum_mn_name(cc1**10,i_cc10m)
        if (i_gcc1m/=0) call sum_mn_name(gcc1   ,i_gcc1m)
        if (i_gcc2m/=0) call sum_mn_name(gcc1**2,i_gcc2m)
        if (i_gcc3m/=0) call sum_mn_name(gcc1**3,i_gcc3m)
        if (i_gcc4m/=0) call sum_mn_name(gcc1**4,i_gcc4m)
        if (i_gcc5m/=0) call sum_mn_name(gcc1**5,i_gcc5m)
        if (i_gcc6m/=0) call sum_mn_name(gcc1**6,i_gcc6m)
        if (i_gcc7m/=0) call sum_mn_name(gcc1**7,i_gcc7m)
        if (i_gcc8m/=0) call sum_mn_name(gcc1**8,i_gcc8m)
        if (i_gcc9m/=0) call sum_mn_name(gcc1**9,i_gcc9m)
        if (i_gcc10m/=0)call sum_mn_name(gcc1**10,i_gcc10m)
      endif
!
    endsubroutine dlncc_dt
!***********************************************************************
    subroutine rprint_pscalar(lreset,lwrite)
!
!  reads and registers print parameters relevant for magnetic fields
!
!   6-jul-02/axel: coded
!
      use Sub
!
      integer :: iname,inamez
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        i_rhoccm=0; i_ccmax=0; i_lnccm=0; i_lnccmz=0
        i_ucm=0; i_uudcm=0; i_Cz2m=0; i_Cz4m=0; i_Crmsm=0
        i_cc1m=0; i_cc2m=0; i_cc3m=0; i_cc4m=0; i_cc5m=0
        i_cc6m=0; i_cc7m=0; i_cc8m=0; i_cc9m=0; i_cc10m=0
        i_gcc1m=0; i_gcc2m=0; i_gcc3m=0; i_gcc4m=0; i_gcc5m=0
        i_gcc6m=0; i_gcc7m=0; i_gcc8m=0; i_gcc9m=0; i_gcc10m=0
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'rhoccm',i_rhoccm)
        call parse_name(iname,cname(iname),cform(iname),'ccmax',i_ccmax)
        call parse_name(iname,cname(iname),cform(iname),'lnccm',i_lnccm)
        call parse_name(iname,cname(iname),cform(iname),'ucm',i_ucm)
        call parse_name(iname,cname(iname),cform(iname),'uudcm',i_uudcm)
        call parse_name(iname,cname(iname),cform(iname),'Cz2m',i_Cz2m)
        call parse_name(iname,cname(iname),cform(iname),'Cz4m',i_Cz4m)
        call parse_name(iname,cname(iname),cform(iname),'Crmsm',i_Crmsm)
        call parse_name(iname,cname(iname),cform(iname),'cc1m',i_cc1m)
        call parse_name(iname,cname(iname),cform(iname),'cc2m',i_cc2m)
        call parse_name(iname,cname(iname),cform(iname),'cc3m',i_cc3m)
        call parse_name(iname,cname(iname),cform(iname),'cc4m',i_cc4m)
        call parse_name(iname,cname(iname),cform(iname),'cc5m',i_cc5m)
        call parse_name(iname,cname(iname),cform(iname),'cc6m',i_cc6m)
        call parse_name(iname,cname(iname),cform(iname),'cc7m',i_cc7m)
        call parse_name(iname,cname(iname),cform(iname),'cc8m',i_cc8m)
        call parse_name(iname,cname(iname),cform(iname),'cc9m',i_cc9m)
        call parse_name(iname,cname(iname),cform(iname),'cc10m',i_cc10m)
        call parse_name(iname,cname(iname),cform(iname),'gcc1m',i_gcc1m)
        call parse_name(iname,cname(iname),cform(iname),'gcc2m',i_gcc2m)
        call parse_name(iname,cname(iname),cform(iname),'gcc3m',i_gcc3m)
        call parse_name(iname,cname(iname),cform(iname),'gcc4m',i_gcc4m)
        call parse_name(iname,cname(iname),cform(iname),'gcc5m',i_gcc5m)
        call parse_name(iname,cname(iname),cform(iname),'gcc6m',i_gcc6m)
        call parse_name(iname,cname(iname),cform(iname),'gcc7m',i_gcc7m)
        call parse_name(iname,cname(iname),cform(iname),'gcc8m',i_gcc8m)
        call parse_name(iname,cname(iname),cform(iname),'gcc9m',i_gcc9m)
        call parse_name(iname,cname(iname),cform(iname),'gcc10m',i_gcc10m)
      enddo
!
!  check for those quantities for which we want xy-averages
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'lnccmz',i_lnccmz)
      enddo
!
!  write column where which magnetic variable is stored
!
      if (lwr) then
        write(3,*) 'i_rhoccm=',i_rhoccm
        write(3,*) 'i_ccmax=',i_ccmax
        write(3,*) 'i_lnccm=',i_lnccm
        write(3,*) 'i_ucm=',i_ucm
        write(3,*) 'i_uudcm=',i_uudcm
        write(3,*) 'i_lnccmz=',i_lnccmz
        write(3,*) 'i_Cz2m=',i_Cz2m
        write(3,*) 'i_Cz4m=',i_Cz4m
        write(3,*) 'i_Crmsm=',i_Crmsm
        write(3,*) 'i_cc1m=',i_cc1m
        write(3,*) 'i_cc2m=',i_cc2m
        write(3,*) 'i_cc3m=',i_cc3m
        write(3,*) 'i_cc4m=',i_cc4m
        write(3,*) 'i_cc5m=',i_cc5m
        write(3,*) 'i_cc6m=',i_cc6m
        write(3,*) 'i_cc7m=',i_cc7m
        write(3,*) 'i_cc8m=',i_cc8m
        write(3,*) 'i_cc9m=',i_cc9m
        write(3,*) 'i_cc10m=',i_cc10m
        write(3,*) 'i_gcc1m=',i_gcc1m
        write(3,*) 'i_gcc2m=',i_gcc2m
        write(3,*) 'i_gcc3m=',i_gcc3m
        write(3,*) 'i_gcc4m=',i_gcc4m
        write(3,*) 'i_gcc5m=',i_gcc5m
        write(3,*) 'i_gcc6m=',i_gcc6m
        write(3,*) 'i_gcc7m=',i_gcc7m
        write(3,*) 'i_gcc8m=',i_gcc8m
        write(3,*) 'i_gcc9m=',i_gcc9m
        write(3,*) 'i_gcc10m=',i_gcc10m
        write(3,*) 'ilncc=',ilncc
      endif
!
    endsubroutine rprint_pscalar
!***********************************************************************
    subroutine calc_mpscalar
!
!  calculate mean magnetic field from xy- or z-averages
!
!  14-apr-03/axel: adaped from calc_mfield
!
      use Cdata
      use Sub
!
      logical,save :: first=.true.
      real :: lnccm
!
!  Magnetic energy in horizontally averaged field
!  The bxmz and bymz must have been calculated,
!  so they are present on the root processor.
!
      if (i_lnccm/=0) then
        if (i_lnccmz==0) then
          if(first) print*
          if(first) print*,"NOTE: to get lnccm, lnccmz must also be set in xyaver"
          if(first) print*,"      We proceed, but you'll get lnccm=0"
          lnccm=0.
        else
          lnccm=sqrt(sum(fnamez(:,:,i_lnccmz)**2)/(nz*nprocz))
        endif
        call save_name(lnccm,i_lnccm)
      endif
!
    endsubroutine calc_mpscalar
!***********************************************************************
    subroutine tensor_diff(f,df,tensor_pscalar_diff,gcc)
!
!  reads file
!
!  11-jul-02/axel: coded
!
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, save, dimension (nx,ny,nz,3) :: bunit,hhh
      real, dimension (nx,3,3) :: g
      real, dimension (nx,3) :: gcc
      real, dimension (nx) :: tmp,scr
      real :: tensor_pscalar_diff
      integer :: iy,iz,i,j
      logical, save :: first=.true.
!
!  read H and Bunit arrays and keep them in memory
!
      if(first) then
        open(1,file=trim(directory)//'/bunit.dat',form='unformatted')
        print*,'read bunit.dat with dimension: ',nx,ny,nz,3
        read(1) bunit,hhh
        close(1)
        print*,'read bunit.dat; bunit(1,1,1,1)=',bunit(1,1,1,1)
      endif
!
!  tmp = (Bunit.G)^2 + H.G + Bi*Bj*Gij
!  for details, see tex/mhd/thcond/tensor_der.tex
!
      call dot_mn(bunit,gcc,scr)
      call dot_mn(hhh,gcc,tmp)
      tmp=tmp+scr**2
!
!  calculate Hessian matrix of lncc
!
      call g2ij(f,ilncc,g)
!
!  dot with bi*bj
!
      iy=m-m1+1
      iz=n-n1+1
      do j=1,3
      do i=1,3
        tmp=tmp+bunit(:,iy,iz,i)*bunit(:,iy,iz,j)*g(:,i,j)
      enddo
      enddo
!
!  and add result to the dlncc/dt equation
!
      df(l1:l2,m,n,ilncc)=df(l1:l2,m,n,ilncc)+tensor_pscalar_diff*tmp
!
      first=.false.
    endsubroutine tensor_diff
!***********************************************************************

endmodule Pscalar



