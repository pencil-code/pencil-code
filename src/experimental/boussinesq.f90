! $Id$
!
! 23-mar-2012/dintrans: coded
!
! Solve the Poisson equation for pressure when using the Boussinesq
! approximation (the so-called ``projection method'').
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: ldensity = .false.
! CPARAM logical, parameter :: lanelastic = .false.
! CPARAM logical, parameter :: lboussinesq = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 1
! COMMUNICATED AUXILIARIES 1
!
! PENCILS PROVIDED rho; lnrho; rho1; glnrho(3); del2rho; del2lnrho
! PENCILS PROVIDED hlnrho(3,3); grho(3); glnrho2
! PENCILS PROVIDED del6lnrho; uij5glnrho(3); uglnrho; ugrho; sglnrho(3)
! PENCILS PROVIDED ekin; transprho
!
! PENCILS PROVIDED glnrhos(3)
! PENCILS PROVIDED totenergy_rel
!
!***************************************************************
module Density
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  logical :: lcalc_lnrhomean=.false.,lcalc_glnrhomean=.false.,lupw_lnrho=.false., &
             lremove_mean_temperature=.false.
!
  logical :: lwrite_debug=.false.
!
  real, dimension (nz) :: lnrhomz,glnrhomz
  real :: dx_2, dz_2
  real, pointer :: chi
!
  include '../density.h'
  !integer :: pushpars2c, pushdiags2c  ! should be procedure pointer (F2003)
!
  integer :: iorder_z=4
  real, dimension(3) :: beta_glnrho_global=0.0, beta_glnrho_scaled=0.0
!
  namelist /density_run_pars/ iorder_z, lwrite_debug, lremove_mean_temperature
!
  contains
!***********************************************************************
    subroutine register_density()
!
      use FArrayManager, only: farray_register_auxiliary
!
      if (lroot) call svn_id( &
          "$Id$")
!
      call farray_register_auxiliary('pp',ipp,communicated=.true.)
      if (lsphere_in_a_box) lgravr=.true.
!
    endsubroutine register_density
!***********************************************************************
    subroutine initialize_density(f)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  24-nov-02/tony: coded
!
      use EquationOfState, only: select_eos_variable
      use DensityMethods, only: initialize_density_methods
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!  Boussinesq not implemented for entropy
!
      if (lenergy) then
        if (lentropy) call fatal_error("initialize_density", &
          "Boussinesq not implemented for entropy")
!
!  Boussinesq only implemented for ltemperature_nolog
!
        if (.not.ltemperature_nolog) call fatal_error("initialize_density", &
          "Boussinesq is only implemented for ltemperature_nolog")
!
      endif
!
!  Tell the equation of state that we're here and we don't have a
!  variable => isochoric (constant density).
!
      call select_eos_variable('lnrho',-1)   ! the same for rho?
!
!  For the implicit solver
!
      dx_2=1./dx**2
      dz_2=1./dz**2
!
      call initialize_density_methods
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_density
!***********************************************************************
    subroutine init_lnrho(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      f(:,:,:,ipp)=1.
!
!  Test of the Poisson solver
!
!      do n=n1,n2; do m=m1,m2
!        f(l1:l2,m,n,ipp)=-2.*sin(x(l1:l2))*cos(z(n))
!      enddo; enddo
!      write(10) f(l1:l2,m1:m2,n1:n2,ipp)
!      if (iorder_z==2) then
!        call inverse_laplacian_z_2nd(f(l1:l2,m1:m2,n1:n2,ipp))
!      else
!        call inverse_laplacian_z(f(l1:l2,m1:m2,n1:n2,ipp))
!      endif
!      write(11) f(l1:l2,m1:m2,n1:n2,ipp)
!
    endsubroutine init_lnrho
!***********************************************************************
    subroutine read_density_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=density_run_pars, IOSTAT=iostat)
!
    endsubroutine read_density_run_pars
!***********************************************************************
    subroutine read_density_init_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_density_init_pars
!***********************************************************************
    subroutine write_density_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_density_init_pars
!***********************************************************************
    subroutine write_density_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=density_run_pars)
!
    endsubroutine write_density_run_pars
!***********************************************************************
    subroutine density_after_boundary(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
  endsubroutine density_after_boundary
!***********************************************************************
    subroutine pencil_criteria_density()
!
!  All pencils that the Density module depends on are specified here.
!
!  20-11-04/anders: coded
!
    endsubroutine pencil_criteria_density
!***********************************************************************
    subroutine pencil_interdep_density(lpencil_in)
!
!  Interdependency among pencils from the Density module is specified here.
!
!  20-11-04/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      if (lpencil_in(i_ekin)) lpencil_in(i_u2)=.true.
!
    endsubroutine pencil_interdep_density
!***********************************************************************
    subroutine calc_pencils_density(f,p)
!
!  Calculate Density pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  20-11-04/anders: coded
!
      use EquationOfState, only: lnrho0, rho0
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
! rho
      if (lpencil(i_rho)) p%rho=rho0
! lnrho
      if (lpencil(i_lnrho)) p%lnrho=lnrho0
! rho1
      if (lpencil(i_rho1)) p%rho1=1/rho0
! glnrho
      if (lpencil(i_glnrho)) p%glnrho=0.0
! grho
      if (lpencil(i_grho)) p%grho=0.0
! del6lnrho
      if (lpencil(i_del6lnrho)) p%del6lnrho=0.0
! hlnrho
      if (lpencil(i_hlnrho)) p%hlnrho=0.0
! sglnrho
      if (lpencil(i_sglnrho)) p%sglnrho=0.0
! uglnrho
      if (lpencil(i_uglnrho)) p%uglnrho=0.0
! ugrho
      if (lpencil(i_ugrho)) p%ugrho=0.0
! uij5glnrho
      if (lpencil(i_uij5glnrho)) p%uij5glnrho=0.0
! ekin
      if (lpencil(i_ekin)) p%ekin=0.5*p%u2
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_pencils_density
!***********************************************************************
    subroutine density_before_boundary(f)
!
      use Sub, only: remove_mean
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!     
      if (lremove_mean_temperature) call remove_mean(f,iTT)
!
    endsubroutine density_before_boundary
!***********************************************************************
    subroutine dlnrho_dt(f,df,p)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in) :: f,df,p
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine dlnrho_dt
!***********************************************************************
    subroutine density_after_timestep(f,df,dtsub)
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real :: dtsub
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(dtsub)
!
    endsubroutine density_after_timestep
!***********************************************************************
    subroutine split_update_density(f)
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine split_update_density
!***********************************************************************
    subroutine impose_density_floor(f)
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine impose_density_floor
!***********************************************************************
    subroutine rprint_density(lreset,lwrite)
!
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
      call keep_compiler_quiet(lreset)
!
    endsubroutine rprint_density
!***********************************************************************
    subroutine get_slices_density(f,slices)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_density
!***********************************************************************
    subroutine get_slices_pressure(f,slices)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_pressure
!***********************************************************************
    subroutine get_init_average_density(f,init_average_density)
!
!  10-dec-09/piyali: added to pass initial average density
!
    real, dimension (mx,my,mz,mfarray):: f
    real:: init_average_density
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(init_average_density)
!
    endsubroutine get_init_average_density
!***********************************************************************
    subroutine anelastic_after_mn(f, p, df, mass_per_proc)
!
!  14-dec-09/dintrans: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension(1) :: mass_per_proc
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
      call keep_compiler_quiet(mass_per_proc)
!
    endsubroutine anelastic_after_mn
!***********************************************************************
    subroutine dynamical_diffusion(uc)
!   
!  dummy routine
!  
      real, intent(in) :: uc
!  
      call keep_compiler_quiet(uc)
!
    endsubroutine dynamical_diffusion
!***********************************************************************
    subroutine boussinesq(f)
!
!  12-may-12/MR: factors dt removed; updating of ghosts zones 
!                for non-periodicity in z direction added
!  15-may-12/dintrans: the updating of ghost zones before inverting the laplacian 
!                is not needed as vertical BCs are hard-coded in the linear solver
!
      use Poisson, only: inverse_laplacian
      use Sub, only: div, grad
      use Boundcond, only: update_ghosts
      use SharedVariables, only: get_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3) :: gpp
      real, dimension (nx) :: phi_rhs_pencil
      integer :: j, ju, ierr
      real, pointer, save :: Pr
      logical, save :: l1st=.true.
!
      if (lviscosity) then
        call update_ghosts(f,iuu,iuu+2)
      else
        call update_ghosts(f)
        if (l1st) then
          call get_shared_variable('Pr',Pr,ierr)
          if (ierr/=0) call fatal_error('implicit_diffusion',&
              'pb to get Pr')
          print*, 'get Pr=', Pr
          l1st=.false.
        endif
!
!  Implicit advance of both viscous and radiative diffusion terms
!
        do j=1,3
          ju=j+iuu-1
          if (nprocz>1) then
            call implicit_diffusion_MPI(f,ju,Pr)
          else
            call implicit_diffusion(f,ju,Pr)
          endif
        enddo
        if (nprocz>1) then
          call implicit_diffusion_MPI(f,iTT,1.)
        else
          call implicit_diffusion(f,iTT,1.)
        endif
      endif
!
!  Find the divergence of uu
!
      do n=n1,n2; do m=m1,m2
        call div(f,iuu,phi_rhs_pencil)
        f(l1:l2,m,n,ipp)=phi_rhs_pencil
      enddo; enddo
      if (lwrite_debug) write(31) f(l1:l2,m1:m2,n1:n2,ipp)
!
      if (lperi(3)) then
        call inverse_laplacian(f(l1:l2,m1:m2,n1:n2,ipp))
      else
        if (iorder_z==2) then
          call inverse_laplacian_z_2nd(f(l1:l2,m1:m2,n1:n2,ipp))
        else
          call inverse_laplacian_z(f(l1:l2,m1:m2,n1:n2,ipp))
        endif
      endif
      if (lwrite_debug) write(32) f(l1:l2,4,n1:n2,ipp)
!
!  refresh the ghost zones for the new pressure
!
      call update_ghosts(f,ipp)
!
!  Correct the velocity field with the gradient term
!
      do n=n1,n2; do m=m1,m2
        call grad(f,ipp,gpp)
        do j=1,3
          ju=j+iuu-1
          f(l1:l2,m,n,ju)=f(l1:l2,m,n,ju)-gpp(:,j)
        enddo
      enddo; enddo
!      f(:,:,:,ipp)=f(:,:,:,ipp)/dt
!
    endsubroutine boussinesq
!***********************************************************************
    subroutine inverse_laplacian_z(phi,bctype,bval_bot,bval_top)
!
!  10-avr-2012/dintrans: coded
!         2017/MR: optional parameters for various BCs added
!                  (not yet used)
!  Fourth-order in the vertical direction that uses pendag().
!  Note: the (kx=0,ky=0) mode is computed using a Green function.
!
      use Fourier, only: fourier_transform_xy, kx_fft, ky_fft
      use Mpicomm, only: transp_xz, transp_zx
      use General, only: pendag
!
      real,                   dimension(nx,ny,nz),      intent(INOUT):: phi
      character(LEN=labellen),dimension(2),    optional,intent(IN)   :: bctype
      real,                   dimension(nx,ny),optional,intent(IN)   :: bval_bot,bval_top

      real, dimension(nx,ny,nz) :: b1
      !real, dimension(nx,ny) :: bval_bot_hat, bval_top_hat

      integer, parameter :: nxt = nx / nprocz
      real, dimension(nzgrid,nxt) :: phit, b1t
      real, dimension (nzgrid) :: a_tri, b_tri, c_tri, d_tri, e_tri, r_tri, u_tri
!
      logical, save :: l1st = .true.
      real,    save :: dz2h, dz_2
      integer, save :: ikx0, iky0
!
      complex, dimension(nzgrid) :: cz
!
      integer :: ikx, iky, iz, iz1
      real    :: ky2, k2
!
!  Initialize the array wsave and other constants for future use.
!
      if (l1st) then
        dz_2 = 1./dz**2
        dz2h = 0.5 * dz * dz
        ikx0 = ipz * nxt
        iky0 = ipy * ny
        l1st = .false.
      endif
!
!  Forward transform in xy
!
      b1 = 0.
      call fourier_transform_xy(phi, b1)
!
      !bval_bot_hat = 0.; bval_top_hat = 0.
      !call fourier_transform_xy(bval_bot, bval_bot_hat)
      !call fourier_transform_xy(bval_top, bval_top_hat)
!
!  Convolution in z
!
      do iky = 1, ny
        call transp_xz(phi(:,iky,:),  phit)
        call transp_xz(b1(:,iky,:), b1t)
        ky2 = ky_fft(iky0+iky)**2
        do ikx = 1, nxt
          k2 = kx_fft(ikx0+ikx)**2+ky2
          if (k2/=0.0) then
            a_tri(:)=-1./12.*dz_2
            b_tri(:)=+4./3.*dz_2
            c_tri(:)=-5./2.*dz_2-k2
            d_tri(:)=+4./3.*dz_2
            e_tri(:)=-1./12.*dz_2
            d_tri(1)=2.*d_tri(1)
            e_tri(1)=2.*e_tri(1)
            e_tri(2)=2.*e_tri(2)
            a_tri(nzgrid)=2.*a_tri(nzgrid)
            b_tri(nzgrid)=2.*b_tri(nzgrid)
            a_tri(nzgrid-1)=2.*a_tri(nzgrid-1)
!
            r_tri=phit(:,ikx)
            call pendag(nzgrid,a_tri,b_tri,c_tri,d_tri,e_tri,r_tri,u_tri)
            phit(:,ikx)=u_tri
!
            r_tri=b1t(:,ikx)
            call pendag(nzgrid,a_tri,b_tri,c_tri,d_tri,e_tri,r_tri,u_tri)
            b1t(:,ikx)=u_tri
          else
            do iz = 1, nzgrid 
              cz(iz) = 0.5*dz2h*(                       &
                cmplx(phit(1,ikx),b1t(1,ikx))*abs(iz-1) &
               +cmplx(phit(nzgrid,ikx),b1t(nzgrid,ikx))*abs(iz-nzgrid))
              do iz1 = 2, nzgrid-1
                cz(iz) = cz(iz) + cmplx(phit(iz1,ikx), b1t(iz1,ikx)) * &
                         abs(iz - iz1) * dz2h
              enddo
            enddo
            phit(:,ikx) = real(cz(1:nzgrid))
            b1t(:,ikx)  = aimag(cz(1:nzgrid))
          endif
        enddo
        call transp_zx(phit, phi(:,iky,:))
        call transp_zx(b1t, b1(:,iky,:))
      enddo
!
!  Inverse transform in xy
!
      call fourier_transform_xy(phi,b1,linv=.true.)
!
      return
!
    endsubroutine inverse_laplacian_z
!***********************************************************************
    subroutine inverse_laplacian_z_2nd(phi)
!
!  19-mar-2012/dintrans: coded
!  Second-order version in the vertical direction that uses tridag().
!  Note: the (kx=0,ky=0) mode is computed using a Green function.
!
      use Fourier, only: fourier_transform_xy, kx_fft, ky_fft
      use Mpicomm, only: transp_xz, transp_zx
      use General, only: tridag
!
      real, dimension(nx,ny,nz) :: phi, b1
      integer, parameter :: nxt = nx / nprocz
      real, dimension(nzgrid,nxt) :: phit, b1t
      real, dimension (nzgrid) :: a_tri, b_tri, c_tri, r_tri, u_tri
!
      logical, save :: l1st = .true.
      real,    save :: dz2h, dz_2
      integer, save :: ikx0, iky0
!
      complex, dimension(nzgrid) :: cz
!
      integer :: ikx, iky, iz, iz1
      real    :: ky2, k2
!
!  Initialize the array wsave and other constants for future use.
!
      if (l1st) then
        dz_2 = 1./dz**2
        dz2h = 0.5 * dz * dz
        ikx0 = ipz * nxt
        iky0 = ipy * ny
        l1st = .false.
      endif
!
!  Forward transform in xy
!
      b1 = 0.
      call fourier_transform_xy(phi, b1)
!
!  Convolution in z
!
      do iky = 1, ny
        call transp_xz(phi(:,iky,:),  phit)
        call transp_xz(b1(:,iky,:), b1t)
        ky2 = ky_fft(iky0+iky)**2
        do ikx = 1, nxt
          k2 = kx_fft(ikx0+ikx)**2+ky2
          if (k2/=0.0) then
            c_tri(:)=dz_2
            b_tri(:)=-2.*dz_2-k2
            a_tri(:)=dz_2
            c_tri(1)=2.*c_tri(1)
            a_tri(nzgrid)=2.*a_tri(nzgrid)
!
            r_tri=phit(:,ikx)
            call tridag(a_tri,b_tri,c_tri,r_tri,u_tri)
            phit(:,ikx)=u_tri
!
            r_tri=b1t(:,ikx)
            call tridag(a_tri,b_tri,c_tri,r_tri,u_tri)
            b1t(:,ikx)=u_tri
          else
            do iz = 1, nzgrid 
              cz(iz) = 0.5*dz2h*(                       &
                cmplx(phit(1,ikx),b1t(1,ikx))*abs(iz-1) &
               +cmplx(phit(nzgrid,ikx),b1t(nzgrid,ikx))*abs(iz-nzgrid))
              do iz1 = 2, nzgrid-1
                cz(iz) = cz(iz) + cmplx(phit(iz1,ikx), b1t(iz1,ikx)) * &
                         abs(iz - iz1) * dz2h
              enddo
            enddo
            phit(:,ikx) = real(cz(1:nzgrid))
            b1t(:,ikx)  = aimag(cz(1:nzgrid))
          endif
        enddo
        call transp_zx(phit, phi(:,iky,:))
        call transp_zx(b1t, b1(:,iky,:))
      enddo
!
!  Inverse transform in xy
!
      call fourier_transform_xy(phi,b1,linv=.true.)
!
      return
!
    endsubroutine inverse_laplacian_z_2nd
!***********************************************************************
    subroutine implicit_diffusion(f,ivar,cdiff)
!
!  06-June-2012/dintrans: coded
!  2-D ADI scheme for a laplacian-like diffusion term and a constant
!  diffusion coefficient cdiff. The ADI scheme is of Yakonov's form:
!
!    (1-dt/2*Lamba_x)*T^(n+1/2) = T^n + Lambda_x(T^n) + Lambda_z(T^n)
!    (1-dt/2*Lamba_z)*T^(n+1)   = T^(n+1/2)
!
!  where Lambda_x and Lambda_z denote diffusion operators.
!  Note: this form is more adapted for a parallelisation compared the 
!  Peaceman & Rachford one.
!
      use General, only: tridag, cyclic
!
      implicit none
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,mz) :: TT, finter
      real, dimension(nx)    :: ax, bx, cx, rhsx
      real, dimension(nz)    :: az, bz, cz, rhsz
      real    :: aalpha, bbeta
      integer :: l, n, ivar
      real    :: cdiff
!
      TT=f(:,4,:,ivar)
!
!  rows dealt implicitly
!
      ax(:)=-cdiff*dt*dx_2/2.
      bx(:)=1.+cdiff*dt*dx_2
      cx(:)=ax
      aalpha=cx(nx) ; bbeta=ax(1)  ! x-direction periodic
      do n=n1,n2
        rhsx=TT(l1:l2,n)+     &
             cdiff*dt*dz_2/2.*(TT(l1:l2,n+1)-2.*TT(l1:l2,n)+TT(l1:l2,n-1))
        rhsx=rhsx+cdiff*dt*dx_2/2.* &
             (TT(l1+1:l2+1,n)-2.*TT(l1:l2,n)+TT(l1-1:l2-1,n))
        call cyclic(ax,bx,cx,aalpha,bbeta,rhsx,finter(l1:l2,n),nx)
      enddo
!
!  columns dealt implicitly
!
      az(:)=-cdiff*dt*dz_2/2.
      bz(:)=1.+cdiff*dt*dz_2
      cz(:)=az
      if (ivar.eq.iTT .or. ivar.eq.iuz) then
        bz(1)=1.  ; cz(1)=0.  ; rhsz(1)=0.   ! T = uz = 0
        bz(nz)=1. ; az(nz)=0. ; rhsz(nz)=0.  ! T = uz = 0
      else
        cz(1)=2.*cz(1)    ! ux' = 0
        az(nz)=2.*az(nz)  ! ux' = 0
      endif
      do l=l1,l2
        rhsz=finter(l,n1:n2)
        call tridag(az,bz,cz,rhsz,f(l,4,n1:n2,ivar))
      enddo
!
    endsubroutine implicit_diffusion
!***********************************************************************
    subroutine implicit_diffusion_MPI(f,ivar,cdiff)
!
!  06-June-2012/dintrans: coded
!  parallel version of implicit_diffusion
!
      use General, only: tridag, cyclic
      use Mpicomm, only: transp_xz, transp_zx
!
      implicit none
!
      integer, parameter :: nxt=nx/nprocz
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,mz)      :: TT, finter
      real, dimension(nzgrid,nxt) :: fintert, wtmp
      real, dimension(nx)         :: ax, bx, cx, rhsx
      real, dimension(nzgrid)     :: az, bz, cz, rhsz
      real    :: aalpha, bbeta
      integer :: l, n, ivar
      real    :: cdiff
!
      TT=f(:,4,:,ivar)
!
!  rows dealt implicitly
!
      ax(:)=-cdiff*dt*dx_2/2.
      bx(:)=1.+cdiff*dt*dx_2
      cx(:)=ax
      aalpha=cx(nx) ; bbeta=ax(1)  ! x-direction periodic
      do n=n1,n2
        rhsx=TT(l1:l2,n)+     &
             cdiff*dt*dz_2/2.*(TT(l1:l2,n+1)-2.*TT(l1:l2,n)+TT(l1:l2,n-1))
        rhsx=rhsx+cdiff*dt*dx_2/2.* &
             (TT(l1+1:l2+1,n)-2.*TT(l1:l2,n)+TT(l1-1:l2-1,n))
        call cyclic(ax,bx,cx,aalpha,bbeta,rhsx,finter(l1:l2,n),nx)
      enddo
!
!  columns dealt implicitly
!
      az(:)=-cdiff*dt*dz_2/2.
      bz(:)=1.+cdiff*dt*dz_2
      cz(:)=az
      if (ivar.eq.iTT .or. ivar.eq.iuz) then
        bz(1)=1.  ; cz(1)=0.  ; rhsz(1)=0.               ! T = uz = 0
        bz(nzgrid)=1. ; az(nzgrid)=0. ; rhsz(nzgrid)=0.  ! T = uz = 0
      else
        cz(1)=2.*cz(1)            ! ux' = 0
        az(nzgrid)=2.*az(nzgrid)  ! ux' = 0
      endif
!
      call transp_xz(finter(l1:l2,n1:n2), fintert)
      do l=1,nxt
        rhsz=fintert(:,l)
        call tridag(az,bz,cz,rhsz,wtmp(:,l))
      enddo
      call transp_zx(wtmp, f(l1:l2,4,n1:n2,ivar))
!
    endsubroutine implicit_diffusion_MPI
!***********************************************************************
    function mean_density(f)
!
!  Return mean density as rho0 from eos
!
!  1-mar-15/MR: derived from mean_density in density
!
      use EquationOfState, only: rho0
!
      real :: mean_density
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(in) :: f
!
      mean_density=rho0
!
      call keep_compiler_quiet(f)
!
    endfunction mean_density
!***********************************************************************
    subroutine update_char_vel_density(f)
!
!  Updates characteristic veelocity for slope-limited diffusion.
!  Most likely not yet a good method.
!
!  21-oct-15/MR: coded
!
      use EquationOfState, only: rho0
!
      real, dimension(mx,my,mz,mfarray), intent(INOUT) :: f
!
      if (lslope_limit_diff) f(2:mx-2,2:my-2,2:mz-2,iFF_char_c) &
                            =f(2:mx-2,2:my-2,2:mz-2,iFF_char_c) + rho0**2
!
    endsubroutine update_char_vel_density
!***********************************************************************
    subroutine impose_density_ceiling(f)
!
!  Dummy routine.
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f

      call keep_compiler_quiet(f)

    endsubroutine impose_density_ceiling
!***********************************************************************
    subroutine calc_diagnostics_density(f,p)

      real, dimension (mx,my,mz,mfarray) :: f
      type(pencil_case) :: p

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)

    endsubroutine calc_diagnostics_density
!***********************************************************************
    subroutine pushpars2c(p_par)

    integer, parameter :: n_pars=0
    integer(KIND=ikind8), dimension(:) :: p_par

    call keep_compiler_quiet(p_par)

    endsubroutine pushpars2c
!***********************************************************************
    subroutine pushdiags2c(p_diag)

    integer(KIND=ikind8), dimension(:) :: p_diag

    call keep_compiler_quiet(p_diag)

    endsubroutine pushdiags2c
!***********************************************************************
endmodule Density
