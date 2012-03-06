! $Id$
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: ldensity = .false.
! CPARAM logical, parameter :: lanelastic = .true.
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 4
! COMMUNICATED AUXILIARIES 4
!
! PENCILS PROVIDED rho; lnrho; rho1; glnrho(3); del2rho; del2lnrho
! PENCILS PROVIDED hlnrho(3,3); grho(3); glnrho2
! PENCILS PROVIDED del6lnrho; uij5glnrho(3); uglnrho; ugrho; sglnrho(3)
! PENCILS PROVIDED ekin; transprho; gpp; del6TT; del6lnTT
!
!***************************************************************
module Density
!
  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  logical :: lcalc_glnrhomean=.false.,lupw_lnrho=.false.
  real, dimension (nz,3) :: glnrhomz
!
  include 'density.h'
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
      call farray_register_auxiliary('rhs',irhs,vector=3,communicated=.true.)
!
    endsubroutine register_density
!***********************************************************************
    subroutine initialize_density(f,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  24-nov-02/tony: coded
!
      use EquationOfState, only: select_eos_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
!  Tell the equation of state that we're here and we don't have a
!  variable => isochoric (constant density).
!
      call select_eos_variable('lnrho',-1)
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_density
!***********************************************************************
    subroutine init_lnrho(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      f(:,:,:,ipp)=0.
!
    endsubroutine init_lnrho
!***********************************************************************
    subroutine calc_ldensity_pars(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
  endsubroutine calc_ldensity_pars
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
      call keep_compiler_quiet(lpencil_in)
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
      if (lpencil(i_ekin)) p%ekin=0.0
! pb with noidealgas: del6TT, del6lnTT
       if (lpencil(i_del6TT)) p%del6TT=0.0
       if (lpencil(i_del6lnTT)) p%del6lnTT=0.0
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_pencils_density
!***********************************************************************
    subroutine density_before_boundary(f)
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine density_before_boundary
!***********************************************************************
    subroutine dlnrho_dt(f,df,p)
!
      use Sub, only: identify_bcs
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in) :: f,df,p
!
      if (headtt) then
        call identify_bcs('pp',ipp)
        call identify_bcs('rhs',irhs)
        call identify_bcs('rhs',irhs+1)
        call identify_bcs('rhs',irhs+2)
      endif
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine dlnrho_dt
!***********************************************************************
    subroutine impose_density_floor(f)
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine impose_density_floor
!***********************************************************************
    subroutine read_density_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
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
    subroutine read_density_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_density_run_pars
!***********************************************************************
    subroutine write_density_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_density_run_pars
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
    subroutine dynamical_diffusion(umax)
!   
!  dummy routine
!  
      real, intent(in) :: umax
!  
      call keep_compiler_quiet(umax)
!
    endsubroutine dynamical_diffusion
!***********************************************************************
    subroutine anelastic_after_mn(f, p, df, mass_per_proc)
!
      use Poisson, only: inverse_laplacian
      use Mpicomm, only: initiate_isendrcv_bdry, finalize_isendrcv_bdry
      use Boundcond, only: update_ghosts
      use Sub, only: div, grad
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx,3) :: gpp
      real, dimension (nx) :: phi_rhs_pencil
      real, dimension (1)  :: mass_per_proc
      integer :: j, ju !, l, i
!
!  Set first the boundary conditions on rhs
!
      call update_ghosts(f,irhs,irhs+2)
!
!  Find the divergence of rhs
!
      do n=n1,n2; do m=m1,m2
        call div(f,irhs,phi_rhs_pencil)
        f(l1:l2,m,n,ipp)=phi_rhs_pencil
      enddo; enddo
!
!  get pressure from inverting the Laplacian
!
      if (lperi(3)) then
        call inverse_laplacian(f,f(l1:l2,m1:m2,n1:n2,ipp))
      else
        call inverse_laplacian_z(f(l1:l2,m1:m2,n1:n2,ipp))       
      endif
!
!  refresh the ghost zones for pressure
!
      call update_ghosts(f,ipp)
!
!  Add the pressure gradient term to the NS equation
!
      do n=n1,n2; do m=m1,m2
        call grad(f,ipp,gpp)
        do j=1,3
          ju=j+iuu-1
          df(l1:l2,m,n,ju)=df(l1:l2,m,n,ju)-gpp(:,j)
        enddo
      enddo; enddo
!
    endsubroutine anelastic_after_mn
!***********************************************************************
    subroutine inverse_laplacian_z(phi)
!
!  Solve the Poisson equation by Fourier transforming in the xy-plane and
!  solving the discrete matrix equation in the z-direction.
!
!  02-mar-2012/dintrans+dubuffet: coded
!
      use Fourier, only: fourier_transform_xy
      use General, only: tridag
!
      real, dimension (nx,ny,nz) :: phi, b1
      real, dimension (nz) :: a_tri, b_tri, c_tri, r_tri, u_tri
      real :: k2
      integer :: ikx, iky
      logical :: err
!
!  Identify version.
!
      if (lroot .and. ip<10) call svn_id( &
          '$Id$')
!
!  The right-hand-side of the Poisson equation is purely real.
!
      b1 = 0.0
!
!  Forward transform (to k-space).
!
      call fourier_transform_xy(phi,b1)
!
!  Solve for discrete z-direction
!
      do iky=1,ny
        do ikx=1,nx
          if ((kx_fft(ikx)==0.0) .and. (ky_fft(iky)==0.0)) then
            phi(ikx,iky,:) = 1.0
          else
            k2=kx_fft(ikx)**2+ky_fft(iky)**2
            a_tri=1.0/dz**2
            b_tri=-2.0/dz**2-k2
            c_tri=1.0/dz**2
            r_tri=phi(ikx,iky,:)
! P = 0 (useful to test the Poisson solver)
!            b_tri(1)=1.  ; c_tri(1)=0.  ; r_tri(1)=0.
!            b_tri(nz)=1. ; a_tri(nz)=0. ; r_tri(nz)=0.
! dP/dz = 0
            c_tri(1)=c_tri(1)+a_tri(1)
            a_tri(nz)=a_tri(nz)+c_tri(nz)
!
!            call tridag(a_tri,b_tri,c_tri,r_tri,u_tri,err)
            call mytridag(a_tri,b_tri,c_tri,r_tri,u_tri,err)
            phi(ikx,iky,:)=u_tri
          endif
        enddo
      enddo
!
      do iky=1,ny
        do ikx=1,nx
          if ((kx_fft(ikx)==0.0) .and. (ky_fft(iky)==0.0)) then
            b1(ikx,iky,:) = 0.0
          else
            k2=kx_fft(ikx)**2+ky_fft(iky)**2
            a_tri=1.0/dz**2
            b_tri=-2.0/dz**2-k2
            c_tri=1.0/dz**2
            r_tri=b1(ikx,iky,:)
! P = 0 (useful to test the Poisson solver)
!            b_tri(1)=1.  ; c_tri(1)=0.  ; r_tri(1)=0.
!            b_tri(nz)=1. ; a_tri(nz)=0. ; r_tri(nz)=0.
! dP/dz = 0
            c_tri(1)=c_tri(1)+a_tri(1)
            a_tri(nz)=a_tri(nz)+c_tri(nz)
!
!            call tridag(a_tri,b_tri,c_tri,r_tri,u_tri,err)
            call mytridag(a_tri,b_tri,c_tri,r_tri,u_tri,err)
            b1(ikx,iky,:)=u_tri
          endif
        enddo
      enddo
!
!  Inverse transform (to real space).
!
      call fourier_transform_xy(phi,b1,linv=.true.)
!
    endsubroutine inverse_laplacian_z
!***********************************************************************
    subroutine mytridag(a,b,c,r,u,err)
!
!  Solves a tridiagonal system.
!
!  01-apr-03/tobi: from Numerical Recipes (p42-43).
!
!  | b1 c1 0  ...            | | u1   |   | r1   |
!  | a2 b2 c2 ...            | | u2   |   | r2   |
!  | 0  a3 b3 c3             | | u3   | = | r3   |
!  |          ...            | | ...  |   | ...  |
!  |          an-1 bn-1 cn-1 | | un-1 |   | rn-1 |
!  |          0    a_n  b_n  | | un   |   | rn   |
!
      implicit none
      real, dimension(:), intent(in) :: a,b,c,r
      real, dimension(:), intent(out) :: u
      real, dimension(size(b)) :: gam
      logical, intent(out), optional :: err
      integer :: n,j
      real :: bet
!
      if (present(err)) err=.false.
      n=size(b)
      bet=b(1)
      if (bet==0.0) then
        print*,'tridag: Error at code stage 1'
        if (present(err)) err=.true.
        return
      endif
!
      u(1)=r(1)/bet
      do j=2,n
!        print*, 'j=', j
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if (bet==0.0) then
          print*,'tridag: Error at code stage 2'
          if (present(err)) err=.true.
          return
        endif
        u(j)=(r(j)-a(j)*u(j-1))/bet
      enddo
!
      do j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
      enddo
!
    endsubroutine mytridag
!***********************************************************************
endmodule Density
