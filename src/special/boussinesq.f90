! $Id: van_der_pol.f90 17798 2011-10-31 14:52:44Z boris.dintrans $
!
!  Solve the Boussinesq equations
!
!  05-mar-2012/dintrans: coded
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 2
! COMMUNICATED AUXILIARIES 1
!
!***************************************************************

module Special
!
  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include '../special.h'
!
  integer :: idivu
  real :: Ra_=0.0
!
  namelist /special_run_pars/ Ra_
!
  contains
!
!***********************************************************************
    subroutine register_special()
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!  6-oct-03/tony: coded
!
      use FArrayManager, only: farray_register_auxiliary
!
!  Identify CVS/SVN version information.
!
      if (lroot) call svn_id( &
           "$Id: boussinesq.f90 17798 2011-10-31 14:52:44Z boris.dintrans $")
!
      call farray_register_auxiliary('pp',ipp,communicated=.true.)
      call farray_register_auxiliary('divu',idivu,communicated=.false.)
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f,lstarting)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  06-oct-03/tony: coded
!
      use EquationOfState, only: select_eos_variable
      use Poisson, only: inverse_laplacian
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
!  Initialize any module variables which are parameter dependent
!
!      call select_eos_variable('rho',irho)
!
!      do n=n1,n2; do m=m1,m2
!        f(l1:l2,m,n,ipp)=-2.*sin(x(l1:l2))*cos(z(n))
!      enddo; enddo
!      write(10) f(l1:l2,m1:m2,n1:n2,ipp)
!      call inverse_laplacian(f,f(l1:l2,m1:m2,n1:n2,ipp))
!      call inverse_laplacian_z(f,f(l1:l2,m1:m2,n1:n2,ipp),dt_)
!      print*, 'write test results in binary file #11'
!      write(11) f(l1:l2,m1:m2,n1:n2,ipp)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f)
!
!  initialise special condition; called from start.f90
!  06-oct-2003/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(inout) :: f
!
      f(:,:,:,ipp)=1.
!
    endsubroutine init_special
!***********************************************************************
    subroutine read_special_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
!  write name list
!
      if (present(iostat)) then
        read(unit,NML=special_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=special_run_pars,ERR=99)
      endif
!
99  endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
      integer, intent(in) :: unit
!
!  write name list
!
      write(unit,NML=special_run_pars)
!   
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine special_after_timestep(f,df,dt_)
!
      use Poisson, only: inverse_laplacian
      use Boundcond, only: update_ghosts
      use Sub, only: div, grad, gij, del2
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx,3) :: gpp
      real, dimension (nx) :: phi_rhs_pencil
      real, dimension (nx,ny,nz) :: tmp_div
      real, intent(in) :: dt_
      real :: s, err
      integer :: i,j, ju, iter, l 
!
!  Set first the boundary conditions on rhs
!
!      bcx(1:3)='a'; bcz(1:3)='a'
      call update_ghosts(f,iuu,iuu+2)
!
!  Find the divergence of rhs
!
      do n=n1,n2; do m=m1,m2
        call div(f,iuu,phi_rhs_pencil)
        f(l1:l2,m,n,ipp)=phi_rhs_pencil/dt_
      enddo; enddo
!      write(31) f(l1:l2,m1:m2,n1:n2,ipp)
!
      if (lperi(3)) then
        call inverse_laplacian(f,f(l1:l2,m1:m2,n1:n2,ipp))
      else
        call inverse_laplacian_z(f(l1:l2,m1:m2,n1:n2,ipp))
      endif
!      write(32) f(l1:l2,4,n1:n2,ipp)
!
!  refresh the ghost zones for pressure
!
      call update_ghosts(f,ipp)
!
!  verify that del2 P = div(u)/dt_
!
!      do n=n1,n2; do m=m1,m2
!        call del2(f,ipp,phi_rhs_pencil)
!        f(l1:l2,m,n,idivu)=phi_rhs_pencil
!      enddo; enddo
!      write(34) f(l1:l2,m1:m2,n1:n2,idivu)
!
!  Euler-advance of the velocity field with just the pressure gradient term
!
      do n=n1,n2; do m=m1,m2
        call grad(f,ipp,gpp)
        do j=1,3
          ju=j+iuu-1
          f(l1:l2,m,n,ju)=f(l1:l2,m,n,ju)-dt_*gpp(:,j)
        enddo
      enddo; enddo
!
!  fill in the auxiliary array idivu with divergence of new velocity
!
!      bcx(1:3)='p'; bcz(1:2)='s'; bcz(3)='a'
      call update_ghosts(f,iuu,iuu+2)
      do n=n1,n2; do m=m1,m2
        call div(f,iuu,phi_rhs_pencil)
        f(l1:l2,m,n,idivu)=phi_rhs_pencil
      enddo; enddo
!      write(33) f(l1:l2,m1:m2,n1:n2,idivu)
!
    endsubroutine special_after_timestep
!***********************************************************************
    subroutine special_calc_hydro(f,df,p)
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      if (ltemperature) &
        df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz) + Ra_*f(l1:l2,m,n,iTT)
!
    endsubroutine special_calc_hydro
!***********************************************************************
    subroutine special_calc_entropy(f,df,p)
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      if (ltemperature) &
        df(l1:l2,m,n,iTT) = df(l1:l2,m,n,iTT) + f(l1:l2,m,n,iuz)
!
    endsubroutine special_calc_entropy
!***********************************************************************
    subroutine inverse_laplacian_z(phi)
!
!  19-mar-2012/dintrans: coded
!
      use Fourier, only: fourier_transform_xy
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
            b_tri(:)=-2.0*dz_2-k2
            a_tri(:)=dz_2
!
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
            cz = 0.
            do iz = 1, nzgrid; do iz1 = 1, nzgrid
              cz(iz) = cz(iz) + cmplx(phit(iz1,ikx), b1t(iz1,ikx)) * &
                                abs(iz - iz1) * dz2h
            enddo; enddo
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
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include 'special_dummies.inc'
!********************************************************************

endmodule Special

