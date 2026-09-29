! $Id$
!
!  Special module that damps the velocity field outside a user-specified cuboid.
!  Useful in cases where you do not want waves reflected from the boundary to
!  interfere with the interior of the domain.
!
!  22-Oct-2024/Kishore: Added.
!  29-Sep-2026/Kishore: ported to GPU
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of special_dummies.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 1
!
!***************************************************************
!
module Special
!
  use Cdata
  use Sub, only: step
  use Quiet
  use Messages, only: svn_id, fatal_error, not_implemented
!
  implicit none
!
  include '../special.h'
!
  real :: x_1=impossible, x_2=impossible !corners of the cuboid outside which things should be damped
  real :: y_1=impossible, y_2=impossible
  real :: z_1=impossible, z_2=impossible
  real :: tau=1 !timescale over which the velocity should be damped
  real :: w=0 !width of the step function
  logical :: ldamp_rho=.false. !whether to damp the density to its initial profile
  logical :: ldamp_ss=.false. !whether to damp the entropy to its initial profile
  character (len=labellen) :: far_field_type='initial' !how to determine the profiles to which rho and ss are damped.
!
  integer :: itauinv=0
  real, dimension (mz) :: rho_prof, ss_prof
  logical :: lprof_from_initial=.true.
!
! run parameters
  namelist /special_run_pars/ &
    x_1, x_2, y_1, y_2, z_1, z_2, tau, w, ldamp_rho, ldamp_ss, far_field_type
!
  contains
!***********************************************************************
    subroutine register_special
!
      use FArrayManager, only: farray_register_auxiliary
!
      call farray_register_auxiliary('tauinv',itauinv)
!
      if (itauinv==0) call fatal_error('register_special', 'failed to allocate aux slot for tauinv')
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
      real, contiguous, dimension(:,:,:,:) :: f
!
      if (x_1==impossible) x_1=xyz0(1)
      if (y_1==impossible) y_1=xyz0(2)
      if (z_1==impossible) z_1=xyz0(3)
      if (x_2==impossible) x_2=xyz1(1)
      if (y_2==impossible) y_2=xyz1(2)
      if (z_2==impossible) z_2=xyz1(3)
!
      select case (far_field_type)
        case ('initial')
!         The profile from the initial condition will be used as the
!         reference damping profile.
!         KG: The problem with the way it is currently implemented is
!         KG: that restarting the run leads to the used profile changing.
          call update_profiles(f)
        case ('time-dep')
!         Update at every timestep
          lprof_from_initial=.false.
          call update_profiles(f)
        case default
          call fatal_error('initialize_special', 'Unknown far_field_type = '//trim(far_field_type))
      endselect
!
!     Need to call init_special again since it calculates an auxiliary variable.
!
      if (lrun) call init_special(f)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f)
!
!     While the earlier CPU version of this module simply defined tauinv_prof
!     as a global array, that seems to cause problems for transpilation (as of
!     Pencil git commit 1109640). As a workaround, we store tauinv_prof in the
!     f-array.
!
      real, contiguous, dimension(:,:,:,:) :: f
!
      f(:,:,:,itauinv) = spread(spread(step(x,x_1,-w)+step(x,x_2,w),2,my),3,mz) &
                       + spread(spread(step(y,y_1,-w)+step(y,y_2,w),1,mx),3,mz) &
                       + spread(spread(step(z,z_1,-w)+step(z,z_2,w),1,mx),2,my)
!
      where (f(:,:,:,itauinv)>1) f(:,:,:,itauinv) = 1 !avoid damping being too strong in the corners
!
      f(:,:,:,itauinv) = f(:,:,:,itauinv)/tau
!
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special
!
      lpenc_requested(i_uu)=.true.
      if (ldamp_rho) lpenc_requested(i_rho)=.true.
      if (ldamp_ss) lpenc_requested(i_ss)=.true.
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine read_special_run_pars(iomsg)
!
      use File_io, only: parallel_unit
!
      character(len=iomsglen), intent(out) :: iomsg
      integer :: iostat
!
      read(parallel_unit, NML=special_run_pars, IOSTAT=iostat, IOMSG=iomsg)
      if (iostat==0) iomsg=""
!
    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=special_run_pars)
!
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine special_calc_hydro(f,df,p)
!
      real, contiguous, dimension(:,:,:,:) :: f
      real, contiguous, dimension(:,:,:,:) :: df
      type(pencil_case), intent(in) :: p
!
      real, dimension(nx) :: tauinv
!
      tauinv = f(l1:l2,m,n,itauinv)
!
      df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) - tauinv*p%uu(:,1)
      df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) - tauinv*p%uu(:,2)
      df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz) - tauinv*p%uu(:,3)
!
    endsubroutine special_calc_hydro
!***********************************************************************
    subroutine special_calc_density(f,df,p)
!
      real, contiguous, dimension(:,:,:,:) :: f
      real, contiguous, dimension(:,:,:,:) :: df
      type(pencil_case), intent(in) :: p
!
      real, dimension(nx) :: tauinv
!
      if (ldamp_rho) then
        tauinv = f(l1:l2,m,n,itauinv)
!
        df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) - tauinv*(p%rho/rho_prof(n) - 1)
      endif
!
    endsubroutine special_calc_density
!***********************************************************************
    subroutine special_calc_energy(f,df,p)
!
      real, contiguous, dimension(:,:,:,:) :: f
      real, contiguous, dimension(:,:,:,:) :: df
      type(pencil_case), intent(in) :: p
!
      real, dimension(nx) :: tauinv
!
      if (ldamp_ss) then
        tauinv = f(l1:l2,m,n,itauinv)
!
        df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) - tauinv*(p%ss - ss_prof(n))
      endif
!
    endsubroutine special_calc_energy
!***********************************************************************
    subroutine get_from_xyroot(dest, f, ivar)
!
      use Mpicomm, only: mpibcast, MPI_COMM_XYPLANE
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mz), intent(out) :: dest
      integer, intent(in) :: ivar
!
      if (ipx==0.and.ipy==0) dest = f(l1,m1,:,ivar)
!
      call mpibcast(dest,size(dest),comm=MPI_COMM_XYPLANE)
!
    endsubroutine get_from_xyroot
!***********************************************************************
    subroutine special_after_boundary(f)
!
!     Used in case rho_prof and ss_prof need to be updated every timestep
!
      real, contiguous, dimension(:,:,:,:) :: f
!
      if (.not.lprof_from_initial) then
        call update_profiles(f)
      endif
!
    endsubroutine special_after_boundary
!***********************************************************************
    subroutine update_profiles(f)
!
!     Get rho_prof and ss_prof from the f array
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
!
      if (ldamp_rho) then
        if (ilnrho==0) call fatal_error('initialize_special', 'could not find density variable to be damped')
        if (ldensity_nolog) call not_implemented('initialize_special', 'damping rho with ldensity_nolog=T')
        if (lreference_state) call not_implemented('initialize_special', 'damping rho with lreference_state=T')
        call get_from_xyroot(rho_prof, f, ilnrho)
        rho_prof = exp(rho_prof)
      endif
!
      if (ldamp_ss) then
        if (iss==0) call fatal_error('initialize_special', 'could not find entropy variable to be damped')
        if (pretend_lnTT) call not_implemented('initialize_special', 'damping entropy with pretend_lnTT=T')
        if (lreference_state) call not_implemented('initialize_special', 'damping entropy with lreference_state=T')
        call get_from_xyroot(ss_prof, f, iss)
      endif
!
    endsubroutine update_profiles
!***********************************************************************
    subroutine pushpars2c(p_par)
!
      use Syscalls, only: copy_addr
      use General , only: string_to_enum
!
      integer, parameter :: n_pars=100
      integer(KIND=ikind8), dimension(n_pars) :: p_par
!
      call copy_addr(itauinv,p_par(1)) ! int
      call copy_addr(rho_prof,p_par(2)) ! (mz)
      call copy_addr(ss_prof,p_par(3)) ! (mz)
      call copy_addr(ldamp_rho,p_par(4)) ! bool
      call copy_addr(ldamp_ss,p_par(5)) ! bool
      call copy_addr(lprof_from_initial,p_par(6)) ! bool
!
    endsubroutine pushpars2c
!***********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../special_dummies.inc'
!********************************************************************
endmodule Special
