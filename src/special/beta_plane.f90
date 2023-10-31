! $Id$
!
!  Add Coriolis force in the beta plane approximation to the hydro equation.
!  The x-axis is in the direction of increasing colatitude, while the
!  z-axis is the radial direction.
!
!  31-Oct-2023: Kishore G. Added.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of special_dummies.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
!
module Special
!
  use Cparam
  use Cdata, Omega_fplane => Omega
  use General, only: keep_compiler_quiet
  use Messages, only: svn_id, fatal_error
!
  implicit none
!
  include '../special.h'
!
  real :: Omega = 0 !rotational rate
  real :: R = 1 !Radius of the sphere
  real :: theta_0 = pi/2 !Colatitude about which the Coriolis force is linearized
  real :: cth = impossible, sth = impossible, Rinv=impossible
!
! run parameters
  namelist /special_run_pars/ &
    Omega, R, theta_0
!
  contains
!***********************************************************************
    subroutine initialize_special(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
      cth = cos(theta_0)
      sth = sin(theta_0)
!
      if (.not.lcartesian_coords) call fatal_error("initialize_special", &
        "Cartesian coordinates required for beta-plane approximation")
      if (lgrav.and..not.lgravz) call fatal_error("initialize_special", &
        "Gravity needs to be in the z direction")
      if (Omega_fplane /= 0) call fatal_error("initialize_special", &
        "Do not set Omega in hydro_run_pars")
!
      if (R>0) then
        Rinv = 1/R
      else
        !We allow the user to set R<0 to obtain the f-plane approximation.
        Rinv = 0
      endif
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine pencil_criteria_special
!
      lpenc_requested(i_uu)=.true.
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine read_special_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=special_run_pars, IOSTAT=iostat)
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
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      call keep_compiler_quiet(f)
!
      df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) &
                          - 2*Omega*(cth - sth*x(l1:l2)*Rinv )*p%uu(:,2)
!
      df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) &
                          + 2*Omega*(sth + cth*x(l1:l2)*Rinv)*p%uu(:,3) &
                          + 2*Omega*(cth - sth*x(l1:l2)*Rinv)*p%uu(:,1)
!
      df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz) &
                          - 2*Omega*(sth + cth*x(l1:l2)*Rinv)*p%uu(:,2)
!
    endsubroutine special_calc_hydro
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
