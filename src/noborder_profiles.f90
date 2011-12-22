! $Id$
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lborder_profiles = .false.
!
!***************************************************************
module BorderProfiles
!
  use Cdata
  use Cparam
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  private
!
  include 'border_profiles.h'
!
  contains
!***********************************************************************
    subroutine initialize_border_profiles()
!
!  Position-dependent quenching factor that multiplies rhs of pde
!  by a factor that goes gradually to zero near the boundaries.
!  border_frac_[xyz] is a 2-D array, separately for all three directions.
!  border_frac_[xyz]=1 would affect everything between center and border.
!
      use Messages, only: fatal_error
!
      if (border_frac_x(1)/=0.0.or.border_frac_x(2)/=0.0) then
        if (lroot) then
          print*, 'initialize_border_profiles: must use '// &
              'BORDER_PROFILES =   border_profiles'
          print*, '                            for border_frac_x'
        endif
        call fatal_error('initialize_border_profiles','')
      endif
      if (border_frac_y(1)/=0.0.or.border_frac_y(2)/=0.0) then
        if (lroot) then
          print*, 'initialize_border_profiles: must use '// &
              'BORDER_PROFILES =   border_profiles'
          print*, '                            for border_frac_y'
        endif
        call fatal_error('initialize_border_profiles','')
      endif
      if (border_frac_z(1)/=0.0.or.border_frac_z(2)/=0.0) then
        if (lroot) then
          print*, 'initialize_border_profiles: must use '// &
              'BORDER_PROFILES =   border_profiles'
          print*, '                            for border_frac_z'
        endif
        call fatal_error('initialize_border_profiles','')
      endif
!
    endsubroutine initialize_border_profiles
!***********************************************************************
    subroutine request_border_driving(border_var)
!
      character (len=labellen) :: border_var
      call keep_compiler_quiet(border_var)
!
    endsubroutine request_border_driving
!***********************************************************************
    subroutine pencil_criteria_borderprofiles()
!
    endsubroutine pencil_criteria_borderprofiles
!***********************************************************************
    subroutine calc_pencils_borderprofiles(f,p)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_borderprofiles
!***********************************************************************
    subroutine set_border_initcond(f,ivar,tmp)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: tmp
      integer :: ivar
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(ivar)
      call keep_compiler_quiet(tmp)
!
    endsubroutine set_border_initcond
!***********************************************************************
    subroutine border_driving(f,df,p,f_target,j)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: f_target
      integer :: j
!
!  Dummy routine
!
      call keep_compiler_quiet(j)
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(f_target)
!
    endsubroutine border_driving
!***********************************************************************
    subroutine border_quenching(f,df,j,dt_sub)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real :: dt_sub
      integer :: j
!
!  Dummy routine
!
      call keep_compiler_quiet(dt_sub)
      call keep_compiler_quiet(j)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(f)
!
    endsubroutine border_quenching
!***********************************************************************
endmodule BorderProfiles
