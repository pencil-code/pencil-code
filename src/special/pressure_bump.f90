! $Id$
!
!  Thie module addes pressure bumps to your simulations.
!  They are needed e.g. when simulating the effects of a pressure bump
!  on particle drift, particle trapping and planetesimal formation via
!  gravitational fragmentation.
!
!  For this, an additional force on the gas is exerted, similar
!  to what is done in noentropy.f90:
!
!      df(l1:l2,m,n,(iux-1)+j) = df(l1:l2,m,n,(iux-1)+j) &
!              - p%cs2*beta_glnrho_scaled(j)
!
!  but with spacial dependency on.
!
!  So far, only azimuthal and vertical periodic and homogeneuos
!  pressure bumps are included, i.e. only radial dependence of lnrho.
!
! Implemented profiles:
!       - gauss:           gaussian in x
!       - sine:            sinusoidal in x
!
!  The amplitude of the pressure bump scales with beta_glnrho_global, i.e.
!  beta_glnrho_scaled.
!
!  24-feb-17/andreasSchreiber: creation of the module
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
    use Cdata
    use General, only: keep_compiler_quiet
    use Messages
!
    implicit none
!
    include '../special.h'
!
! For user input:
!
    real :: pb_amplitude=1.0
    character (len=labellen) :: pb_type='none'
!
! For calculations
!
    real, dimension (nx) :: pb_profile
!
    namelist /special_init_pars/ &
        pb_type, pb_amplitude
!
    namelist /special_run_pars/ &
        pb_type, pb_amplitude
!
    contains
!***********************************************************************
subroutine pb_special_setup
!
! Setup pb_profile
!
!  Add pressure bump profile, which is later superimposed on
!  ux*pressure_gradient*constants, see blow.
!
    use EquationOfState, only: cs0

    ! if (lroot) print*, 'x(l1:l2)= ', x(l1:l2)
    ! if (lroot) print*, 'l1= ', l1
    ! if (lroot) print*, 'l2= ', l2
!
    if (pb_type/='none') then
        if (lroot) print*, '!!!! Pressure Bump selected: ', pb_type
        select case (pb_type)

        case ('gauss-x')
          pb_profile = -1*pb_amplitude*exp(-x(l1:l2)**2/(Lxyz(1)/2)**2)

        case ('sinwave-x')
          pb_profile = -1*pb_amplitude*sin(2*pi/Lxyz(1)*x(l1:l2))

        case default
          if (lroot) print*, '! ERROR: Couldnt identify selected pressure bump profile pb_type = ', pb_type

!
! Convert pb_profile into "pb_profile_scaled", replace since orig. not needed
!
        pb_profile = pb_profile*Omega/cs0

        end select
        ! if (lroot) print*, '!pb_special_setup! x coordinates: x(l1:l2) = ', x(l1:l2)
        ! if (lroot) print*, '!pb_special_setup! size(x(l1:l2)) = ', size(x(l1:l2))
        if (lroot) print*, '!!!!!!! pb_profile = ', pb_profile
        ! if (lroot) print*, '!pb_special_setup! size(pb_profile) = ', size(pb_profile)
    endif
!
endsubroutine pb_special_setup
!***********************************************************************
subroutine init_special(f)
!
!  initialise special condition; called from start.f90
!  06-oct-2003/tony: coded
!
    use EquationOfState, only: cs20
    use Density, only:  beta_glnrho_global, beta_glnrho_scaled

    real, dimension (mx,my,mz,mfarray) :: f
!
    intent(inout) :: f
    integer :: ix,iy,iz
!
    if (lroot) print*, '**************** init_special ****************'
    ! call pb_special_setup  is UNNEEDED
!
    if (lroot) print*, 'init_special: set pb_profile sub-Keplerian gas velocity'

    if (lroot) print*, '!!!!!!! m1 = ', m1
    if (lroot) print*, '!!!!!!! m2 = ', m2
    if (lroot) print*, '!!!!!!! n1 = ', n1
    if (lroot) print*, '!!!!!!! n2 = ', n2

    ! if (lroot) print*, '!!!!!!! shape(f(l1:l2,m,n,iux)) = ', shape(f(l1:l2,m,n,iux))
    ! if (lroot) print*, '!!!!!!! shape(f(l1:l2,m,n,iuy)) = ', shape(f(l1:l2,m,n,iuy))
    ! if (lroot) print*, '!!!!!!! shape(f(:,:,:,iux)) = ', shape(f(:,:,:,iux))
    ! if (lroot) print*, '!!!!!!! shape(f(:,:,:,iuy)) = ', shape(f(:,:,:,iuy))

    ! f(l1:l2,m,n,iux) = f(l1:l2,m,n,iux) - 1/(2*Omega)*cs20*beta_glnrho_scaled(2)*pb_profile
    ! f(l1:l2,m,n,iuy) = f(l1:l2,m,n,iuy) + 1/(2*Omega)*cs20*beta_glnrho_scaled(1)*pb_profile

    ! print*, 'mx', mx  = 22
    do ix=l1,l2;do iy=m1,m2;do iz=n1,n2

            ! if (lroot) print*, ''
            if (lroot) print*, 'ix, iy, iz: ', ix, ' ', iy, ' ', iz
            if (lroot) print*, 'ix-l1 = ',ix-l1+1
            ! if (lroot) print*, 'BEFORE:', f(ix,iy,iz,iuy)
            ! if (lroot) print*, 'CHANGE:', 1/(2*Omega)*cs20*beta_glnrho_scaled(1)*pb_profile(ix)
            ! if (lroot) print*, 'Omega:', Omega
            ! if (lroot) print*, 'cs20:', cs20
            ! if (lroot) print*, 'beta_glnrho_scaled:', beta_glnrho_scaled
            ! if (lroot) print*, 'beta_glnrho_scaled(0):', beta_glnrho_scaled(0)
            ! if (lroot) print*, 'beta_glnrho_scaled(1):', beta_glnrho_scaled(1)
            ! if (lroot) print*, 'beta_glnrho_scaled(2):', beta_glnrho_scaled(2)
            ! if (lroot) print*, 'pb_profile:', pb_profile
            if (lroot) print*, 'pb_profile(ix-l1+1):', pb_profile(ix-l1+1)

            f(ix,iy,iz,iuy) = f(ix,iy,iz,iuy) + 1/(2*Omega)*cs20*beta_glnrho_scaled(1)*pb_profile(ix-l1+1)

            ! if (lroot) print*, 'AFTER:', f(ix,iy,iz,iuy)
    enddo;enddo;enddo

    ! if (lroot) print*, 'f(l1:l2,m,n,iux)='
    ! if (lroot) print*, f(l1:l2,m,n,iux)
    ! if (lroot) print*, ''
    ! if (lroot) print*, 'f(l1:l2,m,n,iuy)='
    ! if (lroot) print*, f(l1:l2,m,n,iuy)
    ! if (lroot) print*, ''
    ! if (lroot) print*, 'f(:,:,:,iux)='
    ! if (lroot) print*, f(:,:,:,iux)
    ! if (lroot) print*, ''
    ! if (lroot) print*, 'f(:,:,:,iuy)='
    ! if (lroot) print*, f(:,:,:,iuy)
    ! if (lroot) print*, ''
!
    call keep_compiler_quiet(f)
!
endsubroutine init_special
!***********************************************************************
subroutine initialize_special(f)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  06-oct-03/tony: coded
!
    use EquationOfState, only: cs20
!
    integer :: j
!
    real, dimension (mx,my,mz,mfarray) :: f
    real, dimension (mx,my,mz,mvar) :: df
    type (pencil_case) :: p
!
    if (lroot) print*, '**************** initialize_special ****************'
!
    call pb_special_setup
!
endsubroutine initialize_special
!***********************************************************************
subroutine dspecial_dt(f,df,p)
!
!  calculate right hand side of ONE OR MORE extra coupled PDEs
!  along the 'current' Pencil, i.e. f(l1:l2,m,n) where
!  m,n are global variables looped over in equ.f90
!
!  Due to the multi-step Runge Kutta timestepping used one MUST always
!  add to the present contents of the df array.  NEVER reset it to zero.
!
!  Several precalculated Pencils of information are passed for
!  efficiency.
!
!  06-oct-03/tony: coded
!
    use Density, only: beta_glnrho_global, beta_glnrho_scaled
!
    integer :: j
!
    real, dimension (mx,my,mz,mfarray) :: f
    real, dimension (mx,my,mz,mvar) :: df
    type (pencil_case) :: p
!
    intent(in) :: f,p
    intent(inout) :: df
!
    ! if (lroot) print*, ''
    ! if (lroot) print*, '**************** dspecial_dt ****************'
    ! if (lroot) print*, '!!!! t=',t
!
!   Add pressure force from pressure gradient profile.
!
    if (any(beta_glnrho_global /= 0.) .and. pb_type/='none') then
        ! if (headtt) print*, 'dspecial_dt: adding global pressure gradient profile force'
        do j=1,3
            if (beta_glnrho_global(j) /= 0) then
                ! if (lroot) print*, '~~> m:', m
                ! if (lroot) print*, '~~> n:', n
                ! if (lroot) print*, '!special_calc_hydro! for j=',j
                ! if (lroot) print*, '!special_calc_hydro! Adding pb_profile=',pb_profile
                ! if (lroot) print*, '!special_calc_hydro! with p%cs:', p%cs2
                ! if (lroot) print*, '!special_calc_hydro! with beta_glnrho_scaled(j):', beta_glnrho_scaled(j)
                ! if (lroot) print*, '!special_calc_hydro! df before:',  df(l1:l2,m,n,(iux-1)+j)
                ! if (lroot) print*, '!special_calc_hydro! df change:', p%cs2*beta_glnrho_scaled(j)*pb_profile
                df(l1:l2,m,n,(iux-1)+j) = df(l1:l2,m,n,(iux-1)+j) &
                  - p%cs2*beta_glnrho_scaled(j)*pb_profile
                ! if (lroot) print*, '!special_calc_hydro! df after:', df(l1:l2,m,n,(iux-1)+j)
            endif
        enddo
    endif
!
!  Identify module and boundary conditions.
!
    ! if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dspecial_dt'
!
    call keep_compiler_quiet(f,df)
    call keep_compiler_quiet(p)
!
endsubroutine dspecial_dt
!***********************************************************************
subroutine special_calc_hydro(f,df,p)
!
!  Calculate an additional 'special' term on the right hand side of the
!  momentum equation.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array.
!
!  06-oct-03/tony: coded
!
    integer :: j
!
    real, dimension (mx,my,mz,mfarray), intent(in) :: f
    real, dimension (mx,my,mz,mvar), intent(inout) :: df
    type (pencil_case), intent(in) :: p
!
    ! if (lroot) print*, ''
    ! if (lroot) print*, '**************** special_calc_hydro ****************'

!
    call keep_compiler_quiet(f,df)
    call keep_compiler_quiet(p)
!
endsubroutine special_calc_hydro
!***********************************************************************
subroutine read_special_init_pars(iostat)
!
  use File_io, only: parallel_unit
!
  integer, intent(out) :: iostat
!
  read(parallel_unit, NML=special_init_pars, IOSTAT=iostat)
!
endsubroutine read_special_init_pars
!***********************************************************************
subroutine write_special_init_pars(unit)
!
  integer, intent(in) :: unit
!
  write(unit, NML=special_init_pars)
!
endsubroutine write_special_init_pars
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
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../special_dummies.inc'
!********************************************************************
endmodule Special
