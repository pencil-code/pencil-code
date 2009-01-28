! $Id$
!
!  This modules takes care of viscosity of inertial particles.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_viscosity = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 6
!
!***************************************************************

module Particles_viscosity

  use Cparam
  use Cdata
  use Messages
  use Particles_cdata
  use Sub, only: keep_compiler_quiet

  implicit none

  include 'particles_viscosity.h'

  integer, parameter :: nviscp_max = 4
  character (len=labellen), dimension(nviscp_max) :: iviscp=''
  real :: nup=0.0
  logical :: lviscp_simplified=.false.
  logical :: lviscp_rhop_nup_const=.false.
  logical :: lviscp_nup_const=.false.

  namelist /particles_visc_run_pars/ &
      nup, iviscp

  contains

!***********************************************************************
    subroutine register_particles_viscosity()
!
!  07-oct-08/anders: coded
!
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_particles_viscosity called twice')
      first = .false.
!
!
!  Set indices for auxiliary variables
!
      iuup   = mvar + naux + 1 + (maux_com - naux_com); naux = naux + 3
      ipvisc = mvar + naux + 1 + (maux_com - naux_com); naux = naux + 3
!
!  Identify version number.
!
      if (lroot) call cvs_id( &
           "$Id$")
!
!  Check that we aren't registering too many auxilary variables
!
      if (naux > maux) then
        if (lroot) write(0,*) 'naux = ', naux, ', maux = ', maux
        call fatal_error('register_particles_viscosity: naux > maux','')
      endif
!
    endsubroutine register_particles_viscosity
!***********************************************************************
    subroutine initialize_particles_viscosity(lstarting)
!
!  07-oct-08/anders: coded
!
      use FArrayManager
      use Mpicomm, only: stop_it
      use SharedVariables
!
      logical, intent(in) :: lstarting
!
      integer :: i, ierr
!
!  Some viscosity types need the rate-of-strain tensor and grad(lnrho)
!
      lviscp_simplified=.false.
      lviscp_rhop_nup_const=.false.
      lviscp_nup_const=.false.
!
      do i=1,nviscp_max
        select case (iviscp(i))
        case ('simplified', '0')
          if (lroot) print*,'particle viscosity: nup*del2v'
          lviscp_simplified=.true.
        case('rhop_nup-const', '1')
          if (lroot) print*,'particle viscosity: mup/rhop*(del2v+graddivv/3)'
          lviscp_rhop_nup_const=.true.
        case('nup-const')
          if (lroot) &
              print*,'particle viscosity: nup*(del2v+graddivv/3+2S.glnrhop)'
          lviscp_nup_const=.true.
        case ('none','')
          ! do nothing
        case default
          if (lroot) print*, 'No such value for iviscp(',i,'): ', trim(iviscp(i))
          call fatal_error('initialize_particles_viscosity','')
        endselect
      enddo
!
!  If we're timestepping, die or warn if the viscosity coefficient that
!  corresponds to the chosen viscosity type is not set.
!
      if (lrun) then
        if ( (lviscp_simplified.or.lviscp_rhop_nup_const.or.lviscp_nup_const) &
            .and.nup==0.0) &
            call warning('initialize_particles_viscosity', &
            'Viscosity coefficient nup is zero!')
      endif
!
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_particles_viscosity
!***********************************************************************
    subroutine calc_particles_viscosity(f,fp,ineargrid)
!
!
!
      use Particles_sub, only: map_vvp_grid
      use Sub, only: del2v
!
      real, dimension (mx,my,mz,mfarray) :: f 
      real, dimension (mpar_loc,mpvar) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension (nx,3) :: del2uup
!
      call map_vvp_grid(f,fp,ineargrid)
!
      if (lviscp_simplified) then
        do n=n1,n2; do m=m1,m2
          call del2v(f,iuup,del2uup)
          f(l1:l2,m,n,ipvisc:ipvisc+2)=del2uup
        enddo; enddo
      endif
!
    endsubroutine calc_particles_viscosity
!***********************************************************************
    subroutine calc_particles_viscous_force(df,p)
!
!  calculate viscous force term for right hand side of  equation
!
!  20-nov-02/tony: coded
!   9-jul-04/nils: added Smagorinsky viscosity

      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: nu_smag
      real, dimension (nx,3) :: nuD2uxb
      type (pencil_case) :: p
!
      intent (in) :: p
      intent (inout) :: df 
!
!  Add viscosity to equation of motion
!
      df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) + p%fvisc
!
!  Calculate max total diffusion coefficient for timestep calculation etc.
!
      if (lfirst.and.ldt) then
        diffus_nu =p%diffus_total *dxyz_2
        diffus_nu2=p%diffus_total2*dxyz_4
        diffus_nu3=p%diffus_total3*dxyz_6
      endif
!
    endsubroutine calc_particles_viscous_force
!***********************************************************************
    subroutine read_particles_visc_init_pars(unit,iostat)
!
!
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat).and.NO_WARN) print*,iostat
      if (NO_WARN) print*,unit
!
    endsubroutine read_particles_visc_init_pars
!***********************************************************************
    subroutine write_particles_visc_init_pars(unit)
!
!
!
      integer, intent(in) :: unit
!
      if (NO_WARN) print*,unit
!
    endsubroutine write_particles_visc_init_pars
!***********************************************************************
    subroutine read_particles_visc_run_pars(unit,iostat)
!
!
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_visc_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_visc_run_pars,ERR=99)
      endif
!
99    return
    endsubroutine read_particles_visc_run_pars
!***********************************************************************
    subroutine write_particles_visc_run_pars(unit)
!
!
!
      integer, intent(in) :: unit
!
      write(unit,NML=particles_visc_run_pars)
!
    endsubroutine write_particles_visc_run_pars
!***********************************************************************
    subroutine calc_viscosity(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (NO_WARN) print*,f  !(to keep compiler quiet)
!
    endsubroutine calc_viscosity
!*******************************************************************
    subroutine rprint_particles_viscosity(lreset,lwrite)
!
!  Writes ishock to index.pro file
!
!  07-oct-08/anders: adapted
!
      use Sub
!
      logical :: lreset
      logical, optional :: lwrite
      integer :: iname
!
      if (present(lwrite)) then
        if (lwrite) then
          write(3,*) 'iuup=', iuup
          write(3,*) 'ipvisc=', ipvisc
        endif
      endif
!
      if (NO_WARN) print*,lreset  !(to keep compiler quiet)
!
    endsubroutine rprint_particles_viscosity
!***********************************************************************
endmodule Particles_viscosity
