! $Id: particles_adsorbed.f90 21950 2014-07-08 08:53:00Z michiel.lambrechts $
!
!  This module takes care of everything related to reactive particles.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MPVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
! CPARAM logical, parameter :: lparticles_adsorbed=.true.
!
!! PENCILS PROVIDED adsp
!
!***************************************************************
module Particles_adsorbed
!
  use Cdata
  use Cparam
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
  use Particles_map
  use Particles_mpicomm
  use Particles_sub
  use Particles_radius
 !
  implicit none
!
  include 'particles_adsorbed.h'
!
  character (len=labellen), dimension (ninit) :: init_adsorbed='nothing'
  real :: init_surf_frac
  real :: diffusivity=0.0
  real :: init_thCO=0.0
!
  namelist /particles_ads_init_pars/ &
      init_adsorbed,init_surf_frac,init_thCO, diffusivity
!
  namelist /particles_ads_run_pars/ &
      diffusivity
!
  contains
!***********************************************************************
    subroutine register_particles_ads()
!
!  Set up indices for access to the fp and dfp arrays
!
!  29-aug-14/jonas: coded
!
      use FArrayManager, only: farray_register_auxiliary
!
      if (lroot) call svn_id( &
           "$Id: particles_dust.f90 21950 2014-07-08 08:53:00Z michiel.lambrechts $")
!
!  Indices for particle position.
!
      iCOp=npvar+1
      pvarname(npvar+1)='iCOp'
!
!  Increase npvar accordingly.
!
      npvar=npvar+1
!
!  Check that the fp and dfp arrays are big enough.
!
      if (npvar > mpvar) then
        if (lroot) write(0,*) 'npvar = ', npvar, ', mpvar = ', mpvar
        call fatal_error('register_ads','npvar > mpvar')
      endif
!
    endsubroutine register_particles_ads
!***********************************************************************
    subroutine initialize_particles_ads(f,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  29-aug-14/jonas coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      
!
    end subroutine initialize_particles_ads
!***********************************************************************
    subroutine init_particles_ads(f,fp)
!
!  Initial particle surface fractions
!
!  01-sep-14/jonas: coded
!
      use General, only: random_number_wrapper
      use Mpicomm, only: mpireduce_sum, mpibcast_real
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mpvar) :: fp
      integer :: j
!
!
      intent (out) :: f, fp
!
!
! will have to be adapted to loop over all species
      fp(1:npar_loc,iCOp)=0.
      do j=1,ninit
!
        select case (init_adsorbed(j))
!
        case ('nothing')
          if (lroot .and. j==1) print*, 'init_particles_ads: nothing'
        case ('constant')
!loop over species
          if (lroot) print*, 'init_particles_ads: Initial Surface Fraction of CO'
          fp(1:npar_loc,iTp)=fp(1:npar_loc,iCOp)+init_thCO
        case default
          if (lroot) &
              print*, 'init_particles_ads: No such such value for init_ads: ', &
              trim(init_adsorbed(j))
          call fatal_error('init_ads','')
!
        endselect
!
      enddo
!
    endsubroutine init_particles_ads
!***********************************************************************    
subroutine pencil_criteria_par_ads()
!
!  All pencils that the Particles_adsorbed module depends on are specified here.
!
!  01-sep-14/jonas: coded
!
    endsubroutine pencil_criteria_par_ads
!***********************************************************************
    subroutine dpads_dt(f,df,fp,dfp,ineargrid)
!
!  Evolution of particle surface fractions.
!
!  01-sep-14/jonas: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dpads_dt
!***********************************************************************
    subroutine dpads_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of particle surface fractions
!
!  01-sep-14/jonas: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      type (pencil_case) :: p
      integer, dimension (mpar_loc,3) :: ineargrid
!
      intent (in) :: f, df, fp, ineargrid
      intent (inout) :: dfp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(p)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dpads_dt_pencil
!***********************************************************************
    subroutine read_particles_ads_init_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_ads_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_ads_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_ads_init_pars
!***********************************************************************
    subroutine write_particles_ads_init_pars(unit)
!
      integer, intent (in) :: unit
!
      write(unit,NML=particles_ads_init_pars)
!
    endsubroutine write_particles_ads_init_pars
!***********************************************************************
    subroutine read_particles_ads_run_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_ads_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_ads_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_ads_run_pars
!***********************************************************************
    subroutine write_particles_ads_run_pars(unit)
!
      integer, intent (in) :: unit
!
      write(unit,NML=particles_ads_run_pars)
!
    endsubroutine write_particles_ads_run_pars
!***********************************************************************
    subroutine rprint_particles_ads(lreset,lwrite)
!
!  Read and register print parameters relevant for particles coverage fraction.
!
!  29-aug-14/jonas: coded
!
      logical :: lreset
      logical, optional :: lwrite
!
      logical :: lwr
!
!  Write information to index.pro
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
      if (lwr) write(3,*) 'iox=', iox
!
      call keep_compiler_quiet(lreset)
!
    endsubroutine rprint_particles_ads
!***********************************************************************    
subroutine particles_ads_prepencil_calc(f)
!
!  28-aug-14/jonas+nils: coded
!
      real,dimension(mx,my,mz,mfarray),intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine particles_ads_prepencil_calc
!***********************************************************************
  end module Particles_adsorbed
