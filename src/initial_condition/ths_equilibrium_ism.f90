! $Id$
!
!  Initial conditions for 3D equilibrium supernova driven turbulence 
!  simulations. First run 1D simulation using identical z grid with 
!  initial_condition/ths1D_equilibrium_ism.f90.
!  This module will load the converged 1D profiles derived and saved to
!  init_ism.dat (now init_ism.in).
!
!   Description                               | Relevant function call
!  ------------------------------------------------------------------------
!   Initial condition registration            | register_initial_condition
!     (pre parameter read)                    |
!   Initial condition initialization          | initialize_initial_condition
!     (post parameter read)                   |
!                                             |
!   Initial condition for momentum            | initial_condition_uu
!   Initial condition for density             | initial_condition_lnrho
!   Initial condition for entropy             | initial_condition_ss
!   Initial condition for magnetic potential  | initial_condition_aa
!                                             |
!   Initial condition for all in one call     | initial_condition_all
!     (called last)                           |
!
!   And a similar subroutine for each module with an "init_XXX" call.
!   The subroutines are organized IN THE SAME ORDER THAT THEY ARE CALLED.
!   First uu, then lnrho, then ss, then aa, and so on.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: linitial_condition = .true.
!
!***************************************************************
!
! HOW TO USE THIS FILE
! --------------------
!
! Change the line above to
!    linitial_condition = .true.
! to enable use of custom initial conditions.
!
! The rest of this file may be used as a template for your own initial
! conditions. Simply fill out the prototypes for the features you want
! to use.
!
! Save the file with a meaningful name, e.g. mhs_equilibrium.f90, and
! place it in the $PENCIL_HOME/src/initial_condition directory. This
! path has been created to allow users to optionally check their
! contributions in to the Pencil Code SVN repository. This may be
! useful if you are working on/using an initial condition with
! somebody else or may require some assistance from one from the main
! Pencil Code team. HOWEVER, less general initial conditions should
! not go here (see below).
!
! You can also place initial condition files directly in the run
! directory. Simply create the folder 'initial_condition' at the same
! level as the *.in files and place an initial condition file there.
! With pc_setupsrc this file is linked automatically into the local
! src directory. This is the preferred method for initial conditions
! that are not very general.
!
! To use your additional initial condition code, edit the
! Makefile.local in the src directory under the run directory in which
! you wish to use your initial condition. Add a line that says e.g.
!
!    INITIAL_CONDITION =   initial_condition/mhs_equilibrium
!
! Here mhs_equilibrium is replaced by the filename of your new file,
! not including the .f90 extension.
!
! This module is based on Tony's special module.
!
module InitialCondition
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include '../initial_condition.h'
!
!  Observed number density per cubic cm parameters from Dickey & Lockman
!  Includes neutral hydrogen and warm ionized hydrogen plus helium 
!  proportionately. Multiply by m_u_cgs for gas density
!
  real, parameter, dimension(5) :: nfraction_cgs = & ! particles per cm cubed normalized to 1 at midplane 
                                      (/0.6541, 0.1775, 0.1028, 0.0245, 0.0411/)
  real, parameter, dimension(5) :: hscale_cgs = & ! scale height in cm
                       (/3.9188e20, 9.8125e20, 1.2435e21, 2.1600e20, 2.7771e21/)
  real, dimension(5) :: rho_fraction, hscale 
!
!  Heating function, cooling function and mass movement
!  method selection.
!
  character (len=labellen) :: cooling_select  = 'WSW'
  character (len=labellen) :: heating_select  = 'wolfire'
!
!
  real, parameter :: T_init_cgs=7.088e2
  real :: T_init=impossible
  real :: rhox=1 ! column density comparative to Milky Way default
!
!
!  TT & z-dependent uv-heating profile
!
  real, dimension(mz) :: rho, lnTT
!
!  Magnetic profile - the magnetic field is propto sqrt density
!  Diffusivities propto sqrt(T)
!
  real, parameter :: amplaa_cgs=1e-21  ! 1 nano Gauss (G g^{1/2} cm^-{3/2})
  real :: amplaa=impossible
  real :: ybias_aa = 1.0 ! adjust angle of field default 45deg
!
  character (len=labellen), dimension(ninit) :: initaa='nothing'
!
!  start parameters
!
  namelist /initial_condition_pars/ &
      T_init, amplaa, cooling_select, rhox, &
      heating_select, initaa, &
      ybias_aa
!
  contains
!***********************************************************************
    subroutine register_initial_condition()
!
!  Register variables associated with this module; likely none.
!
!  07-may-09/wlad: coded
!
      if (lroot) call svn_id( &
         "$Id$")
!
    endsubroutine register_initial_condition
!***********************************************************************
    subroutine initialize_initial_condition(f)
!
!  Initialize any module variables which are parameter dependent.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_initial_condition
!***********************************************************************
    subroutine initial_condition_all(f,profiles)
!
!  Initializes all the f arrays in one call.  This subroutine is called last.
!
!  21-dec-10/ccyang: coded
!  15-feb-15/MR: optional parameter 'profiles' added
!
      use Messages, only: fatal_error
      use EquationOfState, only: getmu, eoscalc
!
      real, dimension (mx,my,mz,mfarray), optional, intent(inout):: f
      real, dimension (:,:),              optional, intent(out)  :: profiles
! 
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
      if (present(profiles)) then
        call fatal_error('initial_condition_all', &
          'If profiles are asked for, a real initial condition must be specified')
        call keep_compiler_quiet(profiles)
      endif
!
    endsubroutine initial_condition_all
!*****************************************************************************
    subroutine initial_condition_uu(f)
!
!  Initialize the velocity field.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      f(:,:,:,iuu:iuu+2) = 0.
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_uu
!***********************************************************************
    subroutine initial_condition_lnrho(f)
!
!  Initialize logarithmic density. init_lnrho will take care of
!  converting it to linear density if you use ldensity_nolog.
!
!  07-jul-15/fred: coded
!
      use General, only: notanumber
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      integer, parameter :: ntotal=nz*nprocz
      real, dimension(nz*nprocz) :: tmp1,tmp2
      real :: var1, var2
      integer :: stat, q
      logical :: exist
!
!  Set up physical units.
!
      if (unit_system=='cgs') then
        if (T_init == impossible) T_init = T_init_cgs/unit_temperature
      else if (unit_system=='SI') then
        call fatal_error('initial_condition_lnrho','SI unit conversions not inplemented')
      endif
!
!  Read rho and TT and write into an array.
!  If file is not found in run directory, search under trim(directory).
!
      inquire(file='init_ism.in',exist=exist)
      if (exist) then
        open(31,file='init_ism.in')
      else
        inquire(file=trim(directory)//'/init_ism.ascii',exist=exist)
        if (exist) then
          open(31,file=trim(directory)//'/init_ism.ascii')
        else
          call fatal_error('initial_lnrho','error - no init_ism input file')
        endif
      endif
!
!  Read profiles.
!
      do q=1,ntotal
        read(31,*,iostat=stat) var1,var2
        if (stat<0) exit
        if (ip<5) print*,'rho, TT: ',var1,var2
        tmp1(q)=var1/unit_density
        tmp2(q)=var2/unit_temperature
      enddo
      close(31)
!
!  Assuming no ghost zones in init_ism.in.
!
      do n=n1,n2
        rho(n)=tmp1(n-nghost+ipz*nz)
        if (ldensity_nolog) then
          f(l1:l2,m1:m2,n,irho)=log(rho(n))
        else
          f(l1:l2,m1:m2,n,ilnrho)=log(rho(n))
        endif
        lnTT(n)=log(tmp2(n-nghost+ipz*nz))
      enddo
!
      if (notanumber(f(l1:l2,m1:m2,n1:n2,ilnrho))) &
        call error('initial_condition_lnrho', 'Imaginary density values')
!
    endsubroutine initial_condition_lnrho
!***********************************************************************
    subroutine initial_condition_ss(f)
!
!  Initialize entropy.
!
!  07-may-09/wlad: coded
!
      use EquationOfState, only: eoscalc, ilnrho_lnTT 
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real :: lnrho, ss
!
      do n=n1,n2
        if (ldensity_nolog) then
          lnrho=log(f(l1,m1,n,irho))
        else
          lnrho=f(l1,m1,n,ilnrho)
        endif
!
        call eoscalc(ilnrho_lnTT,lnrho,lnTT(n),ss=ss)
!
        f(l1:l2,m1:m2,n,iss)=ss
!
      enddo
!
    endsubroutine initial_condition_ss
!***********************************************************************
    subroutine initial_condition_aa(f)
!
!  Initialize the magnetic vector potential.
!
!  This routine sets up an initial magnetic field y-parallel(azimuthal) 
!  with a magnitude directly proportional to the density. 
!
      use Mpicomm, only: mpireduce_sum, mpibcast_real
      use General, only: random_number_wrapper
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      integer :: i ,j, l, m, n
      real, dimension(nx) :: nxran
!
      if (.not. lmagnetic) then
!
        call keep_compiler_quiet(f)
      else     
        do j=1,ninit
!  
          select case (initaa(j))
          case ('nothing'); if (lroot .and. j==1) print*,'init_aa: nothing'
          case ('zero', '0'); f(:,:,:,iax:iaz) = 0.0
          case ('gaussian-noise')
            f(:,:,:,iax:iay) = 0.0
            do m=m1,m2
            do n=n1,n2
              call random_number_wrapper(nxran)
              f(l1:l2,m,n,iaz) = amplaa * nxran * sqrt(rho(n))
            enddo
            enddo
          case ('uniform-x')
            f(:,:,:,iax:iay) = 0.0
            do l=l1,l2
            do n=n1,n2
              f(l,m1:m2,n,iaz) = amplaa * y(m1:m2) * sqrt(rho(n))
            enddo
            enddo
          case ('uniform-y')
            f(:,:,:,iax:iay) = 0.0
            do m=m1,m2
            do n=n1,n2
              f(l1:l2,m,n,iaz) = -amplaa * x(l1:l2) * sqrt(rho(n))
            enddo
            enddo
          case ('uniform-x+y')
            f(:,:,:,iax:iay) = 0.0
            do l=l1,l2
            do m=m1,m2
            do n=n1,n2
              f(l,m,n,iaz) = amplaa * (-x(l) + ybias_aa*y(m)) * sqrt(rho(n)) / &
                             sqrt(1 + ybias_aa**2)
            enddo
            enddo
            enddo
          case default
!  
!  Catch unknown values.
!  
            call fatal_error('inititial_condition_aa', &
                'init_aa value "' // trim(initaa(j)) // '" not recognised')
!  
          endselect
!  
!  End loop over initial conditions.
!  
        enddo
!
      endif
!
    endsubroutine initial_condition_aa
!***********************************************************************
    subroutine initial_condition_aatest(f)
!
!  Initialize testfield.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_aatest
!***********************************************************************
    subroutine initial_condition_uutest(f)
!
!  Initialize testflow.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_uutest
!***********************************************************************
    subroutine initial_condition_lncc(f)
!
!  Initialize passive scalar.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_lncc
!***********************************************************************
    subroutine initial_condition_chiral(f)
!
!  Initialize chiral.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_chiral
!***********************************************************************
    subroutine initial_condition_chemistry(f)
!
!  Initialize chemistry.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_chemistry
!***********************************************************************
    subroutine initial_condition_uud(f)
!
!  Initialize dust fluid velocity.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_uud
!***********************************************************************
    subroutine initial_condition_nd(f)
!
!  Initialize dust fluid density.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_nd
!***********************************************************************
    subroutine initial_condition_uun(f)
!
!  Initialize neutral fluid velocity.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_uun
!***********************************************************************
    subroutine initial_condition_lnrhon(f)
!
!  Initialize neutral fluid density.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_lnrhon
!***********************************************************************
    subroutine initial_condition_ecr(f)
!
!  Initialize cosmic rays.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_ecr
!***********************************************************************
    subroutine initial_condition_fcr(f)
!
!  Initialize cosmic ray flux.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_fcr
!***********************************************************************
    subroutine initial_condition_solid_cells(f)
!
!  Initialize solid cells.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_solid_cells
!***********************************************************************
    subroutine initial_condition_cctest(f)
!
!  Initialize testscalar.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_cctest
!***********************************************************************
    subroutine initial_condition_xxp(f,fp)
!
!  Initialize particles' positions.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (:,:), intent(inout) :: fp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
!
    endsubroutine initial_condition_xxp
!***********************************************************************
    subroutine initial_condition_vvp(f,fp)
!
!  Initialize particles' velocities.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (:,:), intent(inout) :: fp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
!
    endsubroutine initial_condition_vvp
!***********************************************************************
    subroutine read_initial_condition_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=initial_condition_pars, IOSTAT=iostat)
!
    endsubroutine read_initial_condition_pars
!***********************************************************************
    subroutine write_initial_condition_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=initial_condition_pars)
!
    endsubroutine write_initial_condition_pars
!***********************************************************************
    subroutine initial_condition_clean_up
!
!  04-may-11/dhruba: coded
! dummy
!
    endsubroutine initial_condition_clean_up
!***********************************************************************
!***********************************************************************
!
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from noinitial_condition.f90 for any    **
!**  InitialCondition routines not implemented in this file        **
!**                                                                **
    include '../initial_condition_dummies.inc'
!********************************************************************
endmodule InitialCondition
