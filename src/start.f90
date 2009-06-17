! $Id$
!
!***********************************************************************
!
!  The Pencil Code is a high-order finite-difference code for compressible
!  hydrodynamic flows with magnetic fields and particles. It is highly
!  modular and can easily be adapted to different types of problems.
!
!      MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
!      MMMMMMM7MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
!      MMMMMMMIMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
!      MMMMMMMIIMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
!      MMMMMMM7IMMMMMMMMMMMMMMMMMMM7IMMMMMMMMMMMMMMMMMMMMDMMMMMMM
!      MMMMMMZIIMMMMMMMMMMMMMMMMMMMIIMMMMMMMMMMMMMMMMMMMMIMMMMMMM
!      MMMMMMIIIZMMMMMMMMMMMMMMMMMMIIMMMMMMMMMMMMMMMMMMMMIMMMMMMM
!      MMMMMMIIIIMMMMMMMMMMMMMMMMMNII$MMMMMMMMMMMMMMMMMM$IMMMMMMM
!      MMMMM8IIIIMMMMMMMMMMMMMMMMM$IIIMMMMMMMMMMMMMMMMMMII7MMMMMM
!      MMMMMD7II7MMMMMMMMMMMMMMMMMIIIIMMMMMMMMMMMMMMMMMMIIIMMMMMM
!      MMMMN..:~=ZMMMMMMMMMMMMMMMMIIIIDMMMMMMMMMMMMMMMMDIIIMMMMMM
!      MMMM8.,:~=?MMMMMMMMMMMMMMMMOII7NMMMMMMMMMMMMMMMMZIIIDMMMMM
!      MMMM. ,::=+MMMMMMMMMMMMMMM..,~=?MMMMMMMMMMMMMMMMIIII$MMMMM
!      MMMM..,:~=+MMMMMMMMMMMMMM8 .,~=+DMMMMMMMMMMMMMMM8II78MMMMM
!      MMMM .,:~==?MMMMMMMMMMMMMN .,~=+NMMMMMMMMMMMMMM8 ,~~+MMMMM
!      MMM7  ,:~+=?MMMMMMMMMMMMM  .,~==?MMMMMMMMMMMMMM..,~~+DMMMM
!      MMM.  ,:~==?MMMMMMMMMMMMM  .,~~=?MMMMMMMMMMMMMM. ,~~+?MMMM
!      MMM.  ,:~~=??MMMMMMMMMMMM  ,,~~=?MMMMMMMMMMMMM8 .,~~=?NMMM
!      MMM.  ,:~~=+?MMMMMMMMMMMI  ,,~~=+?MMMMMMMMMMMM. .,~~=?MMMM
!      MM~. .,:~:==?MMMMMMMMMMM   ,,~~==?MMMMMMMMMMMM  .,~~=+?MMM
!      MMN8 .D,D+=M=8MMMMMMMMMM   ,,~~==?MMMMMMMMMMMM. .,~~==?MMM
!      MM==8.   .8===MMMMMMMMMM  .,,~~==?$MMMMMMMMMM=  .,~~==?MMM
!      MM==D.    +===MMMMMMMMMO$ .I?7$=7=NMMMMMMMMMM.  ,,~~==?$MM
!      MM==D.    +===MMMMMMMMM==O?..  ====MMMMMMMMMM.  ,,~~==??MM
!      MM==D.    +===MMMMMMMMM===.    .===MMMMMMMMMM+$.?=,7==?7MM
!      MM==D.    +===MMMMMMMMM===.    .===MMMMMMMMMZ==8    .8==MM
!      MM==D.    +===MMMMMMMMM===.    .===MMMMMMMMMZ==I    .Z==MM
!      MM==D. .  +===MMMMMMMMM===.... .===MMMMMMMMMZ==I    .Z==MM
!      MM==D.    +===MMMMMMMMM===.    .===MMMMMMMMMZ==I    .Z==MM
!      MM==D.    +===MMMMMMMMM===.    .===MMMMMMMMMZ==I    .Z==MM
!      MM==D.    +===MMMMMMMMM===..   .===MMMMMMMMMZ==I    .Z==MM
!      MM==D. .  +===MMMMMMMMM===.... .===MMMMMMMMMZ==I    .Z==MM
!
!  More information can be found in the Pencil Code manual and at the
!  website http://www.nordita.org/software/pencil-code/.
!
!***********************************************************************
program start
!
  use Boundcond,        only: update_ghosts
  use Cdata
  use Chemistry,        only: init_chemistry
  use Chiral,           only: init_chiral
  use Cosmicrayflux,    only: init_fcr
  use Cosmicray,        only: init_ecr
  use Density,          only: init_lnrho
  use Dustdensity,      only: init_nd
  use Dustvelocity,     only: init_uud
  use Entropy,          only: init_ss
  use EquationOfState
  use FArrayManager,    only: farray_clean_up
  use Filter
  use General
  use Gravity,          only: init_gg
  use Grid
  use Hydro,            only: init_uu
  use Hyperresi_strict, only: hyperresistivity_strict
  use Hypervisc_strict, only: hyperviscosity_strict
  use Initcond
  use Interstellar,     only: init_interstellar
  use IO
  use Magnetic,         only: init_aa
  use Lorenz_gauge,     only: init_lorenz_gauge
  use Messages
  use Mpicomm
  use NeutralDensity,   only: init_lnrhon
  use NeutralVelocity,  only: init_uun
  use Param_IO
  use Particles_main
  use Particles_nbody,  only: particles_nbody_write_snapshot, &
                              particles_nbody_write_spdim
  use PScalar,          only: init_lncc
  use Radiation,        only: init_rad, radtransfer
  use Register
  use Selfgravity,      only: calc_selfpotential
  use SharedVariables,  only: sharedvars_clean_up
  use Snapshot
  use Solid_Cells,      only: init_solid_cells
  use Special,          only: init_special
  use Sub
  use Testfield,        only: init_aatest
  use Testflow,         only: init_uutest
!
  implicit none
!
!  Define parameters.
!  The f-array includes auxiliary variables.
!  Although they are not necessary for start.f90, idl may want to read them,
!  so we define therefore the full array and write it out.
!
  integer :: i,ifilter,stat
  logical :: lnoerase=.false.
  real :: x00,y00,z00
!        real, dimension (mx,my,mz,mfarray) :: f
!        real, dimension (mx,my,mz,mvar) :: df
  real, allocatable, dimension (:,:,:,:) :: f,df
!
  lstart = .true.
!
!  Initialize the message subsystem, eg. color setting etc.
!
  call initialize_messages()
!
!  Allocate large arrays. We need to make them allocatable in order to
!  avoid segfaults at 128^3 (7 variables) with Intel compiler on 32-bit
!  Linux boxes. Not clear why they crashed (we _did_ increase stacksize
!  limits), but they did and now don't. Also, the present approach runs
!  up to nx=ny=nz=135, but not for even slightly larger grids.
!
  allocate( f(mx,my,mz,mfarray),STAT=stat)
  if (stat>0) call stop_it("Couldn't allocate memory for f ")

  allocate(df(mx,my,mz,mvar)   ,STAT=stat)
  if (stat>0) call stop_it("Couldn't allocate memory for df")

!
!  Pre-initialize f and df to absurd value (to crash the code should we
!  later use uninitialized slots of thos fields)
!
  f = huge(1.0);
  df = huge(1.0);

!
!  Register variables in the f array.
!
  call register_modules()
  call particles_register_modules()
!
!  The logical headtt is sometimes referred to in start.x, even though it is
!  not yet defined. So we set it simply to lroot here.
!
  headtt=lroot
!
!  Identify version.
!
  if (lroot) call svn_id( &
       "$Id$")
!
!  Set default values: box of size (2pi)^3.
!
  xyz0 = (/       -pi,        -pi,       -pi /) ! first corner
  xyz1 = (/impossible, impossible, impossible/) ! last corner
  Lxyz = (/impossible, impossible, impossible/) ! box lengths
  lperi =(/.true.,.true.,.true. /) ! all directions periodic
  lequidist=(/.true.,.true.,.true. /) ! all directions equidistant grid
  lshift_origin=(/.false.,.false.,.false./) ! don't shift origin
!
!  Read parameters from start.in.
!  Call also rprint_list, because it writes iuu, ilnrho, iss, and iaa to disk.
!
  call read_startpars(FILE=.true.)
  call rprint_list(.false.)
  call particles_rprint_list(.false.)
!
!  Initialize start time.
!
  t=tstart
!
!  Will we write all slots of f?
!
  if (lwrite_aux) then
    mvar_io = mvar+maux
  else
    mvar_io = mvar
  endif
!
!  Print resolution.
!
  if (lroot) print*, 'nxgrid,nygrid,nzgrid=',nxgrid,nygrid,nzgrid
!
!  Postprocess input parameters.
!
  gamma1 = gamma-1.
!
!  Set up directory names and check whether the directories exist.
!
  call directory_names()
!
!  Unfortunately the following test for existence of directory fails under
!  OSF1:
!        inquire(FILE=trim(directory_snap), EXIST=exist)
!        if (.not. exist) &
!             call stop_it('Need directory <' // trim(directory_snap) // '>')
!
!
!  Set box dimensions, make sure Lxyz and xyz1 are in sync.
!  Box defaults to [-pi,pi] for all directions if none of xyz1 or Lxyz are set.
!  If luniform_z_mesh_aspect_ratio=T, the default Lz scales with nzgrid/nxgrid.
!
  do i=1,3
    if (Lxyz(i) == impossible) then
      if (xyz1(i) == impossible) then
        if (i==3.and.luniform_z_mesh_aspect_ratio) then
          Lxyz(i)=2*pi*real(nzgrid)/real(nxgrid)
          xyz0(i)=-pi*real(nzgrid)/real(nxgrid)
        else
          Lxyz(i)=2*pi    ! default value
        endif
      else
        Lxyz(i) = xyz1(i)-xyz0(i)
      endif
    else                  ! Lxyz was set
      if (xyz1(i) /= impossible) then ! both Lxyz and xyz1 are set
        call stop_it('Cannot set Lxyz and xyz1 at the same time')
      endif
    endif
  enddo
  xyz1=xyz0+Lxyz
!
!  Abbreviations
!
  x0=xyz0(1); y0=xyz0(2); z0=xyz0(3)
  Lx=Lxyz(1); Ly=Lxyz(2); Lz=Lxyz(3)
!
!  position of equator (if any)
!
  if (lequatory) yequator=xyz0(2)+0.5*Lxyz(2)
  if (lequatorz) zequator=xyz0(3)+0.5*Lxyz(3)
!
!  Size of box at local processor
!
  Lxyz_loc(1)=Lxyz(1)/nprocx
  Lxyz_loc(2)=Lxyz(2)/nprocy
  Lxyz_loc(3)=Lxyz(3)/nprocz
  xyz0_loc(1)=xyz0(1)+ipx*Lxyz_loc(1)
  xyz0_loc(2)=xyz0(2)+ipy*Lxyz_loc(2)
  xyz0_loc(3)=xyz0(3)+ipz*Lxyz_loc(3)
  xyz1_loc(1)=xyz0_loc(1)+Lxyz_loc(1)
  xyz1_loc(2)=xyz0_loc(2)+Lxyz_loc(2)
  xyz1_loc(3)=xyz0_loc(3)+Lxyz_loc(3)
!
!  Calculate dimensionality of the run.
!
  dimensionality=min(nxgrid-1,1)+min(nygrid-1,1)+min(nzgrid-1,1)
!
!  Check consistency.
!
  if (.not.lperi(1).and.nxgrid<2) &
      call stop_it('for lperi(1)=F: must have nxgrid>1')
  if (.not.lperi(2).and.nygrid<2) &
      call stop_it('for lperi(2)=F: must have nygrid>1')
  if (.not.lperi(3).and.nzgrid<2) &
      call stop_it('for lperi(3)=F: must have nzgrid>1')
!
!  Initialise random number generator in processor-dependent fashion for
!  random initial data.
!  Slightly tricky, since setting seed=(/iproc,0,0,0,0,0,0,0,.../)
!  would produce very short-period random numbers with the Intel compiler;
!  so we need to retain most of the initial entropy of the generator.
!
  call get_nseed(nseed)   ! get state length of random number generator
  call random_seed_wrapper(GET=seed)
  seed(1) = -(10+iproc)    ! different random numbers on different CPUs
  call random_seed_wrapper(PUT=seed)
!
!  Generate grid.
!
  call construct_grid(x,y,z,dx,dy,dz,x00,y00,z00)
!
!  Write grid.dat file.
!
  call wgrid(trim(directory)//'/grid.dat')
  if (lparticles) &
    call wproc_bounds(trim(directory)//'/proc_bounds.dat')
!
!  Write .general file for data explorer.
!
  if (lroot) call write_dx_general( &
                    trim(datadir)//'/var.general', &
                    x0-nghost*dx, y0-nghost*dy, z0-nghost*dz)
!
!  Populate wavenumber arrays for fft and calculate nyquist wavenumber.
!
  if (nxgrid/=1) then
    kx_fft=cshift((/(i-(nxgrid+1)/2,i=0,nxgrid-1)/),+(nxgrid+1)/2)*2*pi/Lx
    kx_ny = nxgrid/2 * 2*pi/Lx
  else
    kx_fft=0.0
    kx_ny = 0.0
  endif
!
  if (nygrid/=1) then
    ky_fft=cshift((/(i-(nygrid+1)/2,i=0,nygrid-1)/),+(nygrid+1)/2)*2*pi/Ly
    ky_ny = nygrid/2 * 2*pi/Ly
  else
    ky_fft=0.0
    ky_ny = 0.0
  endif
!
  if (nzgrid/=1) then
    kz_fft=cshift((/(i-(nzgrid+1)/2,i=0,nzgrid-1)/),+(nzgrid+1)/2)*2*pi/Lz
    ky_ny = nzgrid/2 * 2*pi/Lz
  else
    kz_fft=0.0
    kz_ny = 0.0
  endif
!
!  Parameter dependent initialization of module variables and final
!  pre-timestepping setup (must be done before need_XXXX can be used, for
!  example).
!
  call initialize_modules(f,lstarting=.true.)
  call particles_initialize_modules(f,lstarting=.true.)
!
!  Initial conditions: by default, we put f=0 (ss=lnrho=uu=0, etc).
!  alternatively: read existing snapshot and overwrite only some fields
!  by the various procedures below.
!
  if (lread_oldsnap) then
    call rsnap(trim(directory_snap)//'/var.dat',f, mvar)
  else
!
! We used to have just f=0. here, but with GRAVITY=gravity_r,
! the gravitational acceleration (which is computed in
! initialize_gravity and stored in the f-array), is set to zero
! by the statement f=0. After running start.csh, this can lead
! to some confusion as to whether the gravity module actually
! does anything or not.
!
!   So now we are more specific:
    f(:,:,:,1:mvar)=0.
  endif
!
!  The following init routines do then only need to add to f.
!  wd: also in the case where we have read in an existing snapshot??
!
  if (lroot) print* !(empty line)
  do i=1,init_loops
    if (lroot .and. init_loops/=1) &
        print '(A33,i3,A25)', 'start: -- performing loop number', i, &
        ' of initial conditions --'
    call init_gg        (f)
    call init_uu        (f)
    call init_lnrho     (f)
    call init_ss        (f)
    call init_aa        (f)
    call init_lorenz_gauge (f)
    call init_aatest    (f)
    call init_uutest    (f)
    call init_rad       (f)
    call init_lncc      (f)
    call init_chiral    (f)
    call init_chemistry (f)
    call init_uud       (f)
    call init_nd        (f)
    call init_uun       (f)
    call init_lnrhon    (f)
    call init_ecr       (f)
    call init_fcr       (f)
    call init_interstellar (f)
    call init_solid_cells(f)
    call init_special   (f)
  enddo
!
  if (lparticles) call particles_init(f)
!
!  If requested, write original stratification to file.
!
  if (lwrite_stratification) then
    call update_ghosts(f)
    open(19,file=trim(directory_snap)//'/stratification.dat')
      write(19,*) f(l1,m1,:,ilnrho)
    close(19)
  endif
!
!  Check whether we want ionization.
!
  if (leos_ionization) call ioninit(f)
  if (leos_temperature_ionization) call ioncalc(f)
  if (lradiation_ray) call radtransfer(f)
!
!  Filter initial velocity.
!  NOTE: this procedure is currently not very efficient,
!  because for all variables boundary conditions and communication
!  are done again and again for each variable.
!
  if (nfilter/=0) then
    do ifilter=1,nfilter
      call rmwig(f,df,iux,iuz,awig)
    enddo
    if (lroot) print*,'DONE: filter initial velocity, nfilter=',nfilter
  endif
!
!  Calculate the potential of the self-gravity (mostly for debugging).
!
  call calc_selfpotential(f)
!
!  For sixth order momentum-conserving, symmetric hyperviscosity with positive
!  definite heating rate we need to precalculate the viscosity term. The 
!  restivitity term for sixth order hyperresistivity with positive definite
!  heating rate must also be precalculated.
!
  if (lhyperviscosity_strict)   call hyperviscosity_strict(f)
  if (lhyperresistivity_strict) call hyperresistivity_strict(f)
!
!  Set random seed independent of processor.
!
  seed(1) = 1812
  call random_seed_wrapper(PUT=seed)
!
!  Write to disk.
!
  if (lwrite_ic) then
    call wsnap(trim(directory_snap)//'/VAR0',f, &
        mvar_io,ENUM=.false.,FLIST='varN.list')
    if (lparticles) &
        call particles_write_snapshot(trim(directory_snap)//'/PVAR0', &
        ENUM=.false.,FLIST='pvarN.list')
    if (lparticles_nbody.and.lroot) &
        call particles_nbody_write_snapshot(&
        trim(datadir)//'/proc0/SPVAR0', &
        ENUM=.false.,FLIST='spvarN.list')
  endif
!
!  The option lnowrite writes everything except the actual var.dat file.
!  This can be useful if auxiliary files are outdated, and don't want
!  to overwrite an existing var.dat.
!
  inquire(FILE="NOERASE", EXIST=lnoerase)
  if (.not.lnowrite .and. .not.lnoerase) then
    if (ip<12) print*,'START: writing to ' // trim(directory_snap)//'/var.dat'
    call wsnap(trim(directory_snap)//'/var.dat',f,mvar_io,ENUM=.false.)
    if (lparticles) &
        call particles_write_snapshot(trim(directory_snap)//'/pvar.dat', &
        ENUM=.false.)
    if (lparticles_nbody.and.lroot) &
         call particles_nbody_write_snapshot(&
         trim(datadir)//'/proc0/spvar.dat', ENUM=.false.)
    call wtime(trim(directory)//'/time.dat',t)
  endif
  call wdim(trim(directory)//'/dim.dat')
!
!  Also write full dimensions to data/.
!
  if (lroot) then
    call wdim(trim(datadir)//'/dim.dat', &
      nxgrid+2*nghost,nygrid+2*nghost,nzgrid+2*nghost)
    if (lparticles) &
      call particles_write_pdim(trim(datadir)//'/pdim.dat')
    if (lparticles_nbody) &
      call particles_nbody_write_spdim(trim(datadir)//'/spdim.dat')
  endif
!
!  Write global variables.
!
  if (mglobal/=0)  &
      call output_globals(trim(directory_snap)//'/global.dat', &
      f(:,:,:,mvar+maux+1:mvar+maux+mglobal),mglobal)
!
!  Write input parameters to a parameter file (for run.x and IDL).
!  Do this late enough, so init_entropy etc. can adjust them.
!
  call wparam()
!
!  Write information about pencils to disc.
!
  call write_pencil_info()
!
  call mpifinalize
!
  if (lroot) print*
  if (lroot) print*,'start.x has completed successfully'
  if (lroot) print* ! (finish with an empty line)
!
!  Free any allocated memory.
!
  call farray_clean_up()
  call sharedvars_clean_up()
!
endprogram
