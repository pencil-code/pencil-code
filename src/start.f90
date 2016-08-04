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
  use General
  use Boundcond,        only: update_ghosts
  use Cdata
  use Chemistry,        only: init_chemistry
  use Chiral,           only: init_chiral
  use Cosmicrayflux,    only: init_fcr
  use Cosmicray,        only: init_ecr
  use Density,          only: init_lnrho
  use Diagnostics

  use Dustdensity,      only: init_nd
  use Dustvelocity,     only: init_uud
  use Energy,           only: init_energy
  use EquationOfState
  use FArrayManager,    only: farray_clean_up
  use Filter
  use Gravity,          only: init_gg
  use Grid
  use Hydro,            only: init_uu
  use Hyperresi_strict, only: hyperresistivity_strict
  use Hypervisc_strict, only: hyperviscosity_strict
  use Initcond
  use InitialCondition, only: initial_condition_all, initial_condition_clean_up
  use Interstellar,     only: init_interstellar
  use IO,               only: wgrid, directory_names, wproc_bounds, output_globals
  use Lorenz_gauge,     only: init_lorenz_gauge
  use Magnetic,         only: init_aa
  use Messages
  use Mpicomm
  use NeutralDensity,   only: init_lnrhon
  use NeutralVelocity,  only: init_uun
  use Param_IO
  use Particles_main
  use Polymer,          only: init_poly
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
  real, allocatable, dimension (:,:,:,:) :: f, df
  integer :: i, ifilter, stat
  logical :: lnoerase=.false.
  real :: dang
!
  lstart = .true.
!
!  Check processor layout, get processor numbers and define whether we are root.
!
  call mpicomm_init
!
!  Identify version.
!
  if (lroot) call svn_id( &
      '$Id$')
!
!  Initialize the message subsystem, eg. color setting etc.
!
  call initialize_messages
!
!  Allocate large arrays. We need to make them allocatable in order to
!  avoid segfaults at 128^3 (7 variables) with Intel compiler on 32-bit
!  Linux boxes. Not clear why they crashed (we _did_ increase stacksize
!  limits), but they did and now don't. Also, the present approach runs
!  up to nx=ny=nz=135, but not for even slightly larger grids.
!
  allocate( f(mx,my,mz,mfarray),STAT=stat)
  if (stat>0) call fatal_error('start','Could not allocate memory for f')
  allocate(df(mx,my,mz,mvar)   ,STAT=stat)
  if (stat>0) call fatal_error('start','Could not allocate memory for df')
!
!  Pre-initialize f and df to absurd value (to crash the code should we
!  later use uninitialized slots of those fields).
!
  f =huge(1.0)
  df=huge(1.0)
!
!  Calculate dimensionality of the run.
!
  dimensionality=min(nxgrid-1,1)+min(nygrid-1,1)+min(nzgrid-1,1)
!
!  Define the lenergy logical
!
  lenergy=lentropy.or.ltemperature.or.lthermal_energy
!
!  Read initialization parameters from "start.in".
!
  call read_all_init_pars
!
  call set_coords_switches
!
!  Set up directory names and check whether the directories exist.
!
  call directory_names
!
!  Initialise MPI communication.
!
  call initialize_mpicomm
!
!  Register variables in the f array.
!
  call register_modules
  if (lparticles) call particles_register_modules
!
!  Call rprint_list to initialize diagnostics and write indices to file.
!
  call rprint_list(.false.)
  if (lparticles) call particles_rprint_list(.false.)
!
!  The logical headtt is sometimes referred to in start.x, even though it is
!  not yet defined. So we set it simply to lroot here.
!
  headtt=lroot
!
!  Initialize start time.
!
  t=tstart
!
!  Will we write all slots of f?
!
  if (lwrite_aux) then
    mvar_io=mvar+maux
  else
    mvar_io=mvar
  endif
!
!  Print resolution.
!
  if (lroot) print*, 'nxgrid, nygrid, nzgrid=', nxgrid, nygrid, nzgrid
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
!
!  FG: force theta coordinate to spam 0:pi for periodic across pole 
!
        elseif (lpole(i)) then
          if (lperi(i)) call fatal_error('start',&
            'lperi and lpole cannot be used together in same component')
          if (.not.lspherical_coords) call fatal_error('start',&
            'lpole only implemented for spherical coordinates')

          if (i==2) then
            xyz0(i) = 0.
            Lxyz(i) = pi 
            xyz0(3) = 0.
            Lxyz(3) = 2.*pi
          else
            call fatal_error('start',&
                'origin/pole not included for components or coordinates')
          endif
        else
          Lxyz(i)=2.*pi    ! default value
        endif
      else
        if (lpole(i)) then ! overwirte xyz0 and xyz1 to 0:pi 
          if (lperi(i)) call fatal_error('start',&
            'lperi and lpole cannot be used together in same component')
          if (.not. lspherical_coords) call fatal_error('start',&
            'lpole only implemented for spherical coordinates')
          if (i==2) then
            xyz0(i) = 0.
            Lxyz(i) = pi 
            xyz0(3) = 0.
            xyz1(3) = 2.*pi
          else
            call fatal_error('start',&
                'origin/pole not included for components or coordinates')
          endif
        else
          Lxyz(i)=xyz1(i)-xyz0(i)
        endif
      endif
    else                            ! Lxyz was set
      if (xyz1(i)/=impossible) then ! both Lxyz and xyz1 are set
        call fatal_error('start','Cannot set Lxyz and xyz1 at the same time')
      endif
    endif
  enddo

  if (lyinyang) then
    if (lroot) &
      print*, 'Setting latitude and longitude intervals for Yin-Yang grid, ignoring input'
!
! Min(dy,dz) put between Yin and Yang grid at closest distance to minimize overlap.
!
    dang=.999*min(1./nygrid,3./nzgrid)*0.5*pi      ! only valid for equidistant grid!!
    xyz0(2:3) = (/ 1./4., 1./4. /)*pi+0.5*dang
    Lxyz(2:3) = (/ 1./2., 3./2. /)*pi-dang
  endif

  xyz1=xyz0+Lxyz
!
!  Abbreviations
!
  x0 = xyz0(1)
  y0 = xyz0(2)
  z0 = xyz0(3)
  Lx = Lxyz(1)
  Ly = Lxyz(2)
  Lz = Lxyz(3)
!
!  Position of equator (if any).
!
  if (lequatory) yequator=xyz0(2)+0.5*Lxyz(2)
  if (lequatorz) zequator=xyz0(3)+0.5*Lxyz(3)
!
!  Set up limits of averaging if needed.
!
  if (lav_smallx) call init_xaver
!
!  Check consistency.
!
  if (.not.lperi(1).and.nxgrid<2) &
      call fatal_error('start','for lperi(1)=F: must have nxgrid>1')
  if (.not.lperi(2).and.nygrid<2) &
      call fatal_error('start','for lperi(2)=F: must have nygrid>1')
  if (.not.lperi(3).and.nzgrid<2) &
      call fatal_error('start','for lperi(3)=F: must have nzgrid>1')
!
!  Initialise random number generator in processor-dependent fashion for
!  random initial data.
!  Slightly tricky, since setting seed=(/iproc,0,0,0,0,0,0,0,.../)
!  would produce very short-period random numbers with the Intel compiler;
!  so we need to retain most of the initial entropy of the generator.
!
  call get_nseed(nseed)   ! get state length of random number generator
  call random_seed_wrapper(GET=seed)
!
!  Different initial seed (seed0) and random numbers on different CPUs
!  The default is seed0=1812 for some obscure Napoleonic reason
!
  seed(1)=-((seed0-1812+1)*10+iproc_world)
  call random_seed_wrapper(PUT=seed)
!
!  Generate grid and initialize specific grid variables.
!
  call construct_grid(x,y,z,dx,dy,dz)
!
!  Size of box at local processor. The if-statement is for
!  backward compatibility.
!  MR: the following code doubled in run.f90. Perhaps to be put in construct_grid?
!
  if (lequidist(1)) then
    Lxyz_loc(1) = Lxyz(1)/nprocx
    xyz0_loc(1) = xyz0(1)+ipx*Lxyz_loc(1)
    xyz1_loc(1) = xyz0_loc(1)+Lxyz_loc(1)
  else
    !
    !  In the equidistant grid, the processor boundaries (xyz[01]_loc) do NOT 
    !  coincide with the l[mn]1[2] points. Also, xyz0_loc[ipx+1]=xyz1_loc[ipx], i.e.,
    !  the inner boundary of one is exactly the outer boundary of the other. Reproduce
    !  this behavior also for non-equidistant grids. 
    !
    if (ipx==0) then 
      xyz0_loc(1) = x(l1)
    else
      xyz0_loc(1) = x(l1) - .5/dx_1(l1)
    endif
    if (ipx==nprocx-1) then 
      xyz1_loc(1) = x(l2)
    else
      xyz1_loc(1) = x(l2+1) - .5/dx_1(l2+1)
    endif
    Lxyz_loc(1) = xyz1_loc(1) - xyz0_loc(1)
  endif  
!
  if (lequidist(2)) then
    Lxyz_loc(2) = Lxyz(2)/nprocy
    xyz0_loc(2) = xyz0(2)+ipy*Lxyz_loc(2)
    xyz1_loc(2) = xyz0_loc(2)+Lxyz_loc(2)
  else
    if (ipy==0) then
      xyz0_loc(2) = y(m1) 
    else
      xyz0_loc(2) = y(m1) - .5/dy_1(m1)
    endif
    if (ipy==nprocy-1) then
      xyz1_loc(2) = y(m2)
    else
      xyz1_loc(2) = y(m2+1) - .5/dy_1(m2+1)
    endif
    Lxyz_loc(2) = xyz1_loc(2) - xyz0_loc(2)
  endif
!
  if (lequidist(3)) then 
    Lxyz_loc(3) = Lxyz(3)/nprocz
    xyz0_loc(3) = xyz0(3)+ipz*Lxyz_loc(3)
    xyz1_loc(3) = xyz0_loc(3)+Lxyz_loc(3)
  else
    if (ipz==0) then
      xyz0_loc(3) = z(n1) 
    else
      xyz0_loc(3) = z(n1) - .5/dz_1(n1)
    endif
    if (ipz==nprocz-1) then
      xyz1_loc(3) = z(n2)
    else
      xyz1_loc(3) = z(n2+1) - .5/dz_1(n2+1)
    endif
    Lxyz_loc(3) = xyz1_loc(3) - xyz0_loc(3)
  endif
!
!  Write grid.dat file.
!
  call wgrid('grid.dat')
  if (lparticles) call wproc_bounds(trim(directory)//'/proc_bounds.dat')
!
!  Update the list of neighboring processes.
!
  call update_neighbors
!
!  Write .general file for data explorer.
!
  if (lroot) call write_dx_general(trim(datadir)//'/var.general', &
      x0-nghost*dx,y0-nghost*dy,z0-nghost*dz)
!
!  Populate wavenumber arrays for fft and calculate nyquist wavenumber.
!
  if (nxgrid/=1) then
    kx_fft=cshift((/(i-(nxgrid+1)/2,i=0,nxgrid-1)/),+(nxgrid+1)/2)*2*pi/Lx
    kx_nyq=nxgrid/2 * 2*pi/Lx
  else
    kx_fft=0.0
    kx_nyq=0.0
  endif
!
  if (nygrid/=1) then
    ky_fft=cshift((/(i-(nygrid+1)/2,i=0,nygrid-1)/),+(nygrid+1)/2)*2*pi/Ly
    ky_nyq=nygrid/2 * 2*pi/Ly
  else
    ky_fft=0.0
    ky_nyq=0.0
  endif
!
  if (nzgrid/=1) then
    kz_fft=cshift((/(i-(nzgrid+1)/2,i=0,nzgrid-1)/),+(nzgrid+1)/2)*2*pi/Lz
    kz_nyq=nzgrid/2 * 2*pi/Lz
  else
    kz_fft=0.0
    kz_nyq=0.0
  endif
!
!  Set random seed independent of processor prior to initial conditions.
!  Do this only if seed0 is modified from its original value.
!
  if (seed0/=1812) then
    seed(1)=seed0
    call random_seed_wrapper(PUT=seed)
  endif
!
!  Parameter dependent initialization of module variables and final
!  pre-timestepping setup (must be done before need_XXXX can be used, for
!  example).
!
  call initialize_modules(f)
!
!  Initial conditions: by default, we put f=0 (ss=lnrho=uu=0, etc).
!  alternatively: read existing snapshot and overwrite only some fields
!  by the various procedures below.
!
!  If lmodify is true, read var.dat into auxiliary array df
!  and read modify_filename into f
!
  if (lmodify) then
    call rsnap('var.dat',df,mvar,lread_nogrid)
    call rsnap(modify_filename,f,mvar,lread_nogrid)
  else
    if (lread_oldsnap) then
      call rsnap('var.dat',f,mvar,lread_nogrid)
    else
      f(:,:,:,1:mvar)=0.0
    endif
  endif
!
  if (lyinyang.and.lroot) &
    call warning('start','Any initial condition depending on y or z will not be set correctly on Yin-Yang grid.')
!
!  The following init routines only need to add to f.
!  wd: also in the case where we have read in an existing snapshot??
!
  if (lroot) print*
  do i=1,init_loops
    if (lroot .and. init_loops/=1) &
        print '(A33,i3,A25)', 'start: -- performing loop number', i, &
        ' of initial conditions --'
    call init_gg(f)
!
!  This is a temporary solution for calculating the correct initial
!  condition for anelastic case where we need to set div(rho u)=0
!  The dust-vortex auto-test requires that uu is initialised before
!  lnrho
! DM: it seems to me that the sequence of calls below is important.
! so this is just a caution : Please do not modify the sequence of
! calls to 'init' routines below.
!
    if (lanelastic) then
      call init_lnrho(f)
      call init_uu(f)
    else
      call init_uu(f)
      call init_lnrho(f)
    endif
    call init_energy(f)
    call init_aa(f)
    call init_lorenz_gauge(f)
    call init_poly(f)
    call init_aatest(f)
    call init_uutest(f)
    call init_rad(f)
    call init_lncc(f)
    call init_chiral(f)
    call init_chemistry(f)
    call init_uud(f)
    call init_nd(f)
    call init_uun(f)
    call init_lnrhon(f)
    call init_ecr(f)
    call init_fcr(f)
    call init_interstellar(f)
    call init_solid_cells(f)
    call init_special(f)
  enddo
!
!  If desired, the f array can be initialized in one call.
!
  if (linitial_condition) call initial_condition_all(f)
!
!  Initialize particle modules.
!
  if (lparticles) then 
    call particles_initialize_modules(f)
    call particles_init(f)
  endif
!
!  If requested, write original (z-dependent) stratification to file.
!
  if (lwrite_stratification) then
    call update_ghosts(f,ilnrho)
    open(19,file=trim(directory_dist)//'/stratification.dat')
      write(19,*) f(l1,m1,:,ilnrho)
    close(19)
  endif
!
!  Check whether we want ionization.
!
  if (leos_ionization) call ioninit(f)
  if (leos_temperature_ionization) call ioncalc(f)
  if (lradiation_ray) then
    call update_ghosts(f)
    call radtransfer(f)
  endif
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
!  Set random seed independent of processor after initial conditions.
!  Do this only if seed0 is not already changed from its original value.
!
  if (seed0==1812) then
    seed(1)=seed0
    call random_seed_wrapper(PUT=seed)
  endif
!
!  Final update of ghost cells, after that 'f' must not be altered, anymore.
!
  call update_ghosts(f)
!
!  Write initial condition to disk.
!
  if (lwrite_ic) then
    if (lparticles) &
        call write_snapshot_particles(directory_dist,f,ENUM=.false.,snapnum=0)
!
    call wsnap('VAR0',f,mvar_io,ENUM=.false.,FLIST='varN.list')
  endif
!
!  The option lnowrite writes everything except the actual var.dat file.
!  This can be useful if auxiliary files are outdated, and don't want
!  to overwrite an existing var.dat.
!
  lnoerase = control_file_exists("NOERASE")
  if (.not.lnowrite .and. .not.lnoerase) then
    if (ip<12) print*,'START: writing to '//trim(directory_snap)//'/var.dat'
    if (lparticles) &
        call write_snapshot_particles(directory_dist,f,ENUM=.false.)
    call wsnap('var.dat',f,mvar_io,ENUM=.false.)
  elseif (lmodify) then
    call wsnap(modify_filename,f,mvar_io,ENUM=.false.)
  endif
  call wdim(trim(directory)//'/dim.dat')
!
!  Also write full dimensions to data/.
!
  if (lroot) then
    call wdim(trim(datadir)//'/dim.dat', &
        nxgrid+2*nghost,nygrid+2*nghost,nzgrid+2*nghost,lglobal=.true.)
    if (lparticles) call write_dim_particles(trim(datadir))
  endif
!
!  Write global variables.
!
  if (mglobal/=0) &
      call output_globals('global.dat', &
          f(:,:,:,mvar+maux+1:mvar+maux+mglobal),mglobal)
!
!  Write input parameters to a parameter file (for run.x and IDL).
!  Do this late enough, so init_entropy etc. can adjust them.
!
  call write_all_init_pars('IDL')
!
!  Write information about pencils to disc.
!
  call write_pencil_info
!
!  Gvie all modules the possibility to exit properly.
!
  call finalize_modules(f)
!
!  Stop MPI.
!
  call mpifinalize
!
!  Free any allocated memory.
!  MR: do we need all these cleanups? The process ends immediately afterwards and all memory is released anyway.
!
  call farray_clean_up
  call sharedvars_clean_up
  call initial_condition_clean_up
  if (lparticles) call particles_cleanup
!
!  Deallocate the arrays allocated for
!  1-D and 2-D diagnostics.
!
  call fnames_clean_up
  call xyaverages_clean_up
  call xzaverages_clean_up
  call yzaverages_clean_up
  if (lwrite_yaverages)    call yaverages_clean_up
  if (lwrite_zaverages)    call zaverages_clean_up
  if (lwrite_phiaverages)  call phiaverages_clean_up
  if (lwrite_phizaverages) call phizaverages_clean_up
  if (lwrite_phiaverages)  call phiaverages_clean_up
  if (lwrite_sound)        call sound_clean_up
!
! announce completion.
!
  if (lroot) print*
  if (lroot) print*,'start.x has completed successfully'
  if (lroot) print* ! (finish with an empty line)
!
endprogram
