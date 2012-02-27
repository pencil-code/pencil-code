! $Id$
!
!  IO of init and run parameters. Subroutines here are `at the end of the
!  food chain', i.e. depend on all physics modules plus possibly others.
!  Using this module is also a compact way of referring to all physics
!  modules at once.
!
module Param_IO
!
  use Cdata
  use Chemistry
  use Chiral
  use Cosmicray
  use CosmicrayFlux
  use Cparam
  use Density
  use Dustdensity
  use Dustvelocity
  use Entropy
  use EquationOfState
  use Forcing
  use General
  use Gravity
  use Hydro
  use InitialCondition
  use Interstellar
  use Lorenz_gauge
  use Magnetic
  use Messages
  use NeutralDensity
  use NeutralVelocity
  use NSCBC
  use Particles_main
  use Poisson
  use Polymer
  use Power_spectrum
  use Pscalar
  use Radiation
  use Selfgravity
  use Shear
  use Signal_handling
  use Shock
  use Solid_Cells
  use Special
  use Sub
  use Testfield
  use Testflow
  use TestPerturb
  use Testscalar
  use Timeavg
  use Viscosity
!
  implicit none
!
  private
!
  public :: get_datadir, get_snapdir
  public :: read_startpars, print_startpars
  public :: read_runpars,   print_runpars
  public :: rparam, wparam, wparam2, write_pencil_info
!
! The following fixes namelist problems withi MIPSpro 7.3.1.3m
! under IRIX -- at least for the moment
!
  character (len=labellen) :: mips_is_buggy='system'
!
  namelist /init_pars/ &
      cvsid, ip, xyz0, xyz1, Lxyz, lperi, lshift_origin, coord_system, &
      lequidist, coeff_grid, zeta_grid0, grid_func, xyz_star, lwrite_ic, &
      lnowrite, luniform_z_mesh_aspect_ratio, unit_system, unit_length, &
      lmodify,modify_filename, &
      unit_velocity, unit_density, unit_temperature, unit_magnetic, c_light, &
      G_Newton, hbar, random_gen, seed0, nfilter, lserial_io, der2_type, &
      lread_oldsnap, lread_oldsnap_nomag, lread_oldsnap_nopscalar, &
      lread_oldsnap_notestfield, lread_oldsnap_notestscalar, &
      lread_aux, lwrite_aux, pretend_lnTT, lprocz_slowest, &
      lcopysnapshots_exp, bcx, bcy, bcz, r_int, r_ext, r_ref, rsmooth, &
      r_int_border, r_ext_border, mu0, force_lower_bound, force_upper_bound, &
      tstart, lseparate_persist, &
      fbcx1, fbcx2, fbcx1_2, fbcx2_2, &
      fbcy1, fbcy2, fbcy1_1, fbcy1_2, fbcy2_1, fbcy2_2, &
      fbcz1, fbcz2, fbcz1_1, fbcz1_2, fbcz2_1, fbcz2_2, &
      fbcx_bot, fbcx_top, fbcy_bot, fbcy_top, fbcz_bot, fbcz_top, &
      xyz_step, xi_step_frac, xi_step_width, dxi_fact, trans_width, &
      lcylinder_in_a_box, lsphere_in_a_box, llocal_iso, init_loops, lwrite_2d, &
      lcylindrical_gravity, &
      border_frac_x, border_frac_y, border_frac_z, lborder_hyper_diff, &
      luse_latitude, lshift_datacube_x, lfargo_advection, yequator, lequatory, &
      lequatorz, zequator, lav_smallx, xav_max, niter_poisson, &
      lconst_advection, u0_advec
!
  namelist /run_pars/ &
      cvsid, ip, nt, it1, it1d, dt, cdt, ddt, cdtv, cdtv2, cdtv3, cdts, cdtr, &
      cdtc, isave, itorder, dsnap, d2davg, dvid, dsound, dtmin, dspec, tmax, iwig, &
      awig, ialive, max_walltime, dtmax, ldt_paronly, vel_spec, mag_spec, &
      uxy_spec, bxy_spec, jxbxy_spec, xy_spec, oo_spec, &
      uxj_spec, vec_spec, ou_spec, ab_spec, azbz_spec, ub_spec, &
      vel_phispec, mag_phispec, &
      uxj_phispec, vec_phispec, ou_phispec, ab_phispec, EP_spec, ro_spec, &
      TT_spec, ss_spec, cc_spec, cr_spec, isaveglobal, lr_spec, r2u_spec, &
      r3u_spec, rhocc_pdf, cc_pdf, lncc_pdf, gcc_pdf, lngcc_pdf, kinflow, &
      lkinflow_as_aux, ampl_kinflow_x, ampl_kinflow_y, ampl_kinflow_z, &
      kx_kinflow, ky_kinflow, kz_kinflow, dtphase_kinflow, &
      random_gen, der2_type, lrmwig_rho, lrmwig_full, lrmwig_xyaverage, &
      ltime_integrals, lnowrite, noghost_for_isave, lwrite_yaverages, &
      lwrite_zaverages, lwrite_phiaverages, lwrite_slices, test_nonblocking, &
      lread_oldsnap_nomag, lread_oldsnap_nopscalar, &
      lread_oldsnap_notestfield, lread_oldsnap_notestscalar, &
      lread_aux, comment_char, ix, iy, iz, iz2, iz3, iz4, slice_position, &
      xbot_slice, xtop_slice, ybot_slice, ytop_slice, zbot_slice, ztop_slice, &
      bcx, bcy, bcz, r_int, r_ext, r_int_border, &
      r_ext_border, lfreeze_varsquare, lfreeze_varint, lfreeze_varext, &
      xfreeze_square, yfreeze_square, rfreeze_int, rfreeze_ext, wfreeze, &
      wfreeze_int, wfreeze_ext, wborder, wborder_int, wborder_ext, tborder, &
      fshift_int, fshift_ext, &
      lreset_tstart, tstart, lseparate_persist, &
      fbcx1, fbcx2, fbcx1_2, fbcx2_2, &
      fbcy1, fbcy2, fbcy1_1, fbcy1_2, fbcy2_1, fbcy2_2, &
      fbcz1, fbcz2, fbcz1_1, fbcz1_2, fbcz2_1, fbcz2_2, &
      fbcx_bot, fbcx_top, fbcy_bot, fbcy_top, fbcz_bot, fbcz_top, &
      Udrift_bc, ttransient, tavg, idx_tavg, lserial_io, nr_directions, &
      lsfu, lsfb, lsfz1, lsfz2, lsfflux, lpdfu, lpdfb, lpdfz1, lpdfz2, &
      onedall, pretend_lnTT, old_cdtv, lmaxadvec_sum, save_lastsnap, &
      lwrite_aux, lwrite_dvar, force_lower_bound, force_upper_bound, &
      oned, twod, lpoint, mpoint, npoint, lpoint2, mpoint2, npoint2, &
      border_frac_x, border_frac_y, border_frac_z, &
      lcylinder_in_a_box, lsphere_in_a_box, ipencil_swap, &
      lpencil_requested_swap, lpencil_diagnos_swap, lpencil_check, &
      lpencil_check_small, lrandom_f_pencil_check, lpencil_check_diagnos_opti, &
      lpencil_init, penc0, lwrite_2d, lbidiagonal_derij, lisotropic_advection, &
      crash_file_dtmin_factor, niter_poisson, lADI, ltestperturb, eps_rkf, &
      eps_stiff, timestep_scaling, lequatory, lequatorz, zequator, &
      lini_t_eq_zero, lav_smallx, xav_max, ldt_paronly, lweno_transport, &
      it_timing, har_spec, hav_spec, j_spec, jb_spec, lread_less, lformat, ltec, &
      llsode, lsplit_second, nu_sts,lsuper_time_stepping, lfargo_advection, &
      ldynamical_diffusion, re_mesh, lconst_advection, u0_advec, &
      loutput_varn_at_exact_tsnap, &
      lstop_on_ioerror, mailaddress, lenforce_redundant
!
  contains
!***********************************************************************
    subroutine get_datadir(dir)
!
!  Overwrite datadir from datadir.in, if that exists.
!
!   2-oct-02/wolf: coded
!  25-oct-02/axel: default is taken from cdata.f90 where it's defined
!
      use Mpicomm, only: mpibcast_logical, mpibcast_char
!
      character (len=*) :: dir
      logical :: exists
      integer :: unit=1
!
!  Let root processor check for existence of datadir.in.
!
      if (lroot) then
        inquire(FILE='datadir.in',EXIST=exists)
        if (exists) then
          open(unit,FILE='datadir.in',FORM='formatted')
          read(unit,'(A)') dir
          close(unit)
        endif
      endif
!
!  Tell other processors whether we need to communicate dir (i.e. datadir).
!
      call mpibcast_logical(exists, 1)
!
!  Let root processor communicate dir (i.e. datadir) to all other processors.
!
      if (exists) call mpibcast_char(dir, len(dir))
!
    endsubroutine get_datadir
!***********************************************************************
    subroutine get_snapdir(dir)
!
!  Read directory_snap from data/directory_snap, if that exists
!  wd: I think we should unify these into a subroutine
!      `overwrite_string_from_file(dir,file,label[optional])'
!
!   2-nov-02/axel: adapted from get_datadir
!
      character (len=*) :: dir
      character (len=fnlen) :: tdir
      logical :: exists
      integer :: unit=1
!
!  check for existence of `data/directory_snap'
!
      inquire(FILE=trim(datadir)//'/directory_snap',EXIST=exists)
      if (exists) then
        open(unit,FILE=trim(datadir)//'/directory_snap',FORM='formatted')
        read(unit,'(A)') tdir
        close(unit)
        if (len(trim(tdir)) > 0) call parse_shell(tdir,dir)
      endif
      if (lroot.and.ip<10) print*,'get_snapdir: dir=',trim(dir)
!
    endsubroutine get_snapdir
!***********************************************************************
    subroutine read_startpars(print,file)
!
!  read input parameters (done by each processor)
!
!  Now no warning is produced when somebody sets lcylindrical.
!  read_startpars: lcylindrical=T is now outdated'
!  use instead: lcylinder_in_a_box=T'
!  This renaming became necessary with the development of'
!  cylindrical coordinates which led to very similar names'
!  (coord_system="cylindrical_coords")'
!
!   6-jul-02/axel: in case of error, print sample namelist
!  21-oct-03/tony: moved sample namelist stuff to a separate procedure
!
      use Mpicomm, only: parallel_open, parallel_close
!
      logical, optional :: print,file
!
      integer :: ierr, unit=1
!
!  Set default to shearing sheet if lshear=.true. (even when Sshear==0.).
!
      if (lshear) bcx(:)='she'
!
!  Open namelist file.
!
      call parallel_open(unit,'start.in','formatted')
!
!  Read through all items that *may* be present in the various modules.
!
      read(unit,NML=init_pars,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('init_pars',ierr)
      rewind(unit)
!
      call read_initial_condition_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('initial_condition_pars',ierr)
      rewind(unit)
!
      call read_eos_init_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('eos_init_pars',ierr)
      rewind(unit)
!
      call read_hydro_init_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('hydro_init_pars',ierr)
      rewind(unit)
!
      call read_density_init_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('density_init_pars',ierr)
      rewind(unit)
!
      call read_forcing_init_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('forcing_init_pars',ierr)
      rewind(unit)
!
      call read_gravity_init_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('grav_init_pars',ierr)
      rewind(unit)
!
      call read_selfgravity_init_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('selfgrav_init_pars',ierr)
      rewind(unit)
!
      call read_poisson_init_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('poisson_init_pars',ierr)
      rewind(unit)
!
      call read_entropy_init_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('entropy_init_pars',ierr)
      rewind(unit)
!
      call read_magnetic_init_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('magnetic_init_pars',ierr)
      rewind(unit)
!
      call read_lorenz_gauge_init_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('lorenz_gauge_init_pars',ierr)
      rewind(unit)
!
      call read_testscalar_init_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('testscalar_init_pars',ierr)
      rewind(unit)
!
      call read_testfield_init_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('testfield_init_pars',ierr)
      rewind(unit)
!
      call read_testflow_init_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('testflow_init_pars',ierr)
      rewind(unit)
!
      call read_radiation_init_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('radiation_init_pars',ierr)
      rewind(unit)
!
      call read_pscalar_init_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('pscalar_init_pars',ierr)
      rewind(unit)
!
      call read_chiral_init_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('chiral_init_pars',ierr)
      rewind(unit)
!
      call read_chemistry_init_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('chemistry_init_pars',ierr)
      rewind(unit)
!
      call read_signal_init_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('signal_init_pars',ierr)
      rewind(unit)
!
      call read_dustvelocity_init_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('dustvelocity_init_pars',ierr)
      rewind(unit)
!
      call read_dustdensity_init_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('dustdensity_init_pars',ierr)
      rewind(unit)
!
      call read_neutralvelocity_init_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('neutralvelocity_init_pars',ierr)
      rewind(unit)
!
      call read_neutraldensity_init_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('neutraldensity_init_pars',ierr)
      rewind(unit)
!
      call read_cosmicray_init_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('cosmicray_init_pars',ierr)
      rewind(unit)
!
      call read_cosmicrayflux_init_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('cosmicrayflux_init_pars',ierr)
      rewind(unit)
!
      call read_interstellar_init_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('interstellar_init_pars',ierr)
      rewind(unit)
!
      call read_shear_init_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('shear_init_pars',ierr)
      rewind(unit)
!
      call read_testperturb_init_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('testperturb_init_pars',ierr)
      rewind(unit)
!
      call read_viscosity_init_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('viscosity_init_pars',ierr)
      rewind(unit)
!
      call read_special_init_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('special_init_pars',ierr)
      rewind(unit)
!
      call read_shock_init_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('shock_init_pars',ierr)
      rewind(unit)
!
      call read_solid_cells_init_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('solid_cells_init_pars',ierr)
      rewind(unit)
!
      call read_NSCBC_init_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('NSCBC_init_pars',ierr)
      rewind(unit)
!
      call read_polymer_init_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('polymer_init_pars',ierr)
      rewind(unit)
!
      call particles_read_startpars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('particles_init_pars_wrap',ierr)
      rewind(unit)
!
      call parallel_close(unit)
!
!  Print SVN id from first line.
!
      if (lroot) call svn_id(cvsid)
!
!  Give online feedback if called with the PRINT optional argument.
!  Note: Some compiler's [like Compaq's] code crashes with the more
!  compact `if (present(print) .and. print)'.
!
      if (present(print)) then
        if (print) then
          call print_startpars()
        endif
      endif
!
!  Write parameters to log file.
!
      if (present(file)) then
        if (file) then
          call print_startpars(FILE=trim(datadir)//'/params.log')
        endif
      endif
!
!  Parse boundary conditions; compound conditions of the form `a:s' allow
!  to have different variables at the lower and upper boundaries.
!
      call parse_bc(bcx,bcx1,bcx2)
      call parse_bc(bcy,bcy1,bcy2)
      call parse_bc(bcz,bcz1,bcz2)
!
      if (lroot.and.ip<14) then
        print*, 'bcx1,bcx2= ', bcx1," : ",bcx2
        print*, 'bcy1,bcy2= ', bcy1," : ",bcy2
        print*, 'bcz1,bcz2= ', bcz1," : ",bcz2
        print*, 'lperi= ', lperi
      endif
!
      call check_consistency_of_lperi('read_startpars')
!
    endsubroutine read_startpars
!***********************************************************************
    subroutine sample_startpars(label,iostat)
!
      character (len=*), optional :: label
      integer, optional :: iostat
!
      if (lroot) then
        print*
        print*,'-----BEGIN sample namelist ------'
                                print*,'&init_pars                 /'
        if (linitial_condition) print*,'&initial_condition_pars    /'
        if (leos              ) print*,'&eos_init_pars             /'
        if (lhydro            ) print*,'&hydro_init_pars           /'
        if (ldensity          ) print*,'&density_init_pars         /'
        if (lgrav             ) print*,'&grav_init_pars            /'
        if (lselfgravity      ) print*,'&selfgrav_init_pars        /'
        if (lpoisson          ) print*,'&poisson_init_pars         /'
        if (lentropy          ) print*,'&entropy_init_pars         /'
        if (ltemperature      ) print*,'&entropy_init_pars         /'
        if (lmagnetic         ) print*,'&magnetic_init_pars        /'
        if (lmagn_mf          ) print*,'&magn_mf_init_pars         /'
        if (lmagn_mf_demfdt   ) print*,'&magn_mf_demfdt_init_pars  /'
        if (llorenz_gauge     ) print*,'&lorenz_gauge_init_pars    /'
        if (ltestscalar       ) print*,'&testscalar_init_pars      /'
        if (ltestfield        ) print*,'&testfield_init_pars       /'
        if (ltestflow         ) print*,'&testflow_init_pars        /'
        if (lradiation        ) print*,'&radiation_init_pars       /'
        if (lpscalar          ) print*,'&pscalar_init_pars         /'
        if (lchiral           ) print*,'&chiral_init_pars          /'
        if (lchemistry        ) print*,'&chemistry_init_pars       /'
        if (lsignal           ) print*,'&signal_init_pars          /'
        if (ldustvelocity     ) print*,'&dustvelocity_init_pars    /'
        if (ldustdensity      ) print*,'&dustdensity_init_pars     /'
        if (lneutralvelocity  ) print*,'&neutralvelocity_init_pars /'
        if (lneutraldensity   ) print*,'&neutraldensity_init_pars  /'
        if (lcosmicray        ) print*,'&cosmicray_init_pars       /'
        if (lcosmicrayflux    ) print*,'&cosmicrayflux_init_pars   /'
        if (linterstellar     ) print*,'&interstellar_init_pars    /'
        if (lshear            ) print*,'&shear_init_pars           /'
        if (ltestperturb      ) print*,'&testperturb_init_pars     /'
        if (lspecial          ) print*,'&special_init_pars         /'
        if (lsolid_cells      ) print*,'&solid_cells_init_pars     /'
        if (lnscbc            ) print*,'&NSCBC_init_pars           /'
        if (lpolymer          ) print*,'&polymer_init_pars  /'
        if (lparticles        ) print*,'&particles_init_pars_wrap  /'
        print*,'------END sample namelist -------'
        print*
        if (present(label))  print*, 'Found error in input namelist "' // trim(label) // '"'
        if (present(iostat)) print*, 'iostat = ', iostat
        if (present(iostat).or.present(label)) &
                           print*,  '-- use sample above.'
      endif
!
      call fatal_error('','')
!
    endsubroutine sample_startpars
!***********************************************************************
    subroutine print_startpars(file)
!
!  Print input parameters.
!
!  4-oct02/wolf: adapted
!
      character (len=*), optional :: file
      character (len=datelen) :: date
      integer :: unit=6         ! default unit is 6=stdout
!
      if (lroot) then
        if (present(file)) then
          unit = 1
          call date_time_string(date)
          open(unit,FILE=file)
          write(unit,*) &
               '! -------------------------------------------------------------'
          write(unit,'(A,A)') ' ! ', 'Initializing'
          write(unit,'(A,A)') ' ! Date: ', trim(date)
          write(unit,*) '! t=', t
        endif
!
        write(unit,NML=init_pars          )
!
        call write_initial_condition_pars(unit)
        call write_eos_init_pars(unit)
        call write_hydro_init_pars(unit)
        call write_density_init_pars(unit)
        call write_forcing_init_pars(unit)
        call write_gravity_init_pars(unit)
        call write_selfgravity_init_pars(unit)
        call write_poisson_init_pars(unit)
        call write_entropy_init_pars(unit)
        call write_magnetic_init_pars(unit)
        call write_lorenz_gauge_init_pars(unit)
        call write_testscalar_init_pars(unit)
        call write_testfield_init_pars(unit)
        call write_testflow_init_pars(unit)
        call write_radiation_init_pars(unit)
        call write_pscalar_init_pars(unit)
        call write_chiral_init_pars(unit)
        call write_chemistry_init_pars(unit)
        call write_signal_init_pars(unit)
        call write_dustvelocity_init_pars(unit)
        call write_dustdensity_init_pars(unit)
        call write_neutralvelocity_init_pars(unit)
        call write_neutraldensity_init_pars(unit)
        call write_cosmicray_init_pars(unit)
        call write_cosmicrayflux_init_pars(unit)
        call write_interstellar_init_pars(unit)
        call write_shear_init_pars(unit)
        call write_testperturb_init_pars(unit)
        call write_viscosity_init_pars(unit)
        call write_special_init_pars(unit)
        call write_shock_init_pars(unit)
        call write_solid_cells_init_pars(unit)
        call write_NSCBC_init_pars(unit)
        call write_polymer_init_pars(unit)
        call particles_wparam(unit)
!
        if (present(file)) then
          close(unit)
        endif
      endif
!
    endsubroutine print_startpars
!***********************************************************************
    subroutine read_runpars(file,annotation)
!
!  Read input parameters.
!
!  14-sep-01/axel: inserted from run.f90
!  31-may-02/wolf: renamed from cread to read_runpars
!   6-jul-02/axel: in case of error, print sample namelist
!  21-oct-03/tony: moved sample namelist stuff to a separate procedure
!  12-nov-10/MR: added read and write calls for namelist power_spectrum_run_pars
!
      use Dustvelocity, only: copy_bcs_dust
      use Mpicomm, only: parallel_open, parallel_close
      use Sub, only: parse_bc
!
      logical, optional :: file
      character (len=*), optional :: annotation
!
      integer :: ierr, unit=1
!
!  Open namelist file.
!
      call parallel_open(unit,'run.in','formatted')
!
!  Read through all items that *may* be present in the various modules.
!  AB: at some point the sgi_fix stuff should probably be removed (see sgi bug)
!
      read(unit,NML=run_pars,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('run_pars',ierr)
      rewind(unit)
!
      call read_eos_run_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('eos_run_pars',ierr)
      rewind(unit)
!
      call read_hydro_run_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('hydro_run_pars',ierr)
      rewind(unit)
!
      call read_density_run_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('density_run_pars',ierr)
      rewind(unit)
!
      call read_forcing_run_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('forcing_run_pars',ierr)
      rewind(unit)
!
      call read_gravity_run_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('grav_run_pars',ierr)
      rewind(unit)
!
      call read_selfgravity_run_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('selfgrav_run_pars',ierr)
      rewind(unit)
!
      call read_poisson_run_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('poisson_run_pars',ierr)
      rewind(unit)
!
      call read_entropy_run_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('entropy_run_pars',ierr)
      rewind(unit)
!
      call read_magnetic_run_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('magnetic_run_pars',ierr)
      rewind(unit)
!
      call read_lorenz_gauge_run_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('lorenz_gauge_run_pars',ierr)
      rewind(unit)
!
      call read_testscalar_run_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('testscalar_run_pars',ierr)
      rewind(unit)
!
      call read_testfield_run_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('testfield_run_pars',ierr)
      rewind(unit)
!
      call read_testflow_run_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('testflow_run_pars',ierr)
      rewind(unit)
!
      call read_radiation_run_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('radiation_run_pars',ierr)
      rewind(unit)
!
      call read_pscalar_run_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('pscalar_run_pars',ierr)
      rewind(unit)
!
      call read_chiral_run_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('chiral_run_pars',ierr)
      rewind(unit)
!
      call read_chemistry_run_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('chemistry_run_pars',ierr)
      rewind(unit)
!
      call read_dustvelocity_run_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('dustvelocity_run_pars',ierr)
      rewind(unit)
!
      call read_dustdensity_run_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('dustdensity_run_pars',ierr)
      rewind(unit)
!
      call read_neutralvelocity_run_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('neutralvelocity_run_pars',ierr)
      rewind(unit)
!
      call read_neutraldensity_run_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('neutraldensity_run_pars',ierr)
      rewind(unit)
!
      call read_cosmicray_run_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('cosmicray_run_pars',ierr)
      rewind(unit)
!
      call read_cosmicrayflux_run_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('cosmicrayflux_run_pars',ierr)
      rewind(unit)
!
      call read_interstellar_run_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('interstellar_run_pars',ierr)
      rewind(unit)
!
      call read_shear_run_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('shear_run_pars',ierr)
      rewind(unit)
!
      call read_testperturb_run_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('testperturb_run_pars',ierr)
      rewind(unit)
!
      call read_viscosity_run_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('viscosity_run_pars',ierr)
      rewind(unit)
!
      call read_special_run_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('special_run_pars',ierr)
      rewind(unit)
!
      call read_shock_run_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('shock_run_pars',ierr)
      rewind(unit)
!
      call read_solid_cells_run_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('solid_cells_run_pars',ierr)
      rewind(unit)
!
      call read_NSCBC_run_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('NSCBC_run_pars',ierr)
      rewind(unit)
!
      call read_polymer_run_pars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('polymer_run_pars',ierr)
      rewind(unit)
!
      call particles_read_runpars(unit,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('particles_run_pars_wrap',ierr)
      rewind(unit)
!
      call read_power_spectrum_runpars(unit,IOSTAT=ierr)
      !if (ierr/=0) call sample_runpars('power_spectrum_run_pars_wrap',ierr)
      rewind(unit)
!
      call parallel_close(unit)
!
!  Print SVN id from first line.
!
      if (lroot) call svn_id(cvsid)
!
!  Set debug logical (easier to use than the combination of ip and lroot).
!
      ldebug=lroot.and.(ip<7)
!
!  Write parameters to log file.
!  [No longer used, since at the time when read_runpars() is called, t
!  is not known yet]
!
      if (present(file)) then
        if (file) then
          if (present(annotation)) then
            call print_runpars(FILE=trim(datadir)//'/params.log', &
                               ANNOTATION=annotation)
          else
            call print_runpars(FILE=trim(datadir)//'/params.log')
          endif
        endif
      endif
!
!  Parse boundary conditions; compound conditions of the form `a:s' allow
!  to have different variables at the lower and upper boundaries.
!
      call parse_bc(bcx,bcx1,bcx2)
      call parse_bc(bcy,bcy1,bcy2)
      call parse_bc(bcz,bcz1,bcz2)
!
      if (lroot.and.ip<14) then
        print*, 'bcx1,bcx2= ', bcx1," : ",bcx2
        print*, 'bcy1,bcy2= ', bcy1," : ",bcy2
        print*, 'bcz1,bcz2= ', bcz1," : ",bcz2
      endif
!
      call check_consistency_of_lperi('read_runpars')
!
    endsubroutine read_runpars
!***********************************************************************
    subroutine sample_runpars(label,iostat)
!
      character (len=*), optional :: label
      integer, optional :: iostat
!
      if (lroot) then
        print*
        print*,'-----BEGIN sample namelist ------'
                              print*,'&run_pars                 /'
        if (leos            ) print*,'&eos_run_pars             /'
        if (lhydro          ) print*,'&hydro_run_pars           /'
        if (lhydro_kinematic) print*,'&hydro_run_pars           /'
        if (ldensity        ) print*,'&density_run_pars         /'
        if (lforcing        ) print*,'&forcing_run_pars         /'
        if (lgrav           ) print*,'&grav_run_pars            /'
        if (lselfgravity    ) print*,'&selfgrav_run_pars        /'
        if (lpoisson        ) print*,'&poisson_run_pars         /'
        if (lentropy        ) print*,'&entropy_run_pars         /'
        if (ltemperature    ) print*,'&entropy_run_pars         /'
        if (lmagnetic       ) print*,'&magnetic_run_pars        /'
        if (lmagn_mf        ) print*,'&magn_mf_run_pars         /'
        if (llorenz_gauge   ) print*,'&lorenz_gauge_run_pars    /'
        if (ltestscalar     ) print*,'&testscalar_run_pars      /'
        if (ltestfield      ) print*,'&testfield_run_pars       /'
        if (ltestflow       ) print*,'&testflow_run_pars        /'
        if (lradiation      ) print*,'&radiation_run_pars       /'
        if (lpscalar        ) print*,'&pscalar_run_pars         /'
        if (lchiral         ) print*,'&chiral_run_pars          /'
        if (lchemistry      ) print*,'&chemistry_run_pars       /'
        if (ldustvelocity   ) print*,'&dustvelocity_run_pars    /'
        if (ldustdensity    ) print*,'&dustdensity_run_pars     /'
        if (lneutralvelocity) print*,'&neutralvelocity_run_pars /'
        if (lneutraldensity ) print*,'&neutraldensity_run_pars  /'
        if (lcosmicray      ) print*,'&cosmicray_run_pars       /'
        if (lcosmicrayflux  ) print*,'&cosmicrayflux_run_pars   /'
        if (linterstellar   ) print*,'&interstellar_run_pars    /'
        if (lshear          ) print*,'&shear_run_pars           /'
        if (ltestperturb    ) print*,'&testperturb_run_pars     /'
        if (lviscosity      ) print*,'&viscosity_run_pars       /'
        if (lspecial        ) print*,'&special_run_pars         /'
        if (lshock          ) print*,'&shock_run_pars           /'
        if (lsolid_cells    ) print*,'&solid_cells_run_pars     /'
        if (lnscbc          ) print*,'&NSCBC_run_pars           /'
        if (lpolymer        ) print*,'&polymer_run_pars         /'
        if (lparticles      ) print*,'&particles_run_pars_wrap  /'
        print*,'------END sample namelist -------'
        print*
        if (present(label))  print*, 'Found error in input namelist "' // trim(label) // '"'
        if (present(iostat)) print*, 'iostat = ', iostat
        if (present(iostat).or.present(label)) &
                           print*,  '-- use sample above.'
      endif
!
      call fatal_error('sample_runpars','')
!
    endsubroutine sample_runpars
!***********************************************************************
    subroutine print_runpars(file,annotation)
!
!  Print input parameters.
!
!  14-sep-01/axel: inserted from run.f90
!  31-may-02/wolf: renamed from cprint to print_runpars
!   4-oct-02/wolf: added log file stuff
!
      character (len=*), optional :: file,annotation
      integer :: unit=6         ! default unit is 6=stdout
      character (len=linelen) :: line
      character (len=datelen) :: date
!
      if (lroot) then
        ! get first line from file RELOAD
        line = read_line_from_file('RELOAD')
        if ((line == '') .and. present(annotation)) then
          line = trim(annotation)
        endif
        if (present(file)) then
          unit = 1
          call date_time_string(date)
          open(unit,FILE=file,position='append')
          write(unit,*) &
              '! -------------------------------------------------------------'
!
!  Add comment from `RELOAD' and time.
!
          write(unit,'(A,A)') ' ! ', trim(line)
          write(unit,'(A,A)') ' ! Date: ', trim(date)
          write(unit,*) '! t=', t
        endif
!
        write(unit,NML=run_pars)
!
        call write_eos_run_pars(unit)
        call write_hydro_run_pars(unit)
        call write_density_run_pars(unit)
        call write_forcing_run_pars(unit)
        call write_gravity_run_pars(unit)
        call write_selfgravity_run_pars(unit)
        call write_poisson_run_pars(unit)
        call write_entropy_run_pars(unit)
        call write_magnetic_run_pars(unit)
        call write_lorenz_gauge_run_pars(unit)
        call write_testscalar_run_pars(unit)
        call write_testfield_run_pars(unit)
        call write_testflow_run_pars(unit)
        call write_radiation_run_pars(unit)
        call write_pscalar_run_pars(unit)
        call write_chiral_run_pars(unit)
        call write_chemistry_run_pars(unit)
        call write_dustvelocity_run_pars(unit)
        call write_dustdensity_run_pars(unit)
        call write_neutralvelocity_run_pars(unit)
        call write_neutraldensity_run_pars(unit)
        call write_cosmicray_run_pars(unit)
        call write_cosmicrayflux_run_pars(unit)
        call write_interstellar_run_pars(unit)
        call write_shear_run_pars(unit)
        call write_testperturb_run_pars(unit)
        call write_viscosity_run_pars(unit)
        call write_special_run_pars(unit)
        call write_shock_run_pars(unit)
        call write_solid_cells_run_pars(unit)
        call write_NSCBC_run_pars(unit)
        call write_power_spectrum_runpars(unit)
        call write_polymer_run_pars(unit)
        call particles_wparam2(unit)
!
        if (present(file)) close(unit)
!
      endif
!
    endsubroutine print_runpars
!***********************************************************************
    subroutine check_consistency_of_lperi(label)
!
!  Check consistency of lperi.
!
!  18-jul-03/axel: coded
!
      character (len=*) :: label
      logical :: lwarning=.true.
      integer :: j
!
!  Identifier.
!
      if (lroot.and.ip<5) print*,'check_consistency_of_lperi: called from ',label
!
!  Make the warnings less dramatic looking, if we are only in start
!  and exit this routine altogether if, in addition, ip > 13.
!
      if (label=='read_startpars'.and.ip>13) return
      if (label=='read_startpars') lwarning=.false.
!
      if (nvar > 0) then
!
!  Check x direction.
!
        j=1
        if (any(bcx(1:nvar)=='p'.or. bcx(1:nvar)=='she').and..not.lperi(j).or.&
            any(bcx(1:nvar)/='p'.and.bcx(1:nvar)/='she').and.lperi(j)) &
            call warning_lperi(lwarning,bcx(1:nvar),lperi,j)
!
!  Check y direction.
!
        j=2
        if (any(bcy(1:nvar)=='p').and..not.lperi(j).or.&
            any(bcy(1:nvar)/='p').and.lperi(j)) &
            call warning_lperi(lwarning,bcy(1:nvar),lperi,j)
!
!  Check z direction.
!
        j=3
        if (any(bcz(1:nvar)=='p').and..not.lperi(j).or.&
            any(bcz(1:nvar)/='p').and.lperi(j)) &
            call warning_lperi(lwarning,bcz(1:nvar),lperi,j)
      endif
!
!  Print final warning.
!  Make the warnings less dramatic looking, if we are only in start.
!
      if (lroot .and. (.not. lwarning)) then
        if (label=='read_startpars') then
          print*,'[bad BCs in start.in only affects post-processing' &
               //' of start data, not the run]'
        else
          print*,'check_consistency_of_lperi(run.in): you better stop and check!'
          print*,'------------------------------------------------------'
          print*
        endif
      endif
!
    endsubroutine check_consistency_of_lperi
!***********************************************************************
    subroutine warning_lperi(lwarning,bc,lperi,j)
!
!  Print consistency warning of lperi.
!
!  18-jul-03/axel: coded
!
      character (len=*), dimension(:), intent(in) :: bc
      logical, dimension(3) :: lperi
      logical :: lwarning
      integer :: j
!
      if (lroot) then
        if (lwarning) then
          print*
          print*,'------------------------------------------------------'
          print*,'W A R N I N G'
          lwarning=.false.
        else
          print*
        endif
!
        print*,'warning_lperi: inconsistency, j=', j, ', lperi(j)=',lperi(j)
        print*,'bc=',bc
        print*,"any(bc=='p'.or. bc=='she'), .not.lperi(j) = ", &
          any(bc=='p'.or. bc=='she'), .not.lperi(j)
        print*, "any(bcx/='p'.and.bcx/='she'), lperi(j) = ", &
          any(bc=='p'.or. bc=='she'), .not.lperi(j)
      endif
!
    endsubroutine warning_lperi
!***********************************************************************
    subroutine wparam()
!
!  Write startup parameters
!
!  21-jan-02/wolf: coded
!
      logical :: lhydro           = lhydro_var
      logical :: ldensity         = ldensity_var
      logical :: lentropy         = lentropy_var
      logical :: ltemperature     = ltemperature_var
      logical :: lshock           = lshock_var
      logical :: lmagnetic        = lmagnetic_var
      logical :: lforcing         = lforcing_var
      logical :: llorenz_gauge    = llorenz_gauge_var
      logical :: ldustvelocity    = ldustvelocity_var
      logical :: ldustdensity     = ldustdensity_var
      logical :: ltestscalar      = ltestscalar_var
      logical :: ltestfield       = ltestfield_var
      logical :: ltestflow        = ltestflow_var
      logical :: linterstellar    = linterstellar_var
      logical :: lcosmicray       = lcosmicray_var
      logical :: lcosmicrayflux   = lcosmicrayflux_var
      logical :: lshear           = lshear_var
      logical :: lpscalar         = lpscalar_var
      logical :: lchiral          = lchiral_var
      logical :: leos             = leos_var
      logical :: lradiation       = lradiation_var
      logical :: lneutralvelocity = lneutralvelocity_var
      logical :: lneutraldensity  = lneutraldensity_var
      logical :: lpolymer         = lpolymer_var
      integer :: unit=1
!
      namelist /lphysics/ &
          lhydro, ldensity, lentropy, lmagnetic, lshear, llorenz_gauge,  &
          ltestscalar, ltestfield, ltestflow, lpscalar, lradiation, &
          ldustvelocity, ldustdensity, lforcing, lgravz, lgravr, &
          ltestperturb, linterstellar, lcosmicray, lcosmicrayflux, &
          lshock, lradiation_fld, leos_ionization, leos_fixed_ionization, &
          lvisc_hyper, lchiral, leos, leos_temperature_ionization, &
          lneutralvelocity, lneutraldensity, ltemperature,lpolymer
!
!  Write the param.nml file only from root processor.
!  However, for pacx-MPI (grid-style computations across different platforms)
!  we'd need this on each site separately (not done yet).
!  (In that case we'd need to identify separate master-like processors
!  one at each site.)
!
      if (lroot) then
        open(unit,FILE=trim(datadir)//'/param.nml',DELIM='apostrophe', &
               STATUS='unknown')
!
!  Write init_pars.
!
        write(unit,NML=init_pars)
!
!  Write each namelist separately.
!
        call write_eos_init_pars(unit)
        call write_hydro_init_pars(unit)
        call write_density_init_pars(unit)
        call write_forcing_init_pars(unit)
        call write_gravity_init_pars(unit)
        call write_selfgravity_init_pars(unit)
        call write_poisson_init_pars(unit)
        call write_entropy_init_pars(unit)
        call write_magnetic_init_pars(unit)
        call write_lorenz_gauge_init_pars(unit)
        call write_testscalar_init_pars(unit)
        call write_testfield_init_pars(unit)
        call write_testflow_init_pars(unit)
        call write_radiation_init_pars(unit)
        call write_pscalar_init_pars(unit)
        call write_chiral_init_pars(unit)
        call write_chemistry_init_pars(unit)
        call write_signal_init_pars(unit)
        call write_dustvelocity_init_pars(unit)
        call write_dustdensity_init_pars(unit)
        call write_neutralvelocity_init_pars(unit)
        call write_neutraldensity_init_pars(unit)
        call write_cosmicray_init_pars(unit)
        call write_cosmicrayflux_init_pars(unit)
        call write_interstellar_init_pars(unit)
        call write_shear_init_pars(unit)
        call write_testperturb_init_pars(unit)
        call write_viscosity_init_pars(unit)
        call write_special_init_pars(unit)
        call write_shock_init_pars(unit)
        call write_solid_cells_init_pars(unit)
        call write_NSCBC_init_pars(unit)
        call write_polymer_init_pars(unit)
        call write_initial_condition_pars(unit)
        call particles_wparam(unit)
        ! The following parameters need to be communicated to IDL
        write(unit,NML=lphysics)
        close(unit)
      endif
!
      call keep_compiler_quiet(lhydro)
      call keep_compiler_quiet(ldensity)
      call keep_compiler_quiet(lentropy)
      call keep_compiler_quiet(ltemperature)
      call keep_compiler_quiet(lmagnetic)
      call keep_compiler_quiet(lforcing)
      call keep_compiler_quiet(ltestscalar,ltestfield,ltestflow)
      call keep_compiler_quiet(lpscalar,lradiation,lcosmicray,lcosmicrayflux)
      call keep_compiler_quiet(linterstellar,lshock)
      call keep_compiler_quiet(ldustdensity,ldustvelocity)
      call keep_compiler_quiet(llorenz_gauge)
      call keep_compiler_quiet(lshear)
      call keep_compiler_quiet(leos)
      call keep_compiler_quiet(lchiral)
      call keep_compiler_quiet(lneutralvelocity)
      call keep_compiler_quiet(lneutraldensity)
      call keep_compiler_quiet(lpolymer)
!
    endsubroutine wparam
!***********************************************************************
    subroutine rparam()
!
!  Read startup parameters.
!
!  21-jan-02/wolf: coded
!
      use Mpicomm, only: parallel_open, parallel_close
!
      integer :: unit=1
!
      call parallel_open(unit,trim(datadir)//'/param.nml')
!
      read(unit,NML=init_pars)
      rewind(unit)
      call read_eos_init_pars(unit)
      rewind(unit)
      call read_hydro_init_pars(unit)
      rewind(unit)
      call read_density_init_pars(unit)
      rewind(unit)
      call read_forcing_init_pars(unit)
      rewind(unit)
      call read_gravity_init_pars(unit)
      rewind(unit)
      call read_selfgravity_init_pars(unit)
      rewind(unit)
      call read_poisson_init_pars(unit)
      rewind(unit)
      call read_entropy_init_pars(unit)
      rewind(unit)
      call read_magnetic_init_pars(unit)
      rewind(unit)
      call read_testscalar_init_pars(unit)
      rewind(unit)
      call read_testfield_init_pars(unit)
      rewind(unit)
      call read_testflow_init_pars(unit)
      rewind(unit)
      call read_radiation_init_pars(unit)
      rewind(unit)
      call read_pscalar_init_pars(unit)
      rewind(unit)
      call read_chiral_init_pars(unit)
      rewind(unit)
      call read_chemistry_init_pars(unit)
      rewind(unit)
      call read_signal_init_pars(unit)
      rewind(unit)
      call read_dustvelocity_init_pars(unit)
      rewind(unit)
      call read_dustdensity_init_pars(unit)
      rewind(unit)
      call read_neutralvelocity_init_pars(unit)
      rewind(unit)
      call read_neutraldensity_init_pars(unit)
      rewind(unit)
      call read_cosmicray_init_pars(unit)
      rewind(unit)
      call read_cosmicrayflux_init_pars(unit)
      rewind(unit)
      call read_interstellar_init_pars(unit)
      rewind(unit)
      call read_shear_init_pars(unit)
      rewind(unit)
      call read_testperturb_init_pars(unit)
      rewind(unit)
      call read_viscosity_init_pars(unit)
      rewind(unit)
      call read_special_init_pars(unit)
      rewind(unit)
      call read_shock_init_pars(unit)
      rewind(unit)
      call read_solid_cells_init_pars(unit)
      rewind(unit)
      call read_NSCBC_init_pars(unit)
      rewind(unit)
      call read_polymer_init_pars(unit)
      rewind(unit)
      call read_initial_condition_pars(unit)
      rewind(unit)
      call particles_rparam(unit)
      rewind(unit)
!
      call parallel_close(unit)
!
    endsubroutine rparam
!***********************************************************************
    subroutine wparam2()
!
!  Write runtime parameters for IDL.
!
!  21-jan-02/wolf: coded
!
      integer :: unit=1
!
      if (lroot) then
        open(unit,FILE=trim(datadir)//'/param2.nml',DELIM='apostrophe')
!
        write(unit,NML=run_pars)
!
        call write_eos_run_pars(unit)
        call write_hydro_run_pars(unit)
        call write_density_run_pars(unit)
        call write_forcing_run_pars(unit)
        call write_gravity_run_pars(unit)
        call write_selfgravity_run_pars(unit)
        call write_poisson_run_pars(unit)
        call write_entropy_run_pars(unit)
        call write_magnetic_run_pars(unit)
        call write_lorenz_gauge_run_pars(unit)
        call write_testscalar_run_pars(unit)
        call write_testfield_run_pars(unit)
        call write_testflow_run_pars(unit)
        call write_radiation_run_pars(unit)
        call write_pscalar_run_pars(unit)
        call write_chiral_run_pars(unit)
        call write_chemistry_run_pars(unit)
        call write_dustvelocity_run_pars(unit)
        call write_dustdensity_run_pars(unit)
        call write_neutralvelocity_run_pars(unit)
        call write_neutraldensity_run_pars(unit)
        call write_cosmicray_run_pars(unit)
        call write_cosmicrayflux_run_pars(unit)
        call write_interstellar_run_pars(unit)
        call write_shear_run_pars(unit)
        call write_testperturb_run_pars(unit)
        call write_viscosity_run_pars(unit)
        call write_special_run_pars(unit)
        call write_shock_run_pars(unit)
        call write_solid_cells_run_pars(unit)
        call write_NSCBC_run_pars(unit)
        call write_power_spectrum_runpars(unit)
        call write_polymer_run_pars(unit)
        call particles_wparam2(unit)
!
        close(unit)
      endif
!
    endsubroutine wparam2
!***********************************************************************
    subroutine write_pencil_info()
!
!  Write information about requested and diagnostic pencils.
!  Do this only when on root processor.
!
     integer :: i
     integer :: unit=1
!
     if (lroot) then
       open(unit,FILE=trim(datadir)//'/pencils.list')
       write(unit,*) 'Pencils requested:'
       do i=1,npencils
         if (lpenc_requested(i)) write(unit,*) i, pencil_names(i)
       enddo
       write(unit,*) ''
       write(unit,*) 'Pencils requested for diagnostics:'
       do i=1,npencils
         if (lpenc_diagnos(i)) write(unit,*) i, pencil_names(i)
       enddo
       write(unit,*) ''
       write(unit,*) 'Pencils requested for 2-D diagnostics:'
       do i=1,npencils
         if (lpenc_diagnos2d(i)) write(unit,*) i, pencil_names(i)
       enddo
       write(unit,*) ''
       write(unit,*) 'Pencils requested video output'
       do i=1,npencils
         if (lpenc_video(i)) write(unit,*) i, pencil_names(i)
       enddo
       print*, 'write_pencil_info: pencil information written to the file pencils.list'
       close(unit)
     endif
!
   endsubroutine write_pencil_info
!***********************************************************************
endmodule Param_IO
