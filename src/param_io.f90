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
!  use Conductivity
  use Density
  use Detonate
  use Dustdensity
  use Dustvelocity
  use Energy
  use EquationOfState
  use Forcing
  use General
  use Gravity
  use Heatflux
  use Hydro
  use ImplicitDiffusion
  use InitialCondition
  use Interstellar
  use Lorenz_gauge
  use Magnetic
  use Messages
  use NeutralDensity
  use NeutralVelocity
  use NSCBC
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
  use Streamlines
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
  public :: read_all_init_pars, read_all_run_pars
  public :: write_all_init_pars, write_all_run_pars
  public :: write_pencil_info
!
  logical :: lforce_shear_bc = .true., ltolerate_namelist_errors=.false.
!
! local quantities
!
  real, dimension(mcom) :: fbcx1=0., fbcx2=0., fbcx1_2=0., fbcx2_2=0., &
                           fbcy1=0., fbcy2=0., fbcy1_1=0., fbcy1_2=0., fbcy2_1=0., fbcy2_2=0., &
                           fbcz1=0., fbcz2=0., fbcz1_1=0., fbcz1_2=0., fbcz2_1=0., fbcz2_2=0.

  namelist /init_pars/ &
      cvsid, ip, xyz0, xyz1, Lxyz, lperi, lshift_origin, lshift_origin_lower, coord_system, &
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
      tstart, lseparate_persist, ldistribute_persist, &
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
      lforce_shear_bc,lread_from_other_prec, &
      pipe_func, glnCrossSec0, CrossSec_x1, CrossSec_x2, CrossSec_w,&
      lcorotational_frame, rcorot, lproper_averages, &
      ldirect_access, ltolerate_namelist_errors
!
  namelist /run_pars/ &
      cvsid, ip, xyz0, xyz1, Lxyz, lperi, lshift_origin, lshift_origin_lower, coord_system, &
      nt, it1, it1d, dt, cdt, ddt, cdtv, cdtv2, cdtv3, cdts, cdtr, &
      cdtc, isave, itorder, dsnap, dsnap_down, d2davg, dvid, dsound, dtmin, dspec, tmax, iwig, &
      dtracers, dfixed_points, unit_system, unit_length, &
      unit_velocity, unit_density, unit_temperature, unit_magnetic, &
      awig, ialive, max_walltime, dtmax, ldt_paronly, vel_spec, mag_spec, &
      uxy_spec, bxy_spec, jxbxy_spec, xy_spec, oo_spec, &
      uxj_spec, vec_spec, ou_spec, ab_spec, azbz_spec, ub_spec, Lor_spec, &
      vel_phispec, mag_phispec, &
      uxj_phispec, vec_phispec, ou_phispec, ab_phispec, EP_spec, ro_spec, &
      TT_spec, ss_spec, cc_spec, cr_spec, sp_spec, isaveglobal, lr_spec, r2u_spec, &
      r3u_spec, rhocc_pdf, cc_pdf, lncc_pdf, gcc_pdf, lngcc_pdf, kinflow, &
      lkinflow_as_aux, ampl_kinflow_x, ampl_kinflow_y, ampl_kinflow_z, &
      kx_kinflow, ky_kinflow, kz_kinflow, dtphase_kinflow, &
      random_gen, der2_type, lrmwig_rho, lrmwig_full, lrmwig_xyaverage, &
      ltime_integrals, lnowrite, noghost_for_isave, nghost_read_fewer, &
      lwrite_yaverages, lwrite_zaverages, lwrite_phiaverages, &
      test_nonblocking, lwrite_tracers, lwrite_fixed_points, &
      lread_oldsnap_lnrho2rho, lread_oldsnap_nomag, lread_oldsnap_nopscalar, &
      lread_oldsnap_notestfield, lread_oldsnap_notestscalar, &
      lwrite_dim_again,&
      lread_aux, comment_char, ix, iy, iz, iz2, iz3, iz4, slice_position, &
      xbot_slice, xtop_slice, ybot_slice, ytop_slice, zbot_slice, ztop_slice, &
      bcx, bcy, bcz, r_int, r_ext, r_int_border, &
      r_ext_border, lfreeze_varsquare, lfreeze_varint, lfreeze_varext, &
      xfreeze_square, yfreeze_square, rfreeze_int, rfreeze_ext, wfreeze, &
      wfreeze_int, wfreeze_ext, wborder, wborder_int, wborder_ext, tborder, &
      luse_oldgrid, luse_xyz1, fshift_int, fshift_ext, &
      lreset_tstart, tstart, lseparate_persist, ldistribute_persist, &
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
      lpencil_check_small, lpencil_check_no_zeros, lpencil_check_diagnos_opti, &
      lpencil_init, penc0, lwrite_2d, lbidiagonal_derij, lisotropic_advection, &
      crash_file_dtmin_factor, niter_poisson, ltestperturb, eps_rkf, &
      eps_stiff, timestep_scaling, lequatory, lequatorz, zequator, &
      lini_t_eq_zero, lav_smallx, xav_max, ldt_paronly, lweno_transport, &
      it_timing, har_spec, hav_spec, j_spec, jb_spec, lread_less, lformat, ltec, &
      llsode, lsplit_second, nu_sts, permute_sts, lfargo_advection, &
      ldynamical_diffusion, ldyndiff_urmsmxy, re_mesh, lreset_seed, &
      loutput_varn_at_exact_tsnap, lstop_on_ioerror, mailaddress, &
      theta_lower_border, wborder_theta_lower, theta_upper_border, &
      wborder_theta_upper, fraction_tborder, lmeridional_border_drive, &
      lread_from_other_prec, downsampl, lfullvar_in_slices, lsubstract_reference_state, &
      ldirect_access, lproper_averages, &
      pipe_func, glnCrossSec0, CrossSec_x1, CrossSec_x2, CrossSec_w
!
  contains
!***********************************************************************
    subroutine get_datadir(dir)
!
!  Overwrite datadir from datadir.in, if that exists.
!
!   2-oct-02/wolf: coded
!  25-oct-02/axel: default is taken from cdata.f90 where it's defined
!  14-jan-15/MR  : corrected call of mpibcast_char
!
      use Mpicomm, only: mpibcast_logical, mpibcast_char
!
      character (len=*) :: dir
      logical :: exists
      integer, parameter :: unit=1
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
      call mpibcast_logical(exists)
!
!  Let root processor communicate dir (i.e. datadir) to all other processors.
!
      if (exists) call mpibcast_char(dir)
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
    subroutine read_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=init_pars, IOSTAT=iostat)
!
    endsubroutine read_init_pars
!***********************************************************************
    subroutine read_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=run_pars, IOSTAT=iostat)
!
    endsubroutine read_run_pars
!***********************************************************************
    subroutine read_all_init_pars(print)
!
!  read input parameters (done by each processor)
!
!  Now no warning is produced when somebody sets lcylindrical.
!  "lcylindrical=T" is now outdated, use instead: "lcylinder_in_a_box=T"
!  This renaming became necessary with the development of
!  cylindrical coordinates which led to very similar names
!  (coord_system="cylindrical_coords")
!
!   6-jul-02/axel: in case of error, print sample namelist
!  21-oct-03/tony: moved sample namelist stuff to a separate procedure
!  18-aug-15/PABourdin: reworked to simplify code and display all errors at once
!  19-aug-15/PABourdin: renamed from read_startpars to read_all_init_pars
!
      use File_io, only: parallel_open, parallel_close
      use Particles_main, only: read_all_particles_init_pars
      use Sub, only: read_namelist
!
      logical, optional, intent(IN) :: print
!
      character(LEN=fnlen) :: file

      lnamelist_error = .false.
!
!  Set default to shearing sheet if lshear=.true. (even when Sshear==0.).
!
      if (lshear .and. lforce_shear_bc) bcx(:)='she'
!
!  Open namelist file and read through all items that *may* be present in the various modules.
!
      if (lstart) then
        file='start.in'
        call parallel_open(file,remove_comments=.true.)
      else
        file=trim(datadir)//'/param.nml'
        call parallel_open(file)
        lstart=.true.                      ! necessary to create correct error messages in read_namelist
      endif
!
      call read_namelist(read_init_pars                ,'')
      call read_namelist(read_initial_condition_pars   ,'initial_condition_pars',linitial_condition)
      call read_namelist(read_streamlines_init_pars    ,'streamlines')
      call read_namelist(read_eos_init_pars            ,'eos'            ,leos)
      call read_namelist(read_hydro_init_pars          ,'hydro'          ,lhydro.or.lhydro_kinematic)
      call read_namelist(read_density_init_pars        ,'density'        ,ldensity)
      call read_namelist(read_gravity_init_pars        ,'grav'           ,lgrav)
      call read_namelist(read_selfgravity_init_pars    ,'selfgrav'       ,lselfgravity)
      call read_namelist(read_poisson_init_pars        ,'poisson'        ,lpoisson)
      call read_namelist(read_energy_init_pars         ,'entropy'        ,lenergy)
      call read_namelist(read_magnetic_init_pars       ,'magnetic'       ,lmagnetic)
      call read_namelist(read_lorenz_gauge_init_pars   ,'lorenz_gauge'   ,llorenz_gauge)
      call read_namelist(read_testscalar_init_pars     ,'testscalar'     ,ltestscalar)
      call read_namelist(read_testfield_init_pars      ,'testfield'      ,ltestfield)
      call read_namelist(read_testflow_init_pars       ,'testflow'       ,ltestflow)
      call read_namelist(read_radiation_init_pars      ,'radiation'      ,lradiation)
      call read_namelist(read_pscalar_init_pars        ,'pscalar'        ,lpscalar)
      call read_namelist(read_chiral_init_pars         ,'chiral'         ,lchiral)
      call read_namelist(read_chemistry_init_pars      ,'chemistry'      ,lchemistry)
      call read_namelist(read_signal_init_pars         ,'signal'         ,lsignal)
      call read_namelist(read_dustvelocity_init_pars   ,'dustvelocity'   ,ldustvelocity)
      call read_namelist(read_dustdensity_init_pars    ,'dustdensity'    ,ldustdensity)
      call read_namelist(read_neutralvelocity_init_pars,'neutralvelocity',lneutralvelocity)
      call read_namelist(read_neutraldensity_init_pars ,'neutraldensity' ,lneutraldensity)
      call read_namelist(read_cosmicray_init_pars      ,'cosmicray'      ,lcosmicray)
      call read_namelist(read_cosmicrayflux_init_pars  ,'cosmicrayflux'  ,lcosmicrayflux)
      call read_namelist(read_interstellar_init_pars   ,'interstellar'   ,linterstellar)
      call read_namelist(read_shear_init_pars          ,'shear'          ,lshear)
      call read_namelist(read_special_init_pars        ,'special'        ,lspecial)
      call read_namelist(read_solid_cells_init_pars    ,'solid_cells'    ,lsolid_cells)
      call read_namelist(read_NSCBC_init_pars          ,'NSCBC'          ,lnscbc)
      call read_namelist(read_polymer_init_pars        ,'polymer'        ,lpolymer)
!
      call read_all_particles_init_pars
!
      call parallel_close
!
      if (lnamelist_error .and. .not.ltolerate_namelist_errors) then
        call sample_pars
        call fatal_error ('read_all_init_pars', 'Please fix all above WARNINGs for file "'//trim(file)//'"')
      endif
!
      if (lrun) lstart=.false.

!  Print SVN id from first line.
!
      if (lroot) call svn_id(cvsid)
!
!  Give online feedback if called with the PRINT optional argument.
!
      if (loptest(print)) call write_all_init_pars
!
!  Write parameters to log file.
!
      if (lstart) call write_all_init_pars(FILE=trim(datadir)//'/params.log')
!
!  Parse boundary conditions; compound conditions of the form `a:s' allow
!  to have different variables at the lower and upper boundaries.
!
      call parse_bc(bcx,bcx12)
      call parse_bc(bcy,bcy12)
      call parse_bc(bcz,bcz12)
!
      fbcx  (:,1) = fbcx1;   fbcx  (:,2) = fbcx2
      fbcx_2(:,1) = fbcx1_2; fbcx_2(:,2) = fbcx2_2

      fbcy  (:,1) = fbcy1;   fbcy  (:,2) = fbcy2
      fbcy_1(:,1) = fbcy1_1; fbcy_1(:,2) = fbcy2_1
      fbcy_2(:,1) = fbcy1_2; fbcy_2(:,2) = fbcy2_2

      fbcz  (:,1) = fbcz1;   fbcz  (:,2) = fbcz2
      fbcz_1(:,1) = fbcz1_1; fbcz_1(:,2) = fbcz2_1
      fbcz_2(:,1) = fbcz1_2; fbcz_2(:,2) = fbcz2_2
!
      if (lroot.and.ip<14) then
        print*, 'bcx1,bcx2= ', bcx12(:,1)," : ",bcx12(:,2)
        print*, 'bcy1,bcy2= ', bcy12(:,1)," : ",bcy12(:,2)
        print*, 'bcz1,bcz2= ', bcz12(:,1)," : ",bcz12(:,2)
        print*, 'lperi= ', lperi
      endif
!
      call check_consistency_of_lperi('read_all_init_pars')
!
    endsubroutine read_all_init_pars
!***********************************************************************
    subroutine read_all_run_pars(logging)
!
!  Read input parameters.
!
!  14-sep-01/axel: inserted from run.f90
!  31-may-02/wolf: renamed from 'cread' to 'read_runpars'
!   6-jul-02/axel: in case of error, print sample namelist
!  21-oct-03/tony: moved sample namelist stuff to a separate procedure
!  18-aug-15/PABourdin: reworked to simplify code and display all errors at once
!  19-aug-15/PABourdin: renamed from 'read_runpars' to 'read_all_run_pars'
!
      use Dustvelocity, only: copy_bcs_dust
      use File_io, only: parallel_open, parallel_close
      use Particles_main, only: read_all_particles_run_pars
      use Sub, only: parse_bc, read_namelist
      use General, only: loptest
!
      logical, optional, intent(in) :: logging
!
      lnamelist_error = .false.
      tstart=impossible
!
!  Open namelist file.
!
      call parallel_open('run.in',remove_comments=.true.)
!
!  Read through all items that *may* be present in the various modules.
!  AB: at some point the sgi_fix stuff should probably be removed (see sgi bug)
!
      call read_namelist(read_run_pars                ,'')
      call read_namelist(read_streamlines_run_pars    ,'streamlines')
      call read_namelist(read_eos_run_pars            ,'eos'               ,leos)
      call read_namelist(read_hydro_run_pars          ,'hydro'             ,lhydro.or.lhydro_kinematic)
      call read_namelist(read_density_run_pars        ,'density'           ,ldensity)
      call read_namelist(read_forcing_run_pars        ,'forcing'           ,lforcing)
      call read_namelist(read_gravity_run_pars        ,'grav'              ,lgrav)
      call read_namelist(read_selfgravity_run_pars    ,'selfgrav'          ,lselfgravity)
      call read_namelist(read_poisson_run_pars        ,'poisson'           ,lpoisson)
      call read_namelist(read_energy_run_pars         ,'entropy'           ,lenergy)
!     call read_namelist(read_conductivity_run_pars   ,'conductivity')
      call read_namelist(read_detonate_run_pars       ,'detonate'          ,ldetonate)
      call read_namelist(read_magnetic_run_pars       ,'magnetic'          ,lmagnetic)
      call read_namelist(read_lorenz_gauge_run_pars   ,'lorenz_gauge'      ,llorenz_gauge)
      call read_namelist(read_testscalar_run_pars     ,'testscalar'        ,ltestscalar)
      call read_namelist(read_testfield_run_pars      ,'testfield'         ,ltestfield)
      call read_namelist(read_testflow_run_pars       ,'testflow'          ,ltestflow)
      call read_namelist(read_radiation_run_pars      ,'radiation'         ,lradiation)
      call read_namelist(read_pscalar_run_pars        ,'pscalar'           ,lpscalar)
      call read_namelist(read_chiral_run_pars         ,'chiral'            ,lchiral)
      call read_namelist(read_chemistry_run_pars      ,'chemistry'         ,lchemistry)
      call read_namelist(read_dustvelocity_run_pars   ,'dustvelocity'      ,ldustvelocity)
      call read_namelist(read_dustdensity_run_pars    ,'dustdensity'       ,ldustdensity)
      call read_namelist(read_neutralvelocity_run_pars,'neutralvelocity'   ,lneutralvelocity)
      call read_namelist(read_neutraldensity_run_pars ,'neutraldensity'    ,lneutraldensity)
      call read_namelist(read_cosmicray_run_pars      ,'cosmicray'         ,lcosmicray)
      call read_namelist(read_cosmicrayflux_run_pars  ,'cosmicrayflux'     ,lcosmicrayflux)
      call read_namelist(read_heatflux_run_pars       ,'heatflux'          ,lheatflux)
      call read_namelist(read_interstellar_run_pars   ,'interstellar'      ,linterstellar)
      call read_namelist(read_shear_run_pars          ,'shear'             ,lshear)
      call read_namelist(read_testperturb_run_pars    ,'testperturb'       ,ltestperturb)
      call read_namelist(read_viscosity_run_pars      ,'viscosity'         ,lviscosity)
      call read_namelist(read_special_run_pars        ,'special'           ,lspecial)
      call read_namelist(read_shock_run_pars          ,'shock'             ,lshock)
      call read_namelist(read_solid_cells_run_pars    ,'solid_cells'       ,lsolid_cells)
      call read_namelist(read_NSCBC_run_pars          ,'NSCBC'             ,lnscbc)
      call read_namelist(read_polymer_run_pars        ,'polymer'           ,lpolymer)
      call read_namelist(read_power_spectrum_run_pars ,'power_spectrum'    ,lpower_spectrum)
      call read_namelist(read_implicit_diff_run_pars  ,'implicit_diffusion',limplicit_diffusion)
!
      call read_all_particles_run_pars
!
      call parallel_close
!
      if (lnamelist_error .and. .not.ltolerate_namelist_errors) then
        call sample_pars
        call fatal_error ('read_all_run_pars', 'Please fix all above WARNINGs for file "run.in"')
      endif
!
!  Print SVN id from first line.
!
      if (lroot) call svn_id(cvsid)
!
!  Set debug logical (easier to use than the combination of ip and lroot).
!
      ldebug = lroot .and. (ip < 7)
!
!  Write parameters to log file, if requested.
!
      if (loptest(logging)) call write_all_run_pars
!
!  Parse boundary conditions; compound conditions of the form `a:s' allow
!  to have different variables at the lower and upper boundaries.
!
      call parse_bc(bcx,bcx12)
      call parse_bc(bcy,bcy12)
      call parse_bc(bcz,bcz12)
!
      fbcx  (:,1) = fbcx1;   fbcx  (:,2) = fbcx2
      fbcx_2(:,1) = fbcx1_2; fbcx_2(:,2) = fbcx2_2

      fbcy  (:,1) = fbcy1;   fbcy  (:,2) = fbcy2
      fbcy_1(:,1) = fbcy1_1; fbcy_1(:,2) = fbcy2_1
      fbcy_2(:,1) = fbcy1_2; fbcy_2(:,2) = fbcy2_2

      fbcz  (:,1) = fbcz1;   fbcz  (:,2) = fbcz2
      fbcz_1(:,1) = fbcz1_1; fbcz_1(:,2) = fbcz2_1
      fbcz_2(:,1) = fbcz1_2; fbcz_2(:,2) = fbcz2_2
!
      if (lroot.and.ip<14) then
        print*, 'bcx1,bcx2= ', bcx12(:,1)," : ",bcx12(:,2)
        print*, 'bcy1,bcy2= ', bcy12(:,1)," : ",bcy12(:,2)
        print*, 'bcz1,bcz2= ', bcz12(:,1)," : ",bcz12(:,2)
      endif
!
      call check_consistency_of_lperi('read_all_run_pars')
!
      if (any(downsampl>1)) then
!
!  If downsampling, calculate local start indices and number of data in
!  output for each direction; inner ghost zones are here disregarded
!
        ldownsampl = .true.
        if (dsnap_down<=0.) dsnap_down=dsnap
!
        call get_downpars(1,nx,ipx)
        call get_downpars(2,ny,ipy)
        call get_downpars(3,nz,ipz)
!
      endif

    endsubroutine read_all_run_pars
!***********************************************************************
    subroutine get_downpars(ind,n,ip)
!
! Calculates start indices & lengths for downsampled output
!
! 13-feb-14/MR: coded
! 19-aug-15/PABourdin: moved, please do not use 'contains' in subroutines
!                      MR: Why not?
! Parameters: coordinate direction, number of inner grid points, processor number
!
      integer, intent(IN) :: ind, n, ip
!
      if ( downsampl(ind)>n ) then
        print*, 'get_downpars: Warning - stepsize for downsampling in '// &
        coornames(ind)//'-direction ', downsampl(ind), 'greater than grid size ', n, &
        '! Set to grid size.'
        downsampl(ind)=n
      endif
!
! first index in direction ind in local farray for output
!
      firstind(ind) = downsampl(ind) - modulo(ip*n-1,downsampl(ind))
!
! number of output items in direction ind *without* ghost zones
!
      ndown(ind) = get_range_no((/firstind(ind),n,downsampl(ind)/),1)
!
      firstind(ind) = firstind(ind) + nghost
!
    endsubroutine get_downpars
!***********************************************************************
    subroutine sample_pars
!
!  standardises handling of errors in namelist, both for init and run
!
!  31-oct-13/MR: coded
!
      character(len=labellen):: partype

      if (lroot) then

        if (lstart) then
          partype='init_pars'
        else
          partype='run_pars'
        endif

        print*
        print*,'-----BEGIN sample namelist ------'
        print*,'&'//partype//'                /'
!
        if (lstart) then
          if (lsignal)             print*,'&signal_'//partype//'         /'
          if (linitial_condition)  print*,'&initial_condition_pars   /'
        endif
!
        if (ltracers)              print*,'&streamlines_'//partype//'    /' !! questionable wg. ltracers
        if (leos)                  print*,'&eos_'//partype//'            /'
        if (lhydro .or. lhydro_kinematic) &
                                   print*,'&hydro_'//partype//'          /'
        if (ldensity)              print*,'&density_'//partype//'        /'
        if (lgrav)                 print*,'&grav_'//partype//'           /'
        if (lselfgravity)          print*,'&selfgrav_'//partype//'       /'
        if (lpoisson)              print*,'&poisson_'//partype//'        /'
        if (lentropy .or. ltemperature) &
                                   print*,'&entropy_'//partype//'        /'
        if (lthermal_energy)       print*,'&energy_'//partype//'         /'
!       if (lconductivity)         print*,'&conductivity_'//partype//'   /'
        if (ldetonate)             print*,'&detonate_'//partype//'       /'
        if (lmagnetic)             print*,'&magnetic_'//partype//'       /'
        if (lmagn_mf)              print*,'&magn_mf_'//partype//'        /'
        if (llorenz_gauge)         print*,'&lorenz_gauge_'//partype//'   /'
        if (ltestscalar)           print*,'&testscalar_'//partype//'     /'
        if (ltestfield)            print*,'&testfield_'//partype//'      /'
        if (ltestflow)             print*,'&testflow_'//partype//'       /'
        if (lradiation)            print*,'&radiation_'//partype//'      /'
        if (lpscalar)              print*,'&pscalar_'//partype//'        /'
        if (lchiral)               print*,'&chiral_'//partype//'         /'
        if (lchemistry)            print*,'&chemistry_'//partype//'      /'
        if (ldustvelocity)         print*,'&dustvelocity_'//partype//'   /'
        if (ldustdensity)          print*,'&dustdensity_'//partype//'    /'
        if (lneutralvelocity)      print*,'&neutralvelocity_'//partype//'/'
        if (lneutraldensity)       print*,'&neutraldensity_'//partype//' /'
        if (lcosmicray)            print*,'&cosmicray_'//partype//'      /'
        if (lcosmicrayflux)        print*,'&cosmicrayflux_'//partype//'  /'
        if (linterstellar)         print*,'&interstellar_'//partype//'   /'
        if (lshear)                print*,'&shear_'//partype//'          /'
        if (ltestperturb)          print*,'&testperturb_'//partype//'    /'
        if (lspecial)              print*,'&special_'//partype//'        /'
        if (lsolid_cells)          print*,'&solid_cells_'//partype//'    /'
        if (lnscbc)                print*,'&NSCBC_'//partype//'          /'
        if (lpolymer)              print*,'&polymer_'//partype//'        /'
!
        if (.not.lstart) then
          if (lforcing)            print*,'&forcing_'//partype//'        /'
          if (lshock)              print*,'&shock_'//partype//'          /'
          if (lviscosity)          print*,'&viscosity_'//partype//'      /'
          if (lpower_spectrum)     print*,'&power_spectrum_'//partype//' /'
          if (limplicit_diffusion) print*,'&implicit_diffusion_run_pars  /'
        endif
!
!  Particles section.
!
        if (lstart) then
          if (lparticles_density)       print*,'&particles_dens_'//partype//'          /'
        endif
!
        if (lparticles)                 print*,'&particles_'//partype//'               /'
        if (lparticles_radius)          print*,'&particles_radius_'//partype//'        /'
!       if (lparticles_potential)       print*,'&particles_potential_'//partype//'     /'
        if (lparticles_spin)            print*,'&particles_spin_'//partype//'          /'
        if (lparticles_sink)            print*,'&particles_sink_'//partype//'          /'
        if (lparticles_number)          print*,'&particles_number_'//partype//'        /'
        if (lparticles_selfgravity)     print*,'&particles_selfgrav_'//partype//'      /'
        if (lparticles_nbody)           print*,'&particles_nbody_'//partype//'         /'
        if (lparticles_stalker)         print*,'&particles_stalker_'//partype//'       /'
        if (lparticles_mass)            print*,'&particles_mass_'//partype//'          /'
        if (lparticles_drag)            print*,'&particles_drag_'//partype//'          /'
        if (lparticles_temperature)     print*,'&particles_TT_'//partype//'            /'
        if (lparticles_adsorbed)        print*,'&particles_ads_'//partype//'           /'
        if (lparticles_surfspec)        print*,'&particles_surf_'//partype//'          /'
        if (lparticles_chemistry)       print*,'&particles_chem_'//partype//'          /'
!
        if (.not.lstart) then
          if (lparticles_adaptation)    print*,'&particles_adapt_'//partype//'         /'
          if (lparticles_coagulation)   print*,'&particles_coag_'//partype//'          /'
          if (lparticles_collisions)    print*,'&particles_coll_'//partype//'          /'
          if (lparticles_stirring)      print*,'&particles_stirring_'//partype//'      /'
          if (lparticles_viscosity)     print*,'&particles_visc_'//partype//'          /'
          if (lparticles_diagnos_dv)    print*,'&particles_diagnos_dv_'//partype//'    /'
          if (lparticles_diagnos_state) print*,'&particles_diagnos_state_'//partype//' /'
        endif

!!!        if (lstart.and.linitial_condition.and.present(label)) partype=''

        print*,'------END sample namelist -------'
        print*
!
      endif
!
    endsubroutine sample_pars
!***********************************************************************
    subroutine write_all_init_pars(file)
!
!  Print input parameters.
!
!  4-oct-02/wolf: adapted
!  19-aug-15/PABourdin: renamed from 'print_startpars' to 'write_all_init_pars'
!
      use Particles_main, only: write_all_particles_init_pars
!
      character (len=*), optional, intent(in) :: file
!
      character (len=datelen) :: date
      integer :: unit
      logical :: lidl_output
!
      lidl_output = .false.
      if (lroot) then
        if (present(file)) then
          unit = 1
          if (file == 'IDL') then
            open(unit,FILE=trim(datadir)//'/param.nml',DELIM='apostrophe',STATUS='unknown')
            lidl_output = .true.
          else
            open(unit,FILE=file)
            call date_time_string(date)
            write(unit,*) '! -------------------------------------------------------------'
            write(unit,'(A,A)') ' ! ', 'Initializing'
            write(unit,'(A,A)') ' ! Date: ', trim(date)
            write(unit,*) '! t=', t
          endif
        else
          unit = 6   ! default unit is 6=stdout
        endif
!
! If the boundary condition has been changed by subroutines then
! corresponding bc arrays need to be changed too. At present (june 2012)
! it is done only for the bcx array because parker_wind.f90 is the only
! subroutine that changes boundary conditions and it does it only for
! the density in the x direction.
!
        if (lidl_output .and. lreset_boundary_values) call inverse_parse_bc(bcx,bcx12)
!
        write(unit,NML=init_pars)
!
        call write_initial_condition_pars(unit)
        call write_streamlines_init_pars(unit)
        call write_eos_init_pars(unit)
        call write_hydro_init_pars(unit)
        call write_density_init_pars(unit)
        call write_gravity_init_pars(unit)
        call write_selfgravity_init_pars(unit)
        call write_poisson_init_pars(unit)
        call write_energy_init_pars(unit)
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
        call write_special_init_pars(unit)
        call write_solid_cells_init_pars(unit)
        call write_NSCBC_init_pars(unit)
        call write_polymer_init_pars(unit)
!
        call write_all_particles_init_pars(unit)
!
        if (lidl_output) call write_IDL_logicals(unit)
!
        if (present(file)) close(unit)
      endif
!
    endsubroutine write_all_init_pars
!***********************************************************************
    subroutine write_all_run_pars(file)
!
!  Print input parameters.
!
!  14-sep-01/axel: inserted from run.f90
!  31-may-02/wolf: renamed from 'cprint' to 'print_runpars'
!   4-oct-02/wolf: added log file stuff
!  19-aug-15/PABourdin: renamed from 'print_runpars' to 'write_all_run_pars'
!
      use Particles_main, only: write_all_particles_run_pars
!
      character (len=*), optional :: file
!
      integer :: unit=1
      character (len=linelen) :: line
      character (len=datelen) :: date
      logical :: lidl_output
!
      lidl_output = .false.
      if (lroot) then
        if (lreloading) then
          ! get first line from file RELOAD
          line = read_line_from_file('RELOAD')
          if ((trim(line) == '') .or. (line == char(0))) line = 'Reloading'
        else
          line = 'Running'
        endif
!
        if (present(file)) then
          if (file == 'stdout') then
            unit = 6
          elseif (file == 'IDL') then
            open(unit,FILE=trim(datadir)//'/param2.nml',DELIM='apostrophe')
            lidl_output = .true.
          else
            open(unit,FILE=file,position='append')
          endif
        else
          open(unit,FILE=trim(datadir)//'/params.log',position='append')
        endif
!
        if (.not. lidl_output) then
          ! Add separator, comment, and time.
          call date_time_string(date)
          write(unit,*) '! -------------------------------------------------------------'
          write(unit,'(A,A)') ' ! ', trim(line)
          write(unit,'(A,A)') ' ! Date: ', trim(date)
          write(unit,*) '! t=', t
        endif
!
        write(unit,NML=run_pars)
!
        call write_streamlines_run_pars(unit)
        call write_eos_run_pars(unit)
        call write_hydro_run_pars(unit)
        call write_density_run_pars(unit)
        call write_forcing_run_pars(unit)
        call write_gravity_run_pars(unit)
        call write_selfgravity_run_pars(unit)
        call write_poisson_run_pars(unit)
        call write_energy_run_pars(unit)
!       call write_conductivity_run_pars(unit)
        call write_detonate_run_pars(unit)
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
        call write_power_spectrum_run_pars(unit)
        call write_polymer_run_pars(unit)
        call write_implicit_diff_run_pars(unit)
!
        call write_all_particles_run_pars(unit)
!
        if (unit /= 6) close(unit)
!
      endif
!
    endsubroutine write_all_run_pars
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
      if (label=='check_consistency_of_lperi'.and.ip>13) return
      if (label=='check_consistency_of_lperi') lwarning=.false.
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
        if (label=='check_consistency_of_lperi') then
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
    subroutine write_IDL_logicals(unit)
!
!  Output logicals in namelist format for IDL.
!
      integer, intent(in) :: unit
!
      write(unit,'(A)') "&lphysics"
      write(unit,'(A,L,A)') " lhydro=", lhydro, ","
      write(unit,'(A,L,A)') " ldensity=", ldensity, ","
      write(unit,'(A,L,A)') " lentropy=", lentropy, ","
      write(unit,'(A,L,A)') " ltemperature=", ltemperature, ","
      write(unit,'(A,L,A)') " lshock=", lshock, ","
      write(unit,'(A,L,A)') " lmagnetic=", lmagnetic, ","
      write(unit,'(A,L,A)') " lforcing=", lforcing, ","
      write(unit,'(A,L,A)') " llorenz_gauge=", llorenz_gauge, ","
      write(unit,'(A,L,A)') " ldustvelocity=", ldustvelocity, ","
      write(unit,'(A,L,A)') " ldustdensity=", ldustdensity, ","
      write(unit,'(A,L,A)') " ltestscalar=", ltestscalar, ","
      write(unit,'(A,L,A)') " ltestfield=", ltestfield, ","
      write(unit,'(A,L,A)') " ltestflow=", ltestflow, ","
      write(unit,'(A,L,A)') " linterstellar=", linterstellar, ","
      write(unit,'(A,L,A)') " lcosmicray=", lcosmicray, ","
      write(unit,'(A,L,A)') " lcosmicrayflux=", lcosmicrayflux, ","
      write(unit,'(A,L,A)') " lshear=", lshear, ","
      write(unit,'(A,L,A)') " lpscalar=", lpscalar, ","
      write(unit,'(A,L,A)') " lradiation=", lradiation, ","
      write(unit,'(A,L,A)') " leos=", leos, ","
      write(unit,'(A,L,A)') " lchiral=", lchiral, ","
      write(unit,'(A,L,A)') " lneutralvelocity=", lneutralvelocity, ","
      write(unit,'(A,L,A)') " lneutraldensity=", lneutraldensity, ","
      write(unit,'(A,L,A)') " lpolymer=", lpolymer, ","
      write(unit,'(A,L,A)') " lsolid_cells=", lsolid_cells, ","
      write(unit,'(A,L,A)') " lpower_spectrum=", lpower_spectrum, ","
      write(unit,'(A,L,A)') " lparticles=", lparticles, ","
      write(unit,'(A,L,A)') " lparticles_drag=", lparticles_drag, ","
      write(unit,'(A)') " /"
!
    endsubroutine write_IDL_logicals
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
