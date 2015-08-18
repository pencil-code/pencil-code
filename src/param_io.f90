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
  use Particles
  use Particles_adaptation
  use Particles_adsorbed
  use Particles_cdata
  use Particles_chemistry
  use Particles_coagulation
  use Particles_collisions
  use Particles_density
  use Particles_diagnos_dv
  use Particles_diagnos_state
  use Particles_drag
  use Particles_main
  use Particles_map
  use Particles_mass
  use Particles_mpicomm
  use Particles_nbody
  use Particles_number
  use Particles_radius
  use Particles_sink
  use Particles_spin
  use Particles_selfgravity
  use Particles_stalker
  use Particles_stirring
  use Particles_sub
  use Particles_surfspec
  use Particles_temperature
  use Particles_viscosity
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
  public :: read_startpars, print_startpars
  public :: read_runpars,   print_runpars
  public :: wparam, wparam2, write_pencil_info
!
! The following fixes namelist problems withi MIPSpro 7.3.1.3m
! under IRIX -- at least for the moment
!
  character (len=labellen) :: mips_is_buggy='system'
!
  logical :: lforce_shear_bc = .true.
  logical :: lnamelist_error, lignore_namelist_errors
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
      ldirect_access
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
      use File_io, only: get_unit
!
      integer, intent(out) :: iostat
      include "parallel_unit.h"
!
      read(parallel_unit, NML=init_pars, IOSTAT=iostat)
!
    endsubroutine read_init_pars
!***********************************************************************
    subroutine read_run_pars(iostat)
!
      use File_io, only: get_unit
!
      integer, intent(out) :: iostat
      include "parallel_unit.h"
!
      read(parallel_unit, NML=run_pars, IOSTAT=iostat)
!
    endsubroutine read_run_pars
!***********************************************************************
    subroutine read_startpars(print)
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
!  31-oct-13/MR  : shortened by use of read_pars
!  18-aug-15/PABourdin : reworked to simplify code and display all errors at once
!
      use File_io, only: parallel_open, parallel_rewind, parallel_close
!
      logical, optional, intent(IN) :: print
!
      include "parallel_unit_declaration.h"
!
      lnamelist_error = .false.
      lignore_namelist_errors = .false.
!
!  Set default to shearing sheet if lshear=.true. (even when Sshear==0.).
!
      if (lshear .and. lforce_shear_bc) bcx(:)='she'
!
!  Open namelist file and read through all items that *may* be present in the various modules.
!
      if (lstart) then
        call parallel_open(parallel_unit,'start.in','formatted')
      else
        call parallel_open(parallel_unit,trim(datadir)//'/param.nml')
        lignore_namelist_errors = .true.
      endif
!
      call read_pars(read_init_pars                ,'init')
      call read_pars(read_initial_condition_pars   ,'initial_condition_pars')
      call read_pars(read_streamlines_init_pars    ,'streamlines')
      call read_pars(read_eos_init_pars            ,'eos')
      call read_pars(read_hydro_init_pars          ,'hydro')
      call read_pars(read_density_init_pars        ,'density')
      call read_pars(read_gravity_init_pars        ,'grav')
      call read_pars(read_selfgravity_init_pars    ,'selfgrav')
      call read_pars(read_poisson_init_pars        ,'poisson')
      call read_pars(read_energy_init_pars         ,'entropy')
      call read_pars(read_magnetic_init_pars       ,'magnetic')
      call read_pars(read_lorenz_gauge_init_pars   ,'lorenz_gauge')
      call read_pars(read_testscalar_init_pars     ,'testscalar')
      call read_pars(read_testfield_init_pars      ,'testfield')
      call read_pars(read_testflow_init_pars       ,'testflow')
      call read_pars(read_radiation_init_pars      ,'radiation')
      call read_pars(read_pscalar_init_pars        ,'pscalar')
      call read_pars(read_chiral_init_pars         ,'chiral')
      call read_pars(read_chemistry_init_pars      ,'chemistry')
      call read_pars(read_signal_init_pars         ,'signal')
      call read_pars(read_dustvelocity_init_pars   ,'dustvelocity')
      call read_pars(read_dustdensity_init_pars    ,'dustdensity')
      call read_pars(read_neutralvelocity_init_pars,'neutralvelocity')
      call read_pars(read_neutraldensity_init_pars ,'neutraldensity')
      call read_pars(read_cosmicray_init_pars      ,'cosmicray')
      call read_pars(read_cosmicrayflux_init_pars  ,'cosmicrayflux')
      call read_pars(read_interstellar_init_pars   ,'interstellar')
      call read_pars(read_shear_init_pars          ,'shear')
      call read_pars(read_special_init_pars        ,'special')
      call read_pars(read_solid_cells_init_pars    ,'solid_cells')
      call read_pars(read_NSCBC_init_pars          ,'NSCBC')
      call read_pars(read_polymer_init_pars        ,'polymer')
!
      call read_pars(read_particles_init_pars      ,'particles',lparticles)
      call read_pars(read_particles_rad_init_pars  ,'particles_radius',lparticles_radius)
!     call read_pars(read_particles_pot_init_pars  ,'particles_potential',lparticles_potential)
      call read_pars(read_particles_spin_init_pars ,'particles_spin',lparticles_spin)
      call read_pars(read_particles_sink_init_pars ,'particles_sink',lparticles_sink)
      call read_pars(read_particles_num_init_pars  ,'particles_number',lparticles_number)
      call read_pars(read_particles_dens_init_pars ,'particles_dens',lparticles_density)
      call read_pars(read_particles_selfg_init_pars,'particles_selfgrav',lparticles_selfgravity)
      call read_pars(read_particles_nbody_init_pars,'particles_nbody',lparticles_nbody)
      call read_pars(read_particles_mass_init_pars ,'particles_mass',lparticles_mass)
      call read_pars(read_particles_drag_init_pars ,'particles_drag',lparticles_drag)
      call read_pars(read_particles_TT_init_pars   ,'particles_TT',lparticles_temperature)
      call read_pars(read_particles_ads_init_pars  ,'particles_ads',lparticles_adsorbed)
      call read_pars(read_particles_surf_init_pars ,'particles_surf',lparticles_surfspec)
      call read_pars(read_particles_chem_init_pars ,'particles_chem',lparticles_chemistry)
      call read_pars(read_pstalker_init_pars       ,'particles_stalker',lparticles_stalker)
!
      call parallel_close(parallel_unit)
!
      if (lnamelist_error) then
        call sample_pars()
        call fatal_error ('read_startpars', 'Please fix all above WARNINGs for "start.in".')
      endif
!
!  Print SVN id from first line.
!
      if (lroot) call svn_id(cvsid)
!
!  Give online feedback if called with the PRINT optional argument.
!
      if (present (print)) then
        if (print) call print_startpars()
      endif
!
!  Write parameters to log file.
!
      if (lstart) call print_startpars(FILE=trim(datadir)//'/params.log')
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
      call check_consistency_of_lperi('read_startpars')
!
    endsubroutine read_startpars
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
        call particles_wparam(unit)
!
        if (present(file)) then
          close(unit)
        endif
      endif
!
    endsubroutine print_startpars
!***********************************************************************
    subroutine read_pars(reader,name,lactive)
!
!  encapsulates reading of pars + error handling
!
!  31-oct-13/MR: coded
!  16-dec-13/MR: handling of optional ierr corrected
!  18-dec-13/MR: changed handling of ierr to avoid compiler trouble
!  19-dec-13/MR: changed ierr into logical lierr to avoid compiler trouble
!  11-jul-14/MR: end-of-file caught to avoid program crash when a namelist is missing
!  18-aug-15/PABourdin: reworked to simplify code and display all errors at once
!
      use File_io, only: parallel_rewind, get_unit
!
      interface
        subroutine reader(iostat)
          integer, intent(out) :: iostat
        endsubroutine reader
      endinterface
!
      character(len=*), intent(in) :: name
      logical, intent(in), optional :: lactive
!
      integer :: ierr
      character(len=linelen) :: type
      include "parallel_unit.h"
!
      if (present (lactive)) then
        if (.not. lactive) return
      endif
!
      call reader(ierr)
!
      if (lstart) then
        type = '_init'
        if (name == 'init') type = ''
      else
        type = '_run'
        if (name == 'run') type = ''
      endif
!
      if (.not. lignore_namelist_errors) then
        if (ierr == -1) then
          if (lroot) call warning ('read_pars', 'namelist "'//trim(name)//trim(type)//'_pars" missing!')
          lnamelist_error = .true.
        elseif (ierr /= 0) then
          if (lroot) call warning ('read_pars', 'namelist "'//trim(name)//trim(type)//'_pars" has an error!')
          lnamelist_error = .true.
        endif
      endif
! 
      call parallel_rewind(parallel_unit)
!
    endsubroutine read_pars
!***********************************************************************
    subroutine sample_pars()
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
        if (lrun) then
          if (lforcing)            print*,'&forcing_'//partype//'        /'
          if (lshock)              print*,'&shock_'//partype//'          /'
          if (lviscosity)          print*,'&viscosity_'//partype//'      /'
          if (lpower_spectrum)     print*,'&power_spectrum_'//partype//' /'
          if (limplicit_diffusion) print*,'&implicit_diffusion_pars  /'
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
        if (lrun) then
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
    subroutine read_runpars(logging)
!
!  Read input parameters.
!
!  14-sep-01/axel: inserted from run.f90
!  31-may-02/wolf: renamed from cread to read_runpars
!   6-jul-02/axel: in case of error, print sample namelist
!  21-oct-03/tony: moved sample namelist stuff to a separate procedure
!  19-dec-13/MR  : adapted calls to read_pars
!  18-aug-15/PABourdin : reworked to simplify code and display all errors at once
!
      use Dustvelocity, only: copy_bcs_dust
      use File_io, only: parallel_open, parallel_rewind, parallel_close
      use Sub, only: parse_bc
!
      logical, optional :: logging
      include "parallel_unit_declaration.h"
!
      lnamelist_error = .false.
      lignore_namelist_errors = .false.
      tstart=impossible
!
!  Open namelist file.
!
      call parallel_open(parallel_unit,'run.in','formatted')
!
!  Read through all items that *may* be present in the various modules.
!  AB: at some point the sgi_fix stuff should probably be removed (see sgi bug)
!
      call read_pars(read_run_pars                ,'run')
      call read_pars(read_streamlines_run_pars    ,'streamlines')
      call read_pars(read_eos_run_pars            ,'eos')
      call read_pars(read_hydro_run_pars          ,'hydro')
      call read_pars(read_density_run_pars        ,'density')
      call read_pars(read_forcing_run_pars        ,'forcing')
      call read_pars(read_gravity_run_pars        ,'grav')
      call read_pars(read_selfgravity_run_pars    ,'selfgrav')
      call read_pars(read_poisson_run_pars        ,'poisson')
      call read_pars(read_energy_run_pars         ,'entropy')
!     call read_pars(read_conductivity_run_pars   ,'conductivity')
      call read_pars(read_detonate_run_pars       ,'detonate')
      call read_pars(read_magnetic_run_pars       ,'magnetic')
      call read_pars(read_lorenz_gauge_run_pars   ,'lorenz_gauge')
      call read_pars(read_testscalar_run_pars     ,'testscalar')
      call read_pars(read_testfield_run_pars      ,'testfield')
      call read_pars(read_testflow_run_pars       ,'testflow')
      call read_pars(read_radiation_run_pars      ,'radiation')
      call read_pars(read_pscalar_run_pars        ,'pscalar')
      call read_pars(read_chiral_run_pars         ,'chiral')
      call read_pars(read_chemistry_run_pars      ,'chemistry')
      call read_pars(read_dustvelocity_run_pars   ,'dustvelocity')
      call read_pars(read_dustdensity_run_pars    ,'dustdensity')
      call read_pars(read_neutralvelocity_run_pars,'neutralvelocity')
      call read_pars(read_neutraldensity_run_pars ,'neutraldensity')
      call read_pars(read_cosmicray_run_pars      ,'cosmicray')
      call read_pars(read_cosmicrayflux_run_pars  ,'cosmicrayflux')
      call read_pars(read_heatflux_run_pars       ,'heatflux')
      call read_pars(read_interstellar_run_pars   ,'interstellar')
      call read_pars(read_shear_run_pars          ,'shear')
      call read_pars(read_testperturb_run_pars    ,'testperturb')
      call read_pars(read_viscosity_run_pars      ,'viscosity')
      call read_pars(read_special_run_pars        ,'special')
      call read_pars(read_shock_run_pars          ,'shock')
      call read_pars(read_solid_cells_run_pars    ,'solid_cells')
      call read_pars(read_NSCBC_run_pars          ,'NSCBC')
      call read_pars(read_polymer_run_pars        ,'polymer')
      call read_pars(read_power_spectrum_runpars  ,'power_spectrum')
      call read_pars(read_implicit_diffusion_pars ,'implicit_diffusion')
!
      call read_pars(read_particles_run_pars          ,'particles',lparticles)
      call read_pars(read_particles_adapt_run_pars    ,'particles_adapt',lparticles_adaptation)
      call read_pars(read_particles_rad_run_pars      ,'particles_radius',lparticles_radius)
!     call read_pars(read_particles_potential_run_pars,'particles_potential',lparticles_potential)
      call read_pars(read_particles_spin_run_pars     ,'particles_spin',lparticles_spin)
      call read_pars(read_particles_sink_run_pars     ,'particles_sink',lparticles_sink)
      call read_pars(read_particles_num_run_pars      ,'particles_number',lparticles_number)
      call read_pars(read_particles_selfg_run_pars    ,'particles_selfgrav',lparticles_selfgravity)
      call read_pars(read_particles_nbody_run_pars    ,'particles_nbody',lparticles_nbody)
      call read_pars(read_particles_visc_run_pars     ,'particles_visc',lparticles_viscosity)
      call read_pars(read_particles_coag_run_pars     ,'particles_coag',lparticles_coagulation)
      call read_pars(read_particles_coll_run_pars     ,'particles_coll',lparticles_collisions)
      call read_pars(read_particles_stir_run_pars     ,'particles_stirring',lparticles_stirring)
      call read_pars(read_pstalker_run_pars           ,'particles_stalker',lparticles_stalker)
      call read_pars(read_pars_diagnos_dv_run_pars    ,'particles_diagnos_dv',lparticles_diagnos_dv)
      call read_pars(read_pars_diag_state_run_pars    ,'particles_diagnos_state',lparticles_diagnos_state)
      call read_pars(read_particles_mass_run_pars     ,'particles_mass',lparticles_mass)
      call read_pars(read_particles_drag_run_pars     ,'particles_drag',lparticles_drag)
      call read_pars(read_particles_TT_run_pars       ,'particles_TT',lparticles_temperature)
      call read_pars(read_particles_ads_run_pars      ,'particles_ads',lparticles_adsorbed)
      call read_pars(read_particles_surf_run_pars     ,'particles_surf',lparticles_surfspec)
      call read_pars(read_particles_chem_run_pars     ,'particles_chem',lparticles_chemistry)
!
      call parallel_close(parallel_unit)
!
      if (lnamelist_error) then
        call sample_pars()
        call fatal_error ('read_startpars', 'Please fix all above WARNINGs for "start.in".')
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
      if (present (logging)) then
        if (logging) call print_runpars()
      endif
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
      call check_consistency_of_lperi('read_runpars')
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
!
      contains
!
!---------------------------------------------------------------------------------------
!
      subroutine get_downpars(ind,n,ip)
!
! Calculates start indices & lengths for downsampled output
!
! 13-feb-14/MR: coded
!
        use General, only: get_range_no
!
! coordinate direction, number of inner grid points, processor number
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
!
    endsubroutine read_runpars
!***********************************************************************
    subroutine print_runpars(file)
!
!  Print input parameters.
!
!  14-sep-01/axel: inserted from run.f90
!  31-may-02/wolf: renamed from cprint to print_runpars
!   4-oct-02/wolf: added log file stuff
!
      character (len=*), optional :: file
!
      integer :: unit=1
      character (len=linelen) :: line
      character (len=datelen) :: date
!
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
          else
            open(unit,FILE=file,position='append')
          endif
        else
          open(unit,FILE=trim(datadir)//'/params.log',position='append')
        endif
!
!  Add separator, comment, and time.
!
        call date_time_string(date)
        write(unit,*) &
            '! -------------------------------------------------------------'
        write(unit,'(A,A)') ' ! ', trim(line)
        write(unit,'(A,A)') ' ! Date: ', trim(date)
        write(unit,*) '! t=', t
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
!        call write_conductivity_run_pars(unit)
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
        call write_power_spectrum_runpars(unit)
        call write_polymer_run_pars(unit)
        call particles_wparam2(unit)
        call write_implicit_diffusion_pars(unit)
!
        if (unit /= 6) close(unit)
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
    subroutine wparam()
!
!  Write startup parameters
!
!  21-jan-02/wolf: coded
!
      logical :: lhydro           = lhydro_var
      logical :: ldensity         = ldensity_var
      logical :: lentropy         = lentropy_var
!      logical :: lenergy          = lenergy_var
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
      logical :: lsolid_cells     = lsolid_cells_var
      logical :: lparticles       = lparticles_var
      logical :: lparticles_drag  = lparticles_drag_var
      integer, parameter :: unit=1
!
      namelist /lphysics/ &
          lhydro, ldensity, lentropy, lenergy, lmagnetic, lshear, &
          llorenz_gauge,  ltestscalar, ltestfield, ltestflow, lpscalar, &
          lradiation, ldustvelocity, ldustdensity, lforcing, lgravz, lgravr, &
          ltestperturb, linterstellar, lcosmicray, lcosmicrayflux, lshock, &
          lradiation_fld, leos_ionization, leos_fixed_ionization, &
          lvisc_hyper, lchiral, leos, leos_temperature_ionization, &
          lneutralvelocity, lneutraldensity, ltemperature,lpolymer, &
          lsolid_cells, &
          lparticles, lparticles_drag
!
! If the boundary condition has been changed by subroutines then
! corresponding bc arrays need to be changed too. At present (june 2012)
! it is done only for the bcx array because parker_wind.f90 is the only
! subroutine that changes boundary conditions and it does it only for
! the density in the x direction.
!
        if  (lreset_boundary_values) then
          call inverse_parse_bc(bcx,bcx12)
        endif
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
        call write_streamlines_init_pars(unit)
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
      call keep_compiler_quiet(lenergy)
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
      call keep_compiler_quiet(lsolid_cells)
!
    endsubroutine wparam
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
        call write_streamlines_run_pars(unit)
        call write_eos_run_pars(unit)
        call write_hydro_run_pars(unit)
        call write_density_run_pars(unit)
        call write_forcing_run_pars(unit)
        call write_gravity_run_pars(unit)
        call write_selfgravity_run_pars(unit)
        call write_poisson_run_pars(unit)
        call write_energy_run_pars(unit)
!        call write_conductivity_run_pars(unit)
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
