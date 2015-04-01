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
!
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
      lcorotational_frame, rcorot, &
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
      uxj_spec, vec_spec, ou_spec, ab_spec, azbz_spec, ub_spec, &
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
!  31-oct-13/MR  : changed for use instead of rparam; shortened by use of read_pars
!
      use Mpicomm, only: parallel_open, parallel_close
      use General, only: rewind
!
      logical, optional, intent(IN) :: print
!
      include 'unit_loc.h'
!
      integer :: ierr
      logical :: lierrl
!
!  Set default to shearing sheet if lshear=.true. (even when Sshear==0.).
!
      if (lshear .and. lforce_shear_bc) bcx(:)='she'
!
!  Open namelist file and read through all items that *may* be present in the various modules.
!
      if (lstart) then
        call parallel_open(unit,'start.in','formatted')
        read(unit,NML=init_pars,IOSTAT=ierr)
        if (ierr/=0) call sample_pars(ierr,'')
        lierrl=.true.
      else
        call parallel_open(unit,trim(datadir)//'/param.nml')
        read(unit,NML=init_pars)
        lierrl=.false.
      endif
!
      call rewind(unit)
!
      call read_pars(unit,read_initial_condition_pars   ,'initial_condition_pars',lierrl)
      call read_pars(unit,read_streamlines_init_pars    ,'streamlines',lierrl)
      call read_pars(unit,read_eos_init_pars            ,'eos',lierrl)
      call read_pars(unit,read_hydro_init_pars          ,'hydro',lierrl)
      call read_pars(unit,read_density_init_pars        ,'density',lierrl)
      call read_pars(unit,read_forcing_init_pars        ,'forcing',lierrl)
      call read_pars(unit,read_gravity_init_pars        ,'grav_init_pars',lierrl)
      call read_pars(unit,read_selfgravity_init_pars    ,'selfgrav',lierrl)
      call read_pars(unit,read_poisson_init_pars        ,'poisson',lierrl)
      call read_pars(unit,read_energy_init_pars         ,'entropy',lierrl)
      call read_pars(unit,read_magnetic_init_pars       ,'magnetic',lierrl)
      call read_pars(unit,read_lorenz_gauge_init_pars   ,'lorenz_gauge',lierrl)
      call read_pars(unit,read_testscalar_init_pars     ,'testscalar',lierrl)
      call read_pars(unit,read_testfield_init_pars      ,'testfield',lierrl)
      call read_pars(unit,read_testflow_init_pars       ,'testflow',lierrl)
      call read_pars(unit,read_radiation_init_pars      ,'radiation',lierrl)
      call read_pars(unit,read_pscalar_init_pars        ,'pscalar',lierrl)
      call read_pars(unit,read_chiral_init_pars         ,'chiral',lierrl)
      call read_pars(unit,read_chemistry_init_pars      ,'chemistry',lierrl)
      call read_pars(unit,read_signal_init_pars         ,'signal',lierrl)
      call read_pars(unit,read_dustvelocity_init_pars   ,'dustvelocity',lierrl)
      call read_pars(unit,read_dustdensity_init_pars    ,'dustdensity',lierrl)
      call read_pars(unit,read_neutralvelocity_init_pars,'neutralvelocity',lierrl)
      call read_pars(unit,read_neutraldensity_init_pars ,'neutraldensity',lierrl)
      call read_pars(unit,read_cosmicray_init_pars      ,'cosmicray',lierrl)
      call read_pars(unit,read_cosmicrayflux_init_pars  ,'cosmicrayflux',lierrl)
      call read_pars(unit,read_heatflux_init_pars       ,'heatflux',lierrl)
      call read_pars(unit,read_interstellar_init_pars   ,'interstellar',lierrl)
      call read_pars(unit,read_shear_init_pars          ,'shear',lierrl)
      call read_pars(unit,read_testperturb_init_pars    ,'testperturb',lierrl)
      call read_pars(unit,read_viscosity_init_pars      ,'viscosity',lierrl)
      call read_pars(unit,read_special_init_pars        ,'special',lierrl)
      call read_pars(unit,read_solid_cells_init_pars    ,'solid_cells',lierrl)
      call read_pars(unit,read_NSCBC_init_pars          ,'NSCBC',lierrl)
      call read_pars(unit,read_polymer_init_pars        ,'polymer',lierrl)
      call read_pars(unit,particles_read_startpars      ,'particles',lierrl)
!
      call parallel_close(unit)
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
        call write_streamlines_init_pars(unit)
        call write_eos_init_pars(unit)
        call write_hydro_init_pars(unit)
        call write_density_init_pars(unit)
        call write_forcing_init_pars(unit)
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
        call write_testperturb_init_pars(unit)
        call write_viscosity_init_pars(unit)
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
    subroutine read_pars(unit,reader,name,lierr)
!
!  encapsulates reading of pars + error handling
!
!  31-oct-13/MR: coded
!  16-dec-13/MR: handling of optional ierr corrected
!  18-dec-13/MR: changed handling of ierr to avoid compiler trouble
!  19-dec-13/MR: changed ierr into logical lierr to avoid compiler trouble
!  11-jul-14/MR: end-of-file caught to avoid program crash when a namelist is missing
!
    use General, only: rewind
    include 'unit.h'
!
    interface
      subroutine reader(unit, iostat)
        include 'unit.h'
        integer, intent(inout), optional :: iostat
      endsubroutine reader
    endinterface
!
    character(LEN=*), intent(IN):: name
    logical,          intent(IN):: lierr
!
    integer :: ierr
!
    if (lierr) then
      ierr=0
      call reader(unit,ierr)
      if (ierr/=0.and.ierr/=-1) call sample_pars(ierr,name)
      if (lroot.and.ierr==-1) then
        if (lstart) then
          print*, 'read_pars: Warning - namelist "', trim(name), '_start_pars" missing!'
        else
          print*, 'read_pars: Warning - namelist "', trim(name), '_run_pars" missing!'
        endif
      endif
    else
      ierr=0
      call reader(unit,ierr)
    endif
! 
    call rewind(unit)
!
    endsubroutine read_pars
!***********************************************************************
    subroutine sample_pars(iostat,label)
!
!  standardises handling of errors in namelist, both for init and run
!
!  31-oct-13/MR: coded
!
      integer         , optional :: iostat
      character(len=*), optional :: label
!
      character(len=labellen):: partype
      character(len=linelen):: msg

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
          if (lsignal           ) print*,'&signal_init_pars       /'
          if (linitial_condition) print*,'&initial_condition_pars /'
        endif

        if (ltracers        ) print*,'&streamlines_'//partype//'    /'     !! questionable wg. ltracers
        if (leos            ) print*,'&eos_'//partype//'            /'
        if (lhydro .or.                     &
            lhydro_kinematic) print*,'&hydro_'//partype//'          /'
        if (ldensity        ) print*,'&density_'//partype//'        /'
        if (lgrav           ) print*,'&grav_'//partype//'           /'
        if (lselfgravity    ) print*,'&selfgrav_'//partype//'       /'
        if (lpoisson        ) print*,'&poisson_'//partype//'        /'
        if (lentropy.or.                    &
            ltemperature    ) print*,'&entropy_'//partype//'        /'
        if (lthermal_energy ) print*,'&energy_'//partype//'         /'
!        if (lconductivity   ) print*,'&conductivity_'//partype//'   /'
        if (ldetonate       ) print*,'&detonate_'//partype//'       /'
        if (lmagnetic       ) print*,'&magnetic_'//partype//'       /'
        if (lmagn_mf        ) print*,'&magn_mf_'//partype//'        /'
        if (llorenz_gauge   ) print*,'&lorenz_gauge_'//partype//'   /'
        if (ltestscalar     ) print*,'&testscalar_'//partype//'     /'
        if (ltestfield      ) print*,'&testfield_'//partype//'      /'
        if (ltestflow       ) print*,'&testflow_'//partype//'       /'
        if (lradiation      ) print*,'&radiation_'//partype//'      /'
        if (lpscalar        ) print*,'&pscalar_'//partype//'        /'
        if (lchiral         ) print*,'&chiral_'//partype//'         /'
        if (lchemistry      ) print*,'&chemistry_'//partype//'      /'
        if (ldustvelocity   ) print*,'&dustvelocity_'//partype//'   /'
        if (ldustdensity    ) print*,'&dustdensity_'//partype//'    /'
        if (lneutralvelocity) print*,'&neutralvelocity_'//partype//'/'
        if (lneutraldensity ) print*,'&neutraldensity_'//partype//' /'
        if (lcosmicray      ) print*,'&cosmicray_'//partype//'      /'
        if (lcosmicrayflux  ) print*,'&cosmicrayflux_'//partype//'  /'
        if (linterstellar   ) print*,'&interstellar_'//partype//'   /'
        if (lshear          ) print*,'&shear_'//partype//'          /'
        if (ltestperturb    ) print*,'&testperturb_'//partype//'    /'
        if (lspecial        ) print*,'&special_'//partype//'        /'
        if (lsolid_cells    ) print*,'&solid_cells_'//partype//'    /'
        if (lnscbc          ) print*,'&NSCBC_'//partype//'          /'
        if (lpolymer        ) print*,'&polymer_'//partype//'        /'
!
        if (lrun) then
          if (lforcing      ) print*,'&forcing_run_pars         /'
          if (lshock        ) print*,'&shock_run_pars           /'
          if (lviscosity    ) print*,'&viscosity_run_pars       /'
          if (lpower_spectrum)print*,'&power_spectrum_run_pars  /'
          if (limplicit_diffusion) &
                              print*,'&implicit_diffusion_pars  /'
        endif

!!!        if (lstart.and.linitial_condition.and.present(label)) partype=''

        print*,'------END sample namelist -------'
        print*
!
        if (present(label)) then
          msg = 'Found error in input namelist "'
          if (label=='') then
            msg = trim(msg)//trim(partype)
          else
            msg = trim(msg)//trim(label)//'_'//trim(partype)
          endif
          print*, trim(msg)//'"'
        endif
        if (present(iostat)) print*, 'iostat = ', iostat
        if (present(iostat).or.present(label)) &
                             print*,  '-- use sample above.'
      endif
!
      call fatal_error('sample_'//partype,'')
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
!  12-nov-10/MR  : added read and write calls for namelist power_spectrum_run_pars
!  18-dec-13/MR  : changed handling of ierr to avoid compiler trouble
!  19-dec-13/MR  : adapted calls to read_pars
!  11-feb-14/MR  : changes for downsampled output
!  13-feb-14/MR  : further preparations for downsampled output
!
      use Dustvelocity, only: copy_bcs_dust
      use Mpicomm, only: parallel_open, parallel_close
      use Sub, only: parse_bc
      use General, only: rewind
!
      include 'unit_loc.h'
!
      integer :: ierr
      logical, optional :: logging
!
!  Open namelist file.
!
      call parallel_open(unit,'run.in','formatted')
!
!  Read through all items that *may* be present in the various modules.
!  AB: at some point the sgi_fix stuff should probably be removed (see sgi bug)
!
      read(unit,NML=run_pars,IOSTAT=ierr)
      if (ierr/=0) call sample_pars(ierr,'')
      call rewind(unit)
!
      call read_pars(unit,read_streamlines_run_pars    ,'streamlines',.true.)
      call read_pars(unit,read_eos_run_pars            ,'eos',.true.)
      call read_pars(unit,read_hydro_run_pars          ,'hydro',.true.)
      call read_pars(unit,read_density_run_pars        ,'density',.true.)
      call read_pars(unit,read_forcing_run_pars        ,'forcing',.true.)
      call read_pars(unit,read_gravity_run_pars        ,'grav',.true.)
      call read_pars(unit,read_selfgravity_run_pars    ,'selfgrav',.true.)
      call read_pars(unit,read_poisson_run_pars        ,'poisson',.true.)
      call read_pars(unit,read_energy_run_pars         ,'entropy',.true.)
!      call read_pars(unit,read_conductivity_run_pars   ,'conductivity',.true.)
      call read_pars(unit,read_detonate_run_pars       ,'detonate',.true.)
      call read_pars(unit,read_magnetic_run_pars       ,'magnetic',.true.)
      call read_pars(unit,read_lorenz_gauge_run_pars   ,'lorenz_gauge',.true.)
      call read_pars(unit,read_testscalar_run_pars     ,'testscalar',.true.)
      call read_pars(unit,read_testfield_run_pars      ,'testfield',.true.)
      call read_pars(unit,read_testflow_run_pars       ,'testflow',.true.)
      call read_pars(unit,read_radiation_run_pars      ,'radiation',.true.)
      call read_pars(unit,read_pscalar_run_pars        ,'pscalar',.true.)
      call read_pars(unit,read_chiral_run_pars         ,'chiral',.true.)
      call read_pars(unit,read_chemistry_run_pars      ,'chemistry',.true.)
      call read_pars(unit,read_dustvelocity_run_pars   ,'dustvelocity',.true.)
      call read_pars(unit,read_dustdensity_run_pars    ,'dustdensity',.true.)
      call read_pars(unit,read_neutralvelocity_run_pars,'neutralvelocity',.true.)
      call read_pars(unit,read_neutraldensity_run_pars ,'neutraldensity',.true.)
      call read_pars(unit,read_cosmicray_run_pars      ,'cosmicray',.true.)
      call read_pars(unit,read_cosmicrayflux_run_pars  ,'cosmicrayflux',.true.)
      call read_pars(unit,read_heatflux_run_pars       ,'heatflux',.true.)
      call read_pars(unit,read_interstellar_run_pars   ,'interstellar',.true.)
      call read_pars(unit,read_shear_run_pars          ,'shear',.true.)
      call read_pars(unit,read_testperturb_run_pars    ,'testperturb',.true.)
      call read_pars(unit,read_viscosity_run_pars      ,'viscosity',.true.)
      call read_pars(unit,read_special_run_pars        ,'special',.true.)
      call read_pars(unit,read_shock_run_pars          ,'shock',.true.)
      call read_pars(unit,read_solid_cells_run_pars    ,'solid_cells',.true.)
      call read_pars(unit,read_NSCBC_run_pars          ,'NSCBC',.true.)
      call read_pars(unit,read_polymer_run_pars        ,'polymer',.true.)
      call read_pars(unit,particles_read_runpars       ,'particles',.true.)
      call read_pars(unit,read_power_spectrum_runpars  ,'power_spectrum',.true.)
      call read_pars(unit,read_implicit_diffusion_pars ,'implicit_diffusion',.true.)
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
!  Write parameters to log file, if requested.
!
      if (present (logging)) then
        if (logging) call print_runpars()
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
      if (any(downsampl>1)) then
!
!  If downsampling, calculate local start indices and number of data in
!  output for each direction; inner ghost zones are here disregarded
!
        ldownsampl = .true.
        if (dsnap_down<=0.) dsnap_down=dsnap
!
        call get_downpars(1,nx,ipx,lfirst_proc_x,llast_proc_x)
        call get_downpars(2,ny,ipy,lfirst_proc_y,llast_proc_y)
        call get_downpars(3,nz,ipz,lfirst_proc_z,llast_proc_z)
!
      endif
!
      contains
!
      subroutine get_downpars(ind,n,ip,lfirst_proc,llast_proc)
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
        logical, intent(IN) :: lfirst_proc,llast_proc
!
      if ( downsampl(ind)>n ) then
        print*, 'read_startpars: Warning - stepsize for downsampling in '// &
        coornames(ind)//'-direction ', downsampl(ind), 'greater than grid size ', n, &
        '! Set to grid size.'
        downsampl(ind)=n
      endif
!
! first index in direction ind in local farray for output
!
        firstind(ind) = downsampl(ind) - modulo(ip*n-1,downsampl(ind))
!
! number of output items in direction ind without ghost zones
!
        ndowni(ind) = get_range_no((/firstind(ind),n,downsampl(ind)/),1)
!
        firstind(ind) = firstind(ind) + nghost
!
! number of output items in direction ind with ghost zones;
! target index corresponing to firstind
!
        ndown(ind) = ndowni(ind)
        if (lfirst_proc) then
          ndown(ind) = ndown(ind)+nghost; startind(ind)=nghost+1
        else
          startind(ind)=1
        endif
        if (llast_proc) ndown(ind) = ndown(ind)+nghost
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
          lsolid_cells
!
! If the boundary condition has been changed by subroutines then
! corresponding bc arrays need to be changed too. At present (june 2012)
! it is done only for the bcx array because parker_wind.f90 is the only
! subroutine that changes boundary conditions and it does it only for
! the density in the x direction.
!
        if  (lreset_boundary_values) then
          call inverse_parse_bc(bcx,bcx1,bcx2)
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
        call write_forcing_init_pars(unit)
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
        call write_testperturb_init_pars(unit)
        call write_viscosity_init_pars(unit)
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
