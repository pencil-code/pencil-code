! $Id: param_io.f90 13354 2010-02-26 08:41:01Z AxelBrandenburg $
!
!  IO of init and run parameters. Subroutines here are `at the end of the
!  food chain', i.e. depend on all physics modules plus possibly others.
!  Using this module is also a compact way of referring to all physics
!  modules at once.
!
module Param_IO
!
  use Cdata
  use Cparam
  use Density
  use Entropy
  use EquationOfState
  use General
  use Gravity
  use Hydro
  use InitialCondition
  use Magnetic
  use Messages
  use Shear
  use Shock
  use Signal_handling
  use Sub
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
!AB/15-Mar-07: the lcylindrical is kept to produce a warning: outdated
  logical :: lcylindrical
!
  namelist /init_pars/ &
      cvsid, ip, xyz0, xyz1, Lxyz, lperi, lshift_origin, &
      coeff_grid, zeta_grid0, xyz_star, lwrite_ic, &
      lnowrite, unit_system, unit_length, &
      unit_velocity, unit_density, unit_temperature, unit_magnetic, c_light, &
      G_Newton, hbar, random_gen, seed0, lserial_io, der2_type, &
      lread_oldsnap, lread_oldsnap_nomag, lread_aux, lwrite_aux, &
      lprocz_slowest, &
      lcopysnapshots_exp, bcx, bcy, bcz, r_ref, rsmooth, &
      mu0, &
      tstart, fbcx1, fbcx2, fbcy1, fbcy2, fbcz1, fbcz2, fbcz1_1, fbcz1_2, &
      fbcz2_1, fbcz2_2, fbcx1_2, fbcx2_2, xyz_step, xi_step_frac, &
      xi_step_width, &
      init_loops, lwrite_2d, &
      yequator, lequatory, &
      lequatorz, zequator, lav_smallx, xav_max 
!
  namelist /run_pars/ &
      cvsid, ip, nt, it1, it1d, dt, cdt, ddt, cdtv, cdts, cdtr, &
      cdtc, isave, itorder, dsnap, d2davg, dvid, dtmin, dspec, tmax, &
      ialive, max_walltime, dtmax, ldt_paronly, &
      isaveglobal, &
      random_gen, der2_type, &
      ltime_integrals, lnowrite, noghost_for_isave, lwrite_yaverages, &
      lwrite_zaverages, test_nonblocking, &
      lread_oldsnap_nomag, &
      lread_aux, comment_char, ix, iy, iz, iz2, slice_position, &
      xbot_slice, xtop_slice, &
      ybot_slice, ytop_slice, &
      zbot_slice, ztop_slice, &
      bcx, bcy, bcz, &
      fbcx1, fbcx2, fbcy1, fbcy2, fbcz1, fbcz2, &
      fbcx1_2, fbcx2_2, fbcz1_1, fbcz1_2, fbcz2_1, fbcz2_2, Udrift_bc, &
      ttransient, lserial_io, &
      lwrite_aux, &
      old_cdtv, lmaxadvec_sum, save_lastsnap, &
      lwrite_dvar, &
      lpoint, mpoint, npoint, lpoint2, mpoint2, &
      npoint2, ipencil_swap, &
      lpencil_requested_swap, lpencil_diagnos_swap, lpencil_check, &
      lpencil_check_small, lrandom_f_pencil_check, lpencil_check_diagnos_opti, &
      lpencil_init, penc0, lwrite_2d, lbidiagonal_derij, &
      crash_file_dtmin_factor, eps_rkf, &
      eps_stiff, timestep_scaling, lequatory, lequatorz, zequator, &
      lini_t_eq_zero, lav_smallx, xav_max, ldt_paronly, &
      it_timing
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
      use Mpicomm
!
      character (len=*) :: dir
      logical :: exist
!
!  Let root processor check for existence of datadir.in.
!
      if (lroot) then
        inquire(FILE='datadir.in',EXIST=exist)
        if (exist) then
          open(1,FILE='datadir.in',FORM='formatted')
          read(1,'(A)') dir
          close(1)
        endif
      endif
!
!  Tell other processors whether we need to communicate dir (i.e. datadir).
!
      call mpibcast_logical(exist, 1)
!
!  Let root processor communicate dir (i.e. datadir) to all other processors.
!
      if (exist) call mpibcast_char(dir, len(dir))
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
      character (len=120) :: tdir
      character (len=10) :: a_format
      logical :: exist
!
!  check for existence of `data/directory_snap'
!
      inquire(FILE=trim(datadir)//'/directory_snap',EXIST=exist)
      if (exist) then
        open(1,FILE=trim(datadir)//'/directory_snap',FORM='formatted')
! NB: the following does not work under IRIX (silently misses reading of
! run parameters):
!        read(1,'(a)') dir
! ..so we do it like this:
        a_format = '(a)'
        read(1,a_format) tdir
        close(1)
        if (len(trim(tdir)) .gt. 0) call parse_shell(tdir,dir)
      endif
      if (lroot.and.ip<20) print*,'get_snapdir: dir=',trim(dir)
!
    endsubroutine get_snapdir
!***********************************************************************
    subroutine read_startpars(print,file)
!
!  read input parameters (done by each processor)
!
!   6-jul-02/axel: in case of error, print sample namelist
!  21-oct-03/tony: moved sample namelist stuff to a separate procedure
!
      integer :: ierr
      logical, optional :: print,file
      character (len=30) :: label='[none]'
!
!  set default to shearing sheet if lshear=.true. (even when Sshear==0.)
!
      if (lshear) bcx(:)='she'
!
! find out if we should open and close the file everytime
! to fix the SGI reading problem
      inquire(FILE='SGIFIX',EXIST=lsgifix)
!
!  open namelist file
!
      open(1,FILE='start.in',FORM='formatted',STATUS='old')
!
!  read through all items that *may* be present
!  in the various modules
!
      label='init_pars'
      read(1,NML=init_pars                 ,ERR=99, IOSTAT=ierr)
!
      call sgi_fix(lsgifix,1,'start.in')
      call read_initial_condition_pars(1,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('initial_condition_pars',ierr)
!
      call sgi_fix(lsgifix,1,'start.in')
      call read_eos_init_pars(1,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('eos_init_pars',ierr)
!
      call sgi_fix(lsgifix,1,'start.in')
      call read_hydro_init_pars(1,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('hydro_init_pars',ierr)
!
      call sgi_fix(lsgifix,1,'start.in')
      call read_density_init_pars(1,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('density_init_pars',ierr)
!
      call sgi_fix(lsgifix,1,'start.in')
      call read_gravity_init_pars(1,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('grav_init_pars',ierr)
!
      call sgi_fix(lsgifix,1,'start.in')
      call read_entropy_init_pars(1,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('entropy_init_pars',ierr)
!
      call sgi_fix(lsgifix,1,'start.in')
      call read_magnetic_init_pars(1,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('magnetic_init_pars',ierr)
!
      call sgi_fix(lsgifix,1,'start.in')
      call read_signal_init_pars(1,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('signal_init_pars',ierr)
!
      call sgi_fix(lsgifix,1,'start.in')
      call read_shear_init_pars(1,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('shear_init_pars',ierr)
      if (ierr/=0) call sample_startpars('testperturb_init_pars',ierr)
!
      call sgi_fix(lsgifix,1,'start.in')
      call read_viscosity_init_pars(1,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('viscosity_init_pars',ierr)
!
      call sgi_fix(lsgifix,1,'start.in')
      call read_shock_init_pars(1,IOSTAT=ierr)
      if (ierr/=0) call sample_startpars('shock_init_pars',ierr)
!
      ! no input parameters for viscosity
      label='[none]'
      close(1)
!
!  Print cvs id from first line.
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
!  Produce a warning when somebody still sets lcylindrical.
!
      if (lcylindrical) then
        if (lroot) then
          print*
          print*,'read_startpars: lcylindrical=T is now outdated'
          print*,'use instead: lcylinder_in_a_box=T'
          print*,'This renaming became necessary with the development of'
          print*,'cylindrical coordinates which led to very similar names'
          print*,'(coord_system="cylindrical_coords")'
          print*
        endif
        call fatal_error('read_startpars','')
      endif
!
!  In case of i/o error: print sample input list.
!
      return
!
99  call sample_startpars(label,ierr)
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
        if (lentropy          ) print*,'&entropy_init_pars         /'
        if (lmagnetic         ) print*,'&magnetic_init_pars        /'
        if (lsignal           ) print*,'&signal_init_pars          /'
        if (lshear            ) print*,'&shear_init_pars           /'
        print*,'------END sample namelist -------'
        print*
        if (present(label))  print*, 'Found error in input namelist "' // trim(label)
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
        call write_gravity_init_pars(unit)
        call write_entropy_init_pars(unit)
        call write_magnetic_init_pars(unit)
        call write_signal_init_pars(unit)
        call write_shear_init_pars(unit)
        call write_viscosity_init_pars(unit)
        call write_shock_init_pars(unit)
!
        if (present(file)) then
          close(unit)
        endif
      endif
!
    endsubroutine print_startpars
!***********************************************************************
    subroutine read_runpars(print,file,annotation)
!
!  Read input parameters.
!
!  14-sep-01/axel: inserted from run.f90
!  31-may-02/wolf: renamed from cread to read_runpars
!   6-jul-02/axel: in case of error, print sample namelist
!  21-oct-03/tony: moved sample namelist stuff to a separate procedure
!
      use Sub, only: parse_bc
      use Slices, only: setup_slices
!
      integer :: ierr
      logical, optional :: print,file
      character (len=*), optional :: annotation
      character (len=30) :: label='[none]'
!
!  Reset some parameters, in particular those where we play tricks with
!  `impossible' values
!
!  set default to shearing sheet if lshear=.true. (even when Sshear==0.)
!
      if (lshear) bcx(:)='she'
!
!  Find out if we should open and close the file everytime
!  to fix the SGI reading problem.
!
      inquire(FILE='SGIFIX',EXIST=lsgifix)
!
!  Open namelist file.
!
      open(1,FILE='run.in',FORM='formatted',STATUS='old')
!
!  Read through all items that *may* be present in the various modules.
!  AB: at some point the sgi_fix stuff should probably be removed (see sgi bug)
!
      label='run_pars'
                         read(1,NML=run_pars              ,ERR=99, IOSTAT=ierr)
!
      call sgi_fix(lsgifix,1,'run.in')
      call read_eos_run_pars(1,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('eos_run_pars',ierr)
!
      call sgi_fix(lsgifix,1,'run.in')
      call read_hydro_run_pars(1,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('hydro_run_pars',ierr)
!
      call sgi_fix(lsgifix,1,'run.in')
      call read_density_run_pars(1,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('density_run_pars',ierr)
!
      call sgi_fix(lsgifix,1,'run.in')
      call read_gravity_run_pars(1,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('grav_run_pars',ierr)
!
      call sgi_fix(lsgifix,1,'run.in')
      call read_entropy_run_pars(1,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('entropy_run_pars',ierr)
!
      call sgi_fix(lsgifix,1,'run.in')
      call read_magnetic_run_pars(1,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('magnetic_run_pars',ierr)

      call read_shear_run_pars(1,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('shear_run_pars',ierr)
!
      call sgi_fix(lsgifix,1,'run.in')
      call read_viscosity_run_pars(1,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('viscosity_run_pars',ierr)
!
      call sgi_fix(lsgifix,1,'run.in')
      call read_shock_run_pars(1,IOSTAT=ierr)
      if (ierr/=0) call sample_runpars('shock_run_pars',ierr)
!
      label='[none]'
      close(1)
!
!  Print cvs id from first line.
!
      if (lroot) call svn_id(cvsid)
!
!  Set debug logical (easier to use than the combination of ip and lroot).
!
      ldebug=lroot.and.(ip<7)
      if (lroot) print*,'ldebug,ip=',ldebug,ip
!
!  Give online feedback if called with the PRINT optional argument.
!  Note: Some compiler's [like Compaq's] code crashes with the more
!  compact `if (present(print) .and. print)'.
!
      if (present(print)) then
        if (print) then
          call print_runpars()
        endif
      endif
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
!  In case of i/o error: print sample input list.
!
      return
!
99    call sample_runpars(label,ierr)
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
        if (ldensity        ) print*,'&density_run_pars         /'
        if (lgrav           ) print*,'&grav_run_pars            /'
        if (lentropy        ) print*,'&entropy_run_pars         /'
        if (lmagnetic       ) print*,'&magnetic_run_pars        /'
        if (lshear          ) print*,'&shear_run_pars           /'
        if (lviscosity      ) print*,'&viscosity_run_pars       /'
        if (lshock          ) print*,'&shock_run_pars           /'
        print*,'------END sample namelist -------'
        print*
        if (present(label))  print*, 'Found error in input namelist "' // trim(label)
        if (present(iostat)) print*, 'iostat = ', iostat
        if (present(iostat).or.present(label)) &
                           print*,  '-- use sample above.'
      endif
!
      call fatal_error('sample_runpars','')
!
    endsubroutine sample_runpars
!***********************************************************************
    subroutine sgi_fix(lfix,lun,file)
!
      logical :: lfix
      integer :: lun
      character (LEN=*) :: file
!
      if (lfix) then
        close (lun)
        open(lun, FILE=file,FORM='formatted')
      endif
!
    endsubroutine sgi_fix
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
        line = read_line_from_file('RELOAD') ! get first line from file RELOAD
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
        write(unit,NML=run_pars             )
!
        call write_eos_run_pars(unit)
        call write_hydro_run_pars(unit)
        call write_gravity_run_pars(unit)
        call write_entropy_run_pars(unit)
        call write_magnetic_run_pars(unit)
        call write_shear_run_pars(unit)
        call write_viscosity_run_pars(unit)
        call write_shock_run_pars(unit)
!
        if (present(file)) then
          close(unit)
        endif
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
      character (len=*), dimension(mvar) :: bc
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
    subroutine wparam ()
!
!  Write startup parameters
!
!  21-jan-02/wolf: coded
!
      logical :: lhydro         = lhydro_var
      logical :: ldensity       = ldensity_var
      logical :: lentropy       = lentropy_var
      logical :: lshock         = lshock_var
      logical :: lmagnetic      = lmagnetic_var
!
      namelist /lphysics/ &
          lhydro,ldensity,lentropy,lmagnetic,&
          lgravz,lshear,lshock,leos
!
!  Write the param.nml file only from root processor.
!  However, for pacx-MPI (grid-style computations across different platforms)
!  we'd need this on each site separately (not done yet).
!  (In that case we'd need to identify separate master-like processors
!  one at each site.)
!
      if (lroot) then
        open(1,FILE=trim(datadir)//'/param.nml',DELIM='apostrophe', &
               STATUS='unknown')
!
!  Write init_pars.
!
        write(1,NML=init_pars)
!
!  Write each namelist separately.
!
        call write_eos_init_pars(1)
        call write_hydro_init_pars(1)
        call write_density_init_pars(1)
        call write_gravity_init_pars(1)
        call write_entropy_init_pars(1)
        call write_magnetic_init_pars(1)
        call write_signal_init_pars(1)
        call write_shear_init_pars(1)
        call write_viscosity_init_pars(1)
        call write_shock_init_pars(1)
        call write_initial_condition_pars(1)
        ! The following parameters need to be communicated to IDL
        write(1,NML=lphysics              )
        close(1)
      endif
!
      call keep_compiler_quiet(lhydro)
      call keep_compiler_quiet(ldensity)
      call keep_compiler_quiet(lentropy)
      call keep_compiler_quiet(lmagnetic)
      call keep_compiler_quiet(lshock)
!
    endsubroutine wparam
!***********************************************************************
    subroutine rparam ()
!
!  Read startup parameters.
!
!  21-jan-02/wolf: coded
!
      open(1,FILE=trim(datadir)//'/param.nml')
      read(1,NML=init_pars)
      call read_eos_init_pars(1)
      call read_hydro_init_pars(1)
      call read_density_init_pars(1)
      call read_gravity_init_pars(1)
      call read_entropy_init_pars(1)
      call read_magnetic_init_pars(1)
      call read_signal_init_pars(1)
      call read_shear_init_pars(1)
      call read_viscosity_init_pars(1)
      call read_shock_init_pars(1)
      call read_initial_condition_pars(1)
      close(1)
!
      if (lroot.and.ip<14) then
        print*, "rho0,gamma=", rho0,gamma
      endif
!
    endsubroutine rparam
!***********************************************************************
    subroutine wparam2 ()
!
!  Write runtime parameters for IDL.
!
!  21-jan-02/wolf: coded
!
      if (lroot) then
        open(1,FILE=trim(datadir)//'/param2.nml',DELIM='apostrophe')
        write(1,NML=run_pars)
        call write_eos_run_pars(1)
        call write_hydro_run_pars(1)
        call write_density_run_pars(1)
        call write_gravity_run_pars(1)
        call write_entropy_run_pars(1)
        call write_shear_run_pars(1)
        call write_viscosity_run_pars(1)
        call write_shock_run_pars(1)
        close(1)
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
!
     if (lroot) then
       open(1,FILE=trim(datadir)//'/pencils.list')
       write(1,*) 'Pencils requested:'
       do i=1,npencils
         if (lpenc_requested(i)) write(1,*) i, pencil_names(i)
       enddo
       write(1,*) ''
       write(1,*) 'Pencils requested for diagnostics:'
       do i=1,npencils
         if (lpenc_diagnos(i)) write(1,*) i, pencil_names(i)
       enddo
       write(1,*) ''
       write(1,*) 'Pencils requested for 2-D diagnostics:'
       do i=1,npencils
         if (lpenc_diagnos2d(i)) write(1,*) i, pencil_names(i)
       enddo
       write(1,*) ''
       write(1,*) 'Pencils requested video output'
       do i=1,npencils
         if (lpenc_video(i)) write(1,*) i, pencil_names(i)
       enddo
       print*, 'write_pencil_info: pencil information written to the file pencils.list'
       close(1)
     endif
!
   endsubroutine write_pencil_info
!***********************************************************************
endmodule Param_IO

