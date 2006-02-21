! $Id: param_io.f90,v 1.222 2006-02-21 15:33:09 nbabkovs Exp $ 

module Param_IO

!
!  IO of init and run parameters. Subroutines here are `at the end of the
!  food chain', i.e. depend on all physics modules plus possibly others.
!  Using this module is also a compact way of referring to all physics
!  modules at once.
!
  use Cdata
  use Sub
  use General
  use Hydro
  use Entropy
  use Density
  use Magnetic
  use Testfield
  use Pscalar
  use Chiral
  use CosmicRay
  use CosmicRayFlux
  use Dustvelocity
  use Dustdensity
  use Radiation
  use EquationOfState
  use Forcing
  use Gravity
  use Interstellar
  use Shear
  use Timeavg
  use Viscosity
  use Special
  use Particles_main
  use Shock
  use Messages
  use Planet
 
  implicit none 

  private

  public :: get_datadir, get_snapdir
  public :: read_startpars, print_startpars
  public :: read_runpars,   print_runpars
  public :: rparam, wparam, wparam2, write_pencil_info

  !
  ! The following fixes namelist problems withi MIPSpro 7.3.1.3m 
  ! under IRIX -- at least for the moment
  !
  character (len=labellen) :: mips_is_buggy='system'

  namelist /init_pars/ &
       cvsid,ip,xyz0,xyz1,Lxyz,lperi,lshift_origin, &
       coord_system,lequidist,coeff_grid,zeta_grid0,grid_func,xyz_star, &
       xyz_step,xi_step_frac,xi_step_width, &
       lwrite_ic,lnowrite, &
       unit_system,unit_length,unit_velocity,unit_density,unit_temperature, &
       random_gen,nfilter,lserial_io, &
       lread_oldsnap,lread_oldsnap_nomag,lread_oldsnap_nopscalar, &
       lwrite_aux,lcalc_cp,pretend_lnTT, &
       lprocz_slowest, lcopysnapshots_exp, &
       bcx,bcy,bcz,r_int,r_ext,r_ref, &
       mu0,force_lower_bound,force_upper_bound, &
       fbcx1,fbcx2,fbcy1,fbcy2,fbcz1,fbcz2, fbcz1_1, fbcz1_2, fbcz2_1, fbcz2_2, &
       xyz_step,xi_step_frac,xi_step_width, &
       lcylindrical,init_loops, H_disk, L_disk, R_star
  namelist /run_pars/ &
       cvsid,ip,nt,it1,dt,cdt,ddt,cdtv,cdts,cdtr,isave,itorder, &
       dsnap,d2davg,dvid,dtmin,dspec,tmax,iwig,awig,ialive, max_walltime, &
       dtmax, &
       vel_spec,mag_spec,uxj_spec,vec_spec,ou_spec,ab_spec,fft_switch, &
       ro_spec,ss_spec,cc_spec,cr_spec,isaveglobal, &
       rhocc_pdf,cc_pdf,lncc_pdf,gcc_pdf,lngcc_pdf, &
       kinflow,random_gen, &
       lrmwig_rho,lrmwig_full,lrmwig_xyaverage, &
       lwrite_yaverages,lwrite_zaverages,lwrite_phiaverages,test_nonblocking, &
       lread_oldsnap_nomag,lread_oldsnap_nopscalar, &
       ix,iy,iz,iz2,slice_position, &
       bcx,bcy,bcz,r_int,r_ext, &
       lfreeze_varint,lfreeze_varext,rfreeze_int,rfreeze_ext,&
       wfreezeint,wfreezeext, &
       fbcx1,fbcx2,fbcy1,fbcy2,fbcz1,fbcz2, fbcz1_1, fbcz1_2, fbcz2_1, fbcz2_2,&
       ttransient,tavg,idx_tavg,lserial_io,nr_directions, &
       lsfu,lsfb,lsfz1,lsfz2,lsfflux,lpdfu,lpdfb,lpdfz1,lpdfz2,oned, &
       lwrite_aux,onedall,lcalc_cp,pretend_lnTT,old_cdtv,lmaxadvec_sum, &
       save_lastsnap, &
       force_lower_bound,force_upper_bound,twod, &
       border_frac,border_frac_x,border_frac_y,border_frac_z, &
       lpoint,mpoint,npoint, &
       lrescaling, lcylindrical, &
       ipencil_swap,lpencil_requested_swap,lpencil_diagnos_swap, &
       lpencil_check,lpencil_check_diagnos_opti,lpencil_init, H_disk, L_disk, R_star
  contains

!***********************************************************************
    subroutine get_datadir(dir)
!
!  Overwrite datadir from datadir.in, if that exists
!
!   2-oct-02/wolf: coded
!  25-oct-02/axel: default is taken from cdata.f90 where it's defined
!
      character (len=*) :: dir
      logical :: exist
!
!  check for existence of datadir.in
!
      inquire(FILE='datadir.in',EXIST=exist)
      if (exist) then
        open(1,FILE='datadir.in',FORM='formatted')
        read(1,*) dir
        close(1)
      endif
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
      if(lroot.and.ip<20) print*,'get_snapdir: dir=',trim(dir)
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

      call sgi_fix(lsgifix,1,'start.in')
      call read_eos_init_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_startpars('eos_init_pars',ierr)

      call sgi_fix(lsgifix,1,'start.in')
      call read_hydro_init_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_startpars('hydro_init_pars',ierr)

      call sgi_fix(lsgifix,1,'start.in')
      call read_density_init_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_startpars('density_init_pars',ierr)

      call sgi_fix(lsgifix,1,'start.in')
      call read_forcing_init_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_startpars('forcing_init_pars',ierr)
      ! no input parameters for forcing

      call sgi_fix(lsgifix,1,'start.in')
      call read_gravity_init_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_startpars('grav_init_pars',ierr)

      call sgi_fix(lsgifix,1,'start.in')
      call read_entropy_init_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_startpars('entropy_init_pars',ierr)

      call sgi_fix(lsgifix,1,'start.in')
      call read_magnetic_init_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_startpars('magnetic_init_pars',ierr)

      call sgi_fix(lsgifix,1,'start.in')
      call read_testfield_init_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_startpars('testfield_init_pars',ierr)

      call sgi_fix(lsgifix,1,'start.in')
      call read_radiation_init_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_startpars('radiation_init_pars',ierr)

      call sgi_fix(lsgifix,1,'start.in')
      call read_pscalar_init_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_startpars('pscalar_init_pars',ierr)

      call sgi_fix(lsgifix,1,'start.in')
      call read_chiral_init_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_startpars('chiral_init_pars',ierr)

      call sgi_fix(lsgifix,1,'start.in')
      call read_dustvelocity_init_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_startpars('dustvelocity_init_pars',ierr)

      call sgi_fix(lsgifix,1,'start.in')
      call read_dustdensity_init_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_startpars('dustdensity_init_pars',ierr)

      call sgi_fix(lsgifix,1,'start.in')
      call read_cosmicray_init_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_startpars('cosmicray_init_pars',ierr)

      call sgi_fix(lsgifix,1,'start.in')
      call read_cosmicrayflux_init_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_startpars('cosmicrayflux_init_pars',ierr)

      call sgi_fix(lsgifix,1,'start.in')
      call read_interstellar_init_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_startpars('interstellar_init_pars',ierr)

      call sgi_fix(lsgifix,1,'start.in')
      call read_shear_init_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_startpars('shear_init_pars',ierr)

      call sgi_fix(lsgifix,1,'start.in')
      call read_viscosity_init_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_startpars('viscosity_init_pars',ierr)

      call sgi_fix(lsgifix,1,'start.in')
      call read_special_init_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_startpars('special_init_pars',ierr)

      call sgi_fix(lsgifix,1,'start.in')
      call read_particles_init_pars_wrap(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_startpars('particles_init_pars',ierr)

      call sgi_fix(lsgifix,1,'start.in')
      call read_shock_init_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_startpars('shock_init_pars',ierr)

      call sgi_fix(lsgifix,1,'start.in')
      call read_planet_init_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_startpars('planet_init_pars',ierr)

      ! no input parameters for viscosity
      label='[none]'
      close(1)
!
!  print cvs id from first line
!
      if(lroot) call cvs_id(cvsid)
!
!  Give online feedback if called with the PRINT optional argument
!  Note: Some compiler's [like Compaq's] code crashes with the more
!  compact `if (present(print) .and. print)'
!
      if (present(print)) then
        if (print) then
          call print_startpars()
        endif
      endif
!
!  Write parameters to log file
!
      if (present(file)) then
        if (file) then
          call print_startpars(FILE=trim(datadir)//'/params.log')
        endif
      endif
!
!  parse boundary conditions; compound conditions of the form `a:s' allow
!  to have different variables at the lower and upper boundaries
!
      call parse_bc(bcx,bcx1,bcx2)
      call parse_bc(bcy,bcy1,bcy2)
      call parse_bc(bcz,bcz1,bcz2)
      if (lroot.and.ip<14) then
        print*, 'bcx1,bcx2= ', bcx1," : ",bcx2
        print*, 'bcy1,bcy2= ', bcy1," : ",bcy2
        print*, 'bcz1,bcz2= ', bcz1," : ",bcz2
      endif
      call check_consistency_of_lperi('read_startpars')
!
!  in case of i/o error: print sample input list
!
      return
99  call sample_startpars(label,ierr)
    endsubroutine read_startpars
!***********************************************************************
    subroutine sample_startpars(label,iostat)
      use Mpicomm, only: stop_it
      
      character (len=*), optional :: label
      integer, optional :: iostat

      if (lroot) then
        print*
        print*,'-----BEGIN sample namelist ------'
                             print*,'&init_pars                /'
                             print*,'&eos_init_pars            /'
        if (lhydro        )  print*,'&hydro_init_pars          /'
        if (ldensity      )  print*,'&density_init_pars        /'
        ! no input parameters for forcing
        if (lgrav         )  print*,'&grav_init_pars           /'
        if (lentropy      )  print*,'&entropy_init_pars        /'
        if (lmagnetic     )  print*,'&magnetic_init_pars       /'
        if (ltestfield    )  print*,'&testfield_init_pars      /'
        if (lradiation    )  print*,'&radiation_init_pars      /'
        if (lpscalar      )  print*,'&pscalar_init_pars        /'
        if (lchiral       )  print*,'&chiral_init_pars         /'
        if (ldustvelocity )  print*,'&dustvelocity_init_pars   /'
        if (ldustdensity  )  print*,'&dustdensity_init_pars    /'
        if (lcosmicray    )  print*,'&cosmicray_init_pars      /'
        if (lcosmicrayflux)  print*,'&cosmicrayflux_init_pars  /'
        if (linterstellar )  print*,'&interstellar_init_pars   /'
        ! no input parameters for planet
        if (lshear        )  print*,'&shear_init_pars          /'
        if (lspecial      )  print*,'&special_init_pars        /'
        if (lparticles       ) print*,'&particles_init_pars       /'
        if (lparticles_radius) print*,'&particles_radius_init_pars/'
        if (lparticles_number) print*,'&particles_number_init_pars/'
        !if (lshock       ) print*,'&shock_init_pars          /'
        ! no input parameters for viscosity
        print*,'------END sample namelist -------'
        print*
        if (present(label))  print*, 'Found error in input namelist "' // trim(label)
        if (present(iostat)) print*, 'iostat = ', iostat
        if (present(iostat).or.present(label)) &
                           print*,  '-- use sample above.'
      endif
      call stop_it('')
!
    endsubroutine sample_startpars
!***********************************************************************
    subroutine print_startpars(file)
!
!  print input parameters
!  4-oct02/wolf: adapted
!
      use Cdata
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

        call write_eos_init_pars(unit)
        call write_hydro_init_pars(unit)
        call write_density_init_pars(unit)
        call write_forcing_init_pars(unit)
        call write_gravity_init_pars(unit)
        call write_entropy_init_pars(unit)
        call write_magnetic_init_pars(unit)
        call write_testfield_init_pars(unit)
        call write_radiation_init_pars(unit)
        call write_pscalar_init_pars(unit)
        call write_chiral_init_pars(unit)
        call write_dustvelocity_init_pars(unit)
        call write_dustdensity_init_pars(unit)
        call write_cosmicray_init_pars(unit)
        call write_cosmicrayflux_init_pars(unit)
        call write_interstellar_init_pars(unit)
        call write_shear_init_pars(unit)
        call write_viscosity_init_pars(unit)
        call write_special_init_pars(unit)
        call write_particles_init_pars_wrap(unit)
        call write_shock_init_pars(unit)
        call write_planet_init_pars(unit)
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
!  read input parameters
!
!  14-sep-01/axel: inserted from run.f90
!  31-may-02/wolf: renamed from cread to read_runpars
!   6-jul-02/axel: in case of error, print sample namelist
!  21-oct-03/tony: moved sample namelist stuff to a separate procedure
!
      use Mpicomm, only: stop_it
      use Sub, only: parse_bc
      use Dustvelocity, only: copy_bcs_dust
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
      !  entropy
      !
      Kbot=impossible
      hcond0=impossible
!
! find out if we should open and close the file everytime
! to fix the SGI reading problem
      inquire(FILE='SGIFIX',EXIST=lsgifix)
!
!  open namelist file
!
      open(1,FILE='run.in',FORM='formatted',STATUS='old')
!
!  read through all items that *may* be present
!  in the various modules
!AB: at some point the sgi_fix stuff should probably be removed (see sgi bug)
!
      label='run_pars'
                         read(1,NML=run_pars              ,ERR=99, IOSTAT=ierr)

      call sgi_fix(lsgifix,1,'run.in')
      call read_eos_run_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_runpars('eos_run_pars',ierr)

      call sgi_fix(lsgifix,1,'run.in')
      call read_hydro_run_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_runpars('hydro_run_pars',ierr)

      call sgi_fix(lsgifix,1,'run.in')
      call read_density_run_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_runpars('density_run_pars',ierr)

      call sgi_fix(lsgifix,1,'run.in')
      call read_forcing_run_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_runpars('forcing_run_pars',ierr)

      call sgi_fix(lsgifix,1,'run.in')
      call read_gravity_run_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_runpars('grav_run_pars',ierr)

      call sgi_fix(lsgifix,1,'run.in')
      call read_entropy_run_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_runpars('entropy_run_pars',ierr)

      call sgi_fix(lsgifix,1,'run.in')
      call read_magnetic_run_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_runpars('magnetic_run_pars',ierr)

      call sgi_fix(lsgifix,1,'run.in')
      call read_testfield_run_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_runpars('testfield_run_pars',ierr)

      call sgi_fix(lsgifix,1,'run.in')
      call read_radiation_run_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_runpars('radiation_run_pars',ierr)

      call sgi_fix(lsgifix,1,'run.in')
      call read_pscalar_run_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_runpars('pscalar_run_pars',ierr)

      call sgi_fix(lsgifix,1,'run.in')
      call read_chiral_run_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_runpars('chiral_run_pars',ierr)

      call sgi_fix(lsgifix,1,'run.in')
      call read_dustvelocity_run_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_runpars('dustvelocity_run_pars',ierr)

      call sgi_fix(lsgifix,1,'run.in')
      call read_dustdensity_run_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_runpars('dustdensity_run_pars',ierr)

      call sgi_fix(lsgifix,1,'run.in')
      call read_cosmicray_run_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_runpars('cosmicray_run_pars',ierr)

      call sgi_fix(lsgifix,1,'run.in')
      call read_cosmicrayflux_run_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_runpars('cosmicrayflux_run_pars',ierr)

      call sgi_fix(lsgifix,1,'run.in')
      call read_interstellar_run_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_runpars('interstellar_run_pars',ierr)

      call sgi_fix(lsgifix,1,'run.in')
      call read_shear_run_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_runpars('shear_run_pars',ierr)

      call sgi_fix(lsgifix,1,'run.in')
      call read_viscosity_run_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_runpars('viscosity_run_pars',ierr)

      call sgi_fix(lsgifix,1,'run.in')
      !call read_special_run_pars(1,IOSTAT=ierr)
      call read_special_run_pars(1)
      if (ierr.ne.0) call sample_runpars('special_run_pars',ierr)

      call sgi_fix(lsgifix,1,'run.in')
      call read_particles_run_pars_wrap(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_runpars('particles_run_pars',ierr)

      call sgi_fix(lsgifix,1,'run.in')
      call read_shock_run_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_runpars('shock_run_pars',ierr)

      call sgi_fix(lsgifix,1,'run.in')
      call read_planet_run_pars(1,IOSTAT=ierr)
      if (ierr.ne.0) call sample_runpars('planet_run_pars',ierr)

      label='[none]'
      close(1)
!
!  Copy boundary conditions on first dust species to all others
!
      call copy_bcs_dust 
!
!  print cvs id from first line
!
      if(lroot) call cvs_id(cvsid)
!
!  set debug logical (easier to use than the combination of ip and lroot)
!
      ldebug=lroot.and.(ip<7)
      if (lroot) print*,'ldebug,ip=',ldebug,ip
!      random_gen=random_gen_tmp
!
!  Give online feedback if called with the PRINT optional argument
!  Note: Some compiler's [like Compaq's] code crashes with the more
!  compact `if (present(print) .and. print)' 
!
      if (present(print)) then
        if (print) then
          call print_runpars()
        endif
      endif
!
!  Write parameters to log file
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
!  set slice position. The default for slice_position is 'p' for periphery,
!  although setting ix, iy, iz, iz2 by hand will overwrite this.
!  If slice_position is not 'p', then ix, iy, iz, iz2 are overwritten.
!
      if (slice_position=='p' .or. slice_position=='S') then
        if (ix<0)  ix=l1
        if (iy<0)  iy=m1
        if (iz<0)  iz=n1
        if (iz2<0) iz2=n2
        lwrite_slice_xy2=(ipz==nprocz-1)
        lwrite_slice_xy=(ipz==0)
        lwrite_slice_xz=(ipy==0)
        lwrite_slice_yz=.true.
!
!  slice position when the first meshpoint in z is the equator (sphere)
!  For one z-processor, iz remains n1, but iz2 is set to the middle.
!
      elseif (slice_position=='m') then
        if (ix<0)  ix=(l1+l2)/2
        if (iy<0)  iy=m1
        if (iz<0)  iz=n1
        if (iz2<0) iz2=n2
        if(nprocy==1) then; iy=(m1+m2)/2; endif
        if(nprocz==1) then; iz=(n1+n2)/2; iz2=(iz+n2)/2; endif
        lwrite_slice_xy2=(ipz==nprocz/2)
        lwrite_slice_xy=(ipz==nprocz/2)
        lwrite_slice_xz=(ipy==nprocy/2)
        lwrite_slice_yz=.true.
!
!  slice position when the first meshpoint in z is the equator (sphere)
!  For one z-processor, iz remains n1, but iz2 is set to the middle.
!
      elseif (slice_position=='e') then
        if (ix<0)  ix=(l1+l2)/2
        if (iy<0)  iy=m1
        if (iz<0)  iz=n1
        if (iz2<0) iz2=n2
        if(nprocy==1) then; iy=(m1+m2)/2; endif
        if(nprocz==1) then; iz2=(iz+n2)/2; endif
        lwrite_slice_xy2=(ipz==nprocz/4)
        lwrite_slice_xy=(ipz==0)
        lwrite_slice_xz=(ipy==nprocy/2)
        lwrite_slice_yz=.true.
!
!  slice position when the first meshpoint in z is the equator (sphere)
!  For one z-processor, iz remains n1, but iz2 is set to the middle.
!
      elseif (slice_position=='c') then
        if (ix<0)  ix=(l1+l2)/2
        if (iy<0)  iy=m1
        if (iz<0)  iz=n1
        if (iz2<0) iz2=n2
        lwrite_slice_xy2=(ipz==nprocz-1)
        lwrite_slice_xy=(ipz==0)
        lwrite_slice_xz=(ipy==0)
        lwrite_slice_yz=.true.
!
!  periphery of the box, but the other way around
!  For one z-processor, iz remains n1, but iz2 is set to the middle.
!
      elseif (slice_position=='q') then
        if (ix<0)  ix=l2
        if (iy<0)  iy=m2
        if (iz<0)  iz=n2
        if (iz2<0) iz2=n1
        lwrite_slice_xy2=(ipz==0)
        lwrite_slice_xy=(ipz==nprocz-1)
        lwrite_slice_xz=(ipy==nprocy-1)
        lwrite_slice_yz=.true.
      else
        if (lroot) print*, &
             "READ_RUN_PARS: No such value for slice_position: '" &
             // slice_position //"'"
        call stop_it("")
      endif
!
!  write slice position to a file (for convenient post-processing)
!
      if (lroot) then
        open(1,file=trim(datadir)//'/slice_position.dat',STATUS='unknown')
        write(1,'(a)') slice_position
        close(1)
      endif
!  
!  make sure ix,iy,iz,iz2 are not outside the boundaries
!
      ix=min(ix,l2); iy=min(iy,m2); iz=min(iz,n2); iz2=min(iz2,n2)
      ix=max(ix,l1); iy=max(iy,m1); iz=max(iz,n1); iz2=max(iz2,n1)
      if (lroot) write(*,'(1x,a,4i4)') &
        'read_runpars: slice position (video files) ix,iy,iz,iz2 =',ix,iy,iz,iz2
!
!  parse boundary conditions; compound conditions of the form `a:s' allow
!  to have different variables at the lower and upper boundaries
!
      call parse_bc(bcx,bcx1,bcx2)
      call parse_bc(bcy,bcy1,bcy2)
      call parse_bc(bcz,bcz1,bcz2)
      if (lroot.and.ip<14) then
        print*, 'bcx1,bcx2= ', bcx1," : ",bcx2
        print*, 'bcy1,bcy2= ', bcy1," : ",bcy2
        print*, 'bcz1,bcz2= ', bcz1," : ",bcz2
      endif
      call check_consistency_of_lperi('read_runpars')
!
!  in case of i/o error: print sample input list
!
      return
99    call sample_runpars(label,ierr)
    endsubroutine read_runpars
!***********************************************************************
    subroutine sample_runpars(label,iostat)
      use Mpicomm, only: stop_it
      
      character (len=*), optional :: label
      integer, optional :: iostat

      if (lroot) then
        print*
        print*,'-----BEGIN sample namelist ------'
                            print*,'&run_pars                 /'
                            print*,'&eos_run_pars             /'
        if (lhydro        ) print*,'&hydro_run_pars           /'
        if (ldensity      ) print*,'&density_run_pars         /'
        if (lforcing      ) print*,'&forcing_run_pars         /'
        if (lgrav         ) print*,'&grav_run_pars            /'
        if (lentropy      ) print*,'&entropy_run_pars         /'
        if (lmagnetic     ) print*,'&magnetic_run_pars        /'
        if (ltestfield    ) print*,'&testfield_run_pars       /'
        if (lradiation    ) print*,'&radiation_run_pars       /'
        if (lpscalar      ) print*,'&pscalar_run_pars         /'
        if (lchiral       ) print*,'&chiral_run_pars          /'
        if (ldustvelocity ) print*,'&dustvelocity_run_pars    /'
        if (ldustdensity  ) print*,'&dustdensity_run_pars     /'
        if (lcosmicray    ) print*,'&cosmicray_run_pars       /'
        if (lcosmicrayflux) print*,'&cosmicrayflux_run_pars   /'
        if (linterstellar ) print*,'&interstellar_run_pars    /'
        if (lshear        ) print*,'&shear_run_pars           /'
        if (lviscosity    ) print*,'&viscosity_run_pars       /'
        if (lspecial      ) print*,'&special_run_pars         /'
        if (lparticles       ) print*,'&particles_run_pars       /'
        if (lparticles_radius) print*,'&particles_radius_run_pars       /'
        if (lparticles_number) print*,'&particles_number_run_pars       /'
        if (lshock        ) print*,'&shock_run_pars           /'
        if (lplanet       ) print*,'&planet_run_pars          /'
        print*,'------END sample namelist -------'
        print*
        if (present(label))  print*, 'Found error in input namelist "' // trim(label)
        if (present(iostat)) print*, 'iostat = ', iostat
        if (present(iostat).or.present(label)) &
                           print*,  '-- use sample above.'
      endif
      call stop_it('')
!
    endsubroutine sample_runpars
!***********************************************************************
    subroutine sgi_fix(lfix,lun,file)
!
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
!  print input parameters
!  14-sep-01/axel: inserted from run.f90
!  31-may-02/wolf: renamed from cprint to print_runpars
!   4-oct-02/wolf: added log file stuff
!
      use Cdata
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
          ! Add comment from `RELOAD' and time
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
        call write_forcing_run_pars(unit)
        call write_gravity_run_pars(unit)
        call write_entropy_run_pars(unit)
        call write_magnetic_run_pars(unit)
        call write_testfield_run_pars(unit)
        call write_radiation_run_pars(unit)
        call write_pscalar_run_pars(unit)
        call write_chiral_run_pars(unit)
        call write_dustvelocity_run_pars(unit)
        call write_dustdensity_run_pars(unit)
        call write_cosmicray_run_pars(unit)
        call write_cosmicrayflux_run_pars(unit)
        call write_interstellar_run_pars(unit)
        call write_shear_run_pars(unit)
        call write_viscosity_run_pars(unit)
        call write_special_run_pars(unit)
        call write_particles_run_pars_wrap(unit)
        call write_shock_run_pars(unit)
        call write_planet_run_pars(unit)
!
        if (present(file)) then
          close(unit)
        endif
!
      endif
!
    endsubroutine print_runpars
!***********************************************************************
    subroutine check_consistency_of_lperi (label)
!
!  check consistency of lperi
!
!  18-jul-03/axel: coded
!
      use Cdata
!
      character (len=*) :: label
      logical :: lwarning=.true.
      integer :: j
!
!  identifier
!
      if(lroot.and.ip<5) print*,'check_consistency_of_lperi: called from ',label
!
!  make the warnings less dramatic looking, if we are only in start
!  and exit this routine altogether if, in addition, ip > 13.
!
      if(label=='read_startpars'.and.ip>13) return
      if(label=='read_startpars') lwarning=.false.
!
!  check x direction
!
      j=1
      if(any(bcx=='p'.or. bcx=='she').and..not.lperi(j).or.&
         any(bcx/='p'.and.bcx/='she').and.lperi(j)) &
           call warning_lperi(lwarning,bcx,lperi,j)
!
!  check y direction
!
      j=2
      if(any(bcy=='p').and..not.lperi(j).or.&
         any(bcy/='p').and.lperi(j)) &
           call warning_lperi(lwarning,bcy,lperi,j)
!
!  check z direction
!
      j=3
      if(any(bcz=='p').and..not.lperi(j).or.&
         any(bcz/='p').and.lperi(j)) &
           call warning_lperi(lwarning,bcz,lperi,j)
!
!  print final warning
!  make the warnings less dramatic looking, if we are only in start
!
      if (lroot .and. (.not. lwarning)) then
        if(label=='read_startpars') then
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
!  print consistency warning of lperi
!
!  18-jul-03/axel: coded
!
      character (len=*), dimension(mvar) :: bc
      logical, dimension(3) :: lperi
      logical :: lwarning
      integer :: j
!
      if (lroot) then
        if(lwarning) then
          print*
          print*,'------------------------------------------------------'
          print*,'W A R N I N G'
          lwarning=.false.
        endif
!
        print*,'warning_lperi: inconsistency, j=',j
        print*,'lperi(j)=',lperi(j)
        print*,'bc=',bc
      endif
!
    endsubroutine warning_lperi
!***********************************************************************
    subroutine wparam ()
!
!  Write startup parameters
!  21-jan-02/wolf: coded
!
      use Cdata, only: lmagnetic_var,lpscalar,lradiation, &
           lforcing,lgravz,lgravr,lshear, &
           ldustvelocity,ldustdensity,lradiation_fld,  &
           leos_ionization,leos_fixed_ionization,lvisc_hyper,lchiral, &
           leos,lspecial, ltestfield, &
           lhydro_var, lentropy_var, ldensity_var, lshock_var, &
           lcosmicray_var, lcosmicrayflux_var, linterstellar_var, &
           lplanet_var,datadir
!
      logical :: lhydro         = lhydro_var
      logical :: ldensity       = ldensity_var
      logical :: lentropy       = lentropy_var
      logical :: lshock         = lshock_var
      logical :: lmagnetic      = lmagnetic_var
      logical :: linterstellar  = linterstellar_var
      logical :: lcosmicray     = lcosmicray_var
      logical :: lcosmicrayflux = lcosmicrayflux_var
      logical :: lplanet        = lplanet_var

      namelist /lphysics/ &
           lhydro,ldensity,lentropy,lmagnetic,ltestfield,lpscalar,lradiation, &
           lforcing,lgravz,lgravr,lshear,linterstellar,lcosmicray, &
           lcosmicrayflux,ldustvelocity,ldustdensity,lshock,lradiation_fld, &
           leos_ionization,leos_fixed_ionization,lvisc_hyper,lchiral, &
           lplanet,leos
!
!  Write this file from each processor; needed for pacx-MPI (grid-style
!  computations across different platforms), where the data/ directories
!  are different for the different computers involved.
!    Need to do this in a cleverer way if it implies performance
!  penalties, but this is not likely, since the file written is quite
!  small.
!      if (lroot) then
        open(1,FILE=trim(datadir)//'/param.nml',DELIM='apostrophe' )

        write(1,NML=init_pars             )

        call write_eos_init_pars(1)
        call write_hydro_init_pars(1)
        call write_density_init_pars(1)
        call write_forcing_init_pars(1)
        call write_gravity_init_pars(1)
        call write_entropy_init_pars(1)
        call write_magnetic_init_pars(1)
        call write_testfield_init_pars(1)
        call write_radiation_init_pars(1)
        call write_pscalar_init_pars(1)
        call write_chiral_init_pars(1)
        call write_dustvelocity_init_pars(1)
        call write_dustdensity_init_pars(1)
        call write_cosmicray_init_pars(1)
        call write_cosmicrayflux_init_pars(1)
        call write_interstellar_init_pars(1)
        call write_shear_init_pars(1)
        call write_viscosity_init_pars(1)
        call write_special_init_pars(1)
        call write_particles_init_pars_wrap(1)
        call write_shock_init_pars(1)
        call write_planet_init_pars(1)
        ! The following parameters need to be communicated to IDL
        write(1,NML=lphysics              ) 
        close(1)
!      endif
!
    endsubroutine wparam
!***********************************************************************
    subroutine rparam ()
!
!  Read startup parameters
!
!  21-jan-02/wolf: coded
!
      use Cdata
!
        open(1,FILE=trim(datadir)//'/param.nml')
        read(1,NML=init_pars             )
        call read_eos_init_pars(1)
        call read_hydro_init_pars(1)
        call read_density_init_pars(1)
        call read_forcing_init_pars(1)
        call read_gravity_init_pars(1)
        call read_entropy_init_pars(1)
        call read_magnetic_init_pars(1)
        call read_testfield_init_pars(1)
        call read_radiation_init_pars(1)
        call read_pscalar_init_pars(1)
        call read_chiral_init_pars(1)
        call read_dustvelocity_init_pars(1)
        call read_dustdensity_init_pars(1)
        call read_cosmicray_init_pars(1)
        call read_cosmicrayflux_init_pars(1)
        call read_interstellar_init_pars(1)
        call read_shear_init_pars(1)
        call read_viscosity_init_pars(1)
        call read_special_init_pars(1)
        call read_particles_init_pars_wrap(1)
        call read_shock_init_pars(1)
        call read_planet_init_pars(1)
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
!  Write runtime parameters for IDL
!
!  21-jan-02/wolf: coded
!
      use Cdata
!
      if (lroot) then
        open(1,FILE=trim(datadir)//'/param2.nml',DELIM='apostrophe')
                           write(1,NML=run_pars             )
        call write_eos_run_pars(1)
        call write_hydro_run_pars(1)
        call write_density_run_pars(1)
        call write_forcing_run_pars(1)
        call write_gravity_run_pars(1)
        call write_entropy_run_pars(1)
        call write_magnetic_run_pars(1)
        call write_testfield_run_pars(1)
        call write_radiation_run_pars(1)
        call write_pscalar_run_pars(1)
        call write_chiral_run_pars(1)
        call write_dustvelocity_run_pars(1)
        call write_dustdensity_run_pars(1)
        call write_cosmicray_run_pars(1)
        call write_cosmicrayflux_run_pars(1)
        call write_interstellar_run_pars(1)
        call write_shear_run_pars(1)
        call write_viscosity_run_pars(1)
        call write_special_run_pars(1)
        call write_particles_run_pars_wrap(1)
        call write_shock_run_pars(1)
        call write_planet_run_pars(1)
        close(1)
      endif
!
    endsubroutine wparam2
!***********************************************************************
    subroutine write_pencil_info()
!
!  Write information about requested and diagnostic pencils
!
     use Cparam
!
     integer :: i
!
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
!
   endsubroutine write_pencil_info
!   
endmodule Param_IO

