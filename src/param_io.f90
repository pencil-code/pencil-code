! $Id: param_io.f90,v 1.94 2003-02-22 13:39:04 brandenb Exp $ 

module Param_IO

!
!  IO of init and run parameters. Subroutines here are `at the end of the
!  food chain', i.e. depend on all physics modules plus possibly others.
!  Using this module is also a compact way of referring to all physics
!  modules at once.
!
  use Sub
  use General
  use Hydro
  use Entropy
  use Magnetic
  use Pscalar
  use Radiation
  use Ionization
  use Forcing
  use Gravity
  use Interstellar
  use Shear
  use Timeavg
  use Viscosity
 
  implicit none 

  ! start parameters (units)
  ! http://www.astro.wisc.edu/~dolan/constants.html

  real, parameter :: hbar_cgs=1.0545726663d-27 ! [erg*s]
  real, parameter :: k_B_cgs=1.38065812d-16    ! [erg/K]
  real, parameter :: m_p_cgs=1.672623110d-24   ! [g]
  real, parameter :: m_e_cgs=9.109389754d-28   ! [g]
  real, parameter :: eV_cgs=1.602177250d-12    ! [erg]

  ! run parameters
  real :: tmax=1e33,awig=1.
  integer :: isave=100,iwig=0,ialive=0
  logical :: lrmwig_rho=.false.,lrmwig_full=.false.,lrmwig_xyaverage=.false.
  logical :: lwrite_zaverages=.false.
  !
  ! The following fixes namelist problems withi MIPSpro 7.3.1.3m 
  ! under IRIX -- at least for the moment
  !
  character (len=labellen) :: mips_is_buggy='system'

  namelist /init_pars/ &
       cvsid,ip,xyz0,xyz1,Lxyz,lperi,lwrite_ic,lnowrite, &
       unit_system,unit_length,unit_velocity,unit_density, &
       random_gen,lserial_io
  namelist /run_pars/ &
       cvsid,ip,nt,it1,dt,cdt,cdtv,isave,itorder, &
       dsnap,dvid,dtmin,dspec,tmax,iwig,awig,ialive, &
       vel_spec,mag_spec,vec_spec,ou_spec,ab_spec,fft_switch, &
       random_gen, &
       lrmwig_rho,lrmwig_full,lrmwig_xyaverage, &
       lwrite_zaverages,test_nonblocking, &
       ix,iy,iz,iz2, &
       bcx,bcy,bcz, &
       ttransient,tavg,idx_tavg,lserial_io,nr_directions, &
       lsfu,lsfb,lsfz1,lsfz2,lpdfu,lpdfb,lpdfz1,lpdfz2,oned
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
        read(1,a_format) dir
        close(1)
      endif
      if(lroot.and.ip<6) print*,'get_snapdir: dir=',trim(dir)
!
    endsubroutine get_snapdir
!***********************************************************************
    subroutine read_startpars(print,file)
!
!  read input parameters (done by each processor)
!
!   6-jul-02/axel: in case of error, print sample namelist
!
      use Mpicomm, only: stop_it
!
      real :: unit_mass,unit_energy,unit_time
      integer :: ierr
      logical, optional :: print,file
      character (len=30) :: label='[none]'
!
!  open namelist file
!
      open(1,FILE='start.in',FORM='formatted')
!
!  read through all items that *may* be present
!  in the various modules
!
      label='init_pars'
                      read(1,NML=init_pars                 ,ERR=99, IOSTAT=ierr)
      label='hydro_init_pars'
      if (lhydro    ) read(1,NML=hydro_init_pars           ,ERR=99, IOSTAT=ierr)
      label='density_init_pars'
      if (ldensity     ) read(1,NML=density_init_pars      ,ERR=99, IOSTAT=ierr)
      ! no input parameters for forcing
      label='grav_init_pars'
      if (lgrav        ) read(1,NML=grav_init_pars         ,ERR=99, IOSTAT=ierr)
      label='entropy_init_pars'
      if (lentropy     ) read(1,NML=entropy_init_pars      ,ERR=99, IOSTAT=ierr)
      label='magnetic_init_pars'
      if (lmagnetic    ) read(1,NML=magnetic_init_pars     ,ERR=99, IOSTAT=ierr)
      label='radiation_init_pars'
      if (lradiation   ) read(1,NML=radiation_init_pars    ,ERR=99, IOSTAT=ierr)
      label='ionization_init_pars'
      if (lionization  ) read(1,NML=ionization_init_pars   ,ERR=99, IOSTAT=ierr)
      label='pscalar_init_pars'
      if (lpscalar     ) read(1,NML=pscalar_init_pars      ,ERR=99, IOSTAT=ierr)
      label='interstellar_init_pars'
      if (linterstellar) read(1,NML=interstellar_init_pars ,ERR=99, IOSTAT=ierr)
      label='shear_init_pars'
      if (lshear       ) read(1,NML=shear_init_pars        ,ERR=99, IOSTAT=ierr)
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
!  evaluate physical units used in ionization
!  (and perhaps later in the interstellar module)
!
      unit_mass=unit_density*unit_length**3
      unit_energy=unit_mass*unit_velocity**2
      unit_time=unit_length/unit_velocity
!
      if (unit_system=='cgs') then
        hbar=hbar_cgs/(unit_energy*unit_time)
        k_B=k_B_cgs/(unit_energy/unit_temperature)
        m_p=m_p_cgs/unit_mass
        m_e=m_e_cgs/unit_mass
        eV=eV_cgs/unit_energy
      elseif (unit_system=='SI') then
        k_B=1e-7*k_B_cgs/(unit_energy/unit_temperature)
        m_p=m_p_cgs*1e-3/unit_mass
        m_e=m_e_cgs*1e-3/unit_mass
        eV=eV_cgs*1e-7/unit_energy
      endif
!
!  set gamma1, cs20, and lnrho0
!
      gamma1=gamma-1.
      cs20=cs0**2
      lnrho0=alog(rho0)
!
!  calculate shear flow velocity; if Sshear is not given
!  then Sshear=-qshear*Omega is calculated.
!
      if (lshear) then
        if (Sshear==impossible) Sshear=-qshear*Omega
      endif
!
      return
!
!  in case of i/o error: print sample input list
!
99    if (lroot) then
        print*
        print*,'-----BEGIN sample namelist ------'
                           print*,'&init_pars                /'
        if (lhydro       ) print*,'&hydro_init_pars          /'
        if (ldensity     ) print*,'&density_init_pars        /'
        ! no input parameters for forcing
        if (lgrav        ) print*,'&grav_init_pars           /'
        if (lentropy     ) print*,'&entropy_init_pars        /'
        if (lmagnetic    ) print*,'&magnetic_init_pars       /'
        if (lradiation   ) print*,'&radiation_init_pars      /'
        if (lionization  ) print*,'&ionization_init_pars     /'
        if (lpscalar     ) print*,'&pscalar_init_pars        /'
        if (linterstellar) print*,'&interstellar_init_pars   /'
        if (lshear       ) print*,'&shear_init_pars          /'
        ! no input parameters for viscosity
       print*,'------END sample namelist -------'
        print*
      endif
      if (lroot) then
        print*, 'Found error in input namelist "' // trim(label)
        print*, 'iostat = ', ierr
        print*,  '-- use sample above.'
      endif
      call stop_it('')
!
    endsubroutine read_startpars
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
               '# -------------------------------------------------------------'
          write(unit,'(A,A)') ' # ', 'Initializing'
          write(unit,'(A,A)') ' # Date: ', trim(date)
          write(unit,*) '# t=', t
        endif
!
                        write(unit,NML=init_pars          )
        if (lhydro       ) write(unit,NML=hydro_init_pars       )
        if (ldensity     ) write(unit,NML=density_init_pars     )
        ! no input parameters for forcing
        if (lgrav        ) write(unit,NML=grav_init_pars        )
        if (lentropy     ) write(unit,NML=entropy_init_pars     )
        if (lmagnetic    ) write(unit,NML=magnetic_init_pars    )
        if (lradiation   ) write(unit,NML=radiation_init_pars   )
        if (lionization  ) write(unit,NML=ionization_init_pars  )
        if (lpscalar     ) write(unit,NML=pscalar_init_pars     )
        if (linterstellar) write(unit,NML=interstellar_init_pars)
        if (lshear       ) write(unit,NML=shear_init_pars       )
        ! no input parameters for viscosity
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
!
      use Sub, only: parse_bc
      use Mpicomm, only: stop_it
!
      integer :: ierr
      logical, optional :: print,file
      character (len=*), optional :: annotation
      character (len=30) :: label='[none]'
!
!  set default values
!
      bcx(1:nvar)='p'
      bcy(1:nvar)='p'
      bcz(1:nvar)='p'
!
!  set default to shearing sheet if lshear=.true.
!  AB: (even when Sshear==0.)
!
      !! if (lshear .AND. Sshear/=0) bcx(1:nvar)='she'
      if (lshear) bcx(1:nvar)='she'
!
!  open namelist file
!
      open(1,file='run.in',form='formatted')
!
!  read through all items that *may* be present
!  in the various modules
!
      label='run_pars'
                      read(1,NML=run_pars          ,ERR=99, IOSTAT=ierr)
      label='hydro_run_pars'
      if (lhydro       ) read(1,NML=hydro_run_pars        ,ERR=99, IOSTAT=ierr)
      label='density_run_pars'
      if (ldensity     ) read(1,NML=density_run_pars      ,ERR=99, IOSTAT=ierr)
      label='forcing_run_pars'
      if (lforcing     ) read(1,NML=forcing_run_pars      ,ERR=99, IOSTAT=ierr)
      label='grav_run_pars'
      if (lgrav        ) read(1,NML=grav_run_pars         ,ERR=99, IOSTAT=ierr)
      label='entropy_run_pars'
      if (lentropy     ) read(1,NML=entropy_run_pars      ,ERR=99, IOSTAT=ierr)
      label='magnetic_run_pars'
      if (lmagnetic    ) read(1,NML=magnetic_run_pars     ,ERR=99, IOSTAT=ierr)
      label='radiation_run_pars'
      if (lradiation   ) read(1,NML=radiation_run_pars    ,ERR=99, IOSTAT=ierr)
      label='pscalar_run_pars'
      if (lpscalar     ) read(1,NML=pscalar_run_pars      ,ERR=99, IOSTAT=ierr)
      label='interstellar_run_pars'
      if (linterstellar) read(1,NML=interstellar_run_pars ,ERR=99, IOSTAT=ierr)
      label='shear_run_pars'
      if (lshear       ) read(1,NML=shear_run_pars        ,ERR=99, IOSTAT=ierr)
      label='viscosity_run_pars'
      if (lviscosity   ) read(1,NML=viscosity_run_pars    ,ERR=99, IOSTAT=ierr)
      label='[none]'
      close(1)
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
!  make sure ix,iy,iz are not outside the boundaries
!
      ix=min(ix,l2); iy=min(iy,m2); iz=min(iz,n2)
      ix=max(ix,l1); iy=max(iy,m1); iz=max(iz,n1)
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
!
!  in case of i/o error: print sample input list
!
      return
99    if (lroot) then
        print*
        print*,'-----BEGIN sample namelist ------'
                        print*,'&run_pars                /'
        if (lhydro       ) print*,'&hydro_run_pars          /'
        if (ldensity     ) print*,'&density_run_pars        /'
        if (lforcing     ) print*,'&forcing_run_pars        /'
        if (lgrav        ) print*,'&grav_run_pars           /'
        if (lentropy     ) print*,'&entropy_run_pars        /'
        if (lmagnetic    ) print*,'&magnetic_run_pars       /'
        if (lradiation   ) print*,'&radiation_run_pars      /'
        if (lionization  ) print*,'&ionization_run_pars     /'
        if (lpscalar     ) print*,'&pscalar_run_pars        /'
        if (linterstellar) print*,'&interstellar_run_pars   /'
        if (lshear       ) print*,'&shear_run_pars          /'
        if (lviscosity   ) print*,'&viscosity_run_pars      /'
        print*,'------END sample namelist -------'
        print*
      endif
      if (lroot) then
        print*, 'Found error in input namelist "' // trim(label)
        print*, 'iostat = ', ierr
        print*,  '-- use sample above.'
      endif
      call stop_it('')
!
    endsubroutine read_runpars
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
               '# -------------------------------------------------------------'
          !
          ! Add comment from `RELOAD' and time
          !
          write(unit,'(A,A)') ' # ', trim(line)
          write(unit,'(A,A)') ' # Date: ', trim(date)
          write(unit,*) '# t=', t
        endif
!
                        write(unit,NML=run_pars          )
        if (lhydro       ) write(unit,NML=hydro_run_pars       )
        if (lforcing     ) write(unit,NML=forcing_run_pars     )
        if (lgrav        ) write(unit,NML=grav_run_pars        )
        if (lentropy     ) write(unit,NML=entropy_run_pars     )
        if (lmagnetic    ) write(unit,NML=magnetic_run_pars    )
        if (lradiation   ) write(unit,NML=radiation_run_pars   )
        if (lionization  ) write(unit,NML=ionization_run_pars  )
        if (lpscalar     ) write(unit,NML=pscalar_run_pars     )
        if (linterstellar) write(unit,NML=interstellar_run_pars)
        if (lshear       ) write(unit,NML=shear_run_pars       )
        if (lviscosity   ) write(unit,NML=viscosity_run_pars   )
!
        if (present(file)) then
          close(unit)
        endif

      endif
!
    endsubroutine print_runpars
!***********************************************************************
    subroutine wparam ()
!
!  Write startup parameters
!  21-jan-02/wolf: coded
!
      use Cdata
!
      namelist /lphysics/ &
           lhydro,ldensity,lentropy,lmagnetic,lpscalar,lradiation, &
           lforcing,lgravz,lgravr,lshear,linterstellar,lionization
!
!  Write this file from each processor; needed for pacx-MPI (grid-style
!  computations across different platforms), where the data/ directories
!  are different for the different computers involved.
!    Need to do this in a cleverer way if it implies performance
!  penalties, but this is not likely, since the file written is quite
!  small.
!      if (lroot) then
        open(1,FILE=trim(datadir)//'/param.nml',DELIM='apostrophe' )
                        write(1,NML=init_pars          )
        if (lhydro       ) write(1,NML=hydro_init_pars       )
        if (ldensity     ) write(1,NML=density_init_pars     )
        ! no input parameters for forcing
        if (lgrav        ) write(1,NML=grav_init_pars        )
        if (lentropy     ) write(1,NML=entropy_init_pars     )
        if (lmagnetic    ) write(1,NML=magnetic_init_pars    )
        if (lradiation   ) write(1,NML=radiation_init_pars   )
        if (lionization  ) write(1,NML=ionization_init_pars  )
        if (lpscalar     ) write(1,NML=pscalar_init_pars     )
        if (linterstellar) write(1,NML=interstellar_init_pars)
        if (lshear       ) write(1,NML=shear_init_pars       )
        ! no input parameters for viscosity
        ! The following parameters need to be communicated to IDL
        ! Note: logicals will be written as Fortran integers
                       write(1,NML=lphysics         ) 
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
                        read(1,NML=init_pars          )
        if (lhydro       ) read(1,NML=hydro_init_pars       )
        if (ldensity     ) read(1,NML=density_init_pars     )
        if (lgrav        ) read(1,NML=grav_init_pars        )
        if (lentropy     ) read(1,NML=entropy_init_pars     )
        if (lmagnetic    ) read(1,NML=magnetic_init_pars    )
        if (lradiation   ) read(1,NML=radiation_init_pars   )
        if (lionization  ) read(1,NML=ionization_init_pars  )
        if (lpscalar     ) read(1,NML=pscalar_init_pars     )
        if (linterstellar) read(1,NML=interstellar_init_pars)
        if (lshear       ) read(1,NML=shear_init_pars       )
        ! no input parameters for viscosity        
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
                        write(1,NML=run_pars          )
        if (lhydro       ) write(1,NML=hydro_run_pars       )
        if (ldensity     ) write(1,NML=density_run_pars     )
        if (lforcing     ) write(1,NML=forcing_run_pars     )
        if (lgrav        ) write(1,NML=grav_run_pars        )
        if (lentropy     ) write(1,NML=entropy_run_pars     )
        if (lmagnetic    ) write(1,NML=magnetic_run_pars    )
        if (lradiation   ) write(1,NML=radiation_run_pars   )
        if (lionization  ) write(1,NML=ionization_run_pars  )
        if (lpscalar     ) write(1,NML=pscalar_run_pars     )
        if (linterstellar) write(1,NML=interstellar_run_pars)
        if (lshear       ) write(1,NML=shear_run_pars       )
        if (lviscosity   ) write(1,NML=viscosity_run_pars   )
      endif
!
    endsubroutine wparam2
!***********************************************************************

endmodule Param_IO

