! $Id: param_io.f90,v 1.53 2002-09-24 17:35:17 brandenb Exp $ 

module Param_IO

!
!  IO of init and run parameters. Subroutines here are `at the end of the
!  food chain', i.e. depend on all physics modules plus possibly others.
!  Using this module is also a compact way of referring to all physics
!  modules at once.
!
  use Sub
  use Hydro
  use Forcing
  use Gravity
  use Entropy
  use Magnetic
  use radiation
  use Pscalar
  use Shear
 
  implicit none 

  ! run parameters
  real :: tmax=1e33,awig=1.
  integer :: isave=100,iwig=0,ialive=0

  namelist /init_pars/ &
       cvsid,ip,xyz0,Lxyz,lperi,lwrite_ic,lnowrite,random_gen
  namelist /run_pars/ &
       cvsid,ip,nt,it1,dt,cdt,cdtv,isave,itorder, &
       dsnap,dvid,dtmin,dspect,tmax,iwig,awig,ialive, &
       vel_spec,mag_spec,vec_spec, &
       directory_snap, &
       bcx,bcy,bcz, &
       ttransient
 
  contains

!***********************************************************************
    subroutine read_inipars()
!
!  read input parameters (done by each processor)
!
!   6-jul-02/axel: in case of error, print sample namelist
!
      use Mpicomm, only: stop_it
!
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
                      read(1,NML=init_pars          ,err=99)
      label='hydro_init_pars'
      if (lhydro    ) read(1,NML=hydro_init_pars    ,err=99)
      label='density_init_pars'
      if (ldensity  ) read(1,NML=density_init_pars  ,err=99)
      label='grav_init_pars'
      if (lgrav     ) read(1,NML=grav_init_pars     ,err=99)
      label='entropy_init_pars'
      if (lentropy  ) read(1,NML=entropy_init_pars  ,err=99)
      label='magnetic_init_pars'
      if (lmagnetic ) read(1,NML=magnetic_init_pars ,err=99)
      label='radiation_init_pars'
      if (lradiation) read(1,NML=radiation_init_pars,err=99)
      label='pscalar_init_pars'
      if (lpscalar  ) read(1,NML=pscalar_init_pars  ,err=99)
      label='shear_init_pars'
      if (lshear    ) read(1,NML=shear_init_pars    ,err=99)
      label='[none]'
      close(1)
!
!  output on the console, but only when root processor
!
!  print cvs id from first line
      if(lroot) call cvs_id(cvsid)
!
      if (lroot.and.ip<14) then
                        write(*,NML=init_pars          )
        if (lhydro    ) write(*,NML=hydro_init_pars    )
        if (ldensity  ) write(*,NML=density_init_pars  )
        if (lgrav     ) write(*,NML=grav_init_pars     )
        if (lentropy  ) write(*,NML=entropy_init_pars  )
        if (lmagnetic ) write(*,NML=magnetic_init_pars )
        if (lradiation) write(*,NML=radiation_init_pars)
        if (lpscalar  ) write(*,NML=pscalar_init_pars  )
        if (lshear    ) write(*,NML=shear_init_pars    )
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
!  in case of i/o error: print sample input list
!
      return
99    if (lroot) then
        print*
        print*,'-----BEGIN sample namelist ------'
                        print*,'&init_pars                /'
        if (lhydro    ) print*,'&hydro_init_pars          /'
        if (ldensity  ) print*,'&density_init_pars        /'
        if (lgrav     ) print*,'&grav_init_pars           /'
        if (lentropy  ) print*,'&entropy_init_pars        /'
        if (lmagnetic ) print*,'&magnetic_init_pars       /'
        if (lradiation) print*,'&radiation_init_pars      /'
        if (lpscalar  ) print*,'&pscalar_init_pars        /'
        if (lshear    ) print*,'&shear_init_pars          /'
        print*,'------END sample namelist -------'
        print*
      endif
      call stop_it('found error in input namelist "' // trim(label) &
           // '": use sample above')
!
    endsubroutine read_inipars
!***********************************************************************
    subroutine read_runpars(print)
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
      logical, optional :: print
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
                      read(1,NML=run_pars          ,err=99)
      label='hydro_run_pars'
      if (lhydro    ) read(1,NML=hydro_run_pars    ,err=99)
      label='density_run_pars'
      if (ldensity  ) read(1,NML=density_run_pars  ,err=99)
      label='forcing_run_pars'
      if (lforcing  ) read(1,NML=forcing_run_pars  ,err=99)
      label='grav_run_pars'
      if (lgrav     ) read(1,NML=grav_run_pars     ,err=99)
      label='entropy_run_pars'
      if (lentropy  ) read(1,NML=entropy_run_pars  ,err=99)
      label='magnetic_run_pars'
      if (lmagnetic ) read(1,NML=magnetic_run_pars ,err=99)
      label='radiation_run_pars'
      if (lradiation) read(1,NML=radiation_run_pars,err=99)
      label='pscalar_run_pars'
      if (lpscalar  ) read(1,NML=pscalar_run_pars  ,err=99)
      label='shear_run_pars'
      if (lshear    ) read(1,NML=shear_run_pars    ,err=99)
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
!  set gamma1, cs20, and lnrho0
!  general parameter checks (and possible adjustments)
!
      gamma1=gamma-1.
      cs20=cs0**2
      lnrho0=alog(rho0)
      if(lforcing) call param_check_forcing
!
!  calculate shear flow velocity; if Sshear is not given
!  then Sshear=-qshear*Omega is calculated.
!
      if (lshear) then
        if (Sshear==impossible) Sshear=-qshear*Omega
      endif
!
!  timestep: if dt=0 (ie not initialized), ldt=.true.
!
      ldt = (dt==0.)            ! need to calculate dt dynamically?
      if (lroot .and. ip<14) then
        if (ldt) then
          print*,'timestep based on CFL cond; cdt=',cdt
        else
          print*, 'absolute timestep dt=', dt
        endif
      endif
!
!  in case of i/o error: print sample input list
!
      return
99    if (lroot) then
        print*
        print*,'-----BEGIN sample namelist ------'
                        print*,'&run_pars                /'
        if (lhydro    ) print*,'&hydro_run_pars          /'
        if (ldensity  ) print*,'&density_run_pars        /'
        if (lforcing  ) print*,'&forcing_run_pars        /'
        if (lgrav     ) print*,'&grav_run_pars           /'
        if (lentropy  ) print*,'&entropy_run_pars        /'
        if (lmagnetic ) print*,'&magnetic_run_pars       /'
        if (lradiation) print*,'&radiation_run_pars      /'
        if (lpscalar  ) print*,'&pscalar_run_pars        /'
        if (lshear    ) print*,'&shear_run_pars          /'
        print*,'------END sample namelist -------'
        print*
      endif
      call stop_it('found error in input namelist "' // trim(label) &
           // '": use sample above')
!
    endsubroutine read_runpars
!***********************************************************************
    subroutine print_runpars
!
!  print input parameters
!  14-sep-01/axel: inserted from run.f90
!  31-may-02/wolf: renamed from cprint to print_runpars
!
      use Cdata
!
      if (lroot.and.ip<14) then
                        write(*,NML=run_pars          )
        if (lhydro    ) write(*,NML=hydro_run_pars    )
        if (lforcing  ) write(*,NML=forcing_run_pars  )
        if (lgrav     ) write(*,NML=grav_run_pars     )
        if (lentropy  ) write(*,NML=entropy_run_pars  )
        if (lmagnetic ) write(*,NML=magnetic_run_pars )
        if (lradiation) write(*,NML=radiation_run_pars)
        if (lpscalar  ) write(*,NML=pscalar_run_pars  )
        if (lshear    ) write(*,NML=shear_run_pars    )
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
           lhydro,ldensity,lgravz,lgravr,lentropy,lmagnetic,lradiation,lpscalar,lforcing,lshear
!
      if (lroot) then
        open(1,FILE='tmp/param.nml',DELIM='apostrophe' )
                        write(1,NML=init_pars          )
        if (lhydro    ) write(1,NML=hydro_init_pars    )
        if (ldensity  ) write(1,NML=density_init_pars  )
        ! no input parameters for forcing
        if (lgrav     ) write(1,NML=grav_init_pars     )
        if (lentropy  ) write(1,NML=entropy_init_pars  )
        if (lmagnetic ) write(1,NML=magnetic_init_pars )
        if (lradiation) write(1,NML=radiation_init_pars)
        if (lpscalar  ) write(1,NML=pscalar_init_pars  )
        if (lshear    ) write(1,NML=shear_init_pars    )
        ! The following parameters need to be communicated to IDL
        ! Note: logicals will be written as Fortran integers
                       write(1,NML=lphysics         ) 
      endif
!
    endsubroutine wparam
!***********************************************************************
    subroutine rparam ()
!
!  Read startup parameters
!  21-jan-02/wolf: coded
!  ?Is there a good reason to have this routine in sub.f90?
!  ?How about register.f90, for example?
!
      use Cdata
!
        open(1,FILE='tmp/param.nml')
                        read(1,NML=init_pars          )
        if (lhydro    ) read(1,NML=hydro_init_pars    )
        if (ldensity  ) read(1,NML=density_init_pars  )
        if (lgrav     ) read(1,NML=grav_init_pars     )
        if (lentropy  ) read(1,NML=entropy_init_pars  )
        if (lmagnetic ) read(1,NML=magnetic_init_pars )
        if (lradiation) read(1,NML=radiation_init_pars)
        if (lpscalar  ) read(1,NML=pscalar_init_pars  )
        if (lshear    ) read(1,NML=shear_init_pars    )
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
!  21-jan-02/wolf: coded
!
      use Cdata
!
      if (lroot) then
        open(1,FILE='tmp/param2.nml',DELIM='apostrophe')
                        write(1,NML=run_pars          )
        if (lhydro    ) write(1,NML=hydro_run_pars    )
        if (ldensity  ) write(1,NML=density_run_pars  )
        if (lforcing  ) write(1,NML=forcing_run_pars  )
        if (lgrav     ) write(1,NML=grav_run_pars     )
        if (lentropy  ) write(1,NML=entropy_run_pars  )
        if (lmagnetic ) write(1,NML=magnetic_run_pars )
        if (lradiation) write(1,NML=radiation_run_pars)
        if (lpscalar  ) write(1,NML=pscalar_run_pars  )
        if (lshear    ) write(1,NML=shear_run_pars    )
      endif
!
    endsubroutine wparam2
!***********************************************************************

endmodule Param_IO

