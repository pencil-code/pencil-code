! $Id: param_io.f90,v 1.1 2002-05-31 20:43:45 dobler Exp $

module Param_IO

!
!  IO of init and run parameters. Subroutines here are `at the end of the
!  food chain', i.e. depend on all physics modules plus possibly others.
!  Using this module is also a compact way of referring to all physics
!  modules at once.
!
  use Hydro
  use Forcing
  use Gravity
  use Entropy
  use Magnetic
!
  implicit none 
!
  namelist /init_pars/ &
       ip,xyz0,Lxyz,lperi
  namelist /run_pars/ &
       ip,nt,it1,dt,cdtv,isave,itorder, &
       dsnap,dvid,dforce,dtmin, &
       tinit,tdamp,dampu,dampuext,rdamp,wdamp, &
       bcx,bcy,bcz
!
  contains

!***********************************************************************
    subroutine read_inipars()
!
!  read input parameters (done by each processor)
!
      open(1,FILE='start.in',FORM='formatted')
                     read(1,NML=init_pars         )
      if (lhydro)    read(1,NML=hydro_init_pars   )
      ! forcing needs no init parameters
      if (lgravz)    read(1,NML=grav_z_init_pars  )
      if (lentropy)  read(1,NML=entropy_init_pars )
      if (lmagnetic) read(1,NML=magnetic_init_pars)
      close(1)
!
!  output on the console, but only when root processor
!
      if (lroot) then
                       write(*,NML=init_pars         )
        if (lhydro   ) write(*,NML=hydro_init_pars   )
        ! forcing needs no init parameters
        if (lgravz)    write(*,NML=grav_z_init_pars  )
        if (lentropy ) write(*,NML=entropy_init_pars )
        if (lmagnetic) write(*,NML=magnetic_init_pars)
      endif
!
    endsubroutine read_inipars
!***********************************************************************
    subroutine read_runpars(print)
!
!  read input parameters
!  14-sep-01/axel: inserted from run.f90
!  31-may-02/wolf: renamed from cread to read_runpars
!
      use Cdata
      use Sub, only: parse_bc
!
      logical, optional :: print
!
      open(1,file='run.in',form='formatted')
                     read(1,NML=run_pars         )
      if (lhydro   ) read(1,NML=hydro_run_pars )
      if (lforcing ) read(1,NML=forcing_run_pars )
      if (lgravz   ) read(1,NML=grav_z_run_pars  )
      if (lentropy ) read(1,NML=entropy_run_pars )
      if (lmagnetic) read(1,NML=magnetic_run_pars)
      close(1)
      cs20=cs0**2 !(goes into cdata module)
      ss0 = (alog(cs20) - gamma1*alog(rho0)-alog(gamma))/gamma   !!AB: this looks like it belongs to entropy
!
!  Write data to file for IDL
!
      call wparam2()
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
      ix=min(ix,mx); iy=min(iy,my); iz=min(iz,mz)
      ix=max(ix,0); iy=max(iy,0); iz=max(iz,0)
!
!  parse boundary conditions; compound conditions of the form `a:s' allow
!  to have different variables at the lower and upper boundaries
!
      call parse_bc(bcx,bcx1,bcx2)
      call parse_bc(bcy,bcy1,bcy2)
      call parse_bc(bcz,bcz1,bcz2)
      if (lroot) then
        print*, 'bcx1,bcx2= ', bcx1," : ",bcx2
        print*, 'bcy1,bcy2= ', bcy1," : ",bcy2
        print*, 'bcz1,bcz2= ', bcz1," : ",bcz2
      endif
!
!  timestep
!
      ldt=dt.lt.0.
      if (ldt) cdt=abs(dt)
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
      if (lroot) then
                       write(*,NML=run_pars         )
        if (lhydro   ) write(*,NML=hydro_run_pars   )
        if (lforcing ) write(*,NML=forcing_run_pars )
        if (lgravz   ) write(*,NML=grav_z_run_pars  )
        if (lentropy ) write(*,NML=entropy_run_pars )
        if (lmagnetic) write(*,NML=magnetic_run_pars)
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
      use Mpicomm
!
      namelist /lphysics/ &
           lhydro,lgravz,lgravr,lentropy,lmagnetic,lforcing
!
      if (lroot) then
        open(1,FILE='tmp/param.nml',DELIM='apostrophe')
                       write(1,NML=init_pars         )
        if (lhydro   ) write(1,NML=hydro_init_pars   )
        ! no input parameters for forcing
        if (lgravz   ) write(1,NML=grav_z_init_pars  )
        if (lentropy ) write(1,NML=entropy_init_pars )
        if (lmagnetic) write(1,NML=magnetic_init_pars)
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
!     use Mpicomm
! ?AB Mpicomm is no longer used, because lroot is now in cdata
!
        open(1,FILE='tmp/param.nml')
                       read(1,NML=init_pars         )
        if (lhydro   ) read(1,NML=hydro_init_pars   )
        if (lforcing ) read(1,NML=forcing_init_pars )
        if (lgravz   ) read(1,NML=grav_z_init_pars  )
        if (lentropy ) read(1,NML=entropy_init_pars )
        if (lmagnetic) read(1,NML=magnetic_init_pars)
!
      if (lroot) then
        print*, "Lx,Ly,Lz=", Lx,Ly,Lz
        print*, "rho0,gamma,gamma1=", rho0,gamma,gamma1
      endif
!
!  read the print parameter list
!
!     call rprint_list
!  
    endsubroutine rparam
!***********************************************************************
    subroutine wparam2 ()
!
!  Write runtime parameters for IDL
!  21-jan-02/wolf: coded
!
      use Cdata
      use Mpicomm
!
      if (lroot) then
        open(1,FILE='tmp/param2.nml',DELIM='apostrophe')
                       write(1,NML=run_pars         )
        if (lhydro   ) write(1,NML=hydro_run_pars   )
        if (lforcing ) write(1,NML=forcing_run_pars )
        if (lgravz   ) write(1,NML=grav_z_run_pars  )
        if (lentropy ) write(1,NML=entropy_run_pars )
        if (lmagnetic) write(1,NML=magnetic_run_pars)
      endif
!
    endsubroutine wparam2
!***********************************************************************

endmodule Param_IO

