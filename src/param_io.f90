! $Id: param_io.f90,v 1.35 2002-07-04 10:10:55 nilshau Exp $ 

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
  use Shear
 
  implicit none 

  ! run parameters
  real :: tmax=1e33
  integer :: iwig=0

  namelist /init_pars/ &
       cvsid,ip,xyz0,Lxyz,lperi,lwrite_ic,lnowrite
  namelist /run_pars/ &
       cvsid,ip,nt,it1,dt,cdt,cdtv,isave,itorder, &
       dsnap,dvid,dtmin,tmax,iwig, &
       bcx,bcy,bcz
 
  contains

!***********************************************************************
    subroutine read_inipars()
!
!  read input parameters (done by each processor)
!
!  open namelist file
!
      open(1,FILE='start.in',FORM='formatted')
!
!  read through all items that *may* be present
!  in the various modules
!
                     read(1,NML=init_pars         )
      if (lhydro)    read(1,NML=hydro_init_pars   )
      if (ldensity)  read(1,NML=density_init_pars )
      ! forcing needs no init parameters
      if (lgrav)     read(1,NML=grav_init_pars    )
      if (lentropy)  read(1,NML=entropy_init_pars )
      if (lmagnetic) read(1,NML=magnetic_init_pars)
      if (lshear) read(1,NML=shear_init_pars)
      close(1)
!
!  output on the console, but only when root processor
!
!  print cvs id from first line
      if(lroot) call cvs_id(cvsid)
!
      if (lroot.and.ip<14) then
                       write(*,NML=init_pars         )
        if (lhydro   ) write(*,NML=hydro_init_pars   )
        if (ldensity ) write(*,NML=density_init_pars   )
        ! forcing needs no init parameters
        if (lgrav)     write(*,NML=grav_init_pars  )
        if (lentropy ) write(*,NML=entropy_init_pars )
        if (lmagnetic) write(*,NML=magnetic_init_pars)
        if (lshear) read(1,NML=shear_init_pars)
      endif
!
!  set gamma1, cs20, and lnrho0
!
      gamma1=gamma-1.
      cs20=cs0**2
      lnrho0=alog(rho0)
!
    endsubroutine read_inipars
!***********************************************************************
    subroutine read_runpars(print)
!
!  read input parameters
!
!  14-sep-01/axel: inserted from run.f90
!  31-may-02/wolf: renamed from cread to read_runpars
!
      use Sub, only: parse_bc
!
!--   character(len=80) Id
      logical, optional :: print
      integer :: i
!
!  set default values
!
      do i=1,mvar; bcx(i)='p'; enddo
      do i=1,mvar; bcy(i)='p'; enddo
      do i=1,mvar; bcz(i)='p'; enddo
!
!  open namelist file
!
      open(1,file='run.in',form='formatted')
!
!  read through all items that *may* be present
!  in the various modules
!
                     read(1,NML=run_pars         )
      if (lhydro   ) read(1,NML=hydro_run_pars )
      if (ldensity ) read(1,NML=density_run_pars )
      if (lforcing ) read(1,NML=forcing_run_pars )
      if (lgrav    ) read(1,NML=grav_run_pars  )
      if (lentropy ) read(1,NML=entropy_run_pars )
      if (lmagnetic) read(1,NML=magnetic_run_pars)
      if (lshear)    read(1,NML=shear_run_pars)
      close(1)
!
!  print cvs id from first line
!
      if(lroot) call cvs_id(cvsid)
!
!  set debug logical (easier to use than the combination of ip and lroot)
!
      ldebug=lroot.and.(ip<7)
      print*,'ldebug,ip=',ldebug,ip
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
!
      gamma1=gamma-1.
      cs20=cs0**2
      lnrho0=alog(rho0)
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
                       write(*,NML=run_pars         )
        if (lhydro   ) write(*,NML=hydro_run_pars   )
        if (lforcing ) write(*,NML=forcing_run_pars )
        if (lgrav    ) write(*,NML=grav_run_pars  )
        if (lentropy ) write(*,NML=entropy_run_pars )
        if (lmagnetic) write(*,NML=magnetic_run_pars)
        if (lshear)    write(*,NML=shear_run_pars)
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
           lhydro,ldensity,lgravz,lgravr,lentropy,lmagnetic,lforcing,lshear
!
      if (lroot) then
        open(1,FILE='tmp/param.nml',DELIM='apostrophe')
                       write(1,NML=init_pars         )
        if (lhydro   ) write(1,NML=hydro_init_pars   )
        if (ldensity ) write(1,NML=density_init_pars   )
        ! no input parameters for forcing
        if (lgrav    ) write(1,NML=grav_init_pars  )
        if (lentropy ) write(1,NML=entropy_init_pars )
        if (lmagnetic) write(1,NML=magnetic_init_pars)
        if (lshear)    write(1,NML=shear_init_pars)
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
                       read(1,NML=init_pars         )
        if (lhydro   ) read(1,NML=hydro_init_pars   )
        if (ldensity ) read(1,NML=density_init_pars   )
!??     if (lforcing ) read(1,NML=forcing_init_pars )
        if (lgrav    ) read(1,NML=grav_init_pars  )
        if (lentropy ) read(1,NML=entropy_init_pars )
        if (lmagnetic) read(1,NML=magnetic_init_pars)
        if (lshear) read(1,NML=shear_init_pars)
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
                       write(1,NML=run_pars         )
        if (lhydro   ) write(1,NML=hydro_run_pars   )
        if (ldensity ) write(1,NML=density_run_pars   )
        if (lforcing ) write(1,NML=forcing_run_pars )
        if (lgrav    ) write(1,NML=grav_run_pars  )
        if (lentropy ) write(1,NML=entropy_run_pars )
        if (lmagnetic) write(1,NML=magnetic_run_pars)
        if (lshear)    write(1,NML=shear_run_pars)
      endif
!
    endsubroutine wparam2
!***********************************************************************

endmodule Param_IO

