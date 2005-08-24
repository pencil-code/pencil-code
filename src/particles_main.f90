! $Id: particles_main.f90,v 1.4 2005-08-24 13:06:30 ajohan Exp $
!
!  This module contains all the main structure needed for particles.
!
module Particles_main

  use Cdata
  use Particles_cdata
  use Particles_sub
  use Particles
  use Particles_radius
  use Messages

  implicit none

  include 'particles_main.h'

  real, dimension (mpar_loc,mpvar) :: fp, dfp

  contains

!***********************************************************************
    subroutine particles_register_modules()
!
!  Register particle modules.
!
!  07-jan-05/anders: coded
!
      call register_particles()
      call register_particles_radius()
!
    endsubroutine particles_register_modules
!***********************************************************************
    subroutine particles_rprint_list(lreset)
!
!  Read names of diagnostic particle variables to print out during run.
!
!  07-jan-05/anders: coded
!
      logical :: lreset
!
      if (lroot) open(3, file=trim(datadir)//'/index.pro', &
          STATUS='old', POSITION='append')
      call rprint_particles(lreset,LWRITE=lroot)
      call rprint_particles_radius(lreset,LWRITE=lroot)
      if (lroot) close(3)
!
    endsubroutine particles_rprint_list
!***********************************************************************
    subroutine particles_initialize_modules(lstarting)
!
!  Initialize particle modules.
!
!  07-jan-05/anders: coded
!
      logical :: lstarting
!
!  Check if there is enough total space allocated for particles.
!
      if (mpar_loc-ncpus*npar_mig<npar/ncpus) then
        if (lroot) then
          print*, 'particles_initialize_modules: '// &
          'total number of particle slots available at the processors '// &
          'is smaller than the number of particles!'
          print*, 'particles_initialize_modules: npar/ncpus=', npar/ncpus
          print*, 'particles_initialize_modules: mpar_loc-ncpus*npar_mig=', &
              mpar_loc-ncpus*npar_mig
        endif
        call fatal_error('particles_initialize_modules','')
      endif
!
      call initialize_particles(lstarting)
      if (lparticles_radius) call initialize_particles_radius(lstarting)
!
    endsubroutine particles_initialize_modules
!***********************************************************************
    subroutine particles_init(f)
!
!  Set up initial conditios for particle modules.
!
!  07-jan-05/anders: coded
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      intent (in) :: f
!
      call init_particles(f,fp)
      if (lparticles_radius) call init_particles_radius(f,fp)
!
    endsubroutine particles_init
!***********************************************************************
    subroutine particles_read_snapshot(filename)
!
!  Read particle snapshot from file.
!
!  07-jan-05/anders: coded
!
      character (len=*) :: filename
!
      call input_particles(filename,fp,npar_loc,ipar)
!
    endsubroutine particles_read_snapshot
!***********************************************************************
    subroutine particles_write_snapshot(chsnap,msnap,enum,flist)
!
!  Write particle snapshot to file.
!
!  07-jan-05/anders: coded
!
      integer :: msnap
      logical :: enum
      character (len=*) :: chsnap,flist
      optional :: flist
!
      logical :: lsnap
!
      if (present(flist)) then
        call wsnap_particles(chsnap,fp,msnap,enum,lsnap,dsnap_par_minor, &
            npar_loc,ipar,flist)
      else
        call wsnap_particles(chsnap,fp,msnap,enum,lsnap,dsnap_par_minor, &
            npar_loc,ipar)
      endif
!
    endsubroutine particles_write_snapshot
!***********************************************************************
    subroutine particles_write_pdim(filename)
!   
!  Write npar and mpvar to file.
!
!  09-jan-05/anders: coded
!
      character (len=*) :: filename
!
      open(1,file=filename)
        write(1,'(2i7)') npar, mpvar
      close(1)
!
    endsubroutine particles_write_pdim
!***********************************************************************
    subroutine particles_timestep_first()
!
!  Setup dfp in the beginning of each itsub.
!
!  07-jan-05/anders: coded
!
      integer :: k
!
      if (itsub==1) then
        do k=1,npar_loc
          dfp(k,:)=0.
        enddo
      else
        do k=1,npar_loc
          dfp(k,:)=alpha(itsub)*dfp(k,:)
        enddo
      endif
!
    endsubroutine particles_timestep_first
!***********************************************************************
    subroutine particles_timestep_second()
!
!  Time evolution of particle variables.
!
!  07-jan-05/anders: coded
!
      integer :: k
!
      do k=1,npar_loc
        fp(k,:) = fp(k,:) + dt_beta(itsub)*dfp(k,:)
      enddo
!
    endsubroutine particles_timestep_second
!***********************************************************************
    subroutine particles_pde(f,df)
!
!  Dynamical evolution of particle variables.
!
!  07-jan-05/anders: coded
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
!
      intent (in) :: f
!
      call boundconds_particles(fp,npar_loc,ipar,dfp=dfp)
      call dxxp_dt(f,fp,dfp)
      call dvvp_dt(f,df,fp,dfp)
      if (lparticles_radius) call dap_dt(f,df,fp,dfp)
!
    endsubroutine particles_pde
!***********************************************************************
    subroutine read_particles_init_pars_wrapper(unit,iostat)
!    
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      call read_particles_init_pars(unit,iostat)
      if (lparticles_radius) call read_particles_radius_init_pars(unit,iostat)
!
    endsubroutine read_particles_init_pars_wrapper
!***********************************************************************
    subroutine write_particles_init_pars_wrapper(unit)
!    
      integer, intent (in) :: unit
!
      call write_particles_init_pars(unit)
      if (lparticles_radius) call write_particles_radius_init_pars(unit)
!
    endsubroutine write_particles_init_pars_wrapper
!***********************************************************************
    subroutine read_particles_run_pars_wrapper(unit,iostat)
!    
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      call read_particles_run_pars(unit,iostat)
      if (lparticles_radius) call read_particles_radius_run_pars(unit,iostat)
!
    endsubroutine read_particles_run_pars_wrapper
!***********************************************************************
    subroutine write_particles_run_pars_wrapper(unit)
!    
      integer, intent (in) :: unit
!
      call write_particles_run_pars(unit)
      if (lparticles_radius) call write_particles_radius_run_pars(unit)
!
    endsubroutine write_particles_run_pars_wrapper
!***********************************************************************


endmodule Particles_main
