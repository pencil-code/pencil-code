!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!   write_default_namelists.f90   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! $Id: write_default_namelists.f90,v 1.1 2006-06-22 16:35:18 dobler Exp $

!***********************************************************************
  program Write_default_namelists
!
!  Write startup and run parameters to files to allow a Perl script to
!  validate start.in and run.in.
!
!  21-jul-05/wolf: coded
!
    use Param_IO
!
    implicit none
!
    call print_startpars(FILE='default_start_pars.nml')
    call print_runpars  (FILE='default_run_pars.nml'  )

  endprogram Write_default_namelists

!!! End of file write_defaul_namelists.f90
