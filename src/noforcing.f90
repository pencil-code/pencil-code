! $Id: noforcing.f90,v 1.19 2006-08-23 16:53:32 mee Exp $

module Forcing

!! Dummy module for hydro/mhd without forcing

  use Cdata
  use General
  use Messages

  implicit none

  include 'forcing.h'

  !namelist /forcing_init_pars/ dummy
  !namelist /forcing_run_pars/  dummy

  ! other variables (needs to be consistent with reset list below)
  integer :: idiag_rufm=0

  contains

!***********************************************************************
    subroutine register_forcing()
!
!  add forcing in timestep()
!  11-may-2002/wolf: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_forcing called twice')
      first = .false.
!
      lforcing = .false.
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: noforcing.f90,v 1.19 2006-08-23 16:53:32 mee Exp $")
!
    endsubroutine register_forcing
!***********************************************************************
    subroutine initialize_forcing(lstarting)
!
!  initialize random number generator in processor-dependent fashion
!  see comments in start.f90 for details
!
      use Cdata
!
      logical :: lstarting
!
      if(NO_WARN) print*,'lstarting=',lstarting !(to keep compiler quiet)
    endsubroutine initialize_forcing
!***********************************************************************
    subroutine addforce(f)
!
      use Cdata
!
!  add forcing in timestep()
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if(ip==1) print*,f !(to remove compiler warnings)
    endsubroutine addforce
!***********************************************************************
    subroutine read_forcing_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
                                                                                                   
      if (present(iostat) .and. (NO_WARN)) print*,iostat
      if (NO_WARN) print*,unit
                                                                                                   
    endsubroutine read_forcing_init_pars
!***********************************************************************
    subroutine write_forcing_init_pars(unit)
      integer, intent(in) :: unit
                                                                                                   
      if (NO_WARN) print*,unit
                                                                                                   
    endsubroutine write_forcing_init_pars
!***********************************************************************
    subroutine read_forcing_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
                                                                                                   
      if (present(iostat) .and. (NO_WARN)) print*,iostat
      if (NO_WARN) print*,unit
                                                                                                   
    endsubroutine read_forcing_run_pars
!***********************************************************************
    subroutine write_forcing_run_pars(unit)
      integer, intent(in) :: unit
                                                                                                   
      if (NO_WARN) print*,unit
    endsubroutine write_forcing_run_pars
!***********************************************************************
    subroutine input_persistent_forcing(id,lun,done)
!
!  Read in the stored time of the next SNI
!
      integer :: id,lun
      logical :: done
!
      if (NO_WARN) print*,id,lun,done
!
    endsubroutine input_persistent_forcing
!***********************************************************************
    subroutine output_persistent_forcing(lun)
!
!  Writes out the time of the next SNI
!
      integer :: lun
!
      if (NO_WARN) print*,lun
!
    endsubroutine output_persistent_forcing
!***********************************************************************
    subroutine rprint_forcing(lreset,lwrite)
!
!  reads and registers print parameters relevant for hydro part
!
!  26-jan-04/axel: coded
!
      use Cdata
      use Sub
!
      integer :: iname
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_rufm=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if(lroot.and.ip<14) print*,'run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'rufm',idiag_rufm)
      enddo
!
!  write column where which magnetic variable is stored
!
      if (lwr) then
        write(3,*) 'i_rufm=',idiag_rufm
      endif
!
    endsubroutine rprint_forcing
!***********************************************************************

endmodule Forcing
