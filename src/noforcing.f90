! $Id: noforcing.f90,v 1.13 2004-01-26 14:46:02 brandenb Exp $

module Forcing

!! Dummy module for hydro/mhd without forcing

  use Cdata
  use General

  implicit none
  integer :: dummy              ! We cannot define empty namelists
 
  namelist /forcing_init_pars/ dummy
  namelist /forcing_run_pars/  dummy

  ! other variables (needs to be consistent with reset list below)
  integer :: i_rufm=0

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
           "$Id: noforcing.f90,v 1.13 2004-01-26 14:46:02 brandenb Exp $")
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
      if(ip==0) print*,'lstarting=',lstarting !(to keep compiler quiet)
    endsubroutine initialize_forcing
!***********************************************************************
    subroutine addforce(df)
!
      use Cdata
!
!  add forcing in timestep()
!
      real, dimension (mx,my,mz,mvar) :: df
!
      if(ip==1) print*,df !(to remove compiler warnings)
    endsubroutine addforce
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
        i_rufm=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if(lroot.and.ip<14) print*,'run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'rufm',i_rufm)
      enddo
!
!  write column where which magnetic variable is stored
!
      if (lwr) then
        write(3,*) 'i_rufm=',i_rufm
      endif
!
    endsubroutine rprint_forcing
!***********************************************************************

endmodule Forcing
