! $Id: register.f90,v 1.52 2002-10-25 09:13:23 dobler Exp $

!!!  A module for setting up the f-array and related variables (`register' the
!!!  entropy, magnetic, etc modules). Didn't know where else to put this:
!!!  Entropy uses Sub and Init must be used by both, Start and Run.


module Register

  implicit none 

  contains

!***********************************************************************
    subroutine initialize()
!
!  Call all initialisation hooks, i.e. initialise MPI and register
!  physics modules.
!
!  6-nov-01/wolf: coded
!
      use Cdata
      use Mpicomm
      use Sub
      use IO
      use Param_io
      use Gravity
      use Hydro
      use Forcing
      use Entropy
      use Magnetic
      use Radiation
      use Pscalar
      use Shear
!
!  initialize all mpi stuff
!
      call mpicomm_init
!
!  initialize nvar; is increased by the following routines
!
      nvar = 0 
!
      call register_io
!
      call register_hydro
      call register_density
      call register_forcing
      call register_ent
      call register_aa
      call register_rad
      call register_lncc
      call register_grav
      call register_shear
!
      if (nvar /= mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('Initialize: nvar /= mvar. Fix mvar in cparam.local')
      endif
!
!  initialize headt for root processor only
!
      if (lroot) headt=.true.
!
!  overwrite datadir from datadir.in, if that exists
!
      call get_datadir(datadir)
!
    endsubroutine initialize
!***********************************************************************
    subroutine rprint_list(lreset)
!
!  read in output times from control file
!
!   3-may-01/axel: coded
!
      use Cdata
      use Hydro
      use Entropy
      use Magnetic
      use Radiation
      use Pscalar
!
      integer :: iname,inamez,ixy
      logical :: lreset,exist
!
!  read in the list of variables to be printed
!
      open(1,file='print.in')
      do iname=1,mname
        read(1,*,end=99) cname(iname)
      enddo
99    nname=iname-1
      if (lroot.and.ip<14) print*,'nname=',nname
      close(1)
!
!  read in the list of variables for xy-averages
!
      inquire(file='xyaver.in',exist=exist)
      if (exist) then
        open(1,file='xyaver.in')
        do inamez=1,mnamez
          read(1,*,end=98) cnamez(inamez)
        enddo
98      nnamez=inamez-1
        close(1)
      endif
      if (lroot.and.ip<14) print*,'nnamez=',nnamez
!
!  read in the list of variables for z-averages
!
      inquire(file='zaver.in',exist=exist)
      if (exist) then
        open(1,file='zaver.in')
        do ixy=1,mnamexy
          read(1,*,end=97) cnamexy(ixy)
        enddo
97      nnamexy=ixy-1
        close(1)
      endif
      if (lroot.and.ip<14) print*,'nnamexy=',nnamexy
!
!  check which variables are set
!  For the convenience of idl users, the indices of variables in
!  the f-array and the time_series.dat files are written to data/index.pro
!
      open(3,file=trim(datadir)//'/index.pro')
      call rprint_general(lreset)
      call rprint_hydro(lreset)
      call rprint_density(lreset)
      call rprint_entropy(lreset)
      call rprint_magnetic(lreset)
      call rprint_radiation(lreset)
      call rprint_pscalar(lreset)
      close(3)
!
    endsubroutine rprint_list
!***********************************************************************
    subroutine rprint_general(lreset)
!
!  reads and registers *general* print parameters
!
!   8-jun-02/axel: adapted from hydro
!
      use Cdata
      use Sub
!
      integer :: iname
      logical :: lreset
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        i_t=0;i_it=0;i_dt=0;i_dtc=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if(lroot.and.ip<14) print*,'run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'t',i_t)
        call parse_name(iname,cname(iname),cform(iname),'it',i_it)
        call parse_name(iname,cname(iname),cform(iname),'dt',i_dt)
        call parse_name(iname,cname(iname),cform(iname),'dtc',i_dtc)
      enddo
!
!  write column where which magnetic variable is stored
!
      write(3,*) 'i_t=',i_t
      write(3,*) 'i_it=',i_it
      write(3,*) 'i_dt=',i_dt
      write(3,*) 'i_dtc=',i_dtc
      write(3,*) 'nname=',nname
!
    endsubroutine rprint_general
!***********************************************************************

endmodule Register

!!! End of file register.f90
