! $Id: snapshot.f90 12535 2009-12-13 07:49:11Z AxelBrandenburg $
!
!  Write snapshot files (variables and power spectra).
!
module Snapshot
!
  use Cdata
  use Cparam
  use Messages
!
  implicit none
!
  private
!
  integer :: lun_output=92
!
  public :: rsnap, wsnap
  public :: output_globals, input_globals
!
  contains
!***********************************************************************
    subroutine wsnap(chsnap,a,msnap,enum,flist,noghost)
!
!  Write snapshot file, labelled consecutively if enum==.true.
!  Otherwise just write a snapshot without label (used for var.dat).
!
!  30-sep-97/axel: coded
!  08-oct-02/tony: expanded file to handle 120 character datadir // '/tsnap.dat'
!   5-apr-03/axel: possibility for additional (hard-to-get) output
!  31-may-03/axel: wsnap can write either w/ or w/o auxiliary variables
!
      use Mpicomm
      use Boundcond, only: update_ghosts
      use General, only: safe_character_assign
      use Sub, only: read_snaptime,update_snaptime
      use IO, only: log_filename_to_file
!
!  The dimension msnap can either be mfarray (for f-array in run.f90)
!  or just mvar (for f-array in start.f90 or df-array in run.f90
!
      character (len=*) :: chsnap, flist
      real, dimension (mx,my,mz,msnap) :: a
      integer :: msnap
      logical :: enum, noghost
      optional :: enum, flist, noghost
!
      real, save :: tsnap
      integer, save :: ifirst=0,nsnap
      logical :: lsnap
      character (len=fnlen) :: file
      character (len=5) :: ch
!
!  Output snapshot with label in 'tsnap' time intervals.
!  File keeps the information about number and time of last snapshot.
!
      if (enum) then
        call safe_character_assign(file,trim(datadir)//'/tsnap.dat')
!
!  At first call, need to initialize tsnap.
!  tsnap calculated in read_snaptime, but only available to root processor.
!
        if (ifirst==0) then
          call read_snaptime(file,tsnap,nsnap,dsnap,t)
          ifirst=1
        endif
!
!  Check whether we want to output snapshot. If so, then
!  update ghost zones for var.dat (cheap, since done infrequently).
!
        call update_snaptime(file,tsnap,nsnap,dsnap,t,lsnap,ch,ENUM=.true.)
        if (lsnap) then
          call update_ghosts(a)
          if (msnap==mfarray) call update_auxiliaries(a)
          call output_snap(chsnap//ch,a,msnap)
          if (ip<=10.and.lroot) print*,'wsnap: written snapshot ',chsnap//ch
          if (present(flist)) call log_filename_to_file(chsnap//ch,flist)
        endif
!
      else
!
!  Write snapshot without label (typically, var.dat). For dvar.dat we need to
!  make sure that ghost zones are not set on df!
!
        if (present(noghost)) then
          if (.not. noghost) call update_ghosts(a)
        else
          call update_ghosts(a)
        endif
        if (msnap==mfarray) call update_auxiliaries(a) ! Not if e.g. dvar.dat.
        call output_snap(chsnap,a,msnap)
        if (present(flist)) call log_filename_to_file(chsnap,flist)
      endif
!
    endsubroutine wsnap
!***********************************************************************
    subroutine rsnap(chsnap,f,msnap)
!
!  Read snapshot file.
!
!  24-jun-05/tony: coded from snap reading code in run.f90
!
      use Mpicomm
!
!  The dimension msnap can either be mfarray (for f-array in run.f90)
!  or just mvar (for f-array in start.f90 or df-array in run.f90.
!
      integer :: msnap
      real, dimension (mx,my,mz,msnap) :: f
      character (len=*) :: chsnap
      integer :: ivar
!
        if (ip<=6.and.lroot) print*,'reading var files'
!
!  No need to read maux variables as they will be calculated
!  at the first time step -- even if lwrite_aux is set.
!  Allow for the possibility to read in files that didn't
!  have magnetic fields or passive scalar in it.
!  NOTE: for this to work one has to modify *manually* data/param.nml
!  by adding an entry for MAGNETIC_INIT_PARS or PSCALAR_INIT_PARS.
!
!DM: I do not understand why we need to shift the data below. 
! I seem to need to set f(:,:,:,iax:iaz) = 0 . Otherwise 
! the vector potential is initialised as junk data. And then
! init_aa just adds to it, so junk remains junk. Anyhow
! the initialisation to zero cannot do any harm.  
!
        if (lread_oldsnap_nomag) then
          f(:,:,:,iax:iaz)=0.
          print*,'read old snapshot file (but without magnetic field)'
          call input_snap(trim(directory_snap)//'/var.dat',f,msnap-3,1)
          ! shift the rest of the data
          if (iaz<mvar) then
            do ivar=iaz+1,mvar
              f(:,:,:,ivar)=f(:,:,:,ivar-3)
            enddo
            f(:,:,:,iax:iaz)=0.
          endif
        else
          call input_snap(chsnap,f,msnap,1)
        endif
!
    endsubroutine rsnap
!***********************************************************************
    subroutine output_snap(file,a,nv)
!
!  Write snapshot file, always write time and mesh, could add other things
!  version for vector field.
!
!  11-apr-97/axel: coded
!
      use Mpicomm, only: start_serialize,end_serialize
!
      integer :: nv
      real, dimension (mx,my,mz,nv) :: a
      character (len=*) :: file
      real :: t_sp   ! t in single precision for backwards compatibility
! 
      t_sp = t
      if (ip<=8.and.lroot) print*,'output_vect: nv =', nv
!
      if (lserial_io) call start_serialize()
      open(lun_output,FILE=file,FORM='unformatted')
      if (lwrite_2d) then
        if (nx==1) then
          write(lun_output) a(l1,:,:,:)
        elseif (ny==1) then
          write(lun_output) a(:,m1,:,:)
        elseif (nz==1) then
          write(lun_output) a(:,:,n1,:)
        else
          call fatal_error('output_snap','lwrite_2d used for 3-D simulation!')
        endif
      else
        write(lun_output) a
      endif
!
!  Write shear at the end of x,y,z,dx,dy,dz.
!  At some good moment we may want to treat deltay like with
!  other modules and call a corresponding i/o parameter module.
!
      if (lshear) then
        write(lun_output) t_sp,x,y,z,dx,dy,dz,deltay
      else
        write(lun_output) t_sp,x,y,z,dx,dy,dz
      endif
!
      close(lun_output)
      if (lserial_io) call end_serialize()
!
    endsubroutine output_snap
!***********************************************************************
    subroutine input_snap(file,a,nv,mode)
!
!  Read snapshot file, possibly with mesh and time (if mode=1).
!
!  11-apr-97/axel: coded
!
      use Mpicomm, only: start_serialize,end_serialize
!
      character (len=*) :: file
      integer :: nv,mode
      real, dimension (mx,my,mz,nv) :: a
      real :: t_sp   ! t in single precision for backwards compatibility
!
      if (lserial_io) call start_serialize()
      open(1,FILE=file,FORM='unformatted')
!      if (ip<=8) print*,'input_snap: open, mx,my,mz,nv=',mx,my,mz,nv
      if (lwrite_2d) then
        if (nx==1) then
          read(1) a(4,:,:,:)
        elseif (ny==1) then
          read(1) a(:,4,:,:)
        elseif (nz==1) then
          read(1) a(:,:,4,:)
        else
          call fatal_error('input_snap','lwrite_2d used for 3-D simulation!')
        endif
      else
        read(1) a
      endif
      if (ip<=8) print*,'input_snap: read ',file
      if (mode==1) then
!
!  Check whether we want to read deltay from snapshot.
!
        if (lshear) then
          read(1) t_sp,x,y,z,dx,dy,dz,deltay
        else
          read(1) t_sp,x,y,z,dx,dy,dz
        endif
        t = t_sp
!
        if (ip<=3) print*,'input_snap: ip,x=',ip,x
        if (ip<=3) print*,'input_snap: y=',y
        if (ip<=3) print*,'input_snap: z=',z
!
      endif
!
      close(1)
      if (lserial_io) call end_serialize()
!
    endsubroutine input_snap
!***********************************************************************
    subroutine output_globals(file,a,nv)
!
!  Write snapshot file of globals, always write mesh.
!
!  10-nov-06/tony: coded
!
      use Mpicomm, only: start_serialize,end_serialize
!
      integer :: nv
      real, dimension (mx,my,mz,nv) :: a
      character (len=*) :: file
!
      if (ip<=8.and.lroot) print*,'output_vect: nv =', nv
!
      if (lserial_io) call start_serialize()
      open(lun_output,FILE=file,FORM='unformatted')

      if (lwrite_2d) then
        if (nx==1) then
          write(lun_output) a(4,:,:,:)
        elseif (ny==1) then
          write(lun_output) a(:,4,:,:)
        elseif (nz==1) then
          write(lun_output) a(:,:,4,:)
        else
          call fatal_error('output_globals','lwrite_2d used for 3-D simulation!')
        endif
      else
        write(lun_output) a
      endif
!
      close(lun_output)
!
      if (lserial_io) call end_serialize()
!
    endsubroutine output_globals
!***********************************************************************
    subroutine input_globals(filename,a,nv)
!
!  Read globals snapshot file, ignoring mesh.
!
!  10-nov-06/tony: coded
!
      use Mpicomm, only: start_serialize,end_serialize
!
      character (len=*) :: filename
      integer :: nv
      real, dimension (mx,my,mz,nv) :: a
!
      if (lserial_io) call start_serialize()
!
      open(1,FILE=filename,FORM='unformatted')
      if (ip<=8) print*,'input_globals: open, mx,my,mz,nv=',mx,my,mz,nv
      if (lwrite_2d) then
        if (nx==1) then
          read(1) a(4,:,:,:)
        elseif (ny==1) then
          read(1) a(:,4,:,:)
        elseif (nz==1) then
          read(1) a(:,:,4,:)
        else
          call fatal_error('input_globals','lwrite_2d used for 3-D simulation!')
        endif
      else 
        read(1) a
      endif
      if (ip<=8) print*,'input_globals: read ',filename
      close(1)
!
      if (lserial_io) call end_serialize()
!
    endsubroutine input_globals
!***********************************************************************
    subroutine update_auxiliaries(a)
!
      use Shock, only: calc_shock_profile,calc_shock_profile_simple
      use EquationOfState, only: ioncalc
      use Viscosity, only: lvisc_first,calc_viscosity
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: a
!
      if (lshock) then
        call calc_shock_profile(a)
        call calc_shock_profile_simple(a)
      endif
!
    endsubroutine update_auxiliaries
!***********************************************************************
endmodule Snapshot
