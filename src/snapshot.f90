! $Id$
!
!  Write snapshot files (variables and power spectra).
!
!  04-nov-11/MR: I/O error handling generally introduced
!  16-nov-11/MR: reworked
!  11-Dec-2011/Bourdin.KIS: major cleanup, split file access in IO module
!
module Snapshot
!
  use Cdata
  use Cparam
  use Messages
  use Gpu, only: copy_farray_from_GPU
!
  implicit none
!
  private
!
  interface output_form
    module procedure output_form_int_0D
  endinterface
!
  public :: rsnap, wsnap, wsnap_down, powersnap, output_form, powersnap_prepare
!
  contains
!***********************************************************************
    subroutine wsnap_down(a,flist)
!
!  Write downsampled snapshot file VARd*, labelled consecutively
!  timestep can be different from timestep for full snapshots
!
!  13-feb-14/MR: coded
!   4-oct-16/MR: removed par msnap; its role is taken over by the global vars
!                mvar_down, maux_down for selecting subsets of variables.
!   5-oct-16/MR: modified call to wgrid
!
      use General, only: get_range_no, indgen
      use Boundcond, only: boundconds_x, boundconds_y, boundconds_z
      use General, only: safe_character_assign
      use IO, only: output_snap, log_filename_to_file, lun_output, wgrid
      use HDF5_IO, only: wdim
      use Sub, only: read_snaptime, update_snaptime
      use Grid, only: save_grid, coords_aux
      use Messages, only: warning
!
      real, dimension (:,:,:,:) :: a
      character (len=*), optional :: flist
!
      real, save :: tsnap
      integer, save :: nsnap
      logical, save :: lfirst_call=.true.
      character (len=fnlen) :: file
      character (len=intlen) :: ch

      integer :: ndx, ndy, ndz, isx, isy, isz, ifx, ify, ifz, iax, iay, iaz, &
                 iex, iey, iez, l2s, l2is, m2s, m2is, n2s, n2is, nv1, nv2
      real, dimension(ndown(1)+2*nghost,ndown(2)+2*nghost,ndown(3)+2*nghost,mvar_down+maux_down) :: buffer
      integer, dimension(nghost) :: inds
      real, dimension(nghost) :: dxs_ghost, dys_ghost, dzs_ghost

      call safe_character_assign(file,trim(datadir)//'/tsnap_down.dat')
!
!  At first call, need to initialize tsnap.
!  tsnap calculated in read_snaptime, but only available to root processor.
!
      if (lfirst_call) &
        call read_snaptime(file,tsnap,nsnap,dsnap_down,t)
!
!  Check whether we want to output snapshot. If so, then
!  update ghost zones for var.dat (cheap, since done infrequently).
!
      call update_snaptime(file,tsnap,nsnap,dsnap_down,t,lsnap_down,ch)
      if (lsnap_down) then
!
        if (.not.lstart.and.lgpu) call copy_farray_from_GPU(a)
!
!  Set the range for the variable index.
!  Not yet possible: both mvar_down and maux_down>0, but mvar_down<mvar,
!  That is, two disjoint index ranges needed.
!
        if (mvar_down>0) then
          nv1=1
          if (mvar_down==mvar) then
            nv2=mvar+maux_down
          else
            nv2=mvar_down
          endif
        else
          nv1=mvar+1; nv2=mvar+maux_down
        endif

        if (maux_down>0) call update_auxiliaries(a)
!
!  Number of downsampled data points on local processor in each direction
!
        ndx = ndown(1); ndy = ndown(2); ndz = ndown(3)
!
!  Stepsize for downsampling in each direction
!
        isx = downsampl(1); isy = downsampl(2); isz = downsampl(3)
!
!  Position of first sampled element in f array in each direction
!
        ifx = firstind(1); ify = firstind(2); ifz = firstind(3)
!
!  Index ranges in downsampled array (stored in buffer)
!  (ghost zones are present, but contain valid data only at boundary)
!
        iax = nghost+1 ;   iay = nghost+1 ;   iaz = nghost+1
        iex = iax+ndx-1;   iey = iay+ndy-1;   iez = iaz+ndz-1

!       print*, 'iproc=', iproc, ifx, iex, ndx, ify, iey, ndy, ifz, iez, ndz
!
!  Copy downsampled data from *inner* grid points
!
        buffer(iax:iex,iay:iey,iaz:iez,1:nv2-nv1+1) = a(ifx:l2:isx,ify:m2:isy,ifz:n2:isz,nv1:nv2) 
!
!  Generate ghost zone data
!  TBDone: periodic BC
!
        ldownsampling=.true.
!
!  Save those grid-related variables which will be temporarily overwritten
!
        call save_grid
!
!  Downsample grid
!
        if (any(.not.lequidist)) &
          call warning('wsnap_down','BCs not correctly set for non-equidistant grids')

        if (any(lperi)) &
          call warning('wsnap_down','periodic BCs not correctly set')

        x(iax:iex) = x(ifx:l2:isx)
        y(iay:iey) = y(ify:m2:isy)
        z(iaz:iez) = z(ifz:n2:isz)
        dx=isx*dx;  dy=isy*dy; dz=isz*dz   !!!  Valid only for equidistant grids
!
!  Generate downsampled ghost zone coordinates (only necessary at boundaries)
!
        inds = indgen(nghost)
        if (lfirst_proc_x.or.llast_proc_x) then
          dxs_ghost = dx*inds              !!!  Valid only for equidistant grids
          if (lfirst_proc_x) x(iax-1:1:-1      ) = x(iax)-dxs_ghost
          if (llast_proc_x ) x(iex+1:iex+nghost) = x(iex)+dxs_ghost
        endif
!
        if (lfirst_proc_y.or.llast_proc_y) then
          dys_ghost = dy*inds
          if (lfirst_proc_y) y(iay-1:1:-1      ) = y(iay)-dys_ghost
          if (llast_proc_y ) y(iey+1:iey+nghost) = y(iey)+dys_ghost
        endif
!
        if (lfirst_proc_z.or.llast_proc_z) then
          dzs_ghost = dz*inds
          if (lfirst_proc_z) z(iaz-1:1:-1      ) = z(iaz)-dzs_ghost
          if (llast_proc_z ) z(iez+1:iez+nghost) = z(iez)+dzs_ghost
        endif

        dx_1=1./dx; dy_1=1./dy; dz_1=1./dz        !!!preliminary:
        dx_tilde=0.; dz_tilde=0.; dz_tilde=0.     !  better call construct_grid
 
        if (lfirst_call) then
!
!  At first call, write downsampled grid and its global and local dimensions
!
          call wgrid('grid_down.dat',iex+nghost,iey+nghost,iez+nghost,lwrite=.true.)
          call wdim('dim_down.dat', iex+nghost, iey+nghost, iez+nghost, &
              ceiling(float(nxgrid)/isx)+2*nghost, ceiling(float(nygrid)/isy)+2*nghost, ceiling(float(nzgrid)/isz)+2*nghost, &
              mvar_down, maux_down)
          lfirst_call=.false.
        endif
!
!  Calculate auxiliary grid variables for downsampled grid
!
        call coords_aux(x(1:iex+nghost), y(1:iey+nghost),z(1:iez+nghost))
!
!  For each direction at boundary: save upper index limit and inner start index for
!  periodic BCs; then calculate values at downsampled ghost points from BCs
!
!  fred: Temporarily switched off until boundary conditions correctly
!  implemented on coarse mesh. ghost zones not required for time correlations.
!
        if (lfirst_proc_x.or.llast_proc_x) then
          l2s=l2; l2is=l2i; l2=iex; l2i=l2-nghost+1
          call boundconds_x(buffer)
          l2=l2s; l2i=l2is
        endif
        if (lfirst_proc_y.or.llast_proc_y) then
          m2s=m2; m2is=m2i; m2=iey; m2i=m2-nghost+1
          call boundconds_y(buffer)
          m2=m2s; m2i=m2is
        endif
        if (lfirst_proc_z.or.llast_proc_z) then
          n2s=n2; n2is=n2i; n2=iez; n2i=n2-nghost+1
          call boundconds_z(buffer)
          n2=n2s; n2i=n2is
        endif
!
!  Downsampled ouput in VARd<n> (n>0) snapshot
!
        call safe_character_assign(file,'VARd'//ch)
        open (lun_output, FILE=trim(directory_snap)//'/'//file, &
              FORM='unformatted', status='replace')
        call output_snap(buffer,1,nv2-nv1+1)
        close(lun_output)
!
!  Restore grid (including auxiliaries)
!
        call save_grid(lrestore=.true.)
        ldownsampling=.false.
!
        if (present(flist)) call log_filename_to_file(file,flist)
!
        lsnap_down=.false.
      endif
!
    endsubroutine wsnap_down
!***********************************************************************
    subroutine wsnap(chsnap,a,msnap,enum,flist,noghost,nv1)
!
!  Write snapshot file, labelled consecutively if enum==.true.
!  Otherwise just write a snapshot without label (used for var.dat).
!
!  30-sep-97/axel: coded
!  08-oct-02/tony: expanded file to handle 120 character datadir // '/tsnap.dat'
!   5-apr-03/axel: possibility for additional (hard-to-get) output
!  31-may-03/axel: wsnap can write either w/ or w/o auxiliary variables
!  28-jun-10/julien: added different file formats
!   8-mar-13/MR  : made a assumed-size to work properly with calls in run.f90
!  28-may-21/axel: added nv1_capitalvar
!
      use Boundcond, only: update_ghosts
      use General, only: safe_character_assign, loptest
      use IO, only: output_snap, output_snap_finalize, log_filename_to_file
      use Persist, only: output_persistent
      use Sub, only: read_snaptime, update_snaptime
!
!  The dimension msnap can either be mfarray (for f-array in run.f90)
!  or just mvar (for f-array in start.f90 or df-array in run.f90
!
      integer, intent(in) :: msnap
      real, dimension(mx,my,mz,msnap), intent(inout) :: a
      character(len=*), intent(in) :: chsnap
      character(len=*), intent(in), optional :: flist
      logical, intent(in), optional :: enum, noghost
      integer, intent(in), optional :: nv1
!
      real, save :: tsnap
      real :: time1
      integer, save :: nsnap
      logical, save :: lfirst_call=.true.
      character (len=fnlen) :: file
      character (len=intlen) :: ch
      integer :: nv1_capitalvar
!
!  Output snapshot with label in 'tsnap' time intervals.
!  File keeps the information about number and time of last snapshot.
!
      if (loptest(enum)) then
        call safe_character_assign(file,trim(datadir)//'/tsnap.dat')
!
!  At first call, need to initialize tsnap.
!  tsnap calculated in read_snaptime, but only available to root processor.
!
        if (lfirst_call) then
          call read_snaptime(file,tsnap,nsnap,dsnap,t)
          lfirst_call=.false.
        endif
!
!  It is sometimes of interest to output only last part of the data
!  in the capital var files (e.g., VAR1, VAR2, etc).
!
        if (present(nv1)) then
          nv1_capitalvar=nv1
        else
          nv1_capitalvar=1
        endif
!
!  Check whether we want to output snapshot. If so, then
!  update ghost zones for var.dat (cheap, since done infrequently).
!
        call update_snaptime(file,tsnap,nsnap,dsnap,t,lsnap,ch)
        if (lsnap) then
          if (.not.lstart.and.lgpu) call copy_farray_from_GPU(a)
          call update_ghosts(a)
          if (msnap==mfarray) call update_auxiliaries(a)
          call safe_character_assign(file,trim(chsnap)//ch)
          call output_snap(a,nv1=nv1_capitalvar,nv2=msnap,file=file)
          if (lpersist) call output_persistent(file)
          call output_snap_finalize
          if (present(flist)) call log_filename_to_file(file,flist)
          lsnap=.false.
        endif
!
      else
!
!  Write snapshot without label (typically, var.dat). For dvar.dat we need to
!  make sure that ghost zones are not set on df!
!
        if (.not.lstart.and.lgpu) call copy_farray_from_GPU(a)
        if (msnap==mfarray) then
          if (.not. loptest(noghost)) call update_ghosts(a)
          call update_auxiliaries(a) ! Not if e.g. dvar.dat.
        endif
        ! update ghosts, because 'update_auxiliaries' may change the data
        if (.not. loptest(noghost).or.ncoarse>1) call update_ghosts(a)
        call safe_character_assign(file,trim(chsnap))
        call output_snap(a,nv2=msnap,file=file)
        if (lpersist) call output_persistent(file)
        call output_snap_finalize
        if (present(flist)) call log_filename_to_file(file,flist)
      endif
!
      if (lformat) call output_snap_form (file,a,msnap)
      if (ltec) call output_snap_tec (file,a,msnap)
!
    endsubroutine wsnap
!***********************************************************************
    subroutine rsnap(chsnap,f,msnap,lread_nogrid)
!
!  Read snapshot file.
!
!  24-jun-05/tony: coded from snap reading code in run.f90
!   5-jan-13/axel: allowed for lread_oldsnap_lnrho2rho=T
!   8-mar-13/MR  : made f assumed-size to work properly with calls in run.f90
!  12-feb-15/MR  : added substraction of reference state
!
      use IO, only: input_snap, input_snap_finalize
      use Persist, only: input_persistent
      use SharedVariables, only: get_shared_variable
!
!  The dimension msnap can either be mfarray (for f-array in run.f90)
!  or just mvar (for f-array in start.f90 or df-array in run.f90.
!
      logical :: lread_nogrid
      integer :: msnap, mode,ipscalar
      real, dimension (mx,my,mz,msnap) :: f
      real, dimension (:,:,:,:), allocatable :: f_oversize
      character (len=*) :: chsnap
!
      integer :: ivar
      real, dimension(:,:), pointer :: reference_state
      real :: time1
!
      if (ip<=6.and.lroot) print*,'reading var files'
!
!  possibility of not reading the mesh data nor the time
!  of the snapshot. The mesh information is then taken from
!  proc*/mesh.dat
!
      if (lread_nogrid) then
        mode=0
      else
        mode=1
      endif
!
!  No need to read maux variables as they will be calculated
!  at the first time step -- even if lwrite_aux is set.
!  Allow for the possibility to read in files that didn't
!  have magnetic fields or passive scalar in it.
!  NOTE: for this to work one has to modify *manually* data/param.nml
!  by adding an entry for MAGNETIC_INIT_PARS or PSCALAR_INIT_PARS.
!
      if (lread_oldsnap_nomag) then
        if (lroot) print*,'read old snapshot file (but without magnetic field)'
        call input_snap('var.dat',f,msnap-3,mode)
        if (lpersist) call input_persistent
        call input_snap_finalize
        ! shift the rest of the data
        if (iaz<mvar) then
          do ivar=iaz+1,mvar
            f(:,:,:,ivar)=f(:,:,:,ivar-3)
          enddo
        endif
        f(:,:,:,iax:iaz)=0.
!
!  Read data without passive scalar into new run with passive scalar.
!
      elseif (lread_oldsnap_nopscalar) then
        if (lroot) print*,'read old snapshot file (but without passive scalar)'
        call input_snap(chsnap,f,msnap-1,mode)
        if (lpersist) call input_persistent
        call input_snap_finalize
        ! shift the rest of the data
!
!      If we are using the pscalar_nolog module then ilncc is zero 
!      so we must use icc. 
!
        ipscalar=ilncc
        if (ilncc.eq.0) ipscalar = icc
        if (ipscalar<msnap) then
          do ivar=ipscalar+1,msnap
            f(:,:,:,ivar)=f(:,:,:,ivar-1)
          enddo
        endif
        f(:,:,:,ipscalar) = 0.
!
!  Read data without testfield into new run with testfield.
!
      elseif (lread_oldsnap_notestfield.or.lread_oldsnap_notestflow) then

        if (lroot) then
          if (lread_oldsnap_notestfield) print*,'read old snapshot file (but without testfield),iaatest,iaztestpq,mvar,msnap=', &
                            iaatest,iaztestpq,mvar,msnap
          if (lread_oldsnap_notestflow) print*,'read old snapshot file (but without testflow),iuutest,iuztestpq,mvar,msnap=', &
                            iuutest,iuztestpq,mvar,msnap
        endif

        print*,'read old snapshot file (but without testfield/testflow),ntestfield,ntestflow=',ntestfield,ntestflow
        call input_snap(chsnap,f,msnap-ntestfield-ntestflow,mode)
        if (lpersist) call input_persistent
        call input_snap_finalize
        ! shift the rest of the data
        if (iaztestpq>0.and.iaztestpq<msnap .or. iuztestpq>0.and.iuztestpq<msnap) then
          do ivar=max(iaztestpq,iuztestpq)+1,msnap
            f(:,:,:,ivar)=f(:,:,:,ivar-ntestfield-ntestflow)
          enddo
          f(:,:,:,iaatest:iaatest+ntestfield+ntestflow-1)=0.
        endif
!
!  Read data without testscalar into new run with testscalar.
!
      elseif (lread_oldsnap_notestscalar) then
        if (lroot) print*,'read old snapshot file (but without testscalar),icctest,mvar,msnap=',icctest,mvar,msnap
        call input_snap(chsnap,f,msnap-ntestscalar,mode)
        if (lpersist) call input_persistent
        call input_snap_finalize
        ! shift the rest of the data
        if (iaztestpq<msnap) then
          do ivar=iaztestpq+1,msnap
            f(:,:,:,ivar)=f(:,:,:,ivar-ntestscalar)
          enddo
          f(:,:,:,icctest:icctest+ntestscalar-1)=0.
        endif
!
!  Read data without hydro or density
!
      elseif (lread_oldsnap_nohydro) then
        if (lroot) print*,'read old snapshot file nohydro mvar,msnap=',mvar,msnap
        call input_snap(chsnap,f,msnap-4,mode)
        !call input_snap(chsnap,f,msnap-nohydro_but_efield,mode)
        if (lpersist) call input_persistent
        call input_snap_finalize
        ! shift the rest of the data
        do ivar=msnap,5,-1
          f(:,:,:,ivar)=f(:,:,:,ivar-4)
        enddo
        do ivar=1,4
          f(:,:,:,ivar)=0.
        enddo
!
!  Read data without hydro or density, but with chiral chemical potential
!
      elseif (lread_oldsnap_nohydro_nomu5) then
        if (lroot) print*,'read old snapshot file nohydro mvar,msnap=',mvar,msnap
        call input_snap(chsnap,f,msnap-3,mode)
        if (lpersist) call input_persistent
        call input_snap_finalize
        ! shift the rest of the data
        do ivar=msnap,9,-1
          f(:,:,:,ivar)=f(:,:,:,ivar-5)
        enddo
        do ivar=7,5,-1
          f(:,:,:,ivar)=f(:,:,:,ivar-4)
        enddo
        do ivar=1,4
          f(:,:,:,ivar)=0.
        enddo
        ivar=8
          f(:,:,:,ivar)=0.
!
!  Read data only with vector potential A
!
      elseif (lread_oldsnap_onlyA) then
        if (lroot) print*,'read old snapshot file onlyA mvar,msnap=',mvar,msnap
        !call input_snap(chsnap,f,msnap-4,mode)
        call input_snap(chsnap,f,3,mode)
        if (lpersist) call input_persistent
        call input_snap_finalize
        ! shift the rest of the data
        do ivar=msnap,8,-1
          f(:,:,:,ivar)=0.
        enddo
        do ivar=7,5,-1
          f(:,:,:,ivar)=f(:,:,:,ivar-4)
        enddo
        do ivar=1,4
          f(:,:,:,ivar)=0.
        enddo
!
!  Read data without hydro or density, but with electric field
!
      elseif (lread_oldsnap_nohydro_efield) then
        if (lroot) print*,'read old snapshot file nohydro_but_efield mvar,msnap=',mvar,msnap
        call input_snap(chsnap,f,msnap-1,mode)
        !call input_snap(chsnap,f,msnap-nohydro_but_efield,mode)
        if (lpersist) call input_persistent
        call input_snap_finalize
        ! shift the rest of the data
        do ivar=msnap,8,-1
          f(:,:,:,ivar)=f(:,:,:,ivar-1)
        enddo
        do ivar=5,7
          f(:,:,:,ivar)=f(:,:,:,ivar-4)
        enddo
        do ivar=1,4
          f(:,:,:,ivar)=0.
        enddo
!
!  Read data without hydro or density, but with electric field (spectral)
!  A(16-18), GW(22-39) --> u,lnrho(1-4), A(5-7), GW(8-25)
!
      elseif (lread_oldsnap_nohydro_ekfield) then
        if (lroot) print*,'read old snapshot file nohydro_ekfield mvar,msnap=',mvar,msnap
        allocate(f_oversize(mx,my,mz,msnap+14))
        call input_snap(chsnap,f_oversize,msnap+14,mode)
        if (lpersist) call input_persistent
        call input_snap_finalize
        ! shift the rest of the data
        do ivar=1,4
          f(:,:,:,ivar)=0.
        enddo
        do ivar=5,7
          f(:,:,:,ivar)=f_oversize(:,:,:,ivar+11)
        enddo
        do ivar=8,msnap
          f(:,:,:,ivar)=f_oversize(:,:,:,ivar+14)
        enddo
        deallocate(f_oversize)
!
!  Read data without any var data
!  Note: this can also be done with ivar_omit (see io_dist.f90 and register.f90)
!
      elseif (lread_oldsnap_mskipvar) then
        if (lroot) print*,'read old snapshot file with mskipvar mskipvar,msnap=',mskipvar,msnap
        allocate(f_oversize(mx,my,mz,msnap+mskipvar))
        call input_snap(chsnap,f_oversize,msnap+mskipvar,mode)
        if (lpersist) call input_persistent
        call input_snap_finalize
        ! shift the rest of the data
        do ivar=1,msnap
          f(:,:,:,ivar)=f_oversize(:,:,:,ivar+mskipvar)
        enddo
        deallocate(f_oversize)
!
!  Use default input configuration.
!
      else
        call input_snap(chsnap,f,msnap,mode)
        if (lpersist) call input_persistent(chsnap)
        call input_snap_finalize
      endif
!
!  Read data using lnrho, and now convert to rho.
!  This assumes that one is now using ldensity_nolog=T.
!
      if (lread_oldsnap_lnrho2rho) then
        print*,'convert lnrho -> rho',ilnrho,irho
        if (irho>0) &
          f(:,:,:,irho)=exp(f(:,:,:,ilnrho))
      endif
!
!  Read data using rho, and now convert to lnrho.
!  This assumes that one is now using ldensity_nolog=F.
!  NB: require lupw_rho->lupw_lnrho and diagnostics ..grho..->..glnrho..
!
      if (lread_oldsnap_rho2lnrho) then
        print*,'convert rho -> lnrho',irho,ilnrho
        if (ilnrho>0) &
          f(:,:,:,ilnrho)=log(f(:,:,:,ilnrho))
      endif
!
      if (lsubstract_reference_state) then
        call get_shared_variable('reference_state',reference_state,caller='rsnap')
        do n=n1,n2
        do m=m1,m2
          f(l1:l2,m,n,irho)=f(l1:l2,m,n,irho)-reference_state(:,iref_rho)
          f(l1:l2,m,n,iss )=f(l1:l2,m,n,iss )-reference_state(:,iref_s)
        enddo
        enddo
      endif
!
    endsubroutine rsnap
!***********************************************************************
    subroutine powersnap_prepare
!
!  Routine called from run.f90 to determime when spectra are needed
!
!   7-aug-02/axel: added
!
      use Sub, only: read_snaptime, update_snaptime
!
      logical, save :: lfirst_call=.true.
      character (len=fnlen) :: file
      integer, save :: nspec
      real, save :: tspec
!
!  Output snapshot in 'tpower' time intervals.
!  File keeps the information about time of last snapshot.
!
      file=trim(datadir)//'/tspec.dat'
!
!  At first call, need to initialize tspec.
!  tspec calculated in read_snaptime, but only available to root processor.
!
      if (lfirst_call) then
        call read_snaptime(file,tspec,nspec,dspec,t)
        lfirst_call=.false.
      endif
!
      call update_snaptime(file,tspec,nspec,dspec,t,lspec)
!
    endsubroutine powersnap_prepare
!***********************************************************************
    subroutine powersnap(f,lwrite_only)
!
!  Write a snapshot of power spectrum.
!
!  30-sep-97/axel: coded
!  07-oct-02/nils: adapted from wsnap
!  08-oct-02/tony: expanded file to handle 120 character datadir // '/tspec.dat'
!  28-dec-02/axel: call structure from here; allow optional lwrite_only
!  22-apr-11/MR: added possibility to get xy-power-spectrum from xy_specs
!
      use Boundcond, only: update_ghosts
      use Particles_main, only: particles_powersnap
      use Power_spectrum
      use Pscalar, only: cc2m, gcc2m, rhoccm
      use Struct_func, only: structure
!AXEL use Sub, only: update_snaptime, read_snaptime, curli
      use Sub, only: update_snaptime, curli
      use General, only: itoa
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical, optional :: lwrite_only
!
      real, dimension (:,:,:), allocatable :: b_vec
!AXEL character (len=fnlen) :: file
      logical :: llwrite_only=.false.,ldo_all
!AXEL integer, save :: nspec
!AXEL logical, save :: lfirst_call=.true.
!AXEL real, save :: tspec
      integer :: ivec,im,in,stat,ipos,ispec
      real, dimension (nx) :: bb
      character (LEN=40) :: str,sp1,sp2
      logical :: lfirstcall, lfirstcall_powerhel, lsqrt=.true.
!
!  Allocate memory for b_vec at run time.
!
      allocate(b_vec(nx,ny,nz),stat=stat)
      if (stat>0) call fatal_error('powersnap', &
          'Could not allocate memory for b_vec')
!
!  Set llwrite_only.
!
      if (present(lwrite_only)) llwrite_only=lwrite_only
      ldo_all=.not.llwrite_only
!
!  Check whether we want to output power snapshot. If so, then
!  update ghost zones for var.dat (cheap, since done infrequently).
!
      if (lspec.or.llwrite_only) then

        if (.not.lstart.and.lgpu) call copy_farray_from_GPU(f)
        lfirstcall_powerhel=.true.

        if (ldo_all)  call update_ghosts(f)
        if (vel_spec) call power(f,'u')
        if (r2u_spec) call power(f,'r2u')
        if (r3u_spec) call power(f,'r3u')
        if (oo_spec)  call power(f,'o')
        if (mag_spec) call power(f,'b')
        if (vec_spec) call power(f,'a')
        if (j_spec)   call power_vec(f,'j')
        if (jb_spec)  call powerhel(f,'j.b',lfirstcall_powerhel) !(not ready yet) ! ready now
        if (ja_spec)  call powerhel(f,'j.a',lfirstcall_powerhel) !(for now, use this instead) ! now does j.b spectra
        if (Lor_spec) call powerLor(f,'Lor')
        if (EMF_spec) call powerEMF(f,'EMF')
        if (Tra_spec) call powerTra(f,'Tra')
        lfirstcall=.true.
        if (StT_spec) call powerGWs(f,'StT',lfirstcall)
        if (StX_spec) call powerGWs(f,'StX',lfirstcall)
        if (GWs_spec) call powerGWs(f,'GWs',lfirstcall)
        if (GWh_spec) call powerGWs(f,'GWh',lfirstcall)
        if (GWm_spec) call powerGWs(f,'GWm',lfirstcall)
        if (Str_spec) call powerGWs(f,'Str',lfirstcall)
        if (Stg_spec) call powerGWs(f,'Stg',lfirstcall)
        if (SCL_spec) call powerGWs(f,'SCL',lfirstcall)
        if (VCT_spec) call powerGWs(f,'VCT',lfirstcall)
        if (Tpq_spec) call powerGWs(f,'Tpq',lfirstcall)
        if (TGW_spec) call powerGWs(f,'TGW',lfirstcall)
        if (GWd_spec) call powerhel(f,'GWd',lfirstcall_powerhel)
        if (GWe_spec) call powerhel(f,'GWe',lfirstcall_powerhel)
        if (GWf_spec) call powerhel(f,'GWf',lfirstcall_powerhel)
        if (GWg_spec) call powerhel(f,'GWg',lfirstcall_powerhel)
        if (uxj_spec) call powerhel(f,'uxj',lfirstcall_powerhel)
        if (ou_spec)  call powerhel(f,'kin',lfirstcall_powerhel)
        if (oun_spec) call powerhel(f,'neu',lfirstcall_powerhel)
        if (ele_spec) call powerhel(f,'ele',lfirstcall_powerhel)
        if (ab_spec)  call powerhel(f,'mag',lfirstcall_powerhel)
        if (ub_spec)  call powerhel(f,'u.b',lfirstcall_powerhel)
        if (azbz_spec)call powerhel(f,'mgz',lfirstcall_powerhel)
        if (bb2_spec) call powerhel(f,'bb2',lfirstcall_powerhel)
        if (jj2_spec) call powerhel(f,'jj2',lfirstcall_powerhel)
        if (uzs_spec) call powerhel(f,'uzs',lfirstcall_powerhel)
        if (EP_spec)  call powerhel(f,'bEP',lfirstcall_powerhel)
        if (ro_spec)  call powerscl(f,'ro')
        !if (lro_spec) call powerscl(f,'ro',lsqrt)
        if (lr_spec)  call powerscl(f,'lr')
        if (pot_spec) call powerscl(f,'po')
        if (ud_spec) then
          do n=1,ndustspec
            call power(f,'ud',n)
          enddo
        endif
        if (nd_spec) then
          do n=1,ndustspec
            call powerscl(f,'nd',n)
          enddo
        endif
        if (np_spec)  call powerscl(f,'np')
        if (np_ap_spec) then
          do n=1,ndustrad
            call powerscl(f,'na',n)
          enddo
        endif
        if (rhop_spec)call powerscl(f,'rp')
        if (TT_spec)  call powerscl(f,'TT')
        if (ss_spec)  call powerscl(f,'ss')
        if (cc_spec)  call powerscl(f,'cc')
        if (cr_spec)  call powerscl(f,'cr')
        if (sp_spec)  call powerscl(f,'sp')
        if (ssp_spec) call powerscl(f,'sp',lsqrt=lsqrt)
        if (sssp_spec)call powerscl(f,'Ssp')
        if (mu_spec)  call powerscl(f,'mu')
        !if (smu_spec) call powerscl(f,'mu',lsqrt)
        if (har_spec) call powerscl(f,'hr')
        if (hav_spec) call powerscl(f,'ha')
        if (oned) then
          if (vel_spec) call power_1d(f,'u',1)
          if (mag_spec) call power_1d(f,'b',1)
          if (vec_spec) call power_1d(f,'a',1)
          if (vel_spec) call power_1d(f,'u',2)
          if (mag_spec) call power_1d(f,'b',2)
          if (vec_spec) call power_1d(f,'a',2)
          if (vel_spec) call power_1d(f,'u',3)
          if (mag_spec) call power_1d(f,'b',3)
          if (vec_spec) call power_1d(f,'a',3)
        endif
        if (twod) then
          if (vel_spec) call power_2d(f,'u')
          if (mag_spec) call power_2d(f,'b')
          if (vec_spec) call power_2d(f,'a')
        endif
!
!  xy power spectra
!
        if (uxy_spec  ) call power_xy(f,'u')
        if (bxy_spec  ) call power_xy(f,'b')
        if (jxbxy_spec) call power_xy(f,'jxb')
!
        do ispec=1,n_spectra
!
          if ( xy_specs(ispec)/='' ) then
!
            ipos = index(xy_specs(ispec), '.'); sp1=''; sp2=''
!
            if ( ipos==0 ) then
              call power_xy(f,trim(xy_specs(ispec)))
            else
              str = xy_specs(ispec)
              if ( ipos>1 ) sp1 = str(1:ipos-1)
              if ( ipos<=len_trim(xy_specs(ispec))-1 ) sp2=str(ipos+1:)
!
              if ( sp1=='' .or. sp2=='' ) then
                print*, 'powersnap: Warning - '//trim(xy_specs(ispec))//' no valid identifier !'
              else
                call power_xy(f,sp1,sp2)
              endif
            endif
          endif
!
        enddo
!
!  phi power spectra (in spherical or cylindrical coordinates)
!
        if (vel_phispec) call power_phi(f,'u')
        if (mag_phispec) call power_phi(f,'b')
        if (vec_phispec) call power_phi(f,'a')
        if (ab_phispec)  call powerhel_phi(f,'mag')
        if (ou_phispec)  call powerhel_phi(f,'kin')
!
!  Spectra of particle variables.
!
        if (lparticles) call particles_powersnap(f)
!
!  Structure functions.
!
        do ivec=1,3
          if (lsfb .or. lsfz1 .or. lsfz2 .or. lsfflux .or. lpdfb .or. &
              lpdfz1 .or. lpdfz2) then
             do n=n1,n2
               do m=m1,m2
                 call curli(f,iaa,bb,ivec)
                 im=m-nghost
                 in=n-nghost
                 b_vec(:,im,in)=bb
               enddo
            enddo
            b_vec=b_vec/sqrt(exp(f(l1:l2,m1:m2,n1:n2,ilnrho)))
          endif
          if (lsfu)     call structure(f,ivec,b_vec,'u')
          if (lsfb)     call structure(f,ivec,b_vec,'b')
          if (lsfz1)    call structure(f,ivec,b_vec,'z1')
          if (lsfz2)    call structure(f,ivec,b_vec,'z2')
          if (lsfflux)  call structure(f,ivec,b_vec,'flux')
          if (lpdfu)    call structure(f,ivec,b_vec,'pdfu')
          if (lpdfb)    call structure(f,ivec,b_vec,'pdfb')
          if (lpdfz1)   call structure(f,ivec,b_vec,'pdfz1')
          if (lpdfz2)   call structure(f,ivec,b_vec,'pdfz2')
        enddo
!
!  Do pdf of passive scalar field (if present).
!
        if (rhocc_pdf) call pdf(f,'rhocc',rhoccm,sqrt(cc2m))
        if (cc_pdf)    call pdf(f,'cc'   ,rhoccm,sqrt(cc2m))
        if (lncc_pdf)  call pdf(f,'lncc' ,rhoccm,sqrt(cc2m))
        if (gcc_pdf)   call pdf(f,'gcc'  ,0.    ,sqrt(gcc2m))
        if (lngcc_pdf) call pdf(f,'lngcc',0.    ,sqrt(gcc2m))
        if (lnspecial_pdf) call pdf(f,'lnspecial',0.,1.)
        if (special_pdf) call pdf(f,'special',0.,1.)
!
!  Do pdf for the magnetic field
!
!       if (uuBB_pdf) call pdf(f,'uuBB',-1.,1.)
!
!  azimuthally averaged spectra in polar coorinates
!
        if (ou_omega) call polar_spectrum(f,'kin_omega')
        if (uut_polar)call polar_spectrum(f,'uut')
        if (ouout_polar)call polar_spectrum(f,'ouout')
        if (cor_uu)   call polar_spectrum(f,'uucor')
        if (ou_polar) call polar_spectrum(f,'kin')
        if (ab_polar) call polar_spectrum(f,'mag')
        if (jb_polar) call polar_spectrum(f,'j.b')
!
!  power spectra of xy averaged fields
!
        if (ou_kzspec) call power1d_plane(f,'kin')
        if (ab_kzspec) call power1d_plane(f,'mag')
!
!  spectra of two-time correlations
!
        if (uut_spec)   call power_cor(f,'uut')
        if (ouout_spec) call power_cor(f,'ouout')
        if (out_spec)   call power_cor(f,'out')
        if (uot_spec)   call power_cor(f,'uot')
!
        if (saffman_mag) call quadratic_invariants(f,'saffman_mag')
        if (saffman_mag_c) call quadratic_invariants(f,'saffman_mag_c')
!
        lspec=.false.
      endif
!
      deallocate(b_vec)
!
    endsubroutine powersnap
!***********************************************************************
    subroutine update_auxiliaries(a)
!
      use EquationOfState, only: ioncalc
      use Radiation, only: radtransfer
      use Shock, only: calc_shock_profile,calc_shock_profile_simple
      use Viscosity, only: lvisc_first,calc_viscosity
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: a
!
      if (lshock) then
        call calc_shock_profile(a)
        call calc_shock_profile_simple(a)
      endif
      if (leos_ionization.or.leos_temperature_ionization) call ioncalc(a)
      if (lradiation_ray)  call radtransfer(a)
      if (lvisc_hyper.or.lvisc_smagorinsky) then
        if (.not.lvisc_first.or.lfirst) call calc_viscosity(a)
      endif
!
    endsubroutine update_auxiliaries
!***********************************************************************
    subroutine output_form_int_0D(file,data,lappend)
!
!  Write formatted integer data to a file.
!  Set lappend to false to overwrite the file, default is to append the data.
!
!  CANDIDATE FOR REMOVAL:
!  This routine is only used from "run.f90" if 'ialive' is set.
!  This is a completely useless feature, because the MPI library does the job
!  of checking, if all processes are still alive.
!
!  10-Sep-2015/Bourdin.KIS: marked as candidate for removal
!
      character (len=*), intent(in) :: file
      integer, intent(in) :: data
      logical, intent(in), optional :: lappend
!
      logical :: lappend_opt
      integer :: lun_output = 1
!
      lappend_opt = .false.
      if (present(lappend)) lappend_opt = lappend
!
      if (lappend_opt) then
        open (lun_output, file=trim(directory_snap)//'/'//file, form='formatted', position='append', status='old')
      else
        open (lun_output, file=trim(directory_snap)//'/'//file, form='formatted', status='replace')
      endif
      write (lun_output,*) data
      close (lun_output)
!
    endsubroutine output_form_int_0D
!***********************************************************************
    subroutine output_snap_form(file,a,nv)
!
!  Write FORMATTED snapshot file
!
!  28-june-10/julien: coded (copy from output_snap)
!
      use IO, only: lun_output
!
      character (len=*), intent(in) :: file
      integer, intent(in) :: nv
      real, dimension (mx,my,mz,nv), intent(in) :: a
!
      integer :: i, j, k
!
      open(lun_output,FILE=trim(directory_dist)//trim(file)//'.form')
!
      if (lwrite_2d) then
!
        if (nx==1) then
          do i = m1, m2
            do j = n1, n2
              write(lun_output,'(40(f12.5))') x(l1),y(i),z(j),dx,dy,dz,a(l1,i,j,:)
            enddo
          enddo
        elseif (ny==1) then
          do i = l1, l2
            do j = n1, n2
              write(lun_output,'(40(f12.5))') x(i),y(m1),z(j),dx,dy,dz,a(i,m1,j,:)
            enddo
          enddo
        elseif (nz==1) then
          do i = l1, l2
            do j = m1, m2
              write(lun_output,'(40(f12.5))') x(i),y(j),z(n1),dx,dy,dz,a(i,j,n1,:)
            enddo
          enddo
        else
          call fatal_error('output_snap_form','lwrite_2d used for 3-D simulation!')
        endif
!
      else if (ny==1.and.nz==1) then
!
        do i = l1, l2
          write(lun_output,'(40(f12.5))') x(i),a(i,m1,n1,:)
        enddo
!
      else
!
        do i = l1, l2
          do j = m1, m2
            do k = n1, n2
              write(lun_output,'(40(f12.5))') x(i),y(j),z(k),dx,dy,dz,a(i,j,k,:)
            enddo
          enddo
        enddo
!
      endif
!
      close(lun_output)
!
    endsubroutine output_snap_form
!***********************************************************************
    subroutine output_snap_tec(file,a,nv)
!
!  Write TECPLOT output files (binary)
!
!  28-june-10/julien: coded
!
      use IO, only: lun_output
      use General, only: itoa
!
      character (len=*), intent(in) :: file
      integer, intent(in) :: nv
      real, dimension (mx,my,mz,nv), intent(in) :: a
!
      integer :: i, j, k, kk
      real, dimension (nx*ny*nz) :: xx, yy, zz
!
      open(lun_output,FILE=trim(directory_dist)//trim(file)//'.tec')
!
      kk = 0
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            xx(kk+i) = x(i)
            yy(kk+i) = y(j)
            zz(kk+i) = z(k)
          enddo
          kk = kk + nx
        enddo
      enddo
!
!  Write header
!
      write(lun_output,*) 'TITLE     = "output"'
!
      if (lwrite_2d) then
        if (nx == 1) then
          write(lun_output,*) 'VARIABLES = "y"'
          write(lun_output,*) '"z"'
        elseif (ny == 1) then
          write(lun_output,*) 'VARIABLES = "x"'
          write(lun_output,*) '"z"'
        elseif (nz == 1) then
          write(lun_output,*) 'VARIABLES = "x"'
          write(lun_output,*) '"y"'
        endif
      else
        if ((ny == 1) .and. (nz == 1)) then
          write(lun_output,*) 'VARIABLES = "x"'
        else
          write(lun_output,*) 'VARIABLES = "x"'
          write(lun_output,*) '"y"'
          write(lun_output,*) '"z"'
        endif
      endif
      do i = 1, nv
        write(lun_output,*) '"VAR_'//trim(itoa(i))//'"'
      enddo
!
      write(lun_output,*) 'ZONE T="Zone"'
!
      if (lwrite_2d) then
        if (nx == 1) then
          write(lun_output,*) ' I=1, J=',ny, ', K=',nz
        endif
        if (ny == 1) then
          write(lun_output,*) ' I=',nx, ', J=1, K=',nz
        endif
        if (nz == 1) then
          write(lun_output,*) ' I=',nx, ', J=',ny, ', K=1'
        endif
      else
        if ((ny == 1) .and. (nz == 1)) then
          write(lun_output,*) ' I=',nx, ', J=  1, K=  1'
        else
          write(lun_output,*) ' I=',nx, ', J=',ny, ', K=',nz
        endif
      endif
!
      write(lun_output,*) ' DATAPACKING=BLOCK'
!
!  Write data
!
      if (lwrite_2d) then
        if (nx == 1) then
          write(lun_output,*) yy
          write(lun_output,*) zz
          do j = 1, nv
            write(lun_output,*) a(l1,m1:m2,n1:n2,j)
          enddo
        elseif (ny == 1) then
          write(lun_output,*) xx
          write(lun_output,*) zz
          do j = 1, nv
            write(lun_output,*) a(l1:l2,m1,n1:n2,j)
          enddo
        elseif (nz == 1) then
          write(lun_output,*) xx
          write(lun_output,*) yy
          do j = 1, nv
            write(lun_output,*) a(l1:l2,m1:m2,n1,j)
          enddo
        else
          call fatal_error('output_snap_tec','lwrite_2d used for 3-D simulation!')
        endif
      else
        if ((ny == 1) .and. (nz == 1)) then
          write(lun_output,*) xx
          do j = 1, nv
            write(lun_output,*) a(l1:l2,m1,n1,j)
          enddo
        else
          write(lun_output,*) xx
          write(lun_output,*) yy
          write(lun_output,*) zz
          do j = 1, nv
            write(lun_output,*) a(l1:l2,m1:m2,n1:n2,j)
          enddo
        endif
      endif
!
      close(lun_output)
!
    endsubroutine output_snap_tec
!***********************************************************************
endmodule Snapshot
