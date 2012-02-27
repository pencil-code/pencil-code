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
!
  implicit none
!
  private
!
  public :: rsnap, wsnap, powersnap
  public :: shift_dt
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
!  28-jun-10/julien: added different file formats
!
      use Boundcond, only: update_ghosts
      use General, only: safe_character_assign
      use IO, only: output_snap, output_snap_finalize, log_filename_to_file
      use Persist, only: output_persistent
      use Sub, only: read_snaptime, update_snaptime
!
!  The dimension msnap can either be mfarray (for f-array in run.f90)
!  or just mvar (for f-array in start.f90 or df-array in run.f90
!
      character (len=*) :: chsnap, flist
      integer :: msnap
      real, dimension (mx,my,mz,msnap) :: a
      logical :: enum,enum_,noghost
      optional :: enum, flist, noghost
!
      real, save :: tsnap
      integer, save :: nsnap
      logical, save :: lfirst_call=.true.
      logical :: lsnap
      character (len=fnlen) :: file
      character (len=intlen) :: ch
!
      if (present(enum)) then
        enum_=enum
      else
        enum_=.false.
      endif
!
!  Output snapshot with label in 'tsnap' time intervals.
!  File keeps the information about number and time of last snapshot.
!
      if (enum_) then
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
!  Check whether we want to output snapshot. If so, then
!  update ghost zones for var.dat (cheap, since done infrequently).
!
        call update_snaptime(file,tsnap,nsnap,dsnap,t,lsnap,ch)
        if (lsnap) then
          call update_ghosts(a)
          if (msnap==mfarray) call update_auxiliaries(a)
          call safe_character_assign(file,trim(chsnap)//ch)
          call output_snap(file,a,msnap)
          call output_persistent(file)
          call output_snap_finalize()
          if (ip<=10.and.lroot) print*,'wsnap: written snapshot ',file
          if (present(flist)) call log_filename_to_file(file,flist)
        endif
!
      else
!
!  Write snapshot without label (typically, var.dat). For dvar.dat we need to
!  make sure that ghost zones are not set on df!
!
        if (msnap==mfarray) then
          if (present(noghost)) then
            if (.not. noghost) call update_ghosts(a)
          else
            call update_ghosts(a)
          endif
          call update_auxiliaries(a) ! Not if e.g. dvar.dat.
        endif
        if (present(noghost)) then
          if (.not. noghost) call update_ghosts(a)
        else
          call update_ghosts(a)
        endif
        call safe_character_assign(file,trim(chsnap))
        call output_snap(file,a,msnap)
        call output_persistent(file)
        call output_snap_finalize()
        if (present(flist)) call log_filename_to_file(file,flist)
      endif
!
      if (lformat) call output_snap_form (file,a,msnap)
      if (ltec) call output_snap_tec (file,a,msnap)
!
    endsubroutine wsnap
!***********************************************************************
    subroutine rsnap(chsnap,f,msnap)
!
!  Read snapshot file.
!
!  24-jun-05/tony: coded from snap reading code in run.f90
!
      use IO, only: input_snap, input_snap_finalize
      use Persist, only: input_persistent
!
!  The dimension msnap can either be mfarray (for f-array in run.f90)
!  or just mvar (for f-array in start.f90 or df-array in run.f90.
!
      integer :: msnap
      real, dimension (mx,my,mz,msnap) :: f
      character (len=*) :: chsnap
!
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
        call input_snap('var.dat',f,msnap-3,1)
        call input_persistent()
        call input_snap_finalize()
        ! shift the rest of the data
        if (iaz<mvar) then
          do ivar=iaz+1,mvar
            f(:,:,:,ivar)=f(:,:,:,ivar-3)
          enddo
          f(:,:,:,iax:iaz)=0.
        endif
!
!  Read data without passive scalar into new run with passive scalar.
!
      elseif (lread_oldsnap_nopscalar) then
        print*,'read old snapshot file (but without passive scalar)'
        call input_snap(chsnap,f,msnap-1,1)
        call input_persistent()
        call input_snap_finalize()
        ! shift the rest of the data
        if (ilncc<mvar) then
          do ivar=ilncc+1,mvar
            f(:,:,:,ivar)=f(:,:,:,ivar-1)
          enddo
          f(:,:,:,ilncc)=0.
        endif
!
!  Read data without testfield into new run with testfield.
!
      elseif (lread_oldsnap_notestfield) then
        print*,'read old snapshot file (but without testfield),iaatest,iaztestpq,mvar,msnap=',iaatest,iaztestpq,mvar,msnap
        call input_snap(chsnap,f,msnap-ntestfield,1)
        call input_persistent()
        call input_snap_finalize()
        ! shift the rest of the data
        if (iaztestpq<msnap) then
          do ivar=iaztestpq+1,msnap
            f(:,:,:,ivar)=f(:,:,:,ivar-ntestfield)
          enddo
          f(:,:,:,iaatest:iaatest+ntestfield-1)=0.
        endif
!
!  Read data without testscalar into new run with testscalar.
!
      elseif (lread_oldsnap_notestscalar) then
        print*,'read old snapshot file (but without testscalar),icctest,mvar,msnap=',icctest,mvar,msnap
        call input_snap(chsnap,f,msnap-ntestscalar,1)
        call input_persistent()
        call input_snap_finalize()
        ! shift the rest of the data
        if (iaztestpq<msnap) then
          do ivar=iaztestpq+1,msnap
            f(:,:,:,ivar)=f(:,:,:,ivar-ntestscalar)
          enddo
          f(:,:,:,icctest:icctest+ntestscalar-1)=0.
        endif
      else
        call input_snap(chsnap,f,msnap,1)
        call input_persistent(chsnap)
        call input_snap_finalize()
      endif
!
    endsubroutine rsnap
!***********************************************************************
   subroutine powersnap(f,lwrite_only)
!
!  Write a snapshot of power spectrum.
!
!  30-sep-97/axel: coded
!  07-oct-02/nils: adapted from wsnap
!  08-oct-02/tony: expanded file to handle 120 character datadir // '/tspec.dat'
!  28-dec-02/axel: call structure from herel; allow optional lwrite_only
!  22-apr-11/MR: added possibility to get quantity for xy-power-spectrum from xy_specs
!
      use Boundcond, only: update_ghosts
      use Particles_main, only: particles_powersnap
      use Power_spectrum
      use Pscalar, only: cc2m, gcc2m, rhoccm
      use Struct_func, only: structure
      use Sub, only: update_snaptime, read_snaptime, curli
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical, optional :: lwrite_only
!
      real, dimension (:,:,:), allocatable :: b_vec
      character (len=fnlen) :: file
      logical :: lspec,llwrite_only=.false.,ldo_all
      integer, save :: nspec
      logical, save :: lfirst_call=.true.
      real, save :: tspec
      integer :: ivec,im,in,stat,ipos,ispec
      real, dimension (nx) :: bb
      character (LEN=40) :: str,sp1,sp2
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
      lspec=.true.
!
!  Output snapshot in 'tpower' time intervals.
!  File keeps the information about time of last snapshot.
!
      file=trim(datadir)//'/tspec.dat'
!
!  At first call, need to initialize tspec.
!  tspec calculated in read_snaptime, but only available to root processor.
!
      if (ldo_all .and. lfirst_call) then
        call read_snaptime(file,tspec,nspec,dspec,t)
        lfirst_call=.false.
      endif
!
!  Check whether we want to output power snapshot. If so, then
!  update ghost zones for var.dat (cheap, since done infrequently).
!
      if (ldo_all) &
           call update_snaptime(file,tspec,nspec,dspec,t,lspec)
      if (lspec.or.llwrite_only) then
        if (ldo_all)  call update_ghosts(f)
        if (vel_spec) call power(f,'u')
        if (r2u_spec) call power(f,'r2u')
        if (r3u_spec) call power(f,'r3u')
        if (oo_spec)  call power(f,'o')
        if (mag_spec) call power(f,'b')
        if (vec_spec) call power(f,'a')
        if (j_spec)   call power_vec(f,'j')
!         if (jb_spec)   call powerhel(f,'jb')
        if (uxj_spec) call powerhel(f,'uxj')
        if (ou_spec)  call powerhel(f,'kin')
        if (ab_spec)  call powerhel(f,'mag')
        if (azbz_spec)call powerhel(f,'mgz')
        if (ub_spec)  call powerhel(f,'u.b')
        if (EP_spec)  call powerhel(f,'bEP')
        if (ro_spec)  call powerscl(f,'ro')
        if (lr_spec)  call powerscl(f,'lr')
        if (TT_spec)  call powerscl(f,'TT')
        if (ss_spec)  call powerscl(f,'ss')
        if (cc_spec)  call powerscl(f,'cc')
        if (cr_spec)  call powerscl(f,'cr')
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
!
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
    subroutine shift_dt(dt_)
!
!  Hack to make the code output the VARN files at EXACTLY the times
!  defined by dsnap, instead of slightly after it.
!
!  03-aug-11/wlad: coded
!
      use Sub, only: read_snaptime
      use General, only: safe_character_assign
!
      real, intent(inout) :: dt_
      real :: tsnap
      integer :: nsnap
      character (len=fnlen) :: file
!
!  Read the output time defined by dsnap.
!
      call safe_character_assign(file,trim(datadir)//'/tsnap.dat')
      call read_snaptime(file,tsnap,nsnap,dsnap,t)
!
!  Adjust the time-step accordingly, so that the next timestepping
!  lands the simulation at the precise time defined by dsnap.
!
      if ((tsnap-t > dtmin).and.(t+dt_ > tsnap)) then
        dt_=tsnap-t
      else
        dt_=dt_
      endif
!
    endsubroutine shift_dt
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
      integer :: i, j, k, io_err
!
      open(lun_output,FILE=trim(directory_dist)//trim(file)//'.form',IOSTAT=io_err)
      if (outlog(io_err,'open',trim(file)//'.form',dist=lun_output)) return
!
      if (lwrite_2d) then
!
        if (nx==1) then
          do i = m1, m2
            do j = n1, n2
              write(lun_output,'(40(f12.5))',IOSTAT=io_err) x(l1),y(i),z(j),dx,dy,dz,a(l1,i,j,:)
              if (outlog(io_err,'write x, y, z, dx, dy, dz, a(l1,i,j,:)')) return
            enddo
          enddo
        elseif (ny==1) then
          do i = l1, l2
            do j = n1, n2
              write(lun_output,'(40(f12.5))',IOSTAT=io_err) x(i),y(m1),z(j),dx,dy,dz,a(i,m1,j,:)
              if (outlog(io_err,'write x, y, z, dx, dy, dz, a(i,m1,j,:)')) return
            enddo
          enddo
        elseif (nz==1) then
          do i = l1, l2
            do j = m1, m2
              write(lun_output,'(40(f12.5))',IOSTAT=io_err) x(i),y(j),z(n1),dx,dy,dz,a(i,j,n1,:)
              if (outlog(io_err,'write x, y, z, dx, dy, dz, a(i,j,n1,:)')) return
            enddo
          enddo
        else
          call fatal_error('output_snap_form','lwrite_2d used for 3-D simulation!')
        endif
!
      else if (ny==1.and.nz==1) then
!
        do i = l1, l2
          write(lun_output,'(40(f12.5))',IOSTAT=io_err) x(i),a(i,m1,n1,:)
          if (outlog(io_err,'write x, a')) return
        enddo
!
      else
!
        do i = l1, l2
          do j = m1, m2
            do k = n1, n2
              write(lun_output,'(40(f12.5))',IOSTAT=io_err) x(i),y(j),z(k),dx,dy,dz,a(i,j,k,:)
              if (outlog(io_err,'write x, y, z, dx, dy, dz, a')) return
            enddo
          enddo
        enddo
!
      endif
!
      close(lun_output,IOSTAT=io_err)
      if (outlog(io_err,'close')) continue
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
      integer :: i, j, k, kk, io_err
      real, dimension (nx*ny*nz) :: xx, yy, zz
!
      open(lun_output,FILE=trim(directory_dist)//trim(file)//'.tec',IOSTAT=io_err)
      if (outlog(io_err,'open',trim(file)//'.tec',dist=lun_output)) return
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
      write(lun_output,*,IOSTAT=io_err) 'TITLE     = "output"'
      if (outlog(io_err,'write TITLE')) return
!
      if (lwrite_2d) then
!
        if (nx==1) then
          write(lun_output,*,IOSTAT=io_err) 'VARIABLES = "y"'
          if (outlog(io_err,'write "VARIABLES = y"')) return
          write(lun_output,*,IOSTAT=io_err) '"z"'
          if (outlog(io_err,'write "z"')) return
        elseif (ny==1) then
          write(lun_output,*,IOSTAT=io_err) 'VARIABLES = "x"'
          if (outlog(io_err,'write "VARIABLES = x"')) return
          write(lun_output,*,IOSTAT=io_err) '"z"'
          if (outlog(io_err,'write "z"')) return
        elseif (nz==1) then
          write(lun_output,*,IOSTAT=io_err) 'VARIABLES = "x"'
          if (outlog(io_err,'write "VARIABLES = x"')) return
          write(lun_output,*,IOSTAT=io_err) '"y"'
          if (outlog(io_err,'write "y"')) return
        endif
!
      else
!
        if (ny==1.and.nz==1) then
          write(lun_output,*,IOSTAT=io_err) 'VARIABLES = "x"'
          if (outlog(io_err,'write "VARIABLES = x"')) return
        else
          write(lun_output,*,IOSTAT=io_err) 'VARIABLES = "x"'
          if (outlog(io_err,'write "VARIABLES = x"')) return
          write(lun_output,*,IOSTAT=io_err) '"y"'
          if (outlog(io_err,'write "VARIABLES = y"')) return
          write(lun_output,*,IOSTAT=io_err) '"z"'
          if (outlog(io_err,'write "VARIABLES = z"')) return
        endif
!
      endif
      do i = 1, nv
        write(lun_output,*,IOSTAT=io_err) '"VAR_'//trim(itoa(i))//'"'
        if (outlog(io_err,'write filename')) return
      enddo
!
      write(lun_output,*,IOSTAT=io_err) 'ZONE T="Zone"'
      if (outlog(io_err,'write "ZONE T=Zone"')) return
!
      if (lwrite_2d) then
        if (nx==1) then
          write(lun_output,*,IOSTAT=io_err) ' I=1, J=',ny, ', K=',nz
          if (outlog(io_err,'write ny, nz')) return
        endif
        if (ny==1) then
          write(lun_output,*,IOSTAT=io_err) ' I=',nx, ', J=1, K=',nz
          if (outlog(io_err,'write nx, nz')) return
        endif
        if (nz==1) then
          write(lun_output,*,IOSTAT=io_err) ' I=',nx, ', J=',ny, ', K=1'
          if (outlog(io_err,'write nx, ny')) return
        endif
      else
        if (ny==1.and.nz==1) then
          write(lun_output,*,IOSTAT=io_err) ' I=',nx, ', J=  1, K=  1'
          if (outlog(io_err,'write nx')) return
        else
          write(lun_output,*,IOSTAT=io_err) ' I=',nx, ', J=',ny, ', K=',nz
          if (outlog(io_err,'write nx, ny, nz')) return
        endif
      endif
!
      write(lun_output,*,IOSTAT=io_err) ' DATAPACKING=BLOCK'
      if (outlog(io_err,'write "DATAPACKING=BLOCK"')) return
!
!
!  Write data
!
      if (lwrite_2d) then
        if (nx==1) then
!
          write(lun_output,*,IOSTAT=io_err) yy
          if (outlog(io_err,'write yy')) return
          write(lun_output,*,IOSTAT=io_err) zz
          if (outlog(io_err,'write zz')) return
!
          do j = 1, nv
            write(lun_output,*,IOSTAT=io_err) a(l1,m1:m2,n1:n2,j)
            if (outlog(io_err,'write a')) return
          enddo
!
        elseif (ny==1) then
!
          write(lun_output,*,IOSTAT=io_err) xx
          if (outlog(io_err,'write xx')) return
!
          write(lun_output,*,IOSTAT=io_err) zz
          if (outlog(io_err,'write zz')) return
!
          do j = 1, nv
            write(lun_output,*,IOSTAT=io_err) a(l1:l2,m1,n1:n2,j)
            if (outlog(io_err,'write a')) return
          enddo
!
        elseif (nz==1) then
          write(lun_output,*,IOSTAT=io_err) xx
          if (outlog(io_err,'write xx')) return
!
          write(lun_output,*,IOSTAT=io_err) yy
          if (outlog(io_err,'write yy')) return
!
          do j = 1, nv
            write(lun_output,*,IOSTAT=io_err) a(l1:l2,m1:m2,n1,j)
            if (outlog(io_err,'write a')) return
          enddo
!
        else
          call fatal_error('output_snap_tec','lwrite_2d used for 3-D simulation!')
        endif
      else if (ny==1.and.nz==1) then
!
             write(lun_output,*,IOSTAT=io_err) xx
             if (outlog(io_err,'write xx')) return
!
             do j = 1, nv
               write(lun_output,*,IOSTAT=io_err) a(l1:l2,m1,n1,j)
               if (outlog(io_err,'write a')) return
             enddo
!
           else
!
             write(lun_output,*,IOSTAT=io_err) xx
             if (outlog(io_err,'write xx')) return
!
             write(lun_output,*,IOSTAT=io_err) yy
             if (outlog(io_err,'write yy')) return
!
             write(lun_output,*,IOSTAT=io_err) zz
             if (outlog(io_err,'write zz')) return
!
             do j = 1, nv
               write(lun_output,*,IOSTAT=io_err) a(l1:l2,m1:m2,n1:n2,j)
               if (outlog(io_err,'write a')) return
             enddo
!
           endif
!
      close(lun_output,IOSTAT=io_err)
      if (outlog(io_err,'close')) continue
!
    endsubroutine output_snap_tec
!***********************************************************************
endmodule Snapshot
