! $Id$
!
!  This module provide a way for users to specify custom
!  (i.e. not in the standard Pencil Code) physics, diagnostics etc.
!
!  The module provides a set of standard hooks into the Pencil-Code and
!  currently allows the following customizations:
!
!  Description                                     | Relevant function call
!  ---------------------------------------------------------------------------
!  Special variable registration                   | register_special
!    (pre parameter read)                          |
!  Special variable initialization                 | initialize_special
!    (post parameter read)                         |
!  Special variable finalization                   | finalize_special
!    (deallocation, etc.)                          |
!                                                  |
!  Special initial condition                       | init_special
!   this is called last so may be used to modify   |
!   the mvar variables declared by this module     |
!   or optionally modify any of the other f array  |
!   variables.  The latter, however, should be     |
!   avoided where ever possible.                   |
!                                                  |
!  Special term in the mass (density) equation     | special_calc_density
!  Special term in the momentum (hydro) equation   | special_calc_hydro
!  Special term in the energy equation             | special_calc_energy
!  Special term in the induction (magnetic)        | special_calc_magnetic
!     equation                                     |
!                                                  |
!  Special equation                                | dspecial_dt
!    NOT IMPLEMENTED FULLY YET - HOOKS NOT PLACED INTO THE PENCIL-CODE
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of special_dummies.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 6
!
! PENCILS PROVIDED stressL; stressT
!***************************************************************
!
! HOW TO USE THIS FILE
! --------------------
!
! Change the line above to
!   lspecial = .true.
! to enable use of special hooks.
!
! The rest of this file may be used as a template for your own
! special module.  Lines which are double commented are intended
! as examples of code.  Simply fill out the prototypes for the
! features you want to use.
!
! Save the file with a meaningful name, eg. geo_kws.f90 and place
! it in the $PENCIL_HOME/src/special directory.  This path has
! been created to allow users ot optionally check their contributions
! in to the Pencil-Code SVN repository.  This may be useful if you
! are working on/using the additional physics with somebodyelse or
! may require some assistance from one of the main Pencil-Code team.
!
! To use your additional physics code edit the Makefile.local in
! the src directory under the run directory in which you wish to
! use your additional physics.  Add a line with all the module
! selections to say something like:
!
!   SPECIAL=special/geo_kws
!
! Where geo_kws it replaced by the filename of your new module
! upto and not including the .f90
!
module Special
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages, only: svn_id, fatal_error
!
  implicit none
!
  include '../special.h'
!
! Declare index of new variables in f array (if any).
!
  integer :: ihhL,ihhT,iggL,iggT,istressL,istressT
!
  logical :: lno_transverse_part=.false.
!
! input parameters
  namelist /special_init_pars/ &
    lno_transverse_part
!
! run parameters
  namelist /special_run_pars/ &
    lno_transverse_part
!
! Diagnostic variables (needs to be consistent with reset list below).
!
  integer :: idiag_ggLpt=0       ! DIAG_DOC: $g_{\rm L}(x_1,y_1,z_1,t)$
!
  contains
!***********************************************************************
    subroutine register_special()
!
!  Set up indices for variables in special modules.
!
!  6-oct-03/tony: coded
!
      use Sub, only: register_report_aux
      use FArrayManager
!
      if (lroot) call svn_id( &
           "$Id$")
!
!  Set indices for auxiliary variables.
!
      call farray_register_pde('hhL',ihhL)
      call farray_register_pde('hhT',ihhT)
      call farray_register_pde('ggL',iggL)
      call farray_register_pde('ggT',iggT)
      call farray_register_pde('stressL',istressL)
      call farray_register_pde('stressT',istressT)
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  Called after reading parameters, but before the time loop.
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
      if (lfargo_advection) then
        print*,''
        print*,'Switch '
        print*,' SPECIAL = special/fargo'
        print*,'in src/Makefile.local if you want to use the fargo algorithm'
        print*,''
        call fatal_error('nospecial','initialize_special()')
      endif
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine finalize_special(f)
!
!  Called right before exiting.
!
!  14-aug-2011/Bourdin.KIS: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine finalize_special
!***********************************************************************
    subroutine init_special(f)
!
!  initialise special condition; called from start.f90
!  06-oct-2003/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      intent(inout) :: f
!!
!!  SAMPLE IMPLEMENTATION
!!
!!      select case (initspecial)
!!        case ('nothing'); if (lroot) print*,'init_special: nothing'
!!        case ('zero', '0'); f(:,:,:,iSPECIAL_VARIABLE_INDEX) = 0.
!!        case default
!!          call fatal_error("init_special: No such value for initspecial:" &
!!              ,trim(initspecial))
!!      endselect
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special()
!
!  All pencils that this special module depends on are specified here.
!
!  18-07-06/tony: coded
!
  !?  lpenc_requested(i_bb)=.true.
  !?  lpenc_requested(i_b2)=.true.
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine pencil_interdep_special(lpencil_in)
!
!  Interdependency among pencils provided by this module are specified here.
!
!  18-07-06/tony: coded
!
      logical, dimension(npencils), intent(inout) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
!  Calculate Special pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  24-nov-04/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
      if (lno_transverse_part) then
        p%stressL=f(l1:l2,m,n,ibx)**2
        p%stressT=f(l1:l2,m,n,ibx)*f(l1:l2,m,n,iby)
      else
        p%stressL=f(l1:l2,m,n,istressL)
        p%stressT=f(l1:l2,m,n,istressT)
      endif
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine dspecial_dt(f,df,p)
!
!  calculate right hand side of ONE OR MORE extra coupled PDEs
!  along the 'current' Pencil, i.e. f(l1:l2,m,n) where
!  m,n are global variables looped over in equ.f90
!
!  Due to the multi-step Runge Kutta timestepping used one MUST always
!  add to the present contents of the df array.  NEVER reset it to zero.
!
!  Several precalculated Pencils of information are passed for
!  efficiency.
!
!  06-oct-03/tony: coded
!
      use Diagnostics
      use Sub, only: del2
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: del2hhL,del2hhT
      type (pencil_case) :: p
!
      intent(in) :: f,p
      intent(inout) :: df
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dspecial_dt'
!!      if (headtt) call identify_bcs('special',ispecial)
!
!  dh/dt = g, d2h/dt2 = dg/dt = del2h + S
!
      df(l1:l2,m,n,ihhL)=df(l1:l2,m,n,ihhL)+f(l1:l2,m,n,iggL)
      df(l1:l2,m,n,ihhT)=df(l1:l2,m,n,ihhT)+f(l1:l2,m,n,iggT)
!
      call del2(f,ihhL,del2hhL)
      call del2(f,ihhT,del2hhT)
      df(l1:l2,m,n,iggL)=df(l1:l2,m,n,iggL)+del2hhL+p%stressL
      df(l1:l2,m,n,iggT)=df(l1:l2,m,n,iggT)+del2hhT+p%stressT
!
!         df(l1:l2,m,n,igij+1)=df(l1:l2,m,n,igij+1)+del2hii(:,2)+ &
!           p%bb(:,2)**2-onethird*p%b2
!         df(l1:l2,m,n,igij+2)=df(l1:l2,m,n,igij+2)+del2hii(:,3)+ &
!           p%bb(:,3)**2-onethird*p%b2
!         df(l1:l2,m,n,igij+3)=df(l1:l2,m,n,igij+3)+del2hij(:,1)+p%bb(:,1)*p%bb(:,2)
!         df(l1:l2,m,n,igij+4)=df(l1:l2,m,n,igij+4)+del2hij(:,2)+p%bb(:,1)*p%bb(:,3)
!         df(l1:l2,m,n,igij+5)=df(l1:l2,m,n,igij+5)+del2hij(:,3)+p%bb(:,2)*p%bb(:,3)
!
!  diagnostics
!
       if (ldiagnos) then
         if (lroot.and.m==mpoint.and.n==npoint) then
           if (idiag_ggLpt/=0) call save_name(f(lpoint,m,n,iggL),idiag_ggLpt)
         endif
       endif
!
      !call keep_compiler_quiet(f,df)
      !call keep_compiler_quiet(p)
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine read_special_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=special_init_pars, IOSTAT=iostat)
!
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=special_init_pars)
!
    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=special_run_pars, IOSTAT=iostat)
!
    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=special_run_pars)
!
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine special_before_boundary(f)
!
!  Possibility to modify the f array before the boundaries are
!  communicated.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array
!
!  30-mar-17/axel: moved stuff from special_after_boundary to here
!
      !use Boundcond, only: zero_ghosts, update_ghosts
      use Sub, only: gij, curl_mn
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(nx,3,3) :: aij
      real, dimension(nx,3) :: bb
!
      integer :: i
!
!  Compute magnetic stress Mij(x,t)=Bi(x,t)*Bj(x,t) in real space
!
!  Find bb if as communicated auxiliary.
!
      !call zero_ghosts(f, iax, iaz)
      !call update_ghosts(f, iax, iaz)
  !   mn_loop: do imn = 1, ny * nz
  !     m = mm(imn)
  !     n = nn(imn)
  !     call gij(f, iaa, aij, 1)
  !     call curl_mn(aij, bb, f(l1:l2,m,n,iax:iaz))
!
!  Add imposed field, if any
!
  !     bext: if (lB_ext_in_comaux) then
  !       call get_bext(b_ext)
  !       forall(j = 1:3, b_ext(j) /= 0.0) bb(:,j) = bb(:,j) + b_ext(j)
  !       if (headtt .and. imn == 1) &
  !           print *, 'magnetic_before_boundary: B_ext = ', b_ext
  !     endif bext
!
  !     f(l1:l2,m,n,ibx:ibz) = bb
  !   enddo mn_loop
!     endif getbb
!
!     do i=1,n_special_modules
!       call caller(special_sub_handles(i,I_SPECIAL_BEFORE_BOUNDARY),1,f)
!     enddo
!
    endsubroutine special_before_boundary
!***********************************************************************
    subroutine special_after_boundary(f)
!
!  dummy routine
!
!  15-jan-08/axel: coded
!
      use Fourier, only: fourier_transform
!
      real, dimension (:,:,:), allocatable :: S_re, S_im, T_re, T_im, k2
      real, dimension (:), allocatable :: kx, ky, kz
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i,ikx,iky,ikz,stat
      logical :: lscale_tobox1=.true.
      real :: scale_factor, P2
      intent(inout) :: f
!
!  Allocate memory for arrays.
!
      allocate(k2(nx,ny,nz),stat=stat)
      if (stat>0) call fatal_error('special_after_boundary','Could not allocate memory for k2')
!
      allocate(S_re(nx,ny,nz),stat=stat)
      if (stat>0) call fatal_error('special_after_boundary','Could not allocate memory for S_re')
      allocate(S_im(nx,ny,nz),stat=stat)
      if (stat>0) call fatal_error('special_after_boundary','Could not allocate memory for S_im')
!
      allocate(T_re(nx,ny,nz),stat=stat)
      if (stat>0) call fatal_error('special_after_boundary','Could not allocate memory for T_re')
      allocate(T_im(nx,ny,nz),stat=stat)
      if (stat>0) call fatal_error('special_after_boundary','Could not allocate memory for T_im')
!
      allocate(kx(nxgrid),stat=stat)
      if (stat>0) call fatal_error('special_after_boundary', &
          'Could not allocate memory for kx')
      allocate(ky(nygrid),stat=stat)
      if (stat>0) call fatal_error('special_after_boundary', &
          'Could not allocate memory for ky')
      allocate(kz(nzgrid),stat=stat)
      if (stat>0) call fatal_error('special_after_boundary', &
          'Could not allocate memory for kz')
!
!  calculate k^2
!
        scale_factor=1
        if (lscale_tobox1) scale_factor=2*pi/Lx
        kx=cshift((/(i-(nxgrid+1)/2,i=0,nxgrid-1)/),+(nxgrid+1)/2)*scale_factor
!
        scale_factor=1
        if (lscale_tobox1) scale_factor=2*pi/Ly
        ky=cshift((/(i-(nygrid+1)/2,i=0,nygrid-1)/),+(nygrid+1)/2)*scale_factor
!
        scale_factor=1
        if (lscale_tobox1) scale_factor=2*pi/Lz
        kz=cshift((/(i-(nzgrid+1)/2,i=0,nzgrid-1)/),+(nzgrid+1)/2)*scale_factor
!
!  Set k^2 array. Note that in Fourier space, kz is the fastest index and has
!  the full nx extent (which, currently, must be equal to nxgrid).
!
        if (lroot .AND. ip<10) &
             print*,'special_after_boundary:fft ...'
        do iky=1,nz
          do ikx=1,ny
            do ikz=1,nx
              k2(ikz,ikx,iky)=kx(ikx+ipy*ny)**2+ky(iky+ipz*nz)**2+kz(ikz+ipx*nx)**2
            enddo
          enddo
        enddo
!
!  compute 1/k2 for components of unit vector
!
        if (lroot) k2(1,1,1) = 1.  ! Avoid division by zero
        k2=1./k2
!
!  Assemble stress
!
!  Do T11
!  Go into Fourier space
!
        T_im=0.0
        T_re=f(l1:l2,m1:m2,n1:n2,ibx)**2
        call fourier_transform(T_re,T_im)
!
!  projection operator
!
        do iky=1,nz
          do ikx=1,ny
            do ikz=1,nx
!
!  Real part of (ux, uy, uz) -> vx, vy, vz
!  (kk.uu)/k2, vi = ui - ki kj uj
!
!  r = .5*P11^2
!
              P2=.5*(1.-kx(ikx+ipy*ny)**2*k2(ikz,ikx,iky))**2
              S_re(ikz,ikx,iky)=P2*T_re(ikz,ikx,iky)
              S_im(ikz,ikx,iky)=P2*T_im(ikz,ikx,iky)

              !v_im(ikz,ikx,iky,1)=u_im(ikz,ikx,iky,1)-kx(ikx+ipy*ny)*r(ikz,ikx,iky)
              !v_im(ikz,ikx,iky,2)=u_im(ikz,ikx,iky,2)-ky(iky+ipz*nz)*r(ikz,ikx,iky)
              !v_im(ikz,ikx,iky,3)=u_im(ikz,ikx,iky,3)-kz(ikz+ipx*nx)*r(ikz,ikx,iky)

!
!  r = .5*P11*P12
!
              P2=.5*(1.-kx(ikx+ipy*ny)**2*k2(ikz,ikx,iky))*(1.-ky(iky+ipz*nz)**2*k2(ikz,ikx,iky))
              T_re(ikz,ikx,iky)=P2*T_re(ikz,ikx,iky)
              T_im(ikz,ikx,iky)=P2*T_im(ikz,ikx,iky)
!
            enddo
          enddo
        enddo
!
!  back to real space
!
        call fourier_transform(S_re,S_im,linv=.true.)
        call fourier_transform(T_re,T_im,linv=.true.)
!
!  add (or set) corresponding stress
!
        f(l1:l2,m1:m2,n1:n2,istressL)=S_re
        f(l1:l2,m1:m2,n1:n2,istressT)=T_re
!
!  Deallocate arrays.
!
      if (allocated(k2))   deallocate(k2)
      if (allocated(S_re)) deallocate(S_re)
      if (allocated(S_im)) deallocate(S_im)
      if (allocated(T_re)) deallocate(T_re)
      if (allocated(T_im)) deallocate(T_im)
      if (allocated(kx)) deallocate(kx)
      if (allocated(ky)) deallocate(ky)
      if (allocated(kz)) deallocate(kz)
!
    endsubroutine special_after_boundary
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
!  Reads and registers print parameters relevant to special.
!
!  06-oct-03/tony: coded
!
      use Diagnostics
!
      integer :: iname
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!!!
!!!  reset everything in case of reset
!!!  (this needs to be consistent with what is defined above!)
!!!
      if (lreset) then
        idiag_ggLpt=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'ggLpt',idiag_ggLpt)
      enddo
!!
!!!  write column where which magnetic variable is stored
!!      if (lwr) then
!!        write(3,*) 'i_SPECIAL_DIAGNOSTIC=',i_SPECIAL_DIAGNOSTIC
!!      endif
!!
    endsubroutine rprint_special
!***********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../special_dummies.inc'
!********************************************************************
endmodule Special
