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
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 12
! MAUX CONTRIBUTION 3
!
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
!!  namelist /special_init_pars/ dummy
!
!!  namelist /special_run_pars/ dummy
!
! Declare index of new variables in f array (if any).
!
  integer :: ihij,igij
!
!! Diagnostic variables (needs to be consistent with reset list below).
!
  integer :: idiag_g22pt=0       ! DIAG_DOC: $g_{22}(x_1,y_1,z_1,t)$
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
      call farray_register_pde('hij',ihij,vector=6)
      call farray_register_pde('gij',igij,vector=6)
!
!  Set indices for auxiliary variables.
!
      !call farray_register_auxiliary('bb',ibb)
      call register_report_aux('bb', ibb, ibx, iby, ibz)
print*,'AXEL1: registered, ibb, ibx, iby, ibz=',ibb, ibx, iby, ibz
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
      lpenc_requested(i_bb)=.true.
      lpenc_requested(i_b2)=.true.
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
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
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
      use Sub, only: del2v
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: del2hii,del2hij
      type (pencil_case) :: p
!
      integer :: j,jhij,jgij
!
      intent(in) :: f,p
      intent(inout) :: df
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dspecial_dt'
!!      if (headtt) call identify_bcs('special',ispecial)
!
        if (lmagnetic) then
!         print*,'AXEL:',p%bb(:,2)
!
!  g11=1, g22=2, g33=3, g12=4, g13=5, g23=6,
!  g11=0, g22=1, g33=2, g12=3, g13=4, g23=5,
!
          do j=1,6
            jhij=ihij-1+j
            jgij=igij-1+j
            df(l1:l2,m,n,jhij)=df(l1:l2,m,n,jhij)+f(l1:l2,m,n,jgij)
          enddo
          call del2v(f,ihij  ,del2hii)
          call del2v(f,ihij+3,del2hij)
          df(l1:l2,m,n,igij+0)=df(l1:l2,m,n,igij+0)+del2hii(:,1)+ &
            p%bb(:,1)**2-onethird*p%b2
          df(l1:l2,m,n,igij+1)=df(l1:l2,m,n,igij+1)+del2hii(:,2)+ &
            p%bb(:,2)**2-onethird*p%b2
          df(l1:l2,m,n,igij+2)=df(l1:l2,m,n,igij+2)+del2hii(:,3)+ &
            p%bb(:,3)**2-onethird*p%b2
          df(l1:l2,m,n,igij+3)=df(l1:l2,m,n,igij+3)+del2hij(:,1)+p%bb(:,1)*p%bb(:,2)
          df(l1:l2,m,n,igij+4)=df(l1:l2,m,n,igij+4)+del2hij(:,2)+p%bb(:,1)*p%bb(:,3)
          df(l1:l2,m,n,igij+5)=df(l1:l2,m,n,igij+5)+del2hij(:,3)+p%bb(:,2)*p%bb(:,3)
        else
          call fatal_error("dspecial_dt","need magnetic field")
        endif
!
!  diagnostics
!
       if (ldiagnos) then
         if (lroot.and.m==mpoint.and.n==npoint) then
           if (idiag_g22pt/=0) call save_name(p%bb(lpoint-nghost,2),idiag_g22pt)
         endif
       endif
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine calc_lspecial_pars(f)
!
!  dummy routine
!
!  15-jan-08/axel: coded
!
      use Fourier, only: fourier_transform
!
      real, dimension (:,:,:,:), allocatable :: B_re, B_im, v_re, v_im
      real, dimension (:,:,:), allocatable :: k2, r
      real, dimension (:), allocatable :: kx, ky, kz
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i,ikx,iky,ikz,stat
      logical :: lscale_tobox1=.true.
      real :: scale_factor
      intent(inout) :: f
!
!  Allocate memory for arrays.
!
      allocate(k2(nx,ny,nz),stat=stat)
      if (stat>0) call fatal_error('calc_lspecial_pars','Could not allocate memory for k2')
      allocate(r(nx,ny,nz),stat=stat)
      if (stat>0) call fatal_error('calc_lspecial_pars','Could not allocate memory for r')
!
      allocate(B_re(nx,ny,nz,3),stat=stat)
      if (stat>0) call fatal_error('calc_lspecial_pars','Could not allocate memory for B_re')
      allocate(B_im(nx,ny,nz,3),stat=stat)
      if (stat>0) call fatal_error('calc_lspecial_pars','Could not allocate memory for B_im')
!
      allocate(v_re(nx,ny,nz,3),stat=stat)
      if (stat>0) call fatal_error('calc_lspecial_pars','Could not allocate memory for v_re')
      allocate(v_im(nx,ny,nz,3),stat=stat)
      if (stat>0) call fatal_error('calc_lspecial_pars','Could not allocate memory for v_im')
!
      allocate(kx(nxgrid),stat=stat)
      if (stat>0) call fatal_error('calc_lspecial_pars', &
          'Could not allocate memory for kx')
      allocate(ky(nygrid),stat=stat)
      if (stat>0) call fatal_error('calc_lspecial_pars', &
          'Could not allocate memory for ky')
      allocate(kz(nzgrid),stat=stat)
      if (stat>0) call fatal_error('calc_lspecial_pars', &
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
             print*,'calc_lspecial_pars:fft ...'
        do iky=1,nz
          do ikx=1,ny
            do ikz=1,nx
              k2(ikz,ikx,iky)=kx(ikx+ipy*ny)**2+ky(iky+ipz*nz)**2+kz(ikz+ipx*nx)**2
            enddo
          enddo
        enddo
        if (lroot) k2(1,1,1) = 1.  ! Avoid division by zero
!
!  Compute Mij(x,t)=Bi(x,t)*Bj(x,t) in real space
!
!  Find bb if as communicated auxiliary.
!
      call zero_ghosts(f, iax, iaz)
      call update_ghosts(f, iax, iaz)
      mn_loop: do imn = 1, ny * nz
        m = mm(imn)
        n = nn(imn)
        call gij(f, iaa, aij, 1)
        call curl_mn(aij, bb, f(l1:l2,m,n,iax:iaz))
!
!  Add imposed field, if any
!
        bext: if (lB_ext_in_comaux) then
          call get_bext(b_ext)
          forall(j = 1:3, b_ext(j) /= 0.0) bb(:,j) = bb(:,j) + b_ext(j)
          if (headtt .and. imn == 1) print *, 'magnetic_before_boundary: B_ext
= ', b_ext
        endif bext
        f(l1:l2,m,n,ibx:ibz) = bb
        enddo mn_loop
      endif getbb

!
!  Go into Fourier space
!
        v_im=0.0
        v_re=f(l1:l2,m1:m2,n1:n2,iax:iaz)
        do i=1,3
          call fourier_transform(v_re(:,:,:,i),v_im(:,:,:,i))
        enddo !i
!
!  Compute the magnetic field in Fourier space from vector potential.
!  In this definition of the Fourier transform, nabla = +ik.
!
      do ikz=1,nz; do iky=1,ny; do ikx=1,nx
!
!  (vx, vy, vz) -> Bx
!
        B_re(ikz,ikx,iky,1)=+( &
          -ky(iky+ipz*nz)*v_im(ikz,ikx,iky,3) &
          +kz(ikz+ipx*nx)*v_im(ikz,ikx,iky,2) )
        B_im(ikz,ikx,iky,1)=+( &
          +ky(iky+ipz*nz)*v_re(ikz,ikx,iky,3) &
          -kz(ikz+ipx*nx)*v_re(ikz,ikx,iky,2) )
!
!  (vx, vy, vz) -> By
!
        B_re(ikz,ikx,iky,2)=+( &
          -kz(ikz+ipx*nx)*v_im(ikz,ikx,iky,1) &
          +kx(ikx+ipy*ny)*v_im(ikz,ikx,iky,3) )
        B_im(ikz,ikx,iky,2)=+( &
          +kz(ikz+ipx*nx)*v_re(ikz,ikx,iky,1) &
          -kx(ikx+ipy*ny)*v_re(ikz,ikx,iky,3) )
!
!  (vx, vy, vz) -> Bz
!
        B_re(ikz,ikx,iky,3)=+( &
          -kx(ikx+ipy*ny)*v_im(ikz,ikx,iky,2) &
          +ky(iky+ipz*nz)*v_im(ikz,ikx,iky,1) )
        B_im(ikz,ikx,iky,3)=+( &
          +kx(ikx+ipy*ny)*v_re(ikz,ikx,iky,2) &
          -ky(iky+ipz*nz)*v_re(ikz,ikx,iky,1) )
!
      enddo; enddo; enddo
!
!  back to real space
!
print*,'AXEL2: registered, ibb, ibx, iby, ibz=',ibb, ibx, iby, ibz
        do i=1,3
          call fourier_transform(B_re(:,:,:,i),B_im(:,:,:,i),linv=.true.)
          f(l1:l2,m1:m2,n1:n2,ibb+i-1)=B_re(:,:,:,i)
        enddo !i
!
    endsubroutine calc_lspecial_pars
!***********************************************************************
    subroutine read_special_init_pars(unit,iostat)
!
      include '../unit.h'
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      call keep_compiler_quiet(present(iostat))
!
    endsubroutine read_special_init_pars
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
        idiag_g22pt=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'g22pt',idiag_g22pt)
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
