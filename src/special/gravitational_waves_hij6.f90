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
! MVAR CONTRIBUTION 12
! MAUX CONTRIBUTION 2
!
! PENCILS PROVIDED stress_ij(6)
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
  use Initcond
  use General, only: keep_compiler_quiet
  use Messages, only: svn_id, fatal_error
!
  implicit none
!
  include '../special.h'
!
! Declare index of new variables in f array (if any).
!
  character (len=labellen) :: inithij='nothing'
  character (len=labellen) :: initgij='nothing'
  real :: amplhij=0., amplgij=0.
  real :: kx_hij=0., ky_hij=0., kz_hij=0.
  real :: kx_gij=0., ky_gij=0., kz_gij=0.
  real :: dummy=0., diffhh=0., diffgg=0.
  logical :: lno_transverse_part=.false., lsame_diffgg_as_hh=.true.
  logical :: lswitch_sign_e_X=.true., ldebug_print=.false.
  real, dimension(3,3) :: ij_table
  real, parameter :: c2=1
!
! input parameters
  namelist /special_init_pars/ &
    inithij, initgij, amplhij, amplgij, &
    kx_hij, ky_hij, kz_hij, &
    kx_gij, ky_gij, kz_gij
!
! run parameters
  namelist /special_run_pars/ &
    diffhh, diffgg, lsame_diffgg_as_hh, ldebug_print
!
! Diagnostic variables (needs to be consistent with reset list below).
!
  integer :: idiag_g11pt=0       ! DIAG_DOC: $g_{11}(x_1,y_1,z_1,t)$
  integer :: idiag_g22pt=0       ! DIAG_DOC: $g_{22}(x_1,y_1,z_1,t)$
  integer :: idiag_g33pt=0       ! DIAG_DOC: $g_{33}(x_1,y_1,z_1,t)$
  integer :: idiag_g12pt=0       ! DIAG_DOC: $g_{12}(x_1,y_1,z_1,t)$
  integer :: idiag_g23pt=0       ! DIAG_DOC: $g_{23}(x_1,y_1,z_1,t)$
  integer :: idiag_g31pt=0       ! DIAG_DOC: $g_{31}(x_1,y_1,z_1,t)$
  integer :: idiag_ggTpt=0       ! DIAG_DOC: $g_{T}(x_1,y_1,z_1,t)$
  integer :: idiag_ggXpt=0       ! DIAG_DOC: $g_{X}(x_1,y_1,z_1,t)$
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
!  register ggT and ggX as auxiliary arrays
!
      call farray_register_auxiliary('ggT',iggT)
      call farray_register_auxiliary('ggX',iggX)
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
      use EquationOfState, only: cs0
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!  Check whether diffgg=diffhh (which  is the default)
!
      if (lsame_diffgg_as_hh) diffgg=diffhh
!
!  set index table
!
      ij_table(1,1)=1
      ij_table(2,2)=2
      ij_table(3,3)=3
      ij_table(1,2)=4
      ij_table(2,3)=5
      ij_table(3,1)=6
      ij_table(2,1)=4
      ij_table(3,2)=5
      ij_table(1,3)=6
!
!  give a warning if cs0**2=1
!
      if (cs0==1.) call fatal_error('gravitational_waves_hij6', &
          'cs0 should probably not be unity')
!
      call keep_compiler_quiet(f)
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
!
!  initial condition for hij
!
      select case (inithij)
        case ('nothing'); if (lroot) print*,'init_special: nothing'
        case ('coswave-kx'); call coswave(amplhij,f,ihij,kx=kx_hij)
        case default
          call fatal_error("init_special: No such value for inithij:" &
              ,trim(inithij))
      endselect
!
!  initial condition for gij
!
      select case (initgij)
        case ('nothing'); if (lroot) print*,'init_special: nothing'
        case ('coswave-kx'); call coswave(amplgij,f,igij,kx=kx_gij)
        case default
          call fatal_error("init_special: No such value for initgij:" &
              ,trim(initgij))
      endselect
!
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special()
!
!  All pencils that this special module depends on are specified here.
!
!  18-07-06/tony: coded
!
      if (lhydro)    lpenc_requested(i_uu)=.true.
      if (lmagnetic) lpenc_requested(i_bb)=.true.
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
      integer :: i, j, ij
!
!  Construct stress tensor
!
      p%stress_ij=0.0
      do j=1,3
      do i=1,j
        ij=ij_table(i,j)
        if (lhydro)    p%stress_ij(:,ij)=p%stress_ij(:,ij)+p%uu(:,i)*p%uu(:,j)
        if (lmagnetic) p%stress_ij(:,ij)=p%stress_ij(:,ij)+p%bb(:,i)*p%bb(:,j)
      enddo
      enddo
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
      real, dimension (nx,6) :: del2hij, del2gij
      type (pencil_case) :: p
!
      integer :: ij,jhij,jgij
!
      intent(in) :: f,p
      intent(inout) :: df
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dspecial_dt'
!
      do ij=1,6
        jhij=ihij-1+ij
        jgij=igij-1+ij
        call del2(f,jhij,del2hij(:,ij))
        df(l1:l2,m,n,jhij)=df(l1:l2,m,n,jhij)+f(l1:l2,m,n,jgij)
        df(l1:l2,m,n,jgij)=df(l1:l2,m,n,jgij)+c2*del2hij(:,ij)+p%stress_ij(:,ij)
        if (diffhh/=0.) then
          call del2(f,jhij,del2hij(:,ij))
          df(l1:l2,m,n,jhij)=df(l1:l2,m,n,jhij)+diffhh*del2hij(:,ij)
        endif
        if (diffgg/=0.) then
          call del2(f,jgij,del2gij(:,ij))
          df(l1:l2,m,n,jgij)=df(l1:l2,m,n,jgij)+diffgg*del2gij(:,ij)
        endif
      enddo
!
!  timestep constraint
!
      if (lfirst.and.ldt) advec_cs2=max(advec_cs2,c2*dxyz_2)
!
!  diagnostics
!
       if (ldiagnos) then
         if (lroot.and.m==mpoint.and.n==npoint) then
           if (idiag_g11pt/=0) call save_name(f(lpoint-nghost,m,n,igij+1-1),idiag_g11pt)
           if (idiag_g22pt/=0) call save_name(f(lpoint-nghost,m,n,igij+2-1),idiag_g22pt)
           if (idiag_g33pt/=0) call save_name(f(lpoint-nghost,m,n,igij+3-1),idiag_g33pt)
           if (idiag_g12pt/=0) call save_name(f(lpoint-nghost,m,n,igij+4-1),idiag_g12pt)
           if (idiag_g23pt/=0) call save_name(f(lpoint-nghost,m,n,igij+5-1),idiag_g23pt)
           if (idiag_g31pt/=0) call save_name(f(lpoint-nghost,m,n,igij+6-1),idiag_g31pt)
           if (idiag_ggTpt/=0) call save_name(f(lpoint-nghost,m,n,iggT),idiag_ggTpt)
           if (idiag_ggXpt/=0) call save_name(f(lpoint-nghost,m,n,iggX),idiag_ggXpt)
         endif
       endif
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
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
    endsubroutine special_before_boundary
!***********************************************************************
    subroutine special_after_boundary(f)
!
!  Possibility to modify the f array after the boundaries are
!  communicated.
!
!  06-jul-06/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      logical :: always=.true.
!
      if (lspec.and.lfirst) then
        call compute_gT_and_gX_from_gij(f)
      endif
!
    endsubroutine special_after_boundary
!***********************************************************************
    subroutine compute_gT_and_gX_from_gij(f)
!
!  Compute the transverse part of the stress tensor by going into Fourier space.
!
!  15-jan-08/axel: coded
!
      use Fourier, only: fourier_transform, fft_xyz_parallel
!
      real, dimension (:,:,:,:), allocatable :: Tpq_re, Tpq_im
      real, dimension (:,:,:), allocatable :: one_over_k2, S_T_re, S_T_im, S_X_re, S_X_im
      real, dimension (:), allocatable :: kx, ky, kz
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (6) :: Pij, e_T, e_X, Sij_re, Sij_im
      real, dimension (3) :: e1, e2
      integer :: i,j,p,q,ikx,iky,ikz,stat,ij,pq,ip,jq,jgij
      logical :: lscale_tobox1=.true.
      real :: scale_factor, fact, k1, k2, k3
      intent(inout) :: f
!
!  For testing purposes, if lno_transverse_part=T, we would not need to
!  compute the Fourier transform, so we would skip the rest.
!
!  Allocate memory for arrays.
!
      allocate(one_over_k2(nx,ny,nz),stat=stat)
      if (stat>0) call fatal_error('compute_gT_and_gX_from_gij','Could not allocate memory for one_over_k2')
!
      allocate(S_T_re(nx,ny,nz),stat=stat)
      if (stat>0) call fatal_error('compute_gT_and_gX_from_gij','Could not allocate memory for S_T_re')
!
      allocate(S_T_im(nx,ny,nz),stat=stat)
      if (stat>0) call fatal_error('compute_gT_and_gX_from_gij','Could not allocate memory for S_T_im')
!
      allocate(S_X_re(nx,ny,nz),stat=stat)
      if (stat>0) call fatal_error('compute_gT_and_gX_from_gij','Could not allocate memory for S_X_re')
!
      allocate(S_X_im(nx,ny,nz),stat=stat)
      if (stat>0) call fatal_error('compute_gT_and_gX_from_gij','Could not allocate memory for S_X_im')
!
      allocate(Tpq_re(nx,ny,nz,6),stat=stat)
      if (stat>0) call fatal_error('compute_gT_and_gX_from_gij','Could not allocate memory for Tpq_re')
!
      allocate(Tpq_im(nx,ny,nz,6),stat=stat)
      if (stat>0) call fatal_error('compute_gT_and_gX_from_gij','Could not allocate memory for Tpq_im')
!
      allocate(kx(nxgrid),stat=stat)
      if (stat>0) call fatal_error('compute_gT_and_gX_from_gij', &
          'Could not allocate memory for kx')
      allocate(ky(nygrid),stat=stat)
      if (stat>0) call fatal_error('compute_gT_and_gX_from_gij', &
          'Could not allocate memory for ky')
      allocate(kz(nzgrid),stat=stat)
      if (stat>0) call fatal_error('compute_gT_and_gX_from_gij', &
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
!  But call it one_over_k2.
!
      do iky=1,nz
        do ikx=1,ny
          do ikz=1,nx
            one_over_k2(ikz,ikx,iky)=kx(ikx+ipy*ny)**2 &
                                    +ky(iky+ipz*nz)**2 &
                                    +kz(ikz+ipx*nx)**2
          enddo
        enddo
      enddo
!
!  compute 1/k2 for components of unit vector
!  Avoid division by zero
!
      if (lroot) one_over_k2(1,1,1) = 1.  ! Avoid division by zero
      one_over_k2=1./one_over_k2
      if (lroot) one_over_k2(1,1,1) = 0.  ! set origin to zero
!
!  Assemble stress, Tpq
!
      Tpq_re=0.0
      Tpq_im=0.0
      do ij=1,6
        jgij=igij-1+ij
        Tpq_re(:,:,:,ij)=f(l1:l2,m1:m2,n1:n2,jgij)
      enddo
!
!  Fourier transform all 6 components
!
      do ij=1,6
        call fourier_transform(Tpq_re(:,:,:,ij),Tpq_im(:,:,:,ij))
!--     call fft_xyz_parallel(Tpq_re(:,:,:,ij),Tpq_im(:,:,:,ij))
      enddo
!
!  P11, P22, P33, P12, P23, P31
!
      S_T_re=0.
      S_T_im=0.
      S_X_re=0.
      S_X_im=0.
      do iky=1,nz
        do ikx=1,ny
          do ikz=1,nx
!
!  compute e_T and e_X; determine first preferred direction,
!  which is a component with the smallest component by modulus.
!
            k1=kx(ikx+ipy*ny)
            k2=ky(iky+ipz*nz)
            k3=kz(ikz+ipx*nx)
!
!  find two vectors e1 and e2 to compute e_T and e_X
!
            if (lroot.and.ikx==1.and.iky==1.and.ikz==1) then
              e1=0.
              e2=0.
              Pij(1)=0.
              Pij(2)=0.
              Pij(3)=0.
              Pij(4)=0.
              Pij(5)=0.
              Pij(6)=0.
            else
              if(abs(k1)<abs(k2)) then
                if(abs(k1)<abs(k3)) then !(k1 is pref dir)
                  e1=(/0.,-k3,+k2/)
                  e2=(/k2**2+k3**2,-k2*k1,-k3*k1/)
                else !(k3 is pref dir)
                  e1=(/k2,-k1,0./)
                  e2=(/k1*k3,k2*k3,-(k1**2+k2**2)/)
                endif
              else !(k2 smaller than k1)
                if(abs(k2)<abs(k3)) then !(k2 is pref dir)
                  e1=(/-k3,0.,+k1/)
                  e2=(/+k1*k2,-(k1**2+k3**2),+k3*k2/)
                else !(k3 is pref dir)
                  e1=(/k2,-k1,0./)
                  e2=(/k1*k3,k2*k3,-(k1**2+k2**2)/)
                endif
              endif
              e1=e1/sqrt(e1(1)**2+e1(2)**2+e1(3)**2)
              e2=e2/sqrt(e2(1)**2+e2(2)**2+e2(3)**2)
              Pij(1)=1.-kx(ikx+ipy*ny)**2*one_over_k2(ikz,ikx,iky)
              Pij(2)=1.-ky(iky+ipz*nz)**2*one_over_k2(ikz,ikx,iky)
              Pij(3)=1.-kz(ikz+ipx*nx)**2*one_over_k2(ikz,ikx,iky)
              Pij(4)=-kx(ikx+ipy*ny)*ky(iky+ipz*nz)*one_over_k2(ikz,ikx,iky)
              Pij(5)=-ky(iky+ipz*nz)*kz(ikz+ipx*nx)*one_over_k2(ikz,ikx,iky)
              Pij(6)=-kz(ikz+ipx*nx)*kx(ikx+ipy*ny)*one_over_k2(ikz,ikx,iky)
            endif
!
!  compute e_T and e_X
!
            do j=1,3
            do i=1,3
              ij=ij_table(i,j)
              e_T(ij)=e1(i)*e1(j)-e2(i)*e2(j)
              e_X(ij)=e1(i)*e2(j)+e2(i)*e1(j)
            enddo
            enddo
!
!  possibility of swapping the sign
!
            if (lswitch_sign_e_X) then
              if (k3<0.) then
                e_X=-e_X
              elseif (k3==0.) then
                if (k2<0.) then
                  e_X=-e_X
                elseif (k2==0.) then
                  if (k1<0.) then
                    e_X=-e_X
                  endif
                endif
              endif
            endif
!
!  Traceless-tansverse projection:
!  Sij = (Pip*Pjq-.5*Pij*Ppq)*Tpq
!
            Sij_re=0.
            Sij_im=0.
            do j=1,3
            do i=1,j
            do q=1,3
            do p=1,3
              ij=ij_table(i,j)
              pq=ij_table(p,q)
              ip=ij_table(i,p)
              jq=ij_table(j,q)
              Sij_re(ij)=Sij_re(ij)+(Pij(ip)*Pij(jq)-.5*Pij(ij)*Pij(pq))*Tpq_re(ikz,ikx,iky,pq)
              Sij_im(ij)=Sij_im(ij)+(Pij(ip)*Pij(jq)-.5*Pij(ij)*Pij(pq))*Tpq_im(ikz,ikx,iky,pq)
            enddo
            enddo
            enddo
            enddo
!
!  compute S_T and S_X
!  AB: it looks like one could reuse part of Tpq_re and Tpq_im for S_T_re, ...,
!  AB: S_X_im to save memory
!
            do j=1,3
            do i=1,3
              ij=ij_table(i,j)
              S_T_re(ikz,ikx,iky)=S_T_re(ikz,ikx,iky)+.5*e_T(ij)*Sij_re(ij)
              S_T_im(ikz,ikx,iky)=S_T_im(ikz,ikx,iky)+.5*e_T(ij)*Sij_im(ij)
              S_X_re(ikz,ikx,iky)=S_X_re(ikz,ikx,iky)+.5*e_X(ij)*Sij_re(ij)
              S_X_im(ikz,ikx,iky)=S_X_im(ikz,ikx,iky)+.5*e_X(ij)*Sij_im(ij)
            enddo
            enddo
         
 ! Showing results for kz = 0, kz = 2 for testing purpose (Alberto Sayan)
          !if (k1==0..and.k2==0..and.abs(k3)==2.) then
          if (abs(k1)==2..and.abs(k2)==2..and.abs(k3)==0.) then
!
!  debug output (perhaps to be removed)
!
          if (ldebug_print) then
            print*,'PRINTING RESULTS FOR K = (+/-2, 0, 0)'
            print*,'AXEL k1,k2,k3=',k1,k2,k3
            print*,'AXEL e_1=',e1
            print*,'AXEL e_2=',e2
            print*,'AXEL e_T=',e_T
            print*,'AXEL e_X=',e_X
            print*,'AXEL S_X_re=',S_X_re(ikz,ikx,iky)
            print*,'AXEL S_X_im=',S_X_im(ikz,ikx,iky)
            print*,'AXEL S_T_re=',S_T_re(ikz,ikx,iky)
            print*,'AXEL S_T_im=',S_T_im(ikz,ikx,iky)
            print*,'AXEL Sij_re=',Sij_re
            print*,'AXEL Sij_im=',Sij_im
            print*,'AXEL Tpq_re=',Tpq_re(ikz,ikx,iky,:)
            print*,'AXEL Tpq_im=',Tpq_im(ikz,ikx,iky,:)
            print*,'AXEL Pij=',Pij
            if (k1==0..and.k2==0..and.k3==0.) then
              print*,'PRINTING RESULTS FOR K = (0, 0, 0)'
              print*,'AXEL k1,k2,k3=',k1,k2,k3
              print*,'AXEL e_T=',e_T
              print*,'AXEL e_X=',e_X
              print*,'AXEL S_X_re=',S_X_re(ikz,ikx,iky)
              print*,'AXEL S_X_im=',S_X_im(ikz,ikx,iky)
              print*,'AXEL S_T_re=',S_T_re(ikz,ikx,iky)
              print*,'AXEL S_T_im=',S_T_im(ikz,ikx,iky)
              print*,'AXEL Sij_re=',Sij_re
              print*,'AXEL Sij_im=',Sij_im
              print*,'AXEL Tpq_re=',Tpq_re(ikz,ikx,iky,:)
              print*,'AXEL Tpq_im=',Tpq_im(ikz,ikx,iky,:)
              print*,'AXEL Pij=',Pij
            endif
          endif
          endif
!
          enddo
        enddo
      enddo
!
!  debug output (perhaps to be removed)
!
      if (ldebug_print) then
        print*,'AXEL: k3Re=',S_T_re(1,:,1)
        print*,'AXEL: k3Im=',S_T_im(1,:,1)
        print*,'AXEL: k2Re=',S_X_re(1,:,1)
        print*,'AXEL: k2Im=',S_X_im(1,:,1)
      endif
!
!  back to real space
!
      call fourier_transform(S_T_re,S_T_im,linv=.true.)
      call fourier_transform(S_X_re,S_X_im,linv=.true.)
!--   call fft_xyz_parallel(S_T_re,S_T_im,linv=.true.)
!--   call fft_xyz_parallel(S_X_re,S_X_im,linv=.true.)
!
!  debug output (perhaps to be removed)
!
      if (ldebug_print) then
        print*,'AXEL: x3Re=',S_T_re(:,1,1)
        print*,'AXEL: x3Im=',S_T_im(:,1,1)
        print*,'AXEL: x2Re=',S_X_re(:,1,1)
        print*,'AXEL: x2Im=',S_X_im(:,1,1)
      endif
!
!  add (or set) corresponding stress
!
      f(l1:l2,m1:m2,n1:n2,iggT)=S_T_re
      f(l1:l2,m1:m2,n1:n2,iggX)=S_X_re
      !f(l1:l2,m1:m2,n1:n2,iggX)=S_X_im
!
!  For the time being, we keep the lswitch_sign_e_X
!  still as an option. Meaningful results are however
!  only found foro complex values of S_X.
!
      if (.not.lswitch_sign_e_X) then
        !f(l1:l2,m1:m2,n1:n2,iggX)=S_X_im
        f(l1:l2,m1:m2,n1:n2,iggX)=S_X_re
      endif
!
!  debug output (perhaps to be removed)
!
      if (ldebug_print) then
        print*,'PRINTING PHYSICAL SPACE S_T AND S_X'
        print*,'AXEL S_T_re=',S_T_re(:,1,1)
        print*,'AXEL S_X_re=',S_X_re(:,1,1)
        print*,'AXEL S_T_im=',S_T_im(:,1,1)
        print*,'AXEL S_X_im=',S_X_im(:,1,1)
      endif
!
!  Deallocate arrays.
!
      if (allocated(one_over_k2)) deallocate(one_over_k2)
      if (allocated(S_T_re)) deallocate(S_T_re)
      if (allocated(S_X_im)) deallocate(S_X_im)
      if (allocated(Tpq_re)) deallocate(Tpq_re)
      if (allocated(Tpq_im)) deallocate(Tpq_im)
      if (allocated(kx)) deallocate(kx)
      if (allocated(ky)) deallocate(ky)
      if (allocated(kz)) deallocate(kz)
!
    endsubroutine compute_gT_and_gX_from_gij
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
        idiag_g11pt=0; idiag_g22pt=0; idiag_g33pt=0
        idiag_g12pt=0; idiag_g23pt=0; idiag_g31pt=0
        idiag_ggTpt=0; idiag_ggXpt=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'g11pt',idiag_g11pt)
        call parse_name(iname,cname(iname),cform(iname),'g22pt',idiag_g22pt)
        call parse_name(iname,cname(iname),cform(iname),'g33pt',idiag_g33pt)
        call parse_name(iname,cname(iname),cform(iname),'g12pt',idiag_g12pt)
        call parse_name(iname,cname(iname),cform(iname),'g23pt',idiag_g23pt)
        call parse_name(iname,cname(iname),cform(iname),'g31pt',idiag_g31pt)
        call parse_name(iname,cname(iname),cform(iname),'ggTpt',idiag_ggTpt)
        call parse_name(iname,cname(iname),cform(iname),'ggXpt',idiag_ggXpt)
      enddo
!!
!!!  write column where which magnetic variable is stored
!!      if (lwr) then
!!        write(3,*) 'i_SPECIAL_DIAGNOSTIC=',i_SPECIAL_DIAGNOSTIC
!!      endif
!!
    endsubroutine rprint_special
!***********************************************************************
    subroutine get_slices_special(f,slices)
!
!  Write slices for animation of Special variables.
!
!  26-jun-06/tony: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      integer :: inamev
!
!  Loop over slices
!
      select case (trim(slices%name))
!
!  ggT
!
        case ('ggT')
          slices%yz=f(ix_loc,m1:m2,n1:n2,iggT)
          slices%xz=f(l1:l2,iy_loc,n1:n2,iggT)
          slices%xy=f(l1:l2,m1:m2,iz_loc,iggT)
          slices%xy2=f(l1:l2,m1:m2,iz2_loc,iggT)
          slices%ready = .true.
!
!  ggX
!
        case ('ggX')
          slices%yz=f(ix_loc,m1:m2,n1:n2,iggX)
          slices%xz=f(l1:l2,iy_loc,n1:n2,iggX)
          slices%xy=f(l1:l2,m1:m2,iz_loc,iggX)
          slices%xy2=f(l1:l2,m1:m2,iz2_loc,iggX)
          slices%ready = .true.
!
      endselect
!
    endsubroutine get_slices_special
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
