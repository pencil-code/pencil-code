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
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!! PENCILS PROVIDED infl_phi; infl_dphi; gphi(3)
!
!** AUTOMATIC REFERENCE-LINK.TEX GENERATION ********************
! Declare relevant citations from pencil-code/doc/citations/ref.bib for this module.
! The entries are taken from pencil-code/doc/citations/notes.tex
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
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages, only: svn_id, fatal_error
!
  implicit none
!
  include '../special.h'
!
!
! Declare index of new variables in f array (if any).
!
  logical :: luschur2_as_aux=.false., lbschur2_as_aux=.false.
  integer :: iuschur2, iuschur2_SH, iuschur2_RR, iuschur2_EL
  integer :: ibschur2, ibschur2_SH, ibschur2_RR, ibschur2_EL
!
  namelist /special_init_pars/ &
      luschur2_as_aux, lbschur2_as_aux
!
  logical :: luij_schur=.false., lbij_schur=.false.
  real, dimension (nx) :: uSH2, uRR2, uEL2
  real, dimension (nx) :: bSH2, bRR2, bEL2
!
  namelist /special_run_pars/ &
      luij_schur, lbij_schur
!
! Diagnostic variables (needs to be consistent with reset list below).
!
  integer :: idiag_uSH2=0      ! DIAG_DOC: $u^\mathrm{SH}$
  integer :: idiag_uRR2=0      ! DIAG_DOC: $u^\mathrm{RR}$
  integer :: idiag_uEL2=0      ! DIAG_DOC: $u^\mathrm{EL}$
  integer :: idiag_bSH2=0      ! DIAG_DOC: $b^\mathrm{SH}$
  integer :: idiag_bRR2=0      ! DIAG_DOC: $b^\mathrm{RR}$
  integer :: idiag_bEL2=0      ! DIAG_DOC: $b^\mathrm{EL}$
!
  contains
!================== External SELECT function for DGEES ==================
logical function selct(WR, WI)
  implicit none
  real, intent(in) :: WR, WI
  ! Return .TRUE. for eigenvalues to move to the TOP-LEFT block.
  ! We choose REAL eigenvalues => complex pair (if any) goes bottom-right.
  selct = (WI == 0.)
end function selct
!****************************************************************************
    subroutine register_special
!
!  Set up indices for variables in special modules.
!
!  6-oct-03/tony: coded
!
      use FArrayManager
      use Sub, only: register_report_aux
      use SharedVariables, only: put_shared_variable
!
      if (lroot) call svn_id( &
           "$Id$")
!
      if (luschur2_as_aux) &
        call register_report_aux('uschur2', iuschur2, iuschur2_SH, iuschur2_RR, iuschur2_EL)
!
      if (lbschur2_as_aux) &
        call register_report_aux('bschur2', ibschur2, ibschur2_SH, ibschur2_RR, ibschur2_EL)
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  Called after reading parameters, but before the time loop.
!
!  06-oct-03/tony: coded
!
      use SharedVariables, only: get_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f)
!
!  initialise special condition; called from start.f90
!  06-oct-2003/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j
!
      intent(inout) :: f
!
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special
!
!  Dummy
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
!  Calculate Special pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  24-nov-04/tony: coded
!
      use Sub, only: grad
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      real, dimension (3,3) :: matA, matV_SH, matV_RR, matV_EL
      real, allocatable :: matV(:,:), matQ(:,:)
      real :: matA2=0.
      integer :: l, nnn=3
!
      intent(inout) :: f
      intent(inout) :: p
!
!  Possibility of applying triple decomposition of uij
!
      if (luij_schur) then
        do l=1,nx
          matA=p%uij(l,:,:)
          matA2=sum(matA**2)
          allocate(matV(nnn,nnn), matQ(nnn,nnn))
          if (matA2/=0.) then
            call schur_standardized(matA, matV, matQ)
            call schur_decompose(matV, matV_SH, matV_RR, matV_EL)
            uSH2(l)=sum(matV_SH**2)
            uRR2(l)=sum(matV_RR**2)
            uEL2(l)=sum(matV_EL**2)
          else
            uSH2(l)=0.
            uRR2(l)=0.
            uEL2(l)=0.
          endif
!
!  Possibility of bSH2, bRR2, and bEL2 as auxiliary arrays
!
          if (luschur2_as_aux) then
            f(l1:l2,m,n,iuschur2_SH)=uSH2
            f(l1:l2,m,n,iuschur2_RR)=uRR2
            f(l1:l2,m,n,iuschur2_EL)=uEL2
          endif
!
          deallocate(matV, matQ)
        enddo
      endif
!
!  Possibility of applying triple decomposition of bij
!
      if (lbij_schur) then
        do l=1,nx
          matA=p%bij(l,:,:)
          matA2=sum(matA**2)
          allocate(matV(nnn,nnn), matQ(nnn,nnn))
          if (matA2/=0.) then
            call schur_standardized(matA,matV,matQ)
            call schur_decompose(matV, matV_SH, matV_RR, matV_EL)
            bSH2(l)=sum(matV_SH**2)
            bRR2(l)=sum(matV_RR**2)
            bEL2(l)=sum(matV_EL**2)
          else
            bSH2(l)=0.
            bRR2(l)=0.
            bEL2(l)=0.
          endif
!
!  Possibility of bSH2, bRR2, and bEL2 as auxiliary arrays
!
          if (lbschur2_as_aux) then
            f(l1:l2,m,n,ibschur2_SH)=bSH2
            f(l1:l2,m,n,ibschur2_RR)=bRR2
            f(l1:l2,m,n,ibschur2_EL)=bEL2
          endif
!
          deallocate(matV, matQ)
        enddo
      endif
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine schur_decompose(matV, matV_SH, matV_RR, matV_EL)
      real, dimension (3,3) :: matV_SH, matV_RR, matV_EL, matV_rest, matV_rest_trans
      real, allocatable :: matV(:,:)
      integer :: i, j, nnn=3
!
!  Shear:
!
      matV_SH(1,1)=0.
      matV_SH(1,2)=matV(1,2)
      matV_SH(1,3)=matV(1,3)
!
      matV_SH(2,1)=0.
      matV_SH(2,2)=0.
      matV_SH(2,3)=sign(1.,matV(2,3))*max(abs(matV(2,3))-abs(matV(3,2)),0.)
!
      matV_SH(3,1)=0.
      matV_SH(3,2)=sign(1.,matV(3,2))*max(abs(matV(3,2))-abs(matV(2,3)),0.)
      matV_SH(3,3)=0.
!
!  Rest:
!
      matV_rest=matV-matV_SH
!
!  Transposed matrix
!
      do j=1,3
        do i=1,3
          matV_rest_trans(i,j)=matV_rest(j,i)
        enddo
      enddo
!
!  Antisymmetric and symmetric parts
!
      matV_RR=.5*(matV_rest-matV_rest_trans)
      matV_EL=.5*(matV_rest+matV_rest_trans)
!
    endsubroutine schur_decompose
!***********************************************************************
  subroutine print_mat(M)
    real            , intent(in) :: M(:,:)
    integer :: r, c
    do r = 1, size(M,1)
       write(*,'(100(1X,F12.6))') ( M(r,c), c=1,size(M,2) )
    end do
  end subroutine print_mat
!***********************************************************************
    subroutine schur_standardized(A_input, T, VS)
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 307)
  integer            :: nnn, lda, ldvs, info, sdim, lwork
  integer            :: i, j
  real               :: tol
  real            :: A_input(3,3)
  real            , allocatable :: A(:,:), T(:,:), VS(:,:), WR(:), WI(:), WORK(:)
  logical,    allocatable :: BWORK(:)

        real             :: a11, a12, a21, a22
        real             :: rt1r, rt1i, rt2r, rt2i, cs, sn
        real             :: R11, R12, R21, R22
        real            , allocatable :: TMP(:,:), TMP2(:,:)

  external dgees, dlanv2            ! LAPACK externals

  ! ---- Example matrix (3x3). Replace with your data. ----
  nnn = 3
  lda = nnn
  ldvs = nnn
  tol = 1.0d-12

  allocate(A(nnn,nnn), WR(nnn), WI(nnn))
  allocate(BWORK(nnn))

  A = reshape(A_input, shape(A_input), order=(/2,1/) )

  T = A

  ! ---- Workspace query for DGEES ----
  lwork = -1
  allocate(WORK(1))
  call dgees('V','S', selct, nnn, T, lda, sdim, WR, WI, VS, ldvs, WORK, lwork, BWORK, info)
  if (info .ne. 0) then
     print *, 'DGEES(work query) failed, INFO=', info
     stop 1
  end if
  lwork = int(WORK(1))
  deallocate(WORK)
  allocate(WORK(lwork))

  ! ---- Real Schur with sorting: real eigenvalues first (complex pair -> bottom-right) ----
  T = A
  call dgees('V','S', selct, nnn, T, lda, sdim, WR, WI, VS, ldvs, WORK, lwork, BWORK, info)
  if (info .ne. 0) then
     print *, 'DGEES failed, INFO=', info
     stop 1
  end if

  ! ---- Standardize the lower-right 2x2 block if present ----
  if (nnn >= 2) then
     i = nnn-1
     j = nnn
     if (abs(T(j,i)) > tol) then
        ! We have a 2x2 block at (i:i+1, i:i+1)

        a11 = T(i,i);   a12 = T(i,j)
        a21 = T(j,i);   a22 = T(j,j)

        ! DLANV2 computes rotation (cs,sn) that brings the 2x2 to standard real Schur form
        call dlanv2(a11, a12, a21, a22, rt1r, rt1i, rt2r, rt2i, cs, sn)

        ! Build R = [ cs  sn; -sn  cs ]
        R11 = cs;  R12 = sn
        R21 = -sn; R22 = cs

        ! Apply similarity on T: T := (I ⊕ R^T) * T * (I ⊕ R)
        allocate(TMP(nnn,2), TMP2(2,nnn))

        ! Right-multiply columns i:i+1 by R
        TMP(:,1) = T(:,i)*R11 + T(:,j)*R21
        TMP(:,2) = T(:,i)*R12 + T(:,j)*R22
        T(:,i) = TMP(:,1)
        T(:,j) = TMP(:,2)

        ! Left-multiply rows i:i+1 by R^T
        TMP2(1,:) =  R11*T(i,:) + R12*T(j,:)
        TMP2(2,:) =  R21*T(i,:) + R22*T(j,:)
        T(i,:) = TMP2(1,:)
        T(j,:) = TMP2(2,:)

        ! Update VS (Q): VS := VS * (I ⊕ R)
        TMP(:,1) = VS(:,i)*R11 + VS(:,j)*R21
        TMP(:,2) = VS(:,i)*R12 + VS(:,j)*R22
        VS(:,i) = TMP(:,1)
        VS(:,j) = TMP(:,2)

        deallocate(TMP, TMP2)

        ! Enforce upper-right >= 0 (beta >= 0) by a column/row sign flip if needed
        if (T(i,j) < 0.0d0) then
           ! Multiply by F = diag(1, -1): columns
           T(:,j)  = -T(:,j)
           VS(:,j) = -VS(:,j)
           ! and rows
           T(j,:)  = -T(j,:)
        end if

        ! Clean tiny roundoff
        if (abs(T(j,i)) < tol) T(j,i) = 0.0d0
        if (abs(T(i,j) + T(j,i)) < tol) T(i,j) = -T(j,i)
     end if
  end if

  ! ---- Print results ----
  !print *, 'T (standardized Schur form):'
  !call print_mat(T)
  !print *, 'Q (VS from DGEES):'
  !call print_mat(VS)

    endsubroutine schur_standardized
!***********************************************************************
    subroutine dspecial_dt(f,df,p)
!
!  The entire module could be renamed to Klein-Gordon or Scalar field equation.
!  Calculate right hand side of ONE OR MORE extra coupled PDEs
!  along the 'current' Pencil, i.e. f(l1:l2,m,n) where
!  m,n are global variables looped over in equ.f90
!
!  Due to the multi-step Runge Kutta timestepping used one MUST always
!  add to the present contents of the df array.  NEVER reset it to zero.
!
!  Several precalculated Pencils of information are passed for
!  efficiency.
!
!   6-oct-03/tony: coded
!   2-nov-21/axel: first set of equations coded
!
      use Diagnostics, only: sum_mn_name
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in) :: f,p
      intent(inout) :: df
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dspecial_dt'
!
!  Diagnostics
!
      if (ldiagnos) then
        call sum_mn_name(uSH2,idiag_uSH2)
        call sum_mn_name(uRR2,idiag_uRR2)
        call sum_mn_name(uEL2,idiag_uEL2)
        call sum_mn_name(bSH2,idiag_bSH2)
        call sum_mn_name(bRR2,idiag_bRR2)
        call sum_mn_name(bEL2,idiag_bEL2)
      endif
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine dspecial_dt_ode
!
      use Diagnostics, only: save_name
      use SharedVariables, only: get_shared_variable
!
    endsubroutine dspecial_dt_ode
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
    subroutine rprint_special(lreset,lwrite)
!
!  Reads and registers print parameters relevant to special.
!
!  06-oct-03/tony: coded
!
      use Diagnostics, only: parse_name
!
      integer :: iname
      logical :: lreset,lwrite
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_uSH2=0; idiag_uRR2=0; idiag_uEL2=0
        idiag_bSH2=0; idiag_bRR2=0; idiag_bEL2=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'uSH2',idiag_uSH2)
        call parse_name(iname,cname(iname),cform(iname),'uRR2',idiag_uRR2)
        call parse_name(iname,cname(iname),cform(iname),'uEL2',idiag_uEL2)
        call parse_name(iname,cname(iname),cform(iname),'bSH2',idiag_bSH2)
        call parse_name(iname,cname(iname),cform(iname),'bRR2',idiag_bRR2)
        call parse_name(iname,cname(iname),cform(iname),'bEL2',idiag_bEL2)
      enddo
!
    endsubroutine rprint_special
!***********************************************************************
    subroutine special_after_boundary(f)
!
!  Possibility to modify the f array after the boundaries are
!  communicated.
!
!  16-apr-25/axel: coded
!
!     use Mpicomm, only: mpireduce_sum, mpiallreduce_sum, mpibcast_real
!     use Sub, only: dot2_mn, grad, curl, dot_mn
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
!
    endsubroutine special_after_boundary
!********************************************************************
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING        *************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../special_dummies.inc'
!***********************************************************************
endmodule Special
