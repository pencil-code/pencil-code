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
! PENCILS PROVIDED uSH2; uRR2; uEL2; uRRm; uELm; bSH2; bRR2; bEL2; bRRm; bELm
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
  logical :: luschurm_as_aux=.false., lbschurm_as_aux=.false.
  integer :: iuschur2, iuschur2_SH, iuschur2_RR, iuschur2_EL, iuschurm_RR, iuschurm_EL
  integer :: ibschur2, ibschur2_SH, ibschur2_RR, ibschur2_EL, ibschurm_RR, ibschurm_EL
!
  logical :: luschur_as_aux=.false., lbschur_as_aux=.false.
  integer :: iuschur_SH, iuschur_RR, iuschur_EL
  integer :: ibschur_SH, ibschur_RR, ibschur_EL
!
  namelist /special_init_pars/ &
      luschur2_as_aux, lbschur2_as_aux, &
      luschurm_as_aux, lbschurm_as_aux, &
      luschur_as_aux, lbschur_as_aux
!
  logical :: luij_schur=.false., lbij_schur=.false., ldiagnos_always=.false.
  logical :: luschur_unprojected=.false., lbschur_unprojected=.false., lQ_schur_QT=.true.
!
  namelist /special_run_pars/ &
      luij_schur, lbij_schur, ldiagnos_always, &
      luschur_unprojected, lbschur_unprojected, lQ_schur_QT, &
      luschur_as_aux, lbschur_as_aux
!
! Diagnostic variables (needs to be consistent with reset list below).
!
  integer :: idiag_uSH2=0      ! DIAG_DOC: $\left<{u^\mathrm{SH}}^2\right>$
  integer :: idiag_uRR2=0      ! DIAG_DOC: $\left<{u^\mathrm{RR}}^2\right>$
  integer :: idiag_uEL2=0      ! DIAG_DOC: $\left<{u^\mathrm{EL}}^2\right>$
  integer :: idiag_uRRm=0      ! DIAG_DOC: $\left<2u^\mathrm{SH}u^\mathrm{RR}\right>$
  integer :: idiag_uELm=0      ! DIAG_DOC: $\left<2u^\mathrm{SH}u^\mathrm{EL}\right>$
  integer :: idiag_bSH2=0      ! DIAG_DOC: $\left<{b^\mathrm{SH}}^2\right>$
  integer :: idiag_bRR2=0      ! DIAG_DOC: $\left<{b^\mathrm{RR}}^2\right>$
  integer :: idiag_bEL2=0      ! DIAG_DOC: $\left<{b^\mathrm{EL}}^2\right>$
  integer :: idiag_bRRm=0      ! DIAG_DOC: $\left<2b^\mathrm{SH}b^\mathrm{RR}\right>$
  integer :: idiag_bELm=0      ! DIAG_DOC: $\left<2b^\mathrm{SH}b^\mathrm{EL}\right>$
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
      if (luschurm_as_aux) then
        call register_report_aux('uschurm_RR', iuschurm_RR)
        call register_report_aux('uschurm_EL', iuschurm_EL)
      endif
!
      if (lbschurm_as_aux) then
        call register_report_aux('bschurm_RR', ibschurm_RR)
        call register_report_aux('bschurm_EL', ibschurm_EL)
      endif
!
      if (luschur_as_aux) then
        call farray_register_auxiliary('uschur_SH', iuschur_SH, vector=9)
        call farray_register_auxiliary('uschur_RR', iuschur_RR, vector=9)
        call farray_register_auxiliary('uschur_EL', iuschur_EL, vector=9)
      endif
!
      if (lbschur_as_aux) then
        call farray_register_auxiliary('bschur_SH', ibschur_SH, vector=9)
        call farray_register_auxiliary('bschur_RR', ibschur_RR, vector=9)
        call farray_register_auxiliary('bschur_EL', ibschur_EL, vector=9)
      endif
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
      real, dimension (3,3) :: matA, matV_SH, matV_RR, matV_EL, SH, RR, EL
      real, allocatable :: matV(:,:), matQ(:,:)
      real :: matA2=0.
      integer :: i,j,kk,ll, ij, l, nnn=3
!
      intent(inout) :: f
      intent(inout) :: p
!
!  Possibility of applying triple decomposition of uij
!  Do only when output onto command line.
!
      if ((luij_schur .and. lout) .or. ldiagnos_always) then
        do l=1,nx
          matA=p%uij(l,:,:)
          matA2=sum(matA**2)
          allocate(matV(nnn,nnn), matQ(nnn,nnn))
          if (matA2/=0.) then
            call schur_standardized(matA, matV, matQ)
            call schur_decompose(matV, matV_SH, matV_RR, matV_EL)
            p%uSH2(l)=sum(matV_SH**2)
            p%uRR2(l)=sum(matV_RR**2)
            p%uEL2(l)=sum(matV_EL**2)
            p%uRRm(l)=2.*sum(matV_SH*matV_RR)
            p%uELm(l)=2.*sum(matV_SH*matV_EL)
          else
            if (ip<10) print*,'AXEL: warning for u, l=',l
            !p%uSH2(l)=0.
            !p%uRR2(l)=0.
            !p%uEL2(l)=0.
            !p%uRRm(l)=0.
            !p%uELm(l)=0.
          endif
!
!  Possibility of uSH, uRR, and uEL matrices as auxiliary arrays
!  Decomposition in the form GradU = Q (SH+RR+EL) Q^T, so
!  GradU_ij = Q_ik (SH+RR+EL)_kl QT_lj = Q_ik Q_jl (SH+RR+EL)_kl.
!  Setting lQ_schur_QT=.true. applies to GradU = Q^T (SH+RR+EL) Q.
!  Otherwise, we'd have GradU = Q^T (SH+RR+EL) Q.
!
          if (luschur_as_aux) then
            if (luschur_unprojected) then
              do i=1,3
              do j=1,3
                ij=3*(i-1)+(j-1)
                f(l1+l-1,m,n,iuschur_SH+ij)=matV_SH(i,j)
                f(l1+l-1,m,n,iuschur_RR+ij)=matV_RR(i,j)
                f(l1+l-1,m,n,iuschur_EL+ij)=matV_EL(i,j)
              enddo
              enddo
            else
              SH=0.
              RR=0.
              EL=0.
              do i=1,3
              do j=1,3
                do kk=1,3
                do ll=1,3
                  if (lQ_schur_QT) then
                    SH(i,j)=SH(i,j)+matQ(i,kk)*matQ(j,ll)*matV_SH(kk,ll)
                    RR(i,j)=RR(i,j)+matQ(i,kk)*matQ(j,ll)*matV_RR(kk,ll)
                    EL(i,j)=EL(i,j)+matQ(i,kk)*matQ(j,ll)*matV_EL(kk,ll)
                  else
                    SH(i,j)=SH(i,j)+matV_SH(kk,ll)*matQ(kk,i)*matQ(ll,j)
                    RR(i,j)=RR(i,j)+matV_RR(kk,ll)*matQ(kk,i)*matQ(ll,j)
                    EL(i,j)=EL(i,j)+matV_EL(kk,ll)*matQ(kk,i)*matQ(ll,j)
                  endif
                enddo
                enddo
                ij=3*(i-1)+(j-1)
                f(l1+l-1,m,n,iuschur_SH+ij)=SH(i,j)
                f(l1+l-1,m,n,iuschur_RR+ij)=RR(i,j)
                f(l1+l-1,m,n,iuschur_EL+ij)=EL(i,j)
              enddo
              enddo
            endif
          endif
!
          deallocate(matV, matQ)
        enddo
!
!  Possibility of uSH2, uRR2, and uEL2 as auxiliary arrays
!
        if (luschur2_as_aux) then
          f(l1:l2,m,n,iuschur2_SH)=p%uSH2
          f(l1:l2,m,n,iuschur2_RR)=p%uRR2
          f(l1:l2,m,n,iuschur2_EL)=p%uEL2
        endif
!
        if (luschurm_as_aux) then
          f(l1:l2,m,n,iuschurm_RR)=p%uRRm
          f(l1:l2,m,n,iuschurm_EL)=p%uELm
!print*,'AXEL m,p%uRRm(1:4)=',m,p%uRRm(1:4)
        endif
!     else
!       uSH2=0.
!       uRR2=0.
!       uEL2=0.
!       uRRm=0.
!       uELm=0.
      endif
!
!  Possibility of applying triple decomposition of bij
!  Do only when output onto command line.
!
      if ((lbij_schur .and. lout) .or. ldiagnos_always) then
        do l=1,nx
          matA=p%bij(l,:,:)
          matA2=sum(matA**2)
          allocate(matV(nnn,nnn), matQ(nnn,nnn))
          if (matA2/=0.) then
            call schur_standardized(matA,matV,matQ)
            call schur_decompose(matV, matV_SH, matV_RR, matV_EL)
            p%bSH2(l)=sum(matV_SH**2)
            p%bRR2(l)=sum(matV_RR**2)
            p%bEL2(l)=sum(matV_EL**2)
            p%bRRm(l)=2.*sum(matV_SH*matV_RR)
            p%bELm(l)=2.*sum(matV_SH*matV_EL)
          else
            print*,'AXEL: warning for b, l=',l
            !p%bSH2(l)=0.
            !p%bRR2(l)=0.
            !p%bEL2(l)=0.
            !p%bRRm(l)=0.
            !p%bELm(l)=0.
          endif
!
!  Possibility of bSH, bRR, and bEL matrices as auxiliary arrays
!
          if (lbschur_as_aux) then
            if (lbschur_unprojected) then
              do i=1,3
              do j=1,3
                ij=3*(i-1)+(j-1)
                f(l1+l-1,m,n,ibschur_SH+ij)=matV_SH(i,j)
                f(l1+l-1,m,n,ibschur_RR+ij)=matV_RR(i,j)
                f(l1+l-1,m,n,ibschur_EL+ij)=matV_EL(i,j)
              enddo
              enddo
            else
              SH=0.
              RR=0.
              EL=0.
              do i=1,3
              do j=1,3
                do kk=1,3
                do ll=1,3
                  if (lQ_schur_QT) then
                    SH(i,j)=SH(i,j)+matQ(i,kk)*matQ(j,ll)*matV_SH(kk,ll)
                    RR(i,j)=RR(i,j)+matQ(i,kk)*matQ(j,ll)*matV_RR(kk,ll)
                    EL(i,j)=EL(i,j)+matQ(i,kk)*matQ(j,ll)*matV_EL(kk,ll)
                  else
                    SH(i,j)=SH(i,j)+matV_SH(kk,ll)*matQ(kk,i)*matQ(ll,j)
                    RR(i,j)=RR(i,j)+matV_RR(kk,ll)*matQ(kk,i)*matQ(ll,j)
                    EL(i,j)=EL(i,j)+matV_EL(kk,ll)*matQ(kk,i)*matQ(ll,j)
                  endif
                enddo
                enddo
                ij=3*(i-1)+(j-1)
                f(l1+l-1,m,n,ibschur_SH+ij)=SH(i,j)
                f(l1+l-1,m,n,ibschur_RR+ij)=RR(i,j)
                f(l1+l-1,m,n,ibschur_EL+ij)=EL(i,j)
              enddo
              enddo
            endif
          endif
!
          deallocate(matV, matQ)
        enddo
!
!  Possibility of bSH2, bRR2, and bEL2 as auxiliary arrays
!  Use f(1,1,1,ibschur2_EL) to encode the correct update time.
!
        if (lbschur2_as_aux) then
          f(l1:l2,m,n,ibschur2_SH)=p%bSH2
          f(l1:l2,m,n,ibschur2_RR)=p%bRR2
          f(l1:l2,m,n,ibschur2_EL)=p%bEL2
          if (l==1 .and. m==m1 .and. n==n1 .and. ip<10) then
            print*,'AXEL: iproc,t=',iproc,t
            f(1,1,1,ibschur2_EL)=t
          endif
        endif
!
        if (lbschurm_as_aux) then
          f(l1:l2,m,n,ibschurm_RR)=p%bRRm
          f(l1:l2,m,n,ibschurm_EL)=p%bELm
        endif
!     else
!       bSH2=0.
!       bRR2=0.
!       bEL2=0.
!       bRRm=0.
!       bELm=0.
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
      if (ldiagnos .or. ldiagnos_always) then
        call sum_mn_name(p%uSH2,idiag_uSH2)
        call sum_mn_name(p%uRR2,idiag_uRR2)
        call sum_mn_name(p%uEL2,idiag_uEL2)
        call sum_mn_name(p%uRRm,idiag_uRRm)
        call sum_mn_name(p%uELm,idiag_uELm)
        call sum_mn_name(p%bSH2,idiag_bSH2)
        call sum_mn_name(p%bRR2,idiag_bRR2)
        call sum_mn_name(p%bEL2,idiag_bEL2)
        call sum_mn_name(p%bRRm,idiag_bRRm)
        call sum_mn_name(p%bELm,idiag_bELm)
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
        idiag_uSH2=0; idiag_uRR2=0; idiag_uEL2=0; idiag_uRRm=0; idiag_uELm=0
        idiag_bSH2=0; idiag_bRR2=0; idiag_bEL2=0; idiag_bRRm=0; idiag_bELm=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'uSH2',idiag_uSH2)
        call parse_name(iname,cname(iname),cform(iname),'uRR2',idiag_uRR2)
        call parse_name(iname,cname(iname),cform(iname),'uEL2',idiag_uEL2)
        call parse_name(iname,cname(iname),cform(iname),'uRRm',idiag_uRRm)
        call parse_name(iname,cname(iname),cform(iname),'uELm',idiag_uELm)
        call parse_name(iname,cname(iname),cform(iname),'bSH2',idiag_bSH2)
        call parse_name(iname,cname(iname),cform(iname),'bRR2',idiag_bRR2)
        call parse_name(iname,cname(iname),cform(iname),'bEL2',idiag_bEL2)
        call parse_name(iname,cname(iname),cform(iname),'bRRm',idiag_bRRm)
        call parse_name(iname,cname(iname),cform(iname),'bELm',idiag_bELm)
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
