! $Id$
!
!  This modules solves mean-field contributions to both the
!  induction and the momentum equations.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lmagn_mf_demfdt = .true.
!
! MVAR CONTRIBUTION 3
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED
!
!***************************************************************
module Magnetic_meanfield_demfdt
!
  use Cdata
  use Cparam
  use Messages, only: fatal_error,inevitably_fatal_error,warning,svn_id,timing
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'meanfield_demfdt.h'
!
!  Declare indices
!
  integer :: iemf,iemfx,iemfy,iemfz
!
! Input parameters
!
  character (len=labellen), dimension(ninit) :: initemf='nothing'
  real, dimension (ninit) :: amplemf=0.
  real :: tau_emf=0., tau1_emf=0., eta_emf=0.
!
  namelist /magn_mf_demfdt_init_pars/ &
      tau_emf, tau1_emf, eta_emf, &
      initemf
!
! Run parameters
!
  namelist /magn_mf_demfdt_run_pars/ &
      tau_emf, tau1_emf, eta_emf
!
! Diagnostic variables (need to be consistent with reset list below)
! 
  integer :: idiag_EMFrms=0     ! DIAG_DOC: $(\left<{\cal E}\right>)_{\rm rms}$
  integer :: idiag_EMFmax=0     ! DIAG_DOC: $\max(\left<{\cal E}\right>)$
  integer :: idiag_EMFmin=0     ! DIAG_DOC: $\min(\left<{\cal E}\right>)$
!
  contains
!***********************************************************************
    subroutine register_magn_mf_demfdt()
!
!  Initialise additional variables
!
!  6-jan-11/axel: adapted from magnetic
!
      use FArrayManager, only: farray_register_pde
!
      call farray_register_pde('emf',iemf,vector=3)
      iemfx=iemf; iemfy=iemf+1; iemfz=iemf+2
!
!  Identify version number.
!
      if (lroot) call svn_id( &
          "$Id$")
!
!  Writing files for use with IDL
!
      if (lroot) then
        if (maux == 0) then
          if (nvar < mvar) write(4,*) ',emf $'
          if (nvar == mvar) write(4,*) ',emf'
        else
          write(4,*) ',emf $'
        endif
        write(15,*) 'emf = fltarr(mx,my,mz,3)*one'
      endif
!
    endsubroutine register_magn_mf_demfdt
!***********************************************************************
    subroutine initialize_magn_mf_demfdt(f,lstarting)
!
!  Perform any post-parameter-read initialization
!
!   6-jan-11/axel: adapted from meanfield.f90
!
      use BorderProfiles, only: request_border_driving
      use FArrayManager
      use SharedVariables, only: put_shared_variable,get_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
      integer :: i,ierr
!
!  Precalculate tau1_emf if tau_emf is given instead
!
      if (tau_emf/=0.) then
        tau1_emf=1./tau_emf
      endif
      if (lroot) print*,'initialize_magn_mf_demfdt: tau1_emf=',tau1_emf
!
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_magn_mf_demfdt
!***********************************************************************
    subroutine init_aa_mf_demfdt(f)
!
!  Initialise mean-field related magnetic field; called from magnetic.f90
!  At the moment, no own initial conditions are allowed, but we need this
!  to call secondary modules
!
!   6-jan-2011/axel: adapted from magnetic
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      integer :: j
!
!  initial conditions for secondary mean-field modules
!
      do j=1,ninit
        select case (initemf(j))
        case ('nothing'); if (lroot .and. j==1) print*,'initemf: nothing'
        case ('zero', '0'); f(:,:,:,iemfx:iemfz)=0.
        case default
!
!  Catch unknown values.
!
          call fatal_error('initemf', &
              'initemf value "' // trim(initemf(j)) // '" not recognised')
!
        endselect
      enddo
!
    endsubroutine init_aa_mf_demfdt
!***********************************************************************
    subroutine pencil_criteria_magn_mf_demfdt()
!
!  Pencils that are connected with the time advance of the EMF
!
!   6-jan-11/axel: adapted from magn_mf
!
      use Mpicomm, only: stop_it
!
    endsubroutine pencil_criteria_magn_mf_demfdt
!***********************************************************************
    subroutine pencil_interdep_magn_mf_demfdt(lpencil_in)
!
!  Interdependency among pencils from the Magnetic module is specified here.
!
!   6-jan-10/axel: adapted from magn_mf
!
      logical, dimension(npencils) :: lpencil_in
!
    endsubroutine pencil_interdep_magn_mf_demfdt
!***********************************************************************
    subroutine calc_pencils_magn_mf_demfdt(f,p)
!
!  Calculate secondary magnetic mean-field pencils.
!
!   6-jan-11/axel: adapted from magn_mf
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(inout) :: f,p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_magn_mf_demfdt
!***********************************************************************
    subroutine demf_dt_meanfield(f,df,p)
!
!  Solve dEMF/dt = EMF_old/tau + eta_emf*del2(EMF)
!
!   6-jan-11/axel: coded
!
      use Diagnostics
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx,3) :: emf, del2emf
!
      intent(inout)  :: f,p
      intent(inout)  :: df
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print*,'demf_dt_meanfield: SOLVE'
!
!  First part (without spatial diffusion): dEMF/dt = (alp*B-etat*J-EMF)/tau.
!  Note that the first part, alp*B-etat*J, is here obtained as p%mf_EMF.
!
      emf=f(l1:l2,m,n,iemfx:iemfz)
      df(l1:l2,m,n,iemfx:iemfz)=df(l1:l2,m,n,iemfx:iemfz)+tau1_emf*(p%mf_EMF-emf)
!
!  Spatial diffusion part, if eta_emf/=0.
!
      if (eta_emf/=0.) then
        call del2v(f,iemf,del2emf)
        df(l1:l2,m,n,iemfx:iemfz)=df(l1:l2,m,n,iemfx:iemfz)+eta_emf*del2emf
      endif
!
!  Apply EMF to dA/dt (was turned off in meanfield.f90 because lmagn_mf_demfdt=T
!
      df(l1:l2,m,n,iax:iaz)=df(l1:l2,m,n,iax:iaz)+f(l1:l2,m,n,iemfx:iemfz)
!
!  Calculate diagnostic quantities.
!  Diagnostic output for mean field dynamos.
!
      if (ldiagnos) then
        if (idiag_EMFrms/=0) call sum_mn_name(+emf,idiag_EMFrms)
        if (idiag_EMFmax/=0) call max_mn_name(+emf,idiag_EMFmax)
        if (idiag_EMFmin/=0) call max_mn_name(-emf,idiag_EMFmax,lneg=.true.)
      endif
!
    endsubroutine demf_dt_meanfield
!***********************************************************************
    subroutine read_magn_mf_demfdt_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=magn_mf_demfdt_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=magn_mf_demfdt_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_magn_mf_demfdt_init_pars
!***********************************************************************
    subroutine write_magn_mf_demfdt_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=magn_mf_demfdt_init_pars)
!
    endsubroutine write_magn_mf_demfdt_init_pars
!***********************************************************************
    subroutine read_magn_mf_demfdt_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=magn_mf_demfdt_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=magn_mf_demfdt_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_magn_mf_demfdt_run_pars
!***********************************************************************
    subroutine write_magn_mf_demfdt_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=magn_mf_demfdt_run_pars)
!
    endsubroutine write_magn_mf_demfdt_run_pars
!***********************************************************************
    subroutine rprint_magn_mf_demfdt(lreset,lwrite)
!
!  Reads and registers print parameters relevant for magnetic fields.
!
!   6-jan-11/axel: adapted from rprint_magn_mf
!
      use Diagnostics
!
      integer :: iname
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  Reset everything in case of RELOAD.
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_EMFrms=0; idiag_EMFmax=0; idiag_EMFmin=0
      endif
!
!  Check for those quantities that we want to evaluate online.
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'EMFrms',idiag_EMFrms)
        call parse_name(iname,cname(iname),cform(iname),'EMFmax',idiag_EMFmax)
        call parse_name(iname,cname(iname),cform(iname),'EMFmin',idiag_EMFmin)
      enddo
!
    endsubroutine rprint_magn_mf_demfdt
!***********************************************************************
endmodule Magnetic_meanfield_demfdt
