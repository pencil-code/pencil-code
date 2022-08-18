! $Id$
!
!  This module is the dummy for the SGS_hydro module
!  in which e.g the SGS force or heat is
!  computed.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lSGS_hydro = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 6
! COMMUNICATED AUXILIARIES 6
!
!! PENCILS PROVIDED fSGS(3)
! PENCILS PROVIDED SGS_heat
!
!***************************************************************
module SGS_hydro
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages, only: fatal_error
!
  implicit none
!
  include 'SGS_hydro.h'

!
  integer :: itauSGSRey=0, itauSGSMax=0, iSGS_heat=0, iSGS_force=0
!
  real :: cReyStress=1., cMaxStress=1. 
  logical :: lSGS_heat_as_aux=.false., lSGS_forc_as_aux=.false.

  namelist /SGS_hydro_run_pars/ cReyStress, cMaxStress, &
                                lSGS_heat_as_aux, lSGS_forc_as_aux

  integer :: idiag_fSGSm=0, idiag_fSGSrmsx=0

  contains
!***********************************************************************
    subroutine register_SGS_hydro
!
!  19-nov-02/tony: coded
!
      use Messages, only: svn_id
      use FArrayManager, only: farray_register_auxiliary
!
!  Identify version number.
!
      if (lroot) call svn_id( &
           "$Id$")
!
!  Register SGS Reynolds stress as auxilliary variable.
!
      call farray_register_auxiliary('tauSGSRey',itauSGSRey,array=6,communicated=.true.)
print*, 'tauSGSRey=', itauSGSRey
!
!  Register SGS Maxwell stress as auxilliary variable.
!
      if (lmagnetic) &
        call farray_register_auxiliary('tauSGSMax',itauSGSMax,array=6,communicated=.true.)
print*, 'tauSGSMax=', itauSGSMax
!
!  Register an extra aux slot for dissipation rate if requested (so
!  SGS_heat is written to snapshots and can be easily analyzed later).
!
      if (lSGS_heat_as_aux) then
        call farray_register_auxiliary('SGS_heat',iSGS_heat)
        aux_var(aux_count)=',SGS_heat'
        if (naux+naux_com <  maux+maux_com) aux_var(aux_count)=trim(aux_var(aux_count))//' $'
        aux_count=aux_count+1
      endif
!
!  Register an 3 extra aux slot for SGS force (accelaration) if requested (so
!  visc_forc is written to snapshots and can be easily analyzed later).
!
      if (lSGS_forc_as_aux) then
        call farray_register_auxiliary('SGS_forc',iSGS_force,vector=3)
        aux_var(aux_count)=',SGS_forc'
        if (naux+naux_com <  maux+maux_com) aux_var(aux_count)=trim(aux_var(aux_count))//' $'
        aux_count=aux_count+3
      endif
!
    endsubroutine register_SGS_hydro
!***********************************************************************
    subroutine initialize_SGS_hydro
!
      if (.not.ldensity) &
        call fatal_error('initialize_SGS_hydro','density needed')

    endsubroutine initialize_SGS_hydro
!***********************************************************************
    subroutine read_SGS_hydro_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=SGS_hydro_run_pars, IOSTAT=iostat)
!
    endsubroutine read_SGS_hydro_run_pars
!***********************************************************************
    subroutine write_SGS_hydro_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=SGS_hydro_run_pars)
!
    endsubroutine write_SGS_hydro_run_pars
!***********************************************************************
    subroutine rprint_SGS_hydro(lreset,lwrite)
!
!  Writes ishock to index.pro file
!
!  24-nov-03/tony: adapted from rprint_ionization
!
      use Diagnostics, only: parse_name
!
      logical :: lreset
      logical, intent(in), optional :: lwrite
      integer :: iname,inamex,inamez,ixy,irz
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        !idiag_dtnu=0; idiag_nu_LES=0; idiag_Sij2m=0
        !idiag_visc_heatm=0; 
        idiag_fSGSm=0
        !idiag_fSGSmsx=0
        !idiag_fviscmz=0; idiag_fviscmx=0; idiag_fviscmxy=0
        !idiag_fviscymxy=0
        !idiag_fviscsmmz=0; idiag_fviscsmmxy=0; idiag_ufviscm=0
        !idiag_fviscmax=0; idiag_fviscmin=0; idiag_fviscrsphmphi=0
        !idiag_viscforcezmz=0; idiag_viscforcezupmz=0; idiag_viscforcezdownmz=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if (lroot.and.ip<1400) print*,'rprint_SGS_hydro: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'fSGSm',idiag_fSGSm)
      enddo
!     
      call keep_compiler_quiet(lreset,lwrite)
!
    endsubroutine rprint_SGS_hydro
!***********************************************************************
    subroutine pencil_criteria_SGS_hydro
!
!  All pencils that the Viscosity module depends on are specified here.
!
!  20-11-04/anders: coded
!
      lpenc_requested(i_graddivu)=.true.
      lpenc_requested(i_sij2)=.true.
      lpenc_requested(i_ss12)=.true.

      if (lmagnetic) then

      endif

    endsubroutine pencil_criteria_SGS_hydro
!***********************************************************************
    subroutine pencil_interdep_SGS_hydro(lpencil_in)
!
!  Interdependency among pencils from the Viscosity module is specified here.
!
!  20-11-04/anders: coded
!
      logical, dimension (npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_SGS_hydro
!***********************************************************************
    subroutine calc_SGS_hydro_force(f,df,p)
!
!  Calculate SGS force and store in pencil.
!  Most basic pencils should come first, as others may depend on them.
!
!  20-11-04/anders: coded
!
      use Sub, only: div, div_other, dot2
      use Diagnostics, only: sum_mn_name

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df

      type (pencil_case) :: p
!
      intent(inout) :: f,p
      intent(out) :: df

      real, dimension(nx,3) :: tmp
      integer :: j
      real, dimension(nx) :: dtp
!
!  Divergence of tensor tau.
!  Correct only for Cartesian coordinates.
!
      call div(f,itauSGSRey,tmp(:,1))
      call div_other(f(:,:,:,(/itauSGSRey+1,itauSGSRey+3,itauSGSRey+4/)),tmp(:,2))
      call div_other(f(:,:,:,(/itauSGSRey+2,itauSGSRey+4,itauSGSRey+5/)),tmp(:,3))

      if (lSGS_forc_as_aux) f(l1:l2,m,n,iSGS_force:iSGS_force+2)=-tmp

      df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) - tmp
!
!  define SGS_heat
!
      if (lpencil(i_SGS_heat)) p%SGS_heat=0.0
!
      if (idiag_fSGSm/=0) then
        call dot2(tmp,dtp)
        call sum_mn_name(dtp,idiag_fSGSm,lsqrt=.true.)
      endif

    endsubroutine calc_SGS_hydro_force
!***********************************************************************
    subroutine SGS_hydro_after_boundary(f,p)

      use Boundcond, only: update_ghosts

      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p

      integer :: mm,nn
      real, dimension(nx) :: E_SGS,uij2
!
!  Requires communication to be finished!
!
!print*, 'in SGS_hydro_after_boundary'

      do nn=n1,n2; do mm=m1,m2

        E_SGS=cReyStress*max(dx,dy,dz)**2*exp(p%lnrho)*sqrt(2.*p%sij2)
        uij2=sum(sum(p%uij(:,:,:)**2,2),2)

        f(l1:l2,mm,nn,itauSGSRey  ) = E_SGS*(sum(p%uij(:,1,:)**2,2)          /uij2 - 1./3.)   ! tau(1,1)
        f(l1:l2,mm,nn,itauSGSRey+1) = E_SGS*(sum(p%uij(:,1,:)*p%uij(:,2,:),2)/uij2        )   ! tau(1,2)
        f(l1:l2,mm,nn,itauSGSRey+2) = E_SGS*(sum(p%uij(:,1,:)*p%uij(:,3,:),2)/uij2        )   ! tau(1,3)
        f(l1:l2,mm,nn,itauSGSRey+3) = E_SGS*(sum(p%uij(:,2,:)**2,2)          /uij2 - 1./3.)   ! tau(2,2)
        f(l1:l2,mm,nn,itauSGSRey+4) = E_SGS*(sum(p%uij(:,2,:)*p%uij(:,3,:),2)/uij2        )   ! tau(2,3)
        f(l1:l2,mm,nn,itauSGSRey+5) = E_SGS*(sum(p%uij(:,3,:)**2,2)          /uij2 - 1./3.)   ! tau(3,3)

      enddo; enddo
!
      call update_ghosts(f,itauSGSRey,itauSGSRey+5)

    endsubroutine SGS_hydro_after_boundary
!***********************************************************************
    subroutine SGS_hydro_before_boundary(f)

      real, dimension (mx,my,mz,mfarray) :: f

      call keep_compiler_quiet(f)
!
    endsubroutine SGS_hydro_before_boundary
!***********************************************************************
    subroutine calc_SGS_hydro_heat(df,p,Hmax)
!
      real, dimension (mx,my,mz,mvar)    :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: Hmax
!
      intent(in) :: df,p,Hmax
!
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
      call keep_compiler_quiet(Hmax)
!
    endsubroutine calc_SGS_hydro_heat
!***********************************************************************
endmodule SGS_hydro
