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
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of special_dummies.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MAUX CONTRIBUTION 12
!
! PENCILS PROVIDED bb(3)
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
  use Messages, only: svn_id, fatal_error, warning
!
  implicit none
!
  include '../special.h'
!
! Declare index of new variables in f array (if any).
!
  character (len=labellen) :: initAAk='nothing'
  character (len=labellen) :: initEEk='nothing'
  character (len=labellen) :: cc_light='1'
  real :: alpha_inflation=0.
  real :: sigma=0.
  logical :: lbb_as_aux=.true., lEE_as_aux=.true.
  logical :: linflation=.false., ldebug_print=.false.
  real, dimension(3) :: AAkre, AAkim, EEkre, EEkim
  real :: c_light2=1.
  real, dimension (:,:,:,:), allocatable :: Tpq_re, Tpq_im
  real :: alpha2_inflation, kscale_factor
!
! input parameters
  namelist /special_init_pars/ &
    alpha_inflation, &
    initAAk, initEEk, &
    linflation, sigma
!
! run parameters
  namelist /special_run_pars/ &
    alpha_inflation, &
    ldebug_print, &
    cc_light, &
    linflation, sigma
!
! Diagnostic variables (needs to be consistent with reset list below).
!
  integer :: idiag_AA2m=0      ! DIAG_DOC: $\langle A^2\rangle$
  integer :: idiag_Akxpt=0     ! DIAG_DOC: $Akx^{pt}$
  integer :: idiag_Ekxpt=0     ! DIAG_DOC: $Ekx^{pt}$
!
  integer :: iAAk, iAAkim, iEEk, iEEkim
  integer :: iAkx, iAky, iAkz, iAkxim, iAkyim, iAkzim
  integer :: iEkx, iEky, iEkz, iEkxim, iEkyim, iEkzim
  integer, parameter :: nk=nxgrid/2
  type, public :: magpectra
    real, dimension(nk) :: mag   ,ele
    real, dimension(nk) :: maghel,elehel
  endtype magpectra

  type(magpectra) :: spectra

  contains
!***********************************************************************
    subroutine register_special
!
!  Set up indices for variables in special modules.
!  Need 2x2x3+3=5*3=15 chunks.
!
!  30-mar-21/axel: adapted from gravitational_waves_hTXk.f90
!
      use Sub, only: register_report_aux
      use FArrayManager
!
      if (lroot) call svn_id( &
           "$Id$")
!
!  Register AAk and EEk as auxiliary arrays
!  May want to do this only when Fourier transform is enabled.
!
      !call register_report_aux('AAk'  , iAAk  , iAkx  , iAky  , iAkz)
      !call register_report_aux('EEk'  , iEEk  , iEkx  , iEky  , iEkz)
      !call register_report_aux('AAkim', iAAkim, iAkxim, iAkyim, iAkzim)
      !call register_report_aux('EEkim', iEEkim, iEkxim, iEkyim, iEkzim)
      call farray_register_auxiliary('AAk',iAAk,vector=3)
      call farray_register_auxiliary('AAkim',iAAkim,vector=3)
      call farray_register_auxiliary('EEk',iEEk,vector=3)
      call farray_register_auxiliary('EEkim',iEEkim,vector=3)
      iAkx  =iAAk  ; iAky  =iAAk  +1; iAkz  =iAAk  +2
      iAkxim=iAAkim; iAkyim=iAAkim+1; iAkzim=iAAkim+2 
      iEkx  =iEEk  ; iEky  =iEEk  +1; iEkz  =iEEk  +2
      iEkxim=iEEkim; iEkyim=iEEkim+1; iEkzim=iEEkim+2 
!
!  register B array as aux
!
      !if (lbb_as_aux) call register_report_aux('bb', ibb, ibx, iby, ibz)
      !if (lEE_as_aux) call register_report_aux('EE', iEE, iEEx, iEEy, iEEz)
      if (lbb_as_aux) call farray_register_auxiliary('bb',ibb,vector=3)
      if (lEE_as_aux) call farray_register_auxiliary('EE',iEE,vector=3)
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
      integer :: stat, i
!
!  set speed of light
!
      select case (cc_light)
        case ('1'); c_light2=1.
        case ('cgs'); c_light2=c_light_cgs**2
        case default
          call fatal_error("initialize_special: No such value for cc_light:" &
              ,trim(cc_light))
      endselect
      if (headt) print*,'c_light2=',c_light2
!
      if (.not.allocated(Tpq_re)) allocate(Tpq_re(nx,ny,nz,6),stat=stat)
      if (stat>0) call fatal_error('compute_bb_from_AAk_and_EEk','Could not allocate memory for Tpq_re')
!
      if (.not.allocated(Tpq_im)) allocate(Tpq_im(nx,ny,nz,6),stat=stat)
      if (stat>0) call fatal_error('compute_bb_from_AAk_and_EEk','Could not allocate memory for Tpq_im')
!
!  calculate kscale_factor (for later binning)
!
      kscale_factor=2*pi/Lx
!
!  compute alpha*(alpha+1) from Sharma+17 paper
!
      alpha2_inflation=alpha_inflation*(alpha_inflation)
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
!  initial condition for AAk
!
      select case (initAAk)
        case ('nothing')
          if (lroot) print*,'init_special: nothing'
          f(:,:,:,iAAk:iAAkim+2)=0.
        case ('single')
          if (lroot) print*,'init_special for A: single'
          f(:,:,:,iAAk:iAAkim+2)=0.
          if (lroot) f(l1+1,m1+1,n1+1,iAAk)=1.
        case ('powerlow??>')
          !call powerlaw ...
print*,'AXEL: later'
        case default
          call fatal_error("init_special: No such value for initAAk:" &
              ,trim(initAAk))
      endselect
!
!  initial condition for EEk
!
      select case (initEEk)
        case ('nothing')
          if (lroot) print*,'init_special: nothing'
          f(:,:,:,iEEk:iEEkim+2)=0.
        case ('single')
          if (lroot) print*,'init_special for E: single'
          f(:,:,:,iEEk:iEEkim+2)=0.
          if (lroot) f(l1+1,m1+1,n1+1,iEEk)=1.
!
!  dA/deta = -E
!  d^2A/deta^2 = -dE/deta = -(k^2-alpha*(alpha+1)
!
        case ('powerlow??>')
          !call powerlaw ...
print*,'AXEL: later'
        case default
          call fatal_error("init_special: No such value for initEEk:" &
              ,trim(initEEk))
      endselect
print*,'AXEL2',iAAk,iAAkim
print*,'AXEL3',iEEk,iEEkim
print*,'AXEL4',f(5,5,5,iEEkim:iEEkim+2)
!
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special
!
!  All pencils that this special module depends on are specified here.
!
!   1-apr-06/axel: coded
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine pencil_interdep_special(lpencil_in)
!
!  Interdependency among pencils provided by this module are specified here.
!
!  18-jul-06/tony: coded
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
!  30-mar-21/axel: adapted from gravitational_waves_hTXk.f90
!
      use Deriv, only: derij
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: prefactor
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
      integer :: i, j, ij
!
!  Construct stress tensor; notice opposite signs for u and b.
!
!     if (lgamma_factor) then
!       prefactor=fourthird_factor/(1.-p%u2)
!     else
!       prefactor=fourthird_factor
!     endif
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
!  This routine computes various diagnostic quantities.
!
!  06-oct-03/tony: coded
!  07-feb-18/axel: added nscale_factor=0 (no expansion), =.5 (radiation era)
!
      use Diagnostics
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real :: scale_factor, stress_prefactor2, sign_switch=0, fac_stress_comp
      type (pencil_case) :: p
!
      intent(in) :: p
      intent(inout) :: f,df
!
!  Identify module and boundary conditions.
!
      if (lfirst) then
        if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dspecial_dt'
!
!  diagnostics
!
       if (ldiagnos) then
         if (idiag_AA2m/=0) call sum_mn_name(( &
             f(l1:l2,m,n,iAkx  )**2+f(l1:l2,m,n,iAky  )**2+f(l1:l2,m,n,iAkz  )**2+ &
             f(l1:l2,m,n,iAkxim)**2+f(l1:l2,m,n,iAkyim)**2+f(l1:l2,m,n,iAkzim)**2 &
                                               )*nwgrid,idiag_AA2m)
!
         if (lproc_pt.and.m==mpoint.and.n==npoint) then
           if (idiag_Akxpt/=0) call save_name(f(lpoint,m,n,iAkx),idiag_Akxpt)
           if (idiag_Ekxpt/=0) call save_name(f(lpoint,m,n,iEkx),idiag_Ekxpt)
         endif
       endif
      else
        if (headtt.or.ldebug) print*,'dspecial_dt: DONT SOLVE dspecial_dt'
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
      use Sub, only: remove_mean_value
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
!  07-aug-17/axel: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
    endsubroutine special_after_boundary
!***********************************************************************
    subroutine special_after_timestep(f,df,dt_,llast)
!
!  Possibility to modify the f and df after df is updated.
!
!  31-mar-21/axel: adapted from gravitational_waves_hTXk.f90
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      real, intent(in) :: dt_
      logical, intent(in) :: llast
!
!  Compute the transverse part of the stress tensor by going into Fourier space.
!
      if (lfirst) call compute_bb_from_AAk_and_EEk(f)
!
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(dt_)
!
    endsubroutine special_after_timestep
!***********************************************************************
    subroutine make_spectra(f)
!
!  31-mar-21/axel: adapted from gravitational_waves_hTXk.f90
!
      real, dimension (mx,my,mz,mfarray) :: f

      integer :: ikx, iky, ikz, q, p, pq, ik
      real :: k1, k2, k3, ksqr,one_over_k2,one_over_k4,sign_switch
      real :: k1mNy, k2mNy, k3mNy, SCL_re, SCL_im

      spectra%mag=0.; spectra%maghel=0.
      spectra%ele=0.; spectra%elehel=0.
!
!  Loop over all positions in k-space.
!
      do ikz=1,nz
        do iky=1,ny
          do ikx=1,nx
!
            k1=kx_fft(ikx+ipx*nx)
            k2=ky_fft(iky+ipy*ny)
            k3=kz_fft(ikz+ipz*nz)
!
            ksqr=k1**2+k2**2+k3**2
!
            if (lroot.and.ikx==1.and.iky==1.and.ikz==1) then
              one_over_k2=0.
            else
              one_over_k2=1./ksqr
            endif
            one_over_k4=one_over_k2**2
!
            ik=1+nint(sqrt(ksqr)/kscale_factor)
!
!  Debug output
!
            if (ldebug_print) then
              if (ik <= 5) write(*,1000) iproc,ik,k1,k2,k3,f(nghost+ikx,nghost+iky,nghost+ikz,iggX  )
              1000 format(2i5,1p,4e11.2)
            endif
!
            if (ik <= nk) then
!
!  Electromagnetic wave spectrum computed from hdot (=g)
!
  !           if (GWs_spec) then
  !             spectra%GWs(ik)=spectra%GWs(ik) &
  !                +f(nghost+ikx,nghost+iky,nghost+ikz,iggX  )**2 &
  !                +f(nghost+ikx,nghost+iky,nghost+ikz,iggXim)**2 &
  !                +f(nghost+ikx,nghost+iky,nghost+ikz,iggT  )**2 &
  !                +f(nghost+ikx,nghost+iky,nghost+ikz,iggTim)**2
  !             spectra%GWshel(ik)=spectra%GWshel(ik)+2*sign_switch*( &
  !                +f(nghost+ikx,nghost+iky,nghost+ikz,iggXim) &
  !                *f(nghost+ikx,nghost+iky,nghost+ikz,iggT  ) &
  !                -f(nghost+ikx,nghost+iky,nghost+ikz,iggX  ) &
  !                *f(nghost+ikx,nghost+iky,nghost+ikz,iggTim) )
  !           endif
!
  !           if (GWs_spec_complex) then
  !             if (k2==0. .and. k3==0.) then
  !               spectra%complex_Str_T(ikx)=cmplx(f(nghost+ikx,nghost+iky,nghost+ikz,ihhT  ), &
  !                                                f(nghost+ikx,nghost+iky,nghost+ikz,ihhTim))
  !               spectra%complex_Str_X(ikx)=cmplx(f(nghost+ikx,nghost+iky,nghost+ikz,ihhX  ), &
  !                                                f(nghost+ikx,nghost+iky,nghost+ikz,ihhXim))
  !             else
  !               spectra%complex_Str_T(ikx)=0.
  !               spectra%complex_Str_X(ikx)=0.
  !             endif
  !           endif
!
!  Gravitational wave strain spectrum computed from h
!
  !           if (GWh_spec) then
  !             spectra%GWh(ik)=spectra%GWh(ik) &
  !                +f(nghost+ikx,nghost+iky,nghost+ikz,ihhX  )**2 &
  !                +f(nghost+ikx,nghost+iky,nghost+ikz,ihhXim)**2 &
  !                +f(nghost+ikx,nghost+iky,nghost+ikz,ihhT  )**2 &
  !                +f(nghost+ikx,nghost+iky,nghost+ikz,ihhTim)**2
  !             spectra%GWhhel(ik)=spectra%GWhhel(ik)+2*sign_switch*( &
  !                +f(nghost+ikx,nghost+iky,nghost+ikz,ihhXim) &
  !                *f(nghost+ikx,nghost+iky,nghost+ikz,ihhT  ) &
  !                -f(nghost+ikx,nghost+iky,nghost+ikz,ihhX  ) &
  !                *f(nghost+ikx,nghost+iky,nghost+ikz,ihhTim) )
  !           endif
!
!  Gravitational wave mixed spectrum computed from h and g
!
  !           if (GWm_spec) then
  !             spectra%GWm(ik)=spectra%GWm(ik) &
  !                +f(nghost+ikx,nghost+iky,nghost+ikz,ihhX)  *f(nghost+ikx,nghost+iky,nghost+ikz,iggX  ) &
  !                +f(nghost+ikx,nghost+iky,nghost+ikz,ihhXim)*f(nghost+ikx,nghost+iky,nghost+ikz,iggXim) &
  !                +f(nghost+ikx,nghost+iky,nghost+ikz,ihhT)  *f(nghost+ikx,nghost+iky,nghost+ikz,iggT  ) &
  !                +f(nghost+ikx,nghost+iky,nghost+ikz,ihhTim)*f(nghost+ikx,nghost+iky,nghost+ikz,iggTim)
  !             spectra%GWmhel(ik)=spectra%GWmhel(ik)-sign_switch*( &
  !                +f(nghost+ikx,nghost+iky,nghost+ikz,ihhXim) &
  !                *f(nghost+ikx,nghost+iky,nghost+ikz,iggT  ) &
  !                +f(nghost+ikx,nghost+iky,nghost+ikz,iggXim) &
  !                *f(nghost+ikx,nghost+iky,nghost+ikz,ihhT  ) &
  !                -f(nghost+ikx,nghost+iky,nghost+ikz,ihhX  ) &
  !                *f(nghost+ikx,nghost+iky,nghost+ikz,iggTim) &
  !                -f(nghost+ikx,nghost+iky,nghost+ikz,iggX  ) &
  !                *f(nghost+ikx,nghost+iky,nghost+ikz,ihhTim) )
  !           endif
!
            endif

          enddo
        enddo
      enddo

    endsubroutine make_spectra
!***********************************************************************
    subroutine special_calc_spectra(f,spectrum,spectrum_hel,lfirstcall,kind)
!
!  Calculates GW spectra. For use with a single special module.
!
!  16-oct-19/MR: carved out from compute_gT_and_gX_from_gij
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (:) :: spectrum,spectrum_hel
      logical :: lfirstcall
      character(LEN=3) :: kind

      if (lfirstcall) then
        call make_spectra(f)
        lfirstcall=.false.
      endif

      select case(kind)
      case ('mag'); spectrum=spectra%mag; spectrum_hel=spectra%maghel
      case ('ele'); spectrum=spectra%ele; spectrum_hel=spectra%elehel
      case default; if (lroot) call warning('special_calc_spectra', &
                      'kind of spectrum "'//kind//'" not implemented')
      endselect

    endsubroutine special_calc_spectra
!***********************************************************************
    subroutine special_calc_spectra_byte(f,spectrum,spectrum_hel,lfirstcall,kind,len)
!
!  Calculates magnetic and electric spectra. For use with multiple special modules.
!
!  30-mar-21/axel: adapted from gravitational_waves_hTXk.f90
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nk) :: spectrum,spectrum_hel
      logical :: lfirstcall
      integer(KIND=ikind1), dimension(3) :: kind
      integer :: len

      character(LEN=3) :: kindstr

      if (lfirstcall) then
        call make_spectra(f)
        lfirstcall=.false.
      endif

      kindstr=char(kind(1))//char(kind(2))//char(kind(3))

      select case(kindstr)
      case ('mag'); spectrum=spectra%mag; spectrum_hel=spectra%maghel
      case ('ele'); spectrum=spectra%ele; spectrum_hel=spectra%elehel
      case default; if (lroot) call warning('special_calc_spectra', &
                      'kind of spectrum "'//kindstr//'" not implemented')
      endselect

    endsubroutine special_calc_spectra_byte
!***********************************************************************
    subroutine compute_bb_from_AAk_and_EEk(f)
!
!  Compute the transverse part of the stress tensor by going into Fourier space.
!
!  07-aug-17/axel: coded
!
      use Fourier, only: fourier_transform, fft_xyz_parallel
      use SharedVariables, only: put_shared_variable
!
      real, dimension (:,:,:,:), allocatable :: BBkre, BBkim
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (3) :: coefAre, coefAim, coefBre, coefBim
      integer :: i,j,p,q,ik,ikx,iky,ikz,stat,ij,pq,ip,jq,jStress_ij
      complex, dimension (3) :: Acomplex, Ecomplex, Acomplex_new, Ecomplex_new
      complex :: discrim, det1, lam1, lam2, explam1t, explam2t
      complex :: cosotA, cosotE, sinotA, sinotE
      real :: sigmaeff
      real :: fact
      real :: ksqr, k1, k2, k3, k1sqr, k2sqr, k3sqr
!     real :: cosot, sinot, sinot_minus, om12, om, om1, om2
      intent(inout) :: f
!
!  For testing purposes, if lno_transverse_part=T, we would not need to
!  compute the Fourier transform, so we would skip the rest.
!
!  Allocate memory for arrays.
!
      allocate(BBkre(nx,ny,nz,3),stat=stat)
      if (stat>0) call fatal_error('compute_bb_from_AAk_and_EEk','Could not allocate memory for BBkre')
!
      allocate(BBkim(nx,ny,nz,3),stat=stat)
      if (stat>0) call fatal_error('compute_bb_from_AAk_and_EEk','Could not allocate memory for BBkim')
!
!  Set k^2 array. Note that in Fourier space, kz is the fastest index and has
!  the full nx extent (which, currently, must be equal to nxgrid).
!
      do ikz=1,nz
        do iky=1,ny
          do ikx=1,nx
!
!  collect k vector and compute k^2 at each point
!
            k1=kx_fft(ikx+ipx*nx)
            k2=ky_fft(iky+ipy*ny)
            k3=kz_fft(ikz+ipz*nz)
            k1sqr=k1**2
            k2sqr=k2**2
            k3sqr=k3**2
            ksqr=k1sqr+k2sqr+k3sqr
!
!  compute eigenvalues
!
            sigmaeff=sigma
            discrim=sqrt(complex(sigmaeff**2-4.*ksqr,0.))
            lam1=.5*(-sigmaeff+discrim)
            lam2=.5*(-sigmaeff-discrim)
!
!  Compute exact solution for hT, hX, gT, and gX in Fourier space.
!
            AAkre=f(nghost+ikx,nghost+iky,nghost+ikz,iAAk  :iAAk  +2)
            AAkim=f(nghost+ikx,nghost+iky,nghost+ikz,iAAkim:iAAkim+2)
!
            EEkre=f(nghost+ikx,nghost+iky,nghost+ikz,iEEk  :iEEk  +2)
            EEkim=f(nghost+ikx,nghost+iky,nghost+ikz,iEEkim:iEEkim+2)
!
            do j=1,3
              Acomplex(j)=complex(AAkre(j),AAkim(j))
              Ecomplex(j)=complex(EEkre(j),EEkim(j))
            enddo
!
!  compute cos(om*dt) and sin(om*dt) to get from one timestep to the next.
!
            if (ksqr/=0.) then
              explam1t=exp(lam1*dt)
              explam2t=exp(lam2*dt)
              det1=1./discrim
              cosotA=det1*(lam1*explam2t-lam2*explam1t)
              cosotE=det1*(lam1*explam1t-lam2*explam2t)
              sinotA=det1*(     explam2t-     explam1t)
              sinotE=-sinotA*lam1*lam2
!
!  Solve wave equation for hT and gT from one timestep to the next.
!
              Acomplex_new=cosotA*Acomplex+sinotA*Ecomplex
              Ecomplex_new=sinotE*Acomplex+cosotE*Ecomplex
              f(nghost+ikx,nghost+iky,nghost+ikz,iAAk  :iAAk  +2)= real(Acomplex_new)
              f(nghost+ikx,nghost+iky,nghost+ikz,iAAkim:iAAkim+2)=aimag(Acomplex_new)
              f(nghost+ikx,nghost+iky,nghost+ikz,iEEk  :iEEk  +2)= real(Ecomplex_new)
              f(nghost+ikx,nghost+iky,nghost+ikz,iEEkim:iEEkim+2)=aimag(Ecomplex_new)
!
!  Debug output
!
              if (ldebug_print.and.lroot) then
                if (ikx==2.and.iky==1.and.ikz==1) then
                  print*,'AXEL: Acomplex_new=',Acomplex_new
                  print*,'AXEL: Ecomplex_new=',Ecomplex_new
                  print*,'AXEL: cosotA=',cosotA,cos(sqrt(ksqr)*dt)
                  print*,'AXEL: cosotE=',cosotE,cos(sqrt(ksqr)*dt)
                  print*,'AXEL: sinotA=',sinotA,-sin(sqrt(ksqr)*dt)/sqrt(ksqr)
                  print*,'AXEL: sinotE=',sinotE,+sin(sqrt(ksqr)*dt)*sqrt(ksqr)
                endif
              endif
!
            else
!
!  Set origin to zero. It is given by (1,1,1) on root processor.
!
              f(nghost+1,nghost+1,nghost+1,iAAk  :iAAk  +2)=0.
              f(nghost+1,nghost+1,nghost+1,iAAkim:iAAkim+2)=0.
              f(nghost+1,nghost+1,nghost+1,iEEk  :iEEk  +2)=0.
              f(nghost+1,nghost+1,nghost+1,iEEkim:iEEkim+2)=0.
!
            endif
!
!  end of ikx, iky, and ikz loops
!
          enddo
        enddo
      enddo
!
!  back to real space
!  Use AAk instead of BBk for now
!
      if (lbb_as_aux) then
        BBkre=f(l1:l2,m1:m2,n1:n2,iAAk  :iAAk  +2)
        BBkim=f(l1:l2,m1:m2,n1:n2,iAAkim:iAAkim+2)
        call fft_xyz_parallel(BBkre,BBkim,linv=.true.)
        f(l1:l2,m1:m2,n1:n2,ibb:ibb+2)=BBkre
      endif
!
!  EE back to real space, use the names BBkre and BBkim for EE.
!
      if (lEE_as_aux) then
        BBkre=f(l1:l2,m1:m2,n1:n2,iEEk  :iEEk  +2)
        BBkim=f(l1:l2,m1:m2,n1:n2,iEEkim:iEEkim+2)
        call fft_xyz_parallel(BBkre,BBkim,linv=.true.)
        f(l1:l2,m1:m2,n1:n2,iEE:iEE+2)=BBkre
      endif
!
    endsubroutine compute_bb_from_AAk_and_EEk
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
!  Reads and registers print parameters relevant to special.
!
!  06-oct-03/tony: coded
!
      use Diagnostics
!!      use FArrayManager, only: farray_index_append
!
      integer :: iname
      logical :: lreset,lwrite
!!!
!!!  reset everything in case of reset
!!!  (this needs to be consistent with what is defined above!)
!!!
      if (lreset) then
        idiag_AA2m=0
        idiag_Akxpt=0
        idiag_Ekxpt=0
        cformv=''
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'AA2m',idiag_AA2m)
        call parse_name(iname,cname(iname),cform(iname),'Akxpt',idiag_Akxpt)
        call parse_name(iname,cname(iname),cform(iname),'Ekxpt',idiag_Ekxpt)
      enddo
!
!  check for those quantities for which we want video slices
!
      if (lwrite_slices) then
        where(cnamev=='hhT'.or.cnamev=='hhX'.or.cnamev=='ggT'.or.cnamev=='ggX'.or. &
              cnamev=='h22'.or.cnamev=='h33'.or.cnamev=='h23') cformv='DEFINED'
      endif
!
!!!  write column where which magnetic variable is stored
!!      if (lwrite) then
!!        call farray_index_append('i_SPECIAL_DIAGNOSTIC',i_SPECIAL_DIAGNOSTIC)
!!      endif
!!
    endsubroutine rprint_special
!***********************************************************************
    subroutine get_slices_special(f,slices)
!
!  Write slices for animation of Special variables.
!
!  11-apr-21/axel: adapted from gravitational_waves_hTXk.f90
!
      use Slices_methods, only: assign_slices_scal
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
!  Loop over slices
!
      select case (trim(slices%name))
!
!  hhT
!
   !    case ('hhT')
   !      if (lreal_space_hTX_as_aux) then
   !        call assign_slices_scal(slices,f,ihhT_realspace)
   !      else
   !        call assign_slices_scal(slices,f,ihhT)
   !      endif
!
!  hhX
!
   !    case ('hhX')
   !      if (lreal_space_hTX_as_aux) then
   !        call assign_slices_scal(slices,f,ihhX_realspace)
   !      else
   !        call assign_slices_scal(slices,f,ihhX)
   !      endif
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
