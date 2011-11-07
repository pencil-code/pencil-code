! $Id$
!
!  This modules deals with all aspects of polymers.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lpolymer = .true.
!
! MVAR CONTRIBUTION 6
! MAUX CONTRIBUTION 1
!
! PENCILS PROVIDED poly(3,3);trp;fr;frC(3,3)
! PENCILS PROVIDED u_dot_gradC(3,3);Cijk(3,3,3); del2poly(3,3);
! PENCILS PROVIDED div_frC(3); divC(3); grad_fr(3)
!
!***************************************************************
module Polymer
!
  use Cdata
  use Cparam
  use Messages
  use Sub
!
  implicit none
!
  include 'record_types.h'
  include 'polymer.h'
!
  real, dimension(nx,3,3) :: CdotGradu,CdotGraduT
!
!  Start parameters.
!
  character (len=labellen), dimension(ninit) :: initpoly='nothing'
!
  namelist /polymer_init_pars/ &
     initpoly
!
!  Run parameters.
!
  logical :: lpolyback=.true.
  logical :: lpolyadvect=.true.
  logical :: lpoly_diffusion=.false.
  logical :: lupw_poly=.false.
  real :: mu_poly=0.,tau_poly=0.,tau_poly1=1.,eta_poly=0.,fenep_L=0.
  character (len=labellen) :: poly_algo='simple',poly_model='oldroyd-B'
!
   namelist /polymer_run_pars/ &
       lpolyback,lpolyadvect,lpoly_diffusion, &
       mu_poly,tau_poly,tau_poly1,eta_poly,fenep_L,lupw_poly,poly_model
!
! Diagnostic variables
!
  integer :: idiag_polytrm=0     ! DIAG_DOC: $\left\langle Tr[C_{ij}]\right\rangle$
!
  contains
!***********************************************************************
    subroutine register_polymer()
!
!  Initialise variables.
!
!  14-Aug-08/Dhruba: coded
!
      use FArrayManager
!
      call farray_register_pde('poly',ipoly,vector=6)
      ip11 = ipoly     ; ip12 = ipoly+1 ; ip13 = ipoly+2
      ip21 = ip12   ; ip22 = ipoly+3 ; ip23 = ipoly+4
      ip31 = ip13   ; ip32 = ip23   ; ip33 = ipoly+5
      call farray_register_auxiliary('polyfr',ipoly_fr)
!
      if (lroot) call svn_id( &
          "$Id$")
!
    endsubroutine register_polymer
!***********************************************************************
    subroutine initialize_polymer(f,lstarting)
!
!  Perform any post-parameter-read initialization.
!  At present does nothing.
!
!  14-aug-08/dhruba: initialize polymer field
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
      if (tau_poly>tini) then
        tau_poly1=1./tau_poly
      else
        tau_poly1=1.
      endif
!
! give warning if polymeric backreaction is set true but polymer
! density, mu_poly is set to zero.
!
      if (lpolyback.and.(mu_poly==0.)) &
          call warning ('initialize_polymer','lpolyback=T but mu_poly=0!')
!
    endsubroutine initialize_polymer
!***********************************************************************
    subroutine init_poly(f)
!
!  Initialise polymer field.
!
!   14-aug-2008/dhruba: coded
!
      use Initcond
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: rsqr
      integer :: j
!
      do j=1,ninit
!
        select case(initpoly(j))
          case('nothing'); if (lroot .and. j==1) print*,'init_poly: nothing'
          case('zero', '0'); f(:,:,:,ipoly:ipoly+5) = 0.
          case('sphere')
            f(:,:,:,ip11) = 1.
            f(:,:,:,ip12) = 0.
            f(:,:,:,ip13) = 0.
            f(:,:,:,ip22) = 1.
            f(:,:,:,ip23) = 0.
            f(:,:,:,ip33) = 1.
          case default
!
!  Catch unknown values.
!
            call fatal_error('init_poly', &
                'init_poly value "' // trim(initpoly(j)) // '" not recognised')
        endselect
!
!  End loop over initial conditions.
!
      enddo
!
!  Interface for user's own initial condition.
!   Not implemented for polymers yet.
!
!      if (linitial_condition) call initial_condition_pp(f)
!
!  Also set the f(r_p) depending on the model of the polymer.
!
      select case (poly_model)
        case ('oldroyd-B')
          f(:,:,:,ipoly_fr) = 1.
        case ('FENE-P')
          do iz=1,mz; do iy=1,my; do ix=1,mx
            rsqr = f(ix,iy,iz,ip11)+f(ix,iy,iz,ip22)+&
                   f(ix,iy,iz,ip33)
            f(ix,iy,iz,ipoly_fr) = (fenep_L**2-3)/(fenep_L**2-rsqr)
          enddo; enddo; enddo
        case default
          call fatal_error('init_poly','no such polymer model')
      endselect
!
    endsubroutine init_poly
!***********************************************************************
    subroutine pencil_criteria_polymer()
!
!   All pencils that the Polymer module depends on are specified here.
!
      lpenc_requested(i_poly)=.true.
      if (tau_poly/=0) then
        lpenc_requested(i_frC)=.true.
        lpenc_requested(i_fr)=.true.
      endif
!
!  If we consider backreaction from the polymer to the fluid.
!
      if (lhydro.and.lpolyback) then
        lpenc_requested(i_div_frC)=.true.
        lpenc_requested(i_divC)=.true.
        lpenc_requested(i_grad_fr)=.true.
      endif
!
!  If advection by the velocity is turned on.
!
      if ((lhydro.or.lhydro_kinematic).and.(lpolyadvect)) then
        lpenc_requested(i_u_dot_gradC)=.true.
        lpenc_requested(i_uu)=.true.
        lpenc_requested(i_uij)=.true.
      endif
!
!  If a diffusive term in the polymer equation is not included: (not default)
!
      if (eta_poly/=0) lpenc_requested(i_del2poly)=.true.
!
!  Different pencils are chosen depending on different algorithms applied.
!
      select case(poly_algo)
        case('simple')
          write(*,*) 'poly_algo:no more pencils needed now'
        case('cholesky')
          call fatal_error('pencil_criteria_polymer', &
              'poly_algo: cholesky decomposition is not implemented yet')
        case('nothing')
          call fatal_error('pencil_criteria_polymer', &
              'poly_algo: please chosse an algorithm to solve the '// &
              'polymer equations')
      endselect
!
! Diagnostic pencils
!
      if (idiag_polytrm/=0) lpenc_requested(i_trp)=.true.
!
    endsubroutine pencil_criteria_polymer
!***********************************************************************
    subroutine pencil_interdep_polymer(lpencil_in)
!
!  Interdependency among pencils from the Polymer module is specified here.
!
!  18-aug-2008/dhruba: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      if (lpencil_in(i_u_dot_gradC)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_Cijk)=.true.
      endif
      if (lpencil_in(i_div_frC)) then
        lpencil_in(i_divC)=.true.
        lpencil_in(i_grad_fr)=.true.
      endif
      if (lpencil_in(i_divC)) lpencil_in(i_Cijk) = .true.
      if (lpencil_in(i_frC)) then
        lpencil_in(i_fr)=.true.
        lpencil_in(i_poly)=.true.
      endif
      if (lpencil_in(i_trp)) lpencil_in(i_poly)=.true.
!
    endsubroutine pencil_interdep_polymer
!***********************************************************************
    subroutine calc_pencils_polymer(f,p)
!
!  Calculate polymer pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  18-aug-08/dhruba: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      integer :: ipi,ipj,ipk
!
      intent(inout) :: f,p
! poly
      if (lpencil(i_poly)) then
        ipk=0
        do ipi=1,3
          do ipj=ipi,3
            p%poly(:,ipi,ipj) = f(l1:l2,m,n,ipoly+ipk)
            ipk=ipk+1
          enddo
        enddo
      endif
      call symmetrise3x3_ut2lt(p%poly)
!
! polymer diffusion
!
      if (lpencil(i_del2poly)) call del2m3x3_sym(f,ipoly,p%del2poly)
! trp
      if (lpencil(i_trp)) call trace_mn(p%poly,p%trp)
! Cijk
      if (lpencil(i_Cijk)) call gijl_symmetric(f,ipoly,p%Cijk)
! u_dot_gradC
      if (lpencil(i_u_dot_gradC))&
          call u_dot_grad_mat(f,ipoly,p%Cijk,p%uu,p%u_dot_gradC, &
            UPWIND=lupw_poly)
!
      select case (poly_model)
        case ('oldroyd-B')
          call calc_pencils_oldroyd_b(f,p)
        case ('FENE-P')
          call calc_pencils_fene_p(f,p)
        case default
          call fatal_error('calc_pencils_polymer','no such polymer model')
      endselect
!
      if (ldiagnos) call polymer_diagnostic(f,p)
!
! Time step constaint
!
      trelax_poly=tau_poly
!
    endsubroutine calc_pencils_polymer
!***********************************************************************
    subroutine calc_pencils_fene_p(f,p)
!
!  Calculate polymer pencils.
!  Most basic pencils should come first, as others may depend on them.
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      real, dimension (nx,3) :: grad_fr_dotC
      integer :: i,j
!
! f(r_p)
      if (lpencil(i_fr)) p%fr=f(l1:l2,m,n,ipoly_fr)
!
! frC
!
      if (lpencil(i_frC)) then
        do i=1,3; do j=1,3
          p%frC(:,i,j) = p%fr(:)*p%poly(:,i,j)
        enddo; enddo
      endif
! div C
      if (lpencil(i_divC)) call div_mn_2tensor(p%Cijk,p%divC)
! grad f(r_p)
      if (lpencil(i_grad_fr)) then
        call  grad(f,ipoly_fr,p%grad_fr)
        call dot_mn_vm(p%grad_fr,p%poly,grad_fr_dotC)
      endif
! div_frC
      if (lpencil(i_div_frC)) then
        do j=1,3
          p%div_frC(:,j)=p%fr(:)*p%divC(:,j)+grad_fr_dotC(:,j)
        enddo
      endif
!
    endsubroutine calc_pencils_fene_p
!***********************************************************************
    subroutine calc_pencils_oldroyd_b(f,p)
!
!  Calculate polymer pencils.
!  Most basic pencils should come first, as others may depend on them.
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f)
!
! f(r_p)
!
      if (lpencil(i_fr)) p%fr(:) = 1.
!
! frC
!
      if (lpencil(i_frC)) p%frC = p%poly
! div C
      if (lpencil(i_divC)) call div_mn_2tensor(p%Cijk,p%divC)
! grad f(r_p)
      if (lpencil(i_grad_fr)) p%grad_fr=0.
! div_frC
      if (lpencil(i_div_frC)) p%div_frC=p%divC
!
    endsubroutine calc_pencils_oldroyd_b
!***********************************************************************
    subroutine polymer_diagnostic(f,p)
!
!  Calculates the diagnostic quantities for polymer module
!  Most basic pencils should come first, as others may depend on them.
!
      use Diagnostics, only: sum_mn_name
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f)
!
      if (idiag_polytrm/=0)   call sum_mn_name(p%trp,idiag_polytrm)
!
    endsubroutine polymer_diagnostic
!***********************************************************************
    subroutine dpoly_dt(f,df,p)
!
!  Polymer evolution.
!
!  18-aug-08/dhruba: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension(nx,3,3) :: uijT
      integer :: ipi,ipj,ipk
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print*,'dpoly_dt: SOLVE'
      if (headtt) then
        call identify_bcs('P11',ip11)
        call identify_bcs('P22',ip22)
        call identify_bcs('P33',ip33)
        call identify_bcs('P12',ip12)
        call identify_bcs('P13',ip13)
        call identify_bcs('P32',ip32)
      endif
!
!  Add backreaction due to the polymer to momentum equation (default).
!
      if (lhydro.and.lpolyback) then
        if (tau_poly/=0.0) then
          df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)+ &
              mu_poly*tau_poly1*p%div_frC
        endif
      endif
!
!  If we are advecting the polymer (which is default).
!
      if (lpolyadvect)  then
        ipk=0
        do ipi=1,3
          do ipj=ipi,3
            df(l1:l2,m,n,ipoly+ipk)= df(l1:l2,m,n,ipoly+ipk) - &
                p%u_dot_gradC(:,ipi,ipj)
            ipk=ipk+1
          enddo
        enddo
      endif
!
!  CdotGradu and CdotGraduT.
!
      call mult_matrix(p%poly,p%uij,CdotGradu)
      call transpose_mn(p%uij,uijT)
      call mult_matrix(uijT,p%poly,CdotGraduT)
!
!  Select which algorithm we are using.
!
      select case(poly_algo)
      case('simple')
        call simple_dpoly_dt (f,df,p)
      case('cholesky')
        call fatal_error('pencil_criteria_polymer', &
            'poly_algo: cholesky decomposition is not implemented yet')
      case('nothing')
        call fatal_error('pencil_criteria_polymer', &
            'poly_algo: please chosse an algorithm to solve '// &
            'the polymer equations')
      endselect
!
!  polymer diffusion (sometime only for numerical stability)
!
      ipk=0
      do ipi=1,3
        do ipj=ipi,3
          df(l1:l2,m,n,ipoly+ipk)= &
              df(l1:l2,m,n,ipoly+ipk)-eta_poly*p%del2poly(:,ipi,ipj)
          ipk=ipk+1
        enddo
      enddo
!
    endsubroutine dpoly_dt
!***********************************************************************
    subroutine simple_dpoly_dt(f,df,p)
!
! the simplest algorithm.
!
!  24-feb-11/dhruba: moved to a subroutine
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      integer :: ipi,ipj,ipk
!
      call keep_compiler_quiet(f)
!
      ipk=0
      do ipi=1,3
        do ipj=ipi,3
          df(l1:l2,m,n,ipoly+ipk)=df(l1:l2,m,n,ipoly+ipk)+ &
              CdotGradu(:,ipi,ipj)+CdotGraduT(:,ipi,ipj)
          if (tau_poly/=0.) &
            df(l1:l2,m,n,ipoly+ipk)=df(l1:l2,m,n,ipoly+ipk)- &
                tau_poly1*(p%frC(:,ipi,ipj)-kronecker_delta(ipi,ipj))
           ipk=ipk+1
        enddo
      enddo
!
    endsubroutine simple_dpoly_dt
!***********************************************************************
    subroutine calc_polymer_after_boundary(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: rsqr
!
      select case (poly_model)
        case ('oldroyd-B')
!         do nothing
          call keep_compiler_quiet(f)
        case ('FENE-P')
          do iz=1,mz; do iy=1,my; do ix=1,mx
            rsqr = f(ix,iy,iz,ip11)+f(ix,iy,iz,ip22)+&
                f(ix,iy,iz,ip33)
            f(ix,iy,iz,ipoly_fr) = (fenep_L**2-3)/(fenep_L**2-rsqr)
          enddo; enddo; enddo
        case default
          call fatal_error('init_poly','no such polymer model')
        endselect
!
    endsubroutine calc_polymer_after_boundary
!***********************************************************************
    subroutine read_polymer_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=polymer_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=polymer_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_polymer_init_pars
!***********************************************************************
    subroutine write_polymer_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=polymer_init_pars)
!
    endsubroutine write_polymer_init_pars
!***********************************************************************
    subroutine read_polymer_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=polymer_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=polymer_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_polymer_run_pars
!***********************************************************************
    subroutine write_polymer_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=polymer_run_pars)
!
    endsubroutine write_polymer_run_pars
!***********************************************************************
    subroutine get_slices_polymer(f,slices)
!
!  Write slices for animation of polymeric variables.
!
!  26-jul-06/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices)
!
    endsubroutine get_slices_polymer
!***********************************************************************
    subroutine rprint_polymer(lreset,lwrite)
!
!  Reads and registers print parameters relevant for polymer.
!
!   8-aug-10/dhruba: aped from pscalar
!
      use Diagnostics
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname !, inamez, inamey, inamex
      logical :: lwr
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  Reset everything in case of reset.
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_polytrm=0
      endif
!
!  Check for those quantities that we want to evaluate online.
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'polytrm',idiag_polytrm)
      enddo
!
!  Check for those quantities for which we want xy-averages.
!
!      do inamez=1,nnamez
!        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
!            'lnccmz',idiag_lnccmz)
!   enddo
!
!  Check for those quantities for which we want xz-averages.
!
!      do inamey=1,nnamey
!        call parse_name(inamey,cnamey(inamey),cformy(inamey), &
!            'lnccmy',idiag_lnccmy)
!      enddo
!
!  Check for those quantities for which we want yz-averages.
!
!      do inamex=1,nnamex
!        call parse_name(inamex,cnamex(inamex),cformx(inamex), &
!            'lnccmx',idiag_lnccmx)
!      enddo
!
!  Write column where which passive scalar variable is stored.
!
      if (lwr) then
        write(3,*) 'ipoly=',ipoly
        write(3,*) 'ip11=',ip11
        write(3,*) 'ip12=',ip12
        write(3,*) 'ip13=',ip13
        write(3,*) 'ip21=',ip21
        write(3,*) 'ip22=',ip22
        write(3,*) 'ip23=',ip23
        write(3,*) 'ip31=',ip31
        write(3,*) 'ip32=',ip32
        write(3,*) 'ip33=',ip33
        write(3,*) 'ipoly_fr=',ipoly_fr
      endif
!
    endsubroutine rprint_polymer
!***********************************************************************
endmodule Polymer
