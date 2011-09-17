! $Id$
!
!  Mean field dynamo equation
!
!  25-feb-07/axel: adapted from nospecial.f90
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 3
! MAUX CONTRIBUTION 0
!
!***************************************************************
module Special
!
  use Cparam
  use Cdata
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include '../special.h'
!
!  square of wave speed for gauge field
!
  real :: etadyn,alpha_const
!
!  cosine and sine function for setting test fields and analysis
!
  real, dimension(mz) :: cz,sz
  real, dimension(mz,2,2) :: alp_ij,eta_ij
!
  ! input parameters
  real :: cphi=1.,ampl=1e-3,kx=1.,ky=0.,kz=0.
  real :: ktestfield=1., ktestfield1=1.
  integer, parameter :: ntestfield=3*njtest
  character(len=50) :: init='zero'
  logical :: ldebug_meanfield=.false.
  namelist /special_init_pars/ &
    init,ampl,kx,ky,kz

  ! run parameters
  namelist /special_run_pars/ &
    etadyn,alpha_const,ktestfield,ldebug_meanfield
!
! Declare index variables in cdata, because we need them in shear.f90
! other variables (needs to be consistent with reset list below)
!
  integer :: idiag_bmx2m=0,idiag_bmy2m=0
  integer :: idiag_bmxpt=0,idiag_bmypt=0,idiag_bmxp2=0,idiag_bmyp2=0
!
  contains

!***********************************************************************
    subroutine register_special()
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!  6-oct-03/tony: coded
!
      use Cdata
      use FArrayManager
!
      call farray_register_pde('am',iam,vector=3)
      iamx=iam; iamy=iamx+1; iamz=iamx+2
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f,lstarting)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
!  Initialize module variables which are parameter dependent
!  set cosine and sine function for setting test fields and analysis
!
      cz=cos(ktestfield*z)
      sz=sin(ktestfield*z)
!
!  Also calculate its inverse, but only if different from zero
!
      if (ktestfield==0) then
        ktestfield1=1.
      else
        ktestfield1=1./ktestfield
      endif
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f)
!
!  initialise special condition; called from start.f90
!  06-oct-2003/tony: coded
!
      use Cdata
      use Initcond
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      select case (init)
        case ('nothing'); if (lroot) print*,'init_special: nothing'
        case ('zero'); f(:,:,:,iamx:iamz)=0.
        case ('sinwave-x'); call sinwave(ampl,f,iam,kx=kx)
        case ('sinwave-y'); call sinwave(ampl,f,iam,ky=ky)
        case ('sinwave-z'); call sinwave(ampl,f,iam,kz=kz)
        case ('gaussian-noise'); call gaunoise(ampl,f,iamx,iamz)

        case default
          !
          !  Catch unknown values
          !
          if (lroot) print*,'init_special: No such value for init: ', trim(init)
          call stop_it("")
      endselect
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special()
!
!  All pencils that this special module depends on are specified here.
!
!  25-feb-07/axel: adapted
!
      if (cphi/=0.) then
!       lpenc_requested(i_diva)=.true.
      endif
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine pencil_interdep_special(lpencil_in)
!
!  Interdependency among pencils provided by this module are specified here.
!
!  18-07-06/tony: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
!  Calculate Hydro pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!   24-nov-04/tony: coded
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
!  several precalculated Pencils of information are passed if for
!  efficiency.
!
!   06-oct-03/tony: coded
!
      use Cdata
      use Diagnostics
      use Mpicomm
      use Sub
      use Deriv
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: damx,damy,d2amx,d2amy
      real, dimension (nx,2) :: emf,bm,jm
      integer :: i,j
!
      intent(in) :: f,p
      intent(inout) :: df
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dSPECIAL_dt'
!!      if (headtt) call identify_bcs('ss',iss)
!
!  calculate mean current and magnetic field
!
      call der(f,iamx,damx,3)
      call der(f,iamy,damy,3)
      bm(:,1)=-damy
      bm(:,2)=+damx
!
      call der2(f,iamx,d2amx,3)
      call der2(f,iamy,d2amy,3)
      jm(:,1)=-d2amx
      jm(:,2)=-d2amy
!
!  solve dynamo equation
!
      do i=1,2
        emf(:,i)=alpha_const*bm(:,i)
        do j=1,2
          emf(:,i)=emf(:,i)+alp_ij(n,i,j)*bm(:,j)-eta_ij(n,i,j)*jm(:,j)
        enddo
      enddo
!
!  debug output
!
      if (ldebug_meanfield) then
        if (m==m1.and.n==n1) then
          print*,'t=',t
          print*,'alp_11,alp_12=',sum(alp_ij(:,1,1))/nz,sum(alp_ij(:,1,2))/nz
          print*,'alp_21,alp_22=',sum(alp_ij(:,2,1))/nz,sum(alp_ij(:,2,2))/nz
          print*,'eta_11,eta_12=',sum(eta_ij(:,1,1))/nz,sum(eta_ij(:,1,2))/nz
          print*,'eta_21,eta_22=',sum(eta_ij(:,2,1))/nz,sum(eta_ij(:,2,2))/nz
        endif
      endif
!
!  assemble rhs of mean field equation
!
      do j=1,2
        df(l1:l2,m,n,iamx+j-1)=df(l1:l2,m,n,iamx+j-1)+emf(:,j)-etadyn*jm(:,j)
      enddo
!
!  diffusive time step, just take the max of diffus_eta (if existent)
!  and whatever is calculated here
!
      if (lfirst.and.ldt) then
        diffus_eta=max(diffus_eta,etadyn*dxyz_2)
      endif
!
!  diagnostics
!
      if (ldiagnos) then
        if (idiag_bmx2m/=0) call sum_mn_name(bm(:,1)**2,idiag_bmx2m)
        if (idiag_bmy2m/=0) call sum_mn_name(bm(:,2)**2,idiag_bmy2m)
!
!  check for point 1
!
        if (lroot.and.m==mpoint.and.n==npoint) then
          if (idiag_bmxpt/=0) call save_name(bm(lpoint-nghost,1),idiag_bmxpt)
          if (idiag_bmypt/=0) call save_name(bm(lpoint-nghost,2),idiag_bmypt)
        endif
!
!  check for point 2
!
        if (lroot.and.m==mpoint2.and.n==npoint2) then
          if (idiag_bmxp2/=0) call save_name(bm(lpoint2-nghost,1),idiag_bmxp2)
          if (idiag_bmyp2/=0) call save_name(bm(lpoint2-nghost,2),idiag_bmyp2)
        endif
!
      endif
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine read_special_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=special_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=special_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
      integer, intent(in) :: unit
!
!  write name list
!
      write(unit,NML=special_init_pars)
!
    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
!  write name list
!
      if (present(iostat)) then
        read(unit,NML=special_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=special_run_pars,ERR=99)
      endif
!
99  endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
      integer, intent(in) :: unit
!
!  write name list
!
      write(unit,NML=special_run_pars)
!
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
!  reads and registers print parameters relevant to special
!
!   06-oct-03/tony: coded
!
      use Cdata
      use Diagnostics
!!
!!!   SAMPLE IMPLEMENTATION
!!
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
        idiag_bmx2m=0; idiag_bmy2m=0
        idiag_bmxpt=0; idiag_bmypt=0; idiag_bmxp2=0; idiag_bmyp2=0
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'bmx2m',idiag_bmx2m)
        call parse_name(iname,cname(iname),cform(iname),'bmy2m',idiag_bmy2m)
        call parse_name(iname,cname(iname),cform(iname),'bmxpt',idiag_bmxpt)
        call parse_name(iname,cname(iname),cform(iname),'bmypt',idiag_bmypt)
        call parse_name(iname,cname(iname),cform(iname),'bmxp2',idiag_bmxp2)
        call parse_name(iname,cname(iname),cform(iname),'bmyp2',idiag_bmyp2)
      enddo
!
!  write column where which magnetic variable is stored
!
      if (lwr) then
        write(3,*) 'i_bmx2m=',idiag_bmx2m
        write(3,*) 'i_bmy2m=',idiag_bmy2m
        write(3,*) 'i_bmxpt=',idiag_bmxpt
        write(3,*) 'i_bmypt=',idiag_bmypt
        write(3,*) 'i_bmxp2=',idiag_bmxp2
        write(3,*) 'i_bmyp2=',idiag_bmyp2
        write(3,*) 'iam=',iam
      endif
!
    endsubroutine rprint_special
!***********************************************************************
    subroutine get_slices_special(f,slices)
!
!  Write slices for animation of electric potential
!
!  26-feb-07/axel: adapted from gross_pitaevskii
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      type (slice_data) :: slices
!
      integer :: inamev
!
!  Loop over slices
!
      select case (trim(slices%name))
!
!  phi
!
        case ('phi')
          slices%yz=f(ix_loc,m1:m2,n1:n2,iam)
          slices%xz=f(l1:l2,iy_loc,n1:n2,iam)
          slices%xy=f(l1:l2,m1:m2,iz_loc,iam)
          slices%xy2=f(l1:l2,m1:m2,iz2_loc,iam)
          slices%ready = .true.
!
      endselect
!
    endsubroutine get_slices_special
!***********************************************************************
    subroutine calc_lspecial_pars(f)
!
!  calculate alp_ij eta_ij tensors, which are needed when ltestfield=.true.
!
!  15-jan-08/axel: coded
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      real, dimension(mz,2,2) :: eta_ij3
      real, dimension (mz,3,njtest) :: Eipq
      integer :: jtest,j,nxy=nxgrid*nygrid,juxb
      real :: fac
!
      intent(inout) :: f
!
      if (.not.ltestfield) then
        alp_ij(:,:,:)=0.
        eta_ij(:,:,:)=0.
      else
!
!  calculate emf, Eipq
!
        fac=1./nxy
        do jtest=1,njtest
          juxb=iuxb+3*(jtest-1)
          do n=n1,n2
            do j=1,2
              Eipq(n,j,jtest)=fac*sum(f(l1:l2,m1:m2,n,juxb+j-1))
            enddo
          enddo
        enddo
!
!  calculate alp_ij and eta_ij
!
        alp_ij(:,1,1)= +cz*Eipq(:,1,1)+sz*Eipq(:,1,2)
        alp_ij(:,2,1)= +cz*Eipq(:,2,1)+sz*Eipq(:,2,2)
        alp_ij(:,1,2)= +cz*Eipq(:,1,3)+sz*Eipq(:,1,4)
        alp_ij(:,2,2)= +cz*Eipq(:,2,3)+sz*Eipq(:,2,4)
        eta_ij3(:,1,1)=(-sz*Eipq(:,1,1)+cz*Eipq(:,1,2))*ktestfield1
        eta_ij3(:,2,1)=(-sz*Eipq(:,2,1)+cz*Eipq(:,2,2))*ktestfield1
        eta_ij3(:,1,2)=(-sz*Eipq(:,1,3)+cz*Eipq(:,1,4))*ktestfield1
        eta_ij3(:,2,2)=(-sz*Eipq(:,2,3)+cz*Eipq(:,2,4))*ktestfield1
!
!  go from eta_ij3 to eta_ij
!
        eta_ij(:,1,1)=+eta_ij3(:,1,2)
        eta_ij(:,2,1)=-eta_ij3(:,1,1)
        eta_ij(:,1,2)=+eta_ij3(:,2,2)
        eta_ij(:,2,2)=-eta_ij3(:,2,1)
      endif
!
    endsubroutine calc_lspecial_pars
!***********************************************************************
    subroutine special_calc_density(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p

!!
!!  SAMPLE IMPLEMENTATION
!!     (remember one must ALWAYS add to df)
!!
!!
!!  df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + SOME NEW TERM
!!
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_density
!***********************************************************************
    subroutine special_calc_hydro(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p

!!
!!  SAMPLE IMPLEMENTATION
!!     (remember one must ALWAYS add to df)
!!
!!
!!  df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) + SOME NEW TERM
!!  df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) + SOME NEW TERM
!!  df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz) + SOME NEW TERM
!!
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_hydro
!***********************************************************************
    subroutine special_calc_magnetic(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p

!!
!!  SAMPLE IMPLEMENTATION
!!     (remember one must ALWAYS add to df)
!!
!!
!!  df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) + SOME NEW TERM
!!  df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) + SOME NEW TERM
!!  df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz) + SOME NEW TERM
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_magnetic
!!***********************************************************************
    subroutine special_calc_entropy(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p

!!
!!  SAMPLE IMPLEMENTATION
!!     (remember one must ALWAYS add to df)
!!
!!
!!  df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) + SOME NEW TERM
!!
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_entropy
!***********************************************************************
    subroutine special_boundconds(f,bc)
!
!   calculate a additional 'special' term on the right hand side of the
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      type (boundary_condition) :: bc
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(bc)
!
    endsubroutine special_boundconds
!***********************************************************************
    subroutine special_before_boundary(f)
!
!   Possibility to modify the f array before the boundaries are
!   communicated.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-jul-06/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine special_before_boundary
!***********************************************************************

!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../special_dummies.inc'
!********************************************************************

endmodule Special

