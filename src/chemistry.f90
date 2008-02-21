! $Id: chemistry.f90,v 1.9 2008-02-21 06:52:46 brandenb Exp $
!  This modules addes chemical species and reactions.

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lchemistry = .true.
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Chemistry

  use Cparam
  use Cdata
  use Messages
  use Sub, only: keep_compiler_quiet

  implicit none

  include 'chemistry.h'
!
!  parameters related to chemical reactions
!
  logical :: lreactions=.true.
  integer :: nreactions=0
  integer, parameter :: mreactions=2*nchemspec
  integer, dimension(nchemspec,mreactions) :: stoichio,Sijm
  real, dimension(mreactions) :: kreactions
!
!  hydro-related parameters
!
  real, dimension(nchemspec) :: amplchemk=0.
  real :: amplchem=1.,kx_chem=1.,ky_chem=1.,kz_chem=1.,widthchem=1.
  real :: chem_diff=0.
  character (len=labellen), dimension (ninit) :: initchem='nothing'

! input parameters
  namelist /chemistry_init_pars/ &
      initchem, amplchem, kx_chem, ky_chem, kz_chem, widthchem, &
      amplchemk

! run parameters
  namelist /chemistry_run_pars/ &
      chem_diff
!
! diagnostic variables (need to be consistent with reset list below)
!
  integer :: idiag_Y1m=0        ! DIAG_DOC: $\left<Y_1\right>$
  integer :: idiag_Y2m=0        ! DIAG_DOC: $\left<Y_2\right>$
  integer :: idiag_Y3m=0        ! DIAG_DOC: $\left<Y_3\right>$
  integer :: idiag_Y4m=0        ! DIAG_DOC: $\left<Y_4\right>$
  integer :: idiag_Y5m=0        ! DIAG_DOC: $\left<Y_5\right>$
  integer :: idiag_Y6m=0        ! DIAG_DOC: $\left<Y_6\right>$
  integer :: idiag_Y7m=0        ! DIAG_DOC: $\left<Y_7\right>$
  integer :: idiag_Y8m=0        ! DIAG_DOC: $\left<Y_8\right>$
!
  contains

!***********************************************************************
    subroutine register_chemistry()
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!  13-aug-07/steveb: coded
!   8-jan-08/axel: added modifications analogously to dustdensity
!
      use Cdata
      use Mpicomm
      use General, only: chn
!
      logical, save :: first=.true.
      integer :: k
      character (len=5) :: schem
!
! A quick sanity check
!
      if (.not. first) call stop_it('register_chemistry called twice')
      first=.false.
!
!  Set ind to consecutive numbers nvar+1, nvar+2, ..., nvar+nchemspec
!
      do k=1,nchemspec
        ichemspec(k)=nvar+k
      enddo
!
!  Increase nvar accordingly
!
      nvar=nvar+nchemspec
!
!  Print some diagnostics
!
      do k=1,nchemspec
        if ((ip<=8) .and. lroot) then
          print*, 'register_chemistry: k = ', k
          print*, 'register_chemistry: nvar = ', nvar
          print*, 'register_chemistry: ichemspec = ', ichemspec(k)
        endif
!
!  Put variable name in array
!
        call chn(k,schem)
        varname(ichemspec(k))='nd('//trim(schem)//')'
      enddo
!
!  identify CVS version information (if checked in to a CVS repository!)
!  CVS should automatically update everything between $Id: chemistry.f90,v 1.9 2008-02-21 06:52:46 brandenb Exp $
!  when the file in committed to a CVS repository.
!
      if (lroot) call cvs_id( &
           "$Id: chemistry.f90,v 1.9 2008-02-21 06:52:46 brandenb Exp $")
!
!
!  Perform some sanity checks (may be meaningless if certain things haven't
!  been configured in a custom module but they do no harm)
!
      if (naux > maux) then
        if (lroot) write(0,*) 'naux = ', naux, ', maux = ', maux
        call stop_it('register_chemistry: naux > maux')
      endif
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('register_chemistry: nvar > mvar')
      endif
!
    endsubroutine register_chemistry
!***********************************************************************
    subroutine initialize_chemistry(f)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  13-aug-07/steveb: coded
!  19-feb-08/axelb: reads in chemistry.dat file
!
      use Cdata
      use Sub, only: keep_compiler_quiet
!
      character (len=80) :: chemicals=''
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: exist
      integer :: i,j
!
!  read chemistry data, if present
!
      inquire(file='chemistry.dat',exist=exist)
      if(exist) then
        open(19,file='chemistry.dat')
        read(19,*) chemicals
        do j=1,mreactions
          read(19,*,end=990) kreactions(j),(stoichio(i,j),i=1,nchemspec)
        enddo
990     nreactions=j-1
!
!  calculate negative part of stoichiometric matrix
!
        Sijm=-min(stoichio,0)
!
!  print input data for verification
!
        if (lroot) then
          print*,'chemicals=',chemicals
          print*,'kreactions=',kreactions(1:nreactions)
          print*,'stoichio=' ; write(*,100),stoichio(:,1:nreactions)
          print*,'Sijm:' ; write(*,100),Sijm(:,1:nreactions)
100       format(8i4)
        endif
      else 
        if (lroot) print*,'no chemistry.dat file to be read.'
        lreactions=.false.
      endif
      close(19)
!
!  that's it
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_chemistry
!***********************************************************************
    subroutine init_chemistry(f,xx,yy,zz)
!
!  initialise chemistry initial condition; called from start.f90
!  13-aug-07/steveb: coded
!
      use Cdata
      use Initcond
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
      integer :: j,k
      logical :: lnothing
!
      intent(in) :: xx,yy,zz
      intent(inout) :: f
!
!  different initializations of nd (called from start)
!
      lnothing=.false.
      do j=1,ninit
        select case(initchem(j))

        case('nothing')
          if (lroot .and. .not. lnothing) print*, 'init_chem: nothing'
          lnothing=.true.
        case('constant')
          do k=1,nchemspec
            f(:,:,:,ichemspec(k))=amplchemk(k)
          enddo
        case('coswave-x')
          do k=1,nchemspec
            call coswave(amplchem,f,ichemspec(k),kx=kx_chem)
          enddo
        case('hatwave-x')
          do k=1,nchemspec
            call hatwave(amplchem,f,ichemspec(k),kx=kx_chem)
          enddo
        case('hatwave-y')
          do k=1,nchemspec
            call hatwave(amplchem,f,ichemspec(k),ky=ky_chem)
          enddo
        case('hatwave-z')
          do k=1,nchemspec
            call hatwave(amplchem,f,ichemspec(k),kz=kz_chem)
          enddo
        case default
!
!  Catch unknown values
!
          if (lroot) print*, 'initchem: No such value for initchem: ', &
              trim(initchem(j))
          call stop_it('')

        endselect
      enddo
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(xx,yy,zz)
!
    endsubroutine init_chemistry
!***********************************************************************
    subroutine pencil_criteria_chemistry()
!
!  All pencils that this chemistry module depends on are specified here.
!
!  13-aug-07/steveb: coded
!
    endsubroutine pencil_criteria_chemistry
!***********************************************************************
    subroutine pencil_interdep_chemistry(lpencil_in)
!
!  Interdependency among pencils provided by this module are specified here
!
!  13-aug-07/steveb: coded
!
      use Sub, only: keep_compiler_quiet
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_chemistry
!***********************************************************************
    subroutine calc_pencils_chemistry(f,p)
!
!  Calculate Hydro pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!   13-aug-07/steveb: coded
!
      use Cdata
      use Sub, only: keep_compiler_quiet
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
    endsubroutine calc_pencils_chemistry
!***********************************************************************
    subroutine dchemistry_dt(f,df,p)
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
!   13-aug-07/steveb: coded
!    8-jan-08/natalia: included advection/diffusion
!   20-feb-08/axel: included reactions
!
      use Cdata
      use Mpicomm
      use Sub
      use Global
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
     
      real, dimension (nx,3) :: gchemspec
      real, dimension (nx) :: ugchemspec,del2chemspec,diff_op,xdot
      real, dimension (nx,mreactions) :: vreactions
      type (pencil_case) :: p
!
!  indices
!
      integer :: j,k
!
      intent(in) :: f,p
      intent(inout) :: df
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'dchemistry_dt: SOLVE dchemistry_dt'
!!      if (headtt) call identify_bcs('ss',iss)
!
!  if we do reactions, we must calculate the reaction speed vector
!  outside the loop where we multiply it by the stoichiometric matrix
!
      if (lreactions) then
        do j=1,nreactions
          vreactions(:,j)=kreactions(j)
!print*,'1) j,v=',j,vreactions(1,j)
          do k=1,nchemspec
            vreactions(:,j)=vreactions(:,j)*f(l1:l2,m,n,ichemspec(k))**Sijm(k,j)
!print*,'2) x,S=',f(l1:l2,m,n,ichemspec(k)),Sijm(k,j)
!print*,'3) k,v=',j,vreactions(1,j)
          enddo
!print*,'4) j,v=',j,vreactions(1,j)
          enddo
      endif
!
!  loop over all chemicals
!
      do k=1,nchemspec
!
!  advection terms
!
        call grad(f,ichemspec(k),gchemspec) 
        call dot_mn(p%uu,gchemspec,ugchemspec)
        df(l1:l2,m,n,ichemspec(k))=df(l1:l2,m,n,ichemspec(k))-ugchemspec
!
!  diffusion operator
!
        if (chem_diff/=0.) then
          call del2(f,ichemspec(k),del2chemspec) 
          if (headtt) print*,'dchemistry_dt: chem_diff=',chem_diff
          call dot_mn(p%glnrho,gchemspec,diff_op)
          diff_op=diff_op+del2chemspec
          df(l1:l2,m,n,ichemspec(k))=df(l1:l2,m,n,ichemspec(k))+chem_diff*diff_op
        endif
!
!  chemical reactions:
!  multiply with stoichiometric matrix with reaction speed
!  d/dt(x_i) = S_ij v_j
!
        if (lreactions) then
          xdot=0.
          do j=1,nreactions
            xdot=xdot+stoichio(k,j)*vreactions(:,j)
!print*,'5) j,Skj,vj,xdot=',j,stoichio(k,j),vreactions(:,j),xdot
          enddo
          df(l1:l2,m,n,ichemspec(k))=df(l1:l2,m,n,ichemspec(k))+xdot
        endif
!
      enddo 
!
!  For the timestep calculation, need maximum diffusion
!
        if (lfirst.and.ldt) then
          diffus_chem=chem_diff*dxyz_2
        endif
!
!  Calculate diagnostic quantities
!
      if (ldiagnos) then
        if (idiag_Y1m/=0) call sum_mn_name(f(l1:l2,m,n,ichemspec(1)),idiag_Y1m)
        if (idiag_Y2m/=0) call sum_mn_name(f(l1:l2,m,n,ichemspec(2)),idiag_Y2m)
        if (idiag_Y3m/=0) call sum_mn_name(f(l1:l2,m,n,ichemspec(3)),idiag_Y3m)
        if (idiag_Y4m/=0) call sum_mn_name(f(l1:l2,m,n,ichemspec(4)),idiag_Y4m)
        if (idiag_Y5m/=0) call sum_mn_name(f(l1:l2,m,n,ichemspec(5)),idiag_Y5m)
        if (idiag_Y6m/=0) call sum_mn_name(f(l1:l2,m,n,ichemspec(6)),idiag_Y6m)
        if (idiag_Y7m/=0) call sum_mn_name(f(l1:l2,m,n,ichemspec(7)),idiag_Y7m)
        if (idiag_Y8m/=0) call sum_mn_name(f(l1:l2,m,n,ichemspec(8)),idiag_Y8m)
      endif
!
! Keep compiler quiet by ensuring every parameter is used
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)

    endsubroutine dchemistry_dt
!***********************************************************************
    subroutine read_chemistry_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=chemistry_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=chemistry_init_pars,ERR=99)
      endif
!
99    return
    endsubroutine read_chemistry_init_pars
!***********************************************************************
    subroutine write_chemistry_init_pars(unit)
!
      integer, intent(in) :: unit

      write(unit,NML=chemistry_init_pars)

    endsubroutine write_chemistry_init_pars
!***********************************************************************
    subroutine read_chemistry_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=chemistry_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=chemistry_run_pars,ERR=99)
      endif
!
99    return
    endsubroutine read_chemistry_run_pars
!***********************************************************************
    subroutine write_chemistry_run_pars(unit)
!
      integer, intent(in) :: unit

      write(unit,NML=chemistry_run_pars)

    endsubroutine write_chemistry_run_pars
!***********************************************************************
    subroutine rprint_chemistry(lreset,lwrite)
!
!  reads and registers print parameters relevant to chemistry
!
!  13-aug-07/steveb: coded
!
      use Cdata
      use Sub
      use General, only: chn
!
      integer :: iname
      logical :: lreset,lwr
      logical, optional :: lwrite
      character (len=5) :: schem,schemspec,snd1,smd1,smi1
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_Y1m=0; idiag_Y2m=0; idiag_Y3m=0; idiag_Y4m=0
        idiag_Y5m=0; idiag_Y6m=0; idiag_Y7m=0; idiag_Y8m=0
      endif
!
      call chn(nchemspec,schemspec)
!
!  check for those quantities that we want to evaluate online
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'Y1m',idiag_Y1m)
        call parse_name(iname,cname(iname),cform(iname),'Y2m',idiag_Y2m)
        call parse_name(iname,cname(iname),cform(iname),'Y3m',idiag_Y3m)
        call parse_name(iname,cname(iname),cform(iname),'Y4m',idiag_Y4m)
        call parse_name(iname,cname(iname),cform(iname),'Y5m',idiag_Y5m)
        call parse_name(iname,cname(iname),cform(iname),'Y6m',idiag_Y6m)
        call parse_name(iname,cname(iname),cform(iname),'Y7m',idiag_Y7m)
        call parse_name(iname,cname(iname),cform(iname),'Y8m',idiag_Y8m)
      enddo
!
!  Write chemistry index in short notation
!
      call chn(ichemspec(1),snd1)
      if (lwr) then
        write(3,*) 'i_Y1m=',idiag_Y1m
        write(3,*) 'i_Y2m=',idiag_Y2m
        write(3,*) 'i_Y3m=',idiag_Y3m
        write(3,*) 'i_Y4m=',idiag_Y4m
        write(3,*) 'i_Y5m=',idiag_Y5m
        write(3,*) 'i_Y6m=',idiag_Y6m
        write(3,*) 'i_Y7m=',idiag_Y7m
        write(3,*) 'i_Y8m=',idiag_Y8m
        write(3,*) 'ichemspec=indgen('//trim(schemspec)//') + '//trim(snd1)
      endif
!
    endsubroutine rprint_chemistry
!***********************************************************************
    subroutine get_slices_chemistry(f,slices)
!
!  Write slices for animation of chemistry variables.
!
!  13-aug-07/steveb: dummy
!
      use Sub, only: keep_compiler_quiet
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_chemistry
!***********************************************************************
    subroutine special_calc_density(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   13-aug-07/steveb: coded
!
      use Cdata
      use Sub, only: keep_compiler_quiet
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
!   13-aug-07/steveb: coded
!
      use Cdata
      use Sub, only: keep_compiler_quiet

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
!   13-aug-07/steveb: coded
!
      use Cdata
      use Sub, only: keep_compiler_quiet

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
!   13-aug-07/steveb: coded
!
      use Cdata
      use Sub, only: keep_compiler_quiet

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
    subroutine chemistry_boundconds(f,bc)
!
!   calculate a additional 'special' term on the right hand side of the
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   13-aug-07/steveb: coded
!
      use Cdata
      use Sub, only: keep_compiler_quiet
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      type (boundary_condition) :: bc
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(bc)
!
    endsubroutine chemistry_boundconds
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
      use Cdata
      use Sub, only: keep_compiler_quiet
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
    include 'special_dummies.inc'
!********************************************************************

endmodule Chemistry

