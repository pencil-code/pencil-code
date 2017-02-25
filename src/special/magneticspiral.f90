! $Id: stellartide.f90,v 1.1 2014/03/13 16:31:16 wlyra Exp $

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of special_dummies.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************

!-------------------------------------------------------------------
!
! HOW TO USE THIS FILE
! --------------------
!
! The rest of this file may be used as a template for your own
! special module.  Lines which are double commented are intended
! as examples of code.  Simply fill out the prototypes for the 
! features you want to use.
!
! Save the file with a meaningful name, eg. geo_kws.f90 and place
! it in the $PENCIL_HOME/src/special directory.  This path has
! been created to allow users ot optionally check their contributions
! in to the Pencil-Code CVS repository.  This may be useful if you
! are working on/using the additional physics with somebodyelse or 
! may require some assistance from one of the main Pencil-Code team.
!
! To use your additional physics code edit the Makefile.local in
! the src directory under the run directory in which you wish to
! use your additional physics.  Add a line with all the module 
! selections to say something like:
!
!    SPECIAL=special/nstar
!
! Where nstar it replaced by the filename of your new module
! upto and not including the .f90
!
!--------------------------------------------------------------------

!
!  14-feb-17/wlad: coded
!
module Special
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
!
  implicit none
!
  include '../special.h'
!
  real :: B0= 0.0064252940531717906
  real :: etamu0 = 0.6
  real :: Omega0=1.0, r0=1.0
  real, dimension(mx,my,mz) :: brad,bphi
!
  namelist /special_init_pars/ B0,etamu0,Omega0,r0
  namelist /special_run_pars/ B0,etamu0,Omega0,r0
!
  type InternalPencils
     real, dimension(nx)   :: beta
     real, dimension(nx,3)   :: jxbr
  endtype InternalPencils
!
  type (InternalPencils) :: q
!
  integer :: idiag_qbetam=0
  integer :: idiag_qbetamin=0
  integer :: idiag_qbetamax=0
!
contains
!
!***********************************************************************
    subroutine register_special()
!
!  Configure pre-initialised (i.e. before parameter read) variables 
!  which should be know to be able to evaluate
!
!  14-jul-09/wlad: coded
!
      use Cdata
      use FArrayManager, only: farray_register_auxiliary
!
      if (lroot) call svn_id( &
           "$Id: stellartide.f90,v 1.1 2014/03/13 16:31:16 wlyra Exp $")
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  14-jul-09/wlad: coded
!
      use Cdata 
      use Mpicomm
      use EquationOfState, only: cs0
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      do n=1,mz
        do m=1,my
          brad(:,m,n) = B0 * r0/x
          bphi(:,m,n) = -2*B0*Omega0*r0**2 * (1./etamu0) * sqrt(r0/x)
        enddo
      enddo
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_special
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
      integer, intent(in) :: unit

      write(unit,NML=special_init_pars)

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
      integer, intent(in) :: unit
!
      write(unit,NML=special_run_pars)
!
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine pencil_criteria_special()
!
      use Mpicomm
!
      lpenc_requested(i_rho1)=.true.
      if (ldiagnos) lpenc_requested(i_cs2)=.true.
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
!  Add tides and Coriolis force to gas.
!
      use Sub, only: curl_mn,cross_mn,multsv_mn
      use Deriv, only: der
!
      real, dimension(mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      real, dimension (nx,3,3) :: bij
      real, dimension (nx,3) :: bb,jj,jxb
      real, dimension (nx) :: tmp,va2
      integer :: j
!
! jxbr
!
      bb(:,1) = brad(l1:l2,m,n)
      bb(:,2) = bphi(l1:l2,m,n)
      bb(:,3) = 0. 
!      
      do j=1,3
        call der(brad,tmp,j)
        bij(:,1,j)=tmp
        call der(bphi,tmp,j)
        bij(:,2,j)=tmp
        bij(:,3,j)=0. 
      enddo
!           
      call curl_mn(bij,jj,A=bb,LCOVARIANT_DERIVATIVE=.false.)
!
      call cross_mn(jj,bb,jxb)
      call multsv_mn(p%rho1,jxb,q%jxbr)
!
      if (ldiagnos) then 
        va2 = (brad(l1:l2,m,n)**2 + bphi(l1:l2,m,n)**2)*p%rho1/mu0
        q%beta = 2*p%cs2/va2
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine special_calc_hydro(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the 
!   momentum equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
      use Cdata
      use Diagnostics
!      
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
      integer :: j,ju
!
!  Modified momentum equation
!
      do j=1,3 
        ju=j+iuu-1
        df(l1:l2,m,n,ju) = df(l1:l2,m,n,ju) + q%jxbr(:,j)
      enddo
!
!  Diagnostics
!      
      if (ldiagnos) then
        if (idiag_qbetam/=0) call sum_mn_name(q%beta,idiag_qbetam)
        if (idiag_qbetamin/=0) call max_mn_name(q%beta,idiag_qbetamax)
        if (idiag_qbetamax/=0) call max_mn_name(-q%beta,idiag_qbetamin,lneg=.true.)
      endif
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_hydro
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
!  Reads and registers print parameters relevant to special.
!
!  06-oct-03/tony: coded
!
      use Diagnostics, only: parse_name
!
      logical :: lreset
      logical, optional :: lwrite
!
      logical :: lwr
      integer :: iname, inamex
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  Reset everything in case of reset.
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_qbetam=0
        idiag_qbetamin=0
        idiag_qbetamax=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'qbetam',idiag_qbetam)
        call parse_name(iname,cname(iname),cform(iname),'qbetamax',idiag_qbetamax)
        call parse_name(iname,cname(iname),cform(iname),'qbetamin',idiag_qbetamin)
      enddo
!
    endsubroutine rprint_special
!***********************************************************************
    
!***********************************************************************
!***********************************************************************
!********************************************************************
!
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

