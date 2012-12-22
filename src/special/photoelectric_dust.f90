! $Id: baroclinic_run.f90 19193 2012-06-30 12:55:46Z wdobler $

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
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
!   Equation of state from Lyra & Kuchner (2012), due to photoelectric
!   heating of the gas by the dust. A polytropic term is added for
!   completeness and to provide support against gravity in stratified
!   simulations. The pressure gradient itself is
!
!    gradP = C*(p%rho*p%grhop(:,j) + p%rhop*p%grho(:,j))
!    C = cs20/(gamma*rho0)
!
!   Perhaps this whole thing should be moved to special, but then check
!   how to do with all the needed modifications to the equation of state.
!   Or should one use special with noeos instead? Hm...
!
!   17-may-12/wlad: coded
!   20-dec-12/wlad: streamlined
!
module Special
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include '../special.h'
!
  real :: dummy
  real :: mu=1.0, Sentropy=0.0
!
  namelist /special_init_pars/ dummy 
!
  namelist /special_run_pars/ mu, Sentropy
!
!  integer, parameter :: nmode_max = 50
!  real, dimension(nmode_max) :: gauss_ampl, rcenter, phicenter
!
  type InternalPencils
    real, dimension(nx,3) :: fpres
    real, dimension(nx,3) :: fpres_photoelectric
    real, dimension(nx,3) :: fpres_polytropic
  endtype InternalPencils
  type (InternalPencils) :: q

!
  real :: const1,const2
  integer :: idiag_fpresphotom=0,idiag_fpresphotomax=0
  integer :: idiag_fpresphotomin=0,idiag_fprespolym=0
  integer :: idiag_fprespolymax=0,idiag_fprespolymin=0
!
  contains

!***********************************************************************
    subroutine register_special()
!
!  Configure pre-initialised (i.e. before parameter read) variables 
!  which should be know to be able to evaluate
!
!  14-jul-09/wlad: coded
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f,lstarting)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  14-jul-09/wlad: coded
!
      use Cdata 
      use Mpicomm
      use EquationOfState, only: gamma1,cs20,rho0
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
      real :: rho01
!
      rho01=1./rho0
!
      const1=Sentropy*mu
      const2=cs20*gamma1*rho01
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine pencils_criteria_special(f,p)
!
      use Mpicomm
!
      real, dimension(mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      !const1=Sentropy*mu
      !const2=cs20*gamma1*rho01
!
      if (const1 /= 0.0) then 
        lpenc_requested(i_rho)=.true.
        lpenc_requested(i_glnrho)=.true.
      endif
!
      if (const2 /=0.0) then 
        lpenc_requested(i_glnrho)=.true.
        if (lparticles) then 
          lpenc_requested(i_rhop)=.true.
          lpenc_requested(i_glnrhop)=.true.
          lpenc_requested(i_glnrho)=.true.
        else
          lpenc_requested(i_rhodsum)=.true.
          lpenc_requested(i_glnrhodsum)=.true.
        endif
      endif
!
    endsubroutine pencils_criteria_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
      use Mpicomm
!
      real, dimension(mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      integer :: j
!
      do j=1,3
        q%fpres_polytropic(:,j) = -const1 * p%rho**(mu-1) * p%glnrho(:,j)
!
        if (lparticles) then
          q%fpres_photoelectric(:,j) = -const2 * p%rhop * (p%glnrhop(:,j) + p%glnrho(:,j))
        else
          q%fpres_photoelectric(:,j) = -const2 * p%rhodsum * (p%glnrhodsum(:,j) + p%glnrho(:,j))
        endif
      enddo
!
      q%fpres = q%fpres_photoelectric + q%fpres_polytropic
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine special_before_boundary(f,lstarting)
!
!  This subroutine calculates the full potential due to the turbulence.
!
!  03-oct-12/wlad: coded
!
      use Mpicomm
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      logical, optional :: lstarting
!!
      if (present(lstarting)) then 
        lstart = lstarting
      else
        lstart=.false.
      endif
!
    endsubroutine special_before_boundary
!***********************************************************************
    subroutine read_special_init_pars(unit,iostat)
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
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=special_init_pars)

    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=special_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=special_run_pars,ERR=99)
      endif
!
99    return
endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
      integer, intent(in) :: unit
!
      write(unit,NML=special_run_pars)
!
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
      use Diagnostics
!
!  reads and registers print parameters relevant to special
!
!   14-jul-09/wlad: coded
!
      integer :: iname
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  Write information to index.pro
!
      if (lreset) then
        idiag_fpresphotom=0;idiag_fpresphotomax=0;idiag_fpresphotomin=0
        idiag_fprespolym=0;idiag_fprespolymax=0;idiag_fprespolymin=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'fpresphotom',idiag_fpresphotom)
        call parse_name(iname,cname(iname),cform(iname),'fpresphotomax',idiag_fpresphotomax) 
        call parse_name(iname,cname(iname),cform(iname),'fpresphotomin',idiag_fpresphotomin) 
        call parse_name(iname,cname(iname),cform(iname),'fprespolym',idiag_fprespolym)
        call parse_name(iname,cname(iname),cform(iname),'fprespolymax',idiag_fprespolymax)
        call parse_name(iname,cname(iname),cform(iname),'fprespolymin',idiag_fprespolymin)
      enddo
!
      if (lwr) then
        write(3,*) 'i_fpresphotom',idiag_fpresphotom
        write(3,*) 'i_fpresphotomax',idiag_fpresphotomax
        write(3,*) 'i_fpresphotomin',idiag_fpresphotomin
        write(3,*) 'i_fprespolym',idiag_fprespolym
        write(3,*) 'i_fprespolymax',idiag_fprespolymax
        write(3,*) 'i_fprespolymin',idiag_fprespolymin
      endif
!
    endsubroutine rprint_special
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
        df(l1:l2,m,n,ju) = df(l1:l2,m,n,ju) - q%fpres(:,j)
      enddo
!
      if (ldiagnos) then 
        if (idiag_fpresphotom/=0)   call sum_mn_name( q%fpres_photoelectric,idiag_fpresphotom)
        if (idiag_fpresphotomax/=0) call max_mn_name( q%fpres_photoelectric,idiag_fpresphotomax)
        if (idiag_fpresphotomin/=0) call max_mn_name(-q%fpres_photoelectric,idiag_fpresphotomin)
        if (idiag_fprespolym/=0)      call sum_mn_name( q%fpres_polytropic,idiag_fprespolym)
        if (idiag_fprespolymax/=0)    call max_mn_name( q%fpres_polytropic,idiag_fprespolymax)
        if (idiag_fprespolymin/=0)    call max_mn_name(-q%fpres_polytropic,idiag_fprespolymin)
      endif
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_hydro
!***********************************************************************
!
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

