! $Id: meanfield_alpm.f90,v 1.2 2005-06-26 17:34:16 eos_merger_tony Exp $

!  This modules solves the passive scalar advection equation

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Special

  use Cparam
  use Cdata

  implicit none

  include 'special.inc'
 !! include 'meanfield_alpm.inc'

  character (len=labellen) :: initalpm='zero'

  ! input parameters
  real :: amplalpm=.1
  real :: kx_alpm=1.,ky_alpm=1.,kz_alpm=1.

  namelist /alpm_init_pars/ &
       initalpm,amplalpm,kx_alpm,ky_alpm,kz_alpm

  ! run parameters
  real :: etat_alpm=0., Rm_alpm=0., kf_alpm=1.

  namelist /alpm_run_pars/ &
       etat_alpm,Rm_alpm,kf_alpm

  ! other variables (needs to be consistent with reset list below)
  integer :: idiag_ammax=0,idiag_amrms=0,idiag_alpmmz=0

  contains

!***********************************************************************
    subroutine register_special()
!
!  Initialise variables which should know that we solve for passive
!  scalar: ialpm; increase nvar accordingly
!
!  6-jul-02/axel: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_alpm called twice')
      first = .false.
!
      lalpm = .true.
      ialpm = nvar+1            ! index to access alpm
      nvar = nvar+1             ! added 1 variable
!
      if ((ip<=8) .and. lroot) then
        print*, 'Register_alpm:  nvar = ', nvar
        print*, 'ialpm = ', ialpm
      endif
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: meanfield_alpm.f90,v 1.2 2005-06-26 17:34:16 eos_merger_tony Exp $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('Register_alpm: nvar > mvar')
      endif
!
!  Put variable name in array
!
      varname(ialpm) = 'alpm'
!
!  Writing files for use with IDL
!
      if (lroot) then
        if (maux == 0) then
          if (nvar < mvar) write(4,*) ',alpm $'
          if (nvar == mvar) write(4,*) ',alpm'
        else
          write(4,*) ',alpm $'
        endif
        write(15,*) 'alpm = fltarr(mx,my,mz)*one'
      endif
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  Perform any necessary post-parameter read initialization
!  Dummy routine
!
!  24-nov-02/tony: coded
!
      real, dimension (mx,my,mz,mvar+maux) :: f
! 
!  set to zero and then call the same initial condition
!  that was used in start.csh
!   
      if(NO_WARN) print*,'f=',f
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f,xx,yy,zz)
!
!  initialise passive scalar field; called from start.f90
!
!   6-jul-2001/axel: coded
!
      use Cdata
      use Mpicomm
      use Sub
      use Initcond
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz)      :: xx,yy,zz,prof
!
      select case(initalpm)
        case('zero'); f(:,:,:,ialpm)=0.
        case('constant'); f(:,:,:,ialpm)=amplalpm
        case default; call stop_it('init_alpm: bad initalpm='//trim(initalpm))
      endselect
!
      if(NO_WARN) print*,xx,yy,zz !(prevent compiler warnings)
    endsubroutine init_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
!  Calculate Hydro pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!   24-nov-04/tony: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f       
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!     
      if(NO_WARN) print*,f(1,1,1,1),p   !(keep compiler quiet)
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine dspecial_dt(f,df,p)
!
!  dynamical alpha quenching equation
!  dalpm/dt=-2*etat*kf2*(EMF*BB/Beq2+alpm/Rm)
!
!  18-nov-04/axel: coded
!
      use Sub
!     use Magnetic, only: meanfield_EMFdotB
!     use Magnetic, only: eta
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: alpm
      type (pencil_case) :: p
      integer :: j
!
      intent(in)  :: f
      intent(out) :: df
!
!  identify module and boundary conditions
!
      if (headtt) call identify_bcs('alpm',ialpm)
!
!  Abbreviations
!        
      alpm=f(l1:l2,m,n,ialpm)
!
!  dynamical quenching equation
!
      if(lmagnetic) df(l1:l2,m,n,ialpm)=df(l1:l2,m,n,ialpm)&
         -2*etat_alpm*kf_alpm**2*(p%mf_EMFdotB+alpm/Rm_alpm)
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine read_special_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=alpm_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=alpm_init_pars,ERR=99)
      endif

99    return
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=alpm_init_pars)

    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=alpm_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=alpm_run_pars,ERR=99)
      endif

99    return
    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=alpm_run_pars)

    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
!  reads and registers print parameters relevant for magnetic fields
!
!   6-jul-02/axel: coded
!
      use Sub
!
      integer :: iname,inamez
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_ammax=0; idiag_amrms=0; idiag_alpmmz=0
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'ammax',idiag_ammax)
        call parse_name(iname,cname(iname),cform(iname),'amrms',idiag_amrms)
      enddo
!
!  check for those quantities for which we want xy-averages
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),&
            'alpmmz',idiag_alpmmz)
      enddo
!
!  write column where which magnetic variable is stored
!
      if (lwr) then
        write(3,*) 'ialpm=',ialpm
      endif
!
    endsubroutine rprint_special
!***********************************************************************
    subroutine special_calc_density(f,df,uu,glnrho,divu,lnrho)
!
!   calculate a additional 'special' term on the right hand side of the 
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      use Cdata
      
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      real, dimension (nx), intent(in) :: uu,glnrho,divu,lnrho 

!!
!!  SAMPLE IMPLEMENTATION
!!     (remember one must ALWAYS add to df)
!!  
!!
!!  df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + SOME NEW TERM
!!
!!

! Keep compiler quiet by ensuring every parameter is used
      if (NO_WARN) print*,f,df,uu,glnrho,divu,lnrho

    endsubroutine special_calc_density
!***********************************************************************
    subroutine special_calc_hydro(f,df,uu,glnrho,divu,rho1,u2,uij)
!
!   calculate a additional 'special' term on the right hand side of the 
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      use Cdata
      
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      real, dimension (nx), intent(in) :: uu,glnrho,divu,rho1,u2,uij 

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

! Keep compiler quiet by ensuring every parameter is used
      if (NO_WARN) print*,f,df,uu,glnrho,divu,rho1,u2,uij

    endsubroutine special_calc_hydro
!***********************************************************************
    subroutine special_calc_magnetic(f,df,uu,rho1,TT1,uij)
!
!   calculate a additional 'special' term on the right hand side of the 
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      use Cdata
      
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      real, dimension (nx), intent(in) :: uu,rho1,TT1,uij 

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

! Keep compiler quiet by ensuring every parameter is used
      if (NO_WARN) print*,f,df,uu,TT1,uij,rho1

    endsubroutine special_calc_magnetic
!!***********************************************************************
    subroutine special_calc_entropy(f,df,uu,glnrho,divu,rho1,lnrho,cs2,TT1)
!
!   calculate a additional 'special' term on the right hand side of the 
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      use Cdata
      
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      real, dimension (nx), intent(in) :: uu,glnrho,divu,rho1,lnrho,cs2,TT1 

!!
!!  SAMPLE IMPLEMENTATION
!!     (remember one must ALWAYS add to df)
!!  
!!
!!  df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) + SOME NEW TERM
!!
!!

! Keep compiler quiet by ensuring every parameter is used
      if (NO_WARN) print*,f,df,uu,glnrho,divu,rho1,lnrho,cs2,TT1

    endsubroutine special_calc_entropy
!***********************************************************************

endmodule Special

