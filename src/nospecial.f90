! $Id: nospecial.f90,v 1.2 2003-10-09 11:14:27 theine Exp $

!  This module provide a way for users to specify custom (i.e. not in the standard Pencil Code)
!  physics, diagnostics etc. 
!
!  The module provides a set of standard hooks into the Pencil-Code and currently allows the following
!  customizations:                                        
!
!   Description                                         | Relevant function call (hook)
!  ------------------------------------------------------------------------------------------------
!   Special variable registration (pre parameter read)   | register_special
!   Special variable initialization (post parameter read)| initialize_special
!
!   Special initial condition                            | init_special
!     this is called last so may be used to modify
!     the mvar variables declared by this module
!     or optionally modify any of the other f array 
!     variables.  The latter, however, should be avoided
!     where ever possible.
!
!
!   Special term in the mass (density) equation          | special_calc_density
!   Special term in the momentum (hydro) equation        | special_calc_hydro
!   Special term in the entropy equation                 | special_calc_entropy
!   Special term in the induction (magnetic) equation    | special_calc_magnetic
!
!   Special equation                                     | dspecial_dt
!     NOT IMPLEMENTED FULLY YET - HOOKS NOT PLACED INTO THE PENCIL-CODE 

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxilliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Special

  use Cparam
  use Cdata

  implicit none
  
  real :: dummy
!!  character, len(50) :: initcustom

!! input parameters
  namelist /special_init_pars/ &
      dummy 
!eg.    initcustom
     
!! run parameters
  namelist /special_run_pars/ &
      dummy

!!
!! Declare any index variables necessary for main or 
!! 
!!   integer :: iSPECIAL_VARIABLE_INDEX=0
!!  
!! other variables (needs to be consistent with reset list below)
!!
!!   integer :: i_POSSIBLEDIAGNOSTIC=0
!!

  

  contains

!***********************************************************************
    subroutine register_special()
!
!  Configure pre-initialised (i.e. before parameter read) variables 
!  which should be know to be able to evaluate
! 
!
!  6-oct-03/tony: coded
!
      use Cdata
      use Mpicomm
      use Sub, only: cvs_id
!
      logical, save :: first=.true.
!
! A quick sanity check
!
      if (.not. first) call stop_it('register_special called twice')
      first = .false.

!!
!! MUST SET lspecial = .true. to enable use of special hooks in the Pencil-Code 
!!
      lspecial = .false.
!
!
!
!! 
!! Set any required f-array indexes to the next available slot 
!!  
!!
!      iSPECIAL_VARIABLE_INDEX = nvar+1             ! index to access entropy
!      nvar = nvar+1
!
!      iSPECIAL_AUXILLIARY_VARIABLE_INDEX = naux+1             ! index to access entropy
!      naux = naux+1
!
!
!  identify CVS version information (if checked in to a CVS repository!)
!  CVS should automatically update everything between $Id: nospecial.f90,v 1.2 2003-10-09 11:14:27 theine Exp $ 
!  when the file in committed to a CVS repository.
!
      if (lroot) call cvs_id( &
           "$Id: nospecial.f90,v 1.2 2003-10-09 11:14:27 theine Exp $")
!
!
!  Perform some sanity checks (may be meaningless if certain things haven't 
!  been configured in a custom module but they do no harm)
!
      if (naux > maux) then
        if (lroot) write(0,*) 'naux = ', nvar, ', maux = ', mvar
        call stop_it('register_special: naux > maux')
      endif
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('register_special: nvar > mvar')
      endif
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special()
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  06-oct-03/tony: coded
!
      use Cdata
!
!!
!!  Initialize any module variables which are parameter dependant  
!!
!
! DO NOTHING
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f,xx,yy,zz)
!
!  initialise special condition; called from start.f90
!  06-oct-2003/tony: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
!
      intent(in) :: xx,yy,zz
      intent(inout) :: f

!!
!!  SAMPLE IMPLEMENTATION
!!
!!      select case(initspecial)
!!        case('nothing'); if(lroot) print*,'init_special: nothing'
!!        case('zero', '0'); f(:,:,:,iSPECIAL_VARIABLE_INDEX) = 0.
!!        case default
!!          !
!!          !  Catch unknown values
!!          !
!!          if (lroot) print*,'init_special: No such value for initspecial: ', trim(initspecial)
!!          call stop_it("")
!!      endselect
!
      if(ip==0) print*,xx(1,1,1),yy(1,1,1),zz(1,1,1)  !(to keep compiler quiet)
!
    endsubroutine init_special
!***********************************************************************
    subroutine dspecial_dt(f,df,uu,glnrho,divu,rho1,lnrho,cs2,TT1)
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
      use Mpicomm
      use Sub
      use Global
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df     
      real, dimension (nx,3) :: uu,glnrho
      real, dimension (nx) :: divu, lnrho,rho1,cs2,TT1
!
      intent(in) :: f,uu,glnrho,rho1,lnrho,cs2,TT1
      intent(inout) :: df
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dSPECIAL_dt'
!!      if (headtt) call identify_bcs('ss',iss)
!
!!
!! SAMPLE DIAGNOSTIC IMPLEMENTATION
!!
!!      if(ldiagnos) then
!!        if(i_SPECIAL_DIAGNOSTIC/=0) then
!!          call sum_mn_name(SOME MATHEMATICAL EXPRESSION,i_SPECIAL_DIAGNOSTIC)
!!! see also integrate_mn_name
!!        endif
!!      endif
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine rprint_special(lreset)
!
!  reads and registers print parameters relevant to special
!
!   06-oct-03/tony: coded
!
      use Cdata
      use Sub

      logical :: lreset

!!
!!!   SAMPLE IMPLEMENTATION
!!
!!      integer :: iname
!!!
!!!  reset everything in case of reset
!!!  (this needs to be consistent with what is defined above!)
!!!
      if (lreset) then
!!        i_SPECIAL_DIAGNOSTIC=0
      endif
!!
!!      do iname=1,nname
!!        call parse_name(iname,cname(iname),cform(iname),'NAMEOFSPECIALDIAGNOSTIC',i_SPECIAL_DIAGNOSTIC)
!!      enddo
!!
!!!  write column where which magnetic variable is stored
!!
!!      write(3,*) 'i_SPECIAL_DIAGNOSTIC=',i_SPECIAL_DIAGNOSTIC
!!

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
      if (ip==0) print*,f(1,1,1,1),df(1,1,1,1),uu(1),glnrho(1),divu(1), &
                      lnrho(1)
                  
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
      if (ip==0) print*,f(1,1,1,1),df(1,1,1,1),uu(1),glnrho(1),divu(1), &
                      rho1(1), u2(1), uij(1)
                  
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
      if (ip==0) print*,f(1,1,1,1),df(1,1,1,1),uu(1), &
                      TT1(1), uij(1)
                  
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
      if (ip==0) print*,f(1,1,1,1),df(1,1,1,1),uu(1),glnrho(1),divu(1), &
                      rho1(1), lnrho(1), cs2(1), TT1(1)
                  
    endsubroutine special_calc_entropy
!***********************************************************************
!***********************************************************************
endmodule Special

