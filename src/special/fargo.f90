! $Id$

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 1
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

module Special
!
  use Cparam
  use Cdata
  use Messages
!
  implicit none
!
  include '../special.h'
! Global arrays
  real, dimension (nx,nz) :: uu_average
! "pencils" 
  real, dimension (nx,3) :: uuadvec_gu
  real, dimension (nx) :: uuadvec_grho,uuadvec_glnrho
  real, dimension (nx) :: uu_residual
  real :: dummy
!
  namelist /special_init_pars/ dummy
!   
  namelist /special_run_pars/ dummy
!
  integer :: idiag_nshift=0
  integer :: idiag_dtuf=0
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
!
      if (lroot) call cvs_id( &
           "$Id$")
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  06-oct-03/tony: coded
!
      use Cdata 
      use Mpicomm
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!   
      if (.not.lfargo_advection) then
        print*,""
        print*,"Switch"
        print*," lfargo_advection=T"
        print*,"in init_pars of start.in if you"
        print*,"want to use the fargo algorithm"
        print*,""
        call stop_it("")
      endif
!
!  Not implemented for other than cylindrical coordinates
!
      if (coord_system/='cylindric') then
        print*,""
        print*,"Fargo advection is only implemented for"
        print*,"cylindrical coordinates. Switch"
        print*," coord_system='cylindric'"
        print*,"in init_pars of start.in if you"
        print*,"want to use the fargo algorithm"
        print*,""
        call stop_it("")
      endif
!
!  Not implemented for 3D runs
!
      if (nzgrid/=1) then
        print*,""
        print*,"You are using nzgrid/=1."
        print*,"By now the fargo algorithm is only"
        print*,"implemented for 2D runs. Stop and check"
        print*,""
        call stop_it("")
      endif
!
!  Not implemented for the induction equation
!
      if (lmagnetic) call stop_it("fargo advection not implemented "//&
          "for the induction equation")
!
!  Not implemented for the energy equation either
!
      if (lentropy.or.ltemperature) call stop_it("fargo advection not "//&
          "implemented for the energy equation")
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f)
!
!  initialise special condition; called from start.f90
!  06-oct-2003/tony: coded
!
      use Cdata
      use Mpicomm
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      intent(inout) :: f
!
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special()
! 
!  All pencils that this special module depends on are specified here.
! 
!  18-07-06/tony: coded
!
      lpenc_requested(i_uu)=.true.
!
!  For continuity equation
!
      lpenc_requested(i_divu)=.true.
      if (ldensity_nolog) then
        lpenc_requested(i_grho)=.true.
      else
        lpenc_requested(i_glnrho)=.true.
      endif
!      
!  For velocity advection
!
      lpenc_requested(i_uij)=.true.
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
      use Sub, only: h_dot_grad
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(nx,3) :: uu_advec,tmp2
      real, dimension(nx) :: tmp
      type (pencil_case) :: p
      integer :: j,nn
!
      nn=n-nghost
!
! Advect by the relative velocity 
!

      uu_residual=p%uu(:,2)-uu_average(:,nn)
!
! The other velocities are the untouched
!
      uu_advec(:,1)=p%uu(:,1)
      uu_advec(:,2)=uu_residual
      uu_advec(:,3)=p%uu(:,3)
!  
!  For the continuity equation
!   
      if (ldensity_nolog) then
        call h_dot_grad(uu_advec,p%grho,uuadvec_grho)
      else
        call h_dot_grad(uu_advec,p%glnrho,uuadvec_glnrho)
      endif
!
!  For velocity advection
!
!  Note: It is tempting to use 
!
!     call h_dot_grad(uu_advec,p%uij,p%uu,uuadvec_gu)
!
!  instead of the lines coded below, but the line just above
!  would introduce the curvature terms with the residual
!  velocity. Yet the curvature terms do not enter the 
!  non-coordinate basis advection that fargo performs. 
!  These terms have to be added manually. If one uses
!  h_dot_grad_vec, then the curvature with residual would
!  have to be removed manually and then the full speed 
!  added again (for the r-component of uuadvec_gu), like 
!  this:
!
!   uuadvec_gu(:,1)=uuadvec_gu(:,1)-&
!      rcyl_mn1*((p%uu(:,2)-uu_advec(:,2))*p%uu(:,2))
!   uuadvec_gu(:,2)=uuadvec_gu(:,2)+&
!      rcyl_mn1*((p%uu(:,1)-uu_advec(:,1))*p%uu(:,2))
!
!
!   Although working (and more line-economically), the 
!   piece of code below is more readable in my opinion. 
!
      do j=1,3
        call h_dot_grad(uu_advec,p%uij(:,j,:),tmp)
        tmp2(:,j)=tmp
      enddo
      tmp2(:,1)=tmp2(:,1)-rcyl_mn1*p%uu(:,2)*p%uu(:,2)
      tmp2(:,2)=tmp2(:,2)+rcyl_mn1*p%uu(:,1)*p%uu(:,2)
!
      uuadvec_gu=tmp2
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
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df     
      type (pencil_case) :: p
      integer :: nn
      real, dimension (nx) :: nshift
!
      intent(in) :: f,p
      intent(inout) :: df
!
!  Estimate the shift by nshift=uavg*dt/dy
!  This is just an approximation, since uavg
!  changes from one subtimestep to another, and
!  the correct cumulative shift is 
!
!    nshift=0.
!    do itsub=1,3
!       nshift=nshift+uavg*dt_sub/dy
!    enddo
!
!  But it also works fairly well this way, since uavg
!  does not change all that much between subtimesteps.  
!
      if (ldiagnos) then 
        if (idiag_nshift/=0) then
          nn=n-nghost
          nshift=uu_average(:,nn)*dt*dy_1(m) 
          call max_mn_name(nshift,idiag_nshift)
        endif
      endif
!
      if (NO_WARN) print*, f, df, p
!
    endsubroutine dspecial_dt
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
!   06-oct-03/tony: coded
!
      integer :: iname
      logical :: lreset,lwr
      logical, optional :: lwrite
!
!  Write information to index.pro
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
      if (lreset) then 
        idiag_nshift=0
        idiag_dtuf=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'nshift',idiag_nshift)
        call parse_name(iname,cname(iname),cform(iname),'dtuf',idiag_dtuf)
      enddo
!
      if (lwr) then
        write(3,*) 'i_nshift=',idiag_nshift
        write(3,*) 'i_nshift=',idiag_dtuf
      endif
!
    endsubroutine rprint_special
!***********************************************************************
    subroutine special_calc_density(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the 
!   mass equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      use Cdata
      use EquationOfState
    
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!    
!  Modified continuity equation
!
      if (ldensity_nolog) then
        df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) - &
            uuadvec_grho   - p%rho*p%divu
      else
        df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) - &
            uuadvec_glnrho - p%divu
      endif
!
    endsubroutine special_calc_density
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
      use Sub
!      
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
!  Modified momentum equation
!
      df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)-uuadvec_gu
!
!  The lines below are not symmetric. This is on purpose, to better 
!  highlight that fargo advects the azimuthal coordinate ONLY!!
!
      if (lfirst.and.ldt) then 
        advec_uu=abs(p%uu(:,1))  *dx_1(l1:l2)+ &
                 abs(uu_residual)*dy_1(  m  )*rcyl_mn1+ &
                 abs(p%uu(:,3))  *dz_1(  n  )
      endif
!
      if (ldiagnos) then 
        if (idiag_dtuf/=0) &
            call max_mn_name(advec_uu/cdt,idiag_dtuf,l_dt=.true.)
      endif
!
    endsubroutine special_calc_hydro
!***********************************************************************
    subroutine special_calc_magnetic(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the 
!   induction equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      use Cdata
!      
      real, dimension (mx,my,mz,mvar+maux), intent(inout) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
    endsubroutine special_calc_magnetic
!***********************************************************************
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
      use Cdata
      
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      if (NO_WARN) print*,df,p
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
      use Cdata
!      
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      type (boundary_condition) :: bc
!
      if (NO_WARN) print*,f(1,1,1,1),bc%bcname
!
    endsubroutine special_boundconds
!***********************************************************************
    subroutine special_before_boundary(f)
!
!  Possibility to modify the f array before the boundaries are 
!  communicated.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array
!
!  06-jul-06/tony: coded
!
      use Cdata
      use Mpicomm, only: mpiallreduce_sum
!      
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (nx,nz) :: fsum_tmp
      real, dimension (nx) :: uphi
      real :: fac
      integer :: n,m,nnghost
!
!  Pre-calculate the average large scale speed of the flow
!
      fac=1.0/nygrid
      fsum_tmp=0.
!
      do m=m1,m2
      do n=n1,n2
        nnghost=n-n1+1
        uphi=f(l1:l2,m,n,iuy)
        fsum_tmp(:,nnghost)=fsum_tmp(:,nnghost)+fac*uphi
      enddo
      enddo
!
! The sum has to be done processor-wise
! Sum over processors of same ipz, and different ipy
! --only relevant for 3D, but is here for generality
!
      call mpiallreduce_sum(fsum_tmp,uu_average,&
          (/nx,nz/),LSUMY=.true.)
!
    endsubroutine special_before_boundary
!***********************************************************************
    subroutine special_after_timestep(f,df,dt_sub)
!
!  Possibility to modify the f array after the evolution equations 
!  are solved.
!
!  In this case, add the fargo shift to the f and df-array, in 
!  fourier space. 
!
!  06-jul-06/tony: coded
!
      use Sub
      use Fourier, only: fourier_shift_y
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,ny,nz) :: tmp
      real :: dt_sub
      real, dimension (nx) :: uua

      integer :: ivar
!
!  Not yet implemented for nzgrid/=1, so this 
!  hardcoded uua=uu_average(:,1) will do no harm.
!  In the (near) future, it has to be generalized so that 
!  fourier_shift_y accepts (nx,nz) arrays. 
!
      uua=uu_average(:,1)
!
      do ivar=1,mvar
        tmp=f(l1:l2,m1:m2,n1:n2,ivar)
        call fourier_shift_y(tmp,uua*dt_sub)
        f(l1:l2,m1:m2,n1:n2,ivar)=tmp
!
!  Also shift df, unless we are at the last subtimestep
!
        if (.not.llast) then
          tmp=df(l1:l2,m1:m2,n1:n2,ivar)
          call fourier_shift_y(tmp,uua*dt_sub)
          df(l1:l2,m1:m2,n1:n2,ivar)=tmp
        endif
      enddo
!
    endsubroutine special_after_timestep
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

