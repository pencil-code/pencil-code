! $Id$
!
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
!
!---------------------------------------------------------------
!
! HOW TO USE THIS FILE
!
!---------------------------------------------------------------
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
  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include '../special.h'
! Global arrays
  real, dimension (nx,nz) :: uu_average
  real, dimension (nx,nz) :: phidot_average
! "pencils" 
  real, dimension (nx,3) :: uuadvec_guu,uuadvec_gaa
  real, dimension (nx) :: uuadvec_grho,uuadvec_glnrho,uuadvec_gss
  real, dimension (nx) :: uu_residual
  real :: dummy
  logical :: lno_radial_advection=.false.
  logical :: lfargoadvection_as_shift=.false.
!
  real, dimension (nxgrid) :: xgrid1
  real :: nygrid1
!
  namelist /special_init_pars/ dummy
!   
  namelist /special_run_pars/ lno_radial_advection, lfargoadvection_as_shift
!
  integer :: idiag_nshift=0
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
      use Mpicomm, only: stop_it
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      logical :: lstarting
!   
      if (.not.lfargo_advection) then
        if (lroot) then
          print*,""
          print*,"Switch"
          print*," lfargo_advection=T"
          print*,"in init_pars of start.in if you"
          print*,"want to use the fargo algorithm"
          print*,""
        endif
        call fatal_error("initialize_special","")
      endif
!
!  Not implemented for other than cylindrical coordinates
!
      if (coord_system/='cylindric') then
        if (lroot) then
          print*,""
          print*,"Fargo advection is only implemented for"
          print*,"cylindrical coordinates. Switch"
          print*," coord_system='cylindric'"
          print*,"in init_pars of start.in if you"
          print*,"want to use the fargo algorithm"
          print*,""
        endif
        call fatal_error("initialize_special","")
      endif
!
!  Not implemented for the energy equation either
!
      if (pretend_lnTT.or.ltemperature) call stop_it("fargo advection not "//&
          "implemented for the temperature equation")
!
!  Stuff that is only calculated once
!
      xgrid1=1./xgrid
      nygrid1=1./nygrid
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
      use Mpicomm
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      intent(inout) :: f
!
      call keep_compiler_quiet(f)
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
!  For the induction equation
!
      if (lmagnetic) then 
        lpenc_requested(i_aa)=.true.
        lpenc_requested(i_aij)=.true.
      endif
!
!  For the entropy equation
!
      if (lentropy) lpenc_requested(i_gss)=.true.
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
      integer :: j,nnghost
!
      nnghost=n-nghost
!
! Advect by the relative velocity 
!
      uu_residual=p%uu(:,2)-uu_average(:,nnghost)
!
! Advect by the original radial and vertical, but residual azimuthal
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
!     call h_dot_grad(uu_advec,p%uij,p%uu,uuadvec_guu)
!
!  instead of the lines coded below, but the line just above
!  would introduce the curvature terms with the residual
!  velocity. Yet the curvature terms do not enter the 
!  non-coordinate basis advection that fargo performs. 
!  These terms have to be added manually. If one uses
!  h_dot_grad_vec, then the curvature with residual would
!  have to be removed manually and then the full speed 
!  added again (for the r-component of uuadvec_guu), like 
!  this:
!
!   uuadvec_guu(:,1)=uuadvec_guu(:,1)-&
!      rcyl_mn1*((p%uu(:,2)-uu_advec(:,2))*p%uu(:,2))
!   uuadvec_guu(:,2)=uuadvec_guu(:,2)+&
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
      uuadvec_guu=tmp2
!
!  Advection of the magnetic potential
!
      if (lmagnetic) then
        do j=1,3
          call h_dot_grad(uu_advec,p%aij(:,j,:),tmp)
          tmp2(:,j)=tmp
        enddo
        tmp2(:,1)=tmp2(:,1)-rcyl_mn1*p%aa(:,2)*p%uu(:,2)
        tmp2(:,2)=tmp2(:,2)+rcyl_mn1*p%aa(:,1)*p%uu(:,2)
!
        uuadvec_gaa=tmp2
      endif
!
!  Advection of entropy
!
      if (lentropy) &
           call h_dot_grad(uu_advec,p%gss,uuadvec_gss)
!      
      call keep_compiler_quiet(f)
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine dspecial_dt(f,df,p)
!
!  Calculate right hand side of ONE OR MORE extra coupled PDEs
!  along the 'current' Pencil, i.e. f(l1:l2,m,n) where
!  m,n are global variables looped over in equ.f90
!
!  Due to the multi-step Runge Kutta timestepping used one MUST always
!  add to the present contents of the df array.  NEVER reset it to zero.
!
!  Several precalculated Pencils of information are passed if for
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
      integer :: nnghost
      real, dimension (nx) :: nshift,phidot
!
      intent(in) :: f,p
      intent(inout) :: df
!
!  Estimate the shift by nshift=phidot*dt/dy
!  This is just an approximation, since uavg
!  changes from one subtimestep to another, and
!  the correct cumulative shift is 
!
!    nshift=0.
!    do itsub=1,3
!       nshift=nshift+phidot*dt_sub/dy
!    enddo
!
!  But it also works fairly well this way, since uavg
!  does not change all that much between subtimesteps.  
!
      if (ldiagnos) then 
        if (idiag_nshift/=0) then
          nnghost=n-nghost
          phidot=uu_average(:,nnghost)*rcyl_mn1
          nshift=phidot*dt*dy_1(m) 
          call max_mn_name(nshift,idiag_nshift)
        endif
      endif
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
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
!  Reads and registers print parameters relevant to special
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
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'nshift',idiag_nshift)
      enddo
!
      if (lwr) then
        write(3,*) 'i_nshift=',idiag_nshift
      endif
!
    endsubroutine rprint_special
!***********************************************************************
    subroutine special_calc_density(f,df,p)
!
!   Calculate a additional 'special' term on the right hand side of the 
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
      call keep_compiler_quiet(f)
!
    endsubroutine special_calc_density
!***********************************************************************
    subroutine special_calc_hydro(f,df,p)
!
!   Calculate a additional 'special' term on the right hand side of the 
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
!
!  Modified momentum equation
!
      df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)-uuadvec_guu
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
      call keep_compiler_quiet(f)
!
    endsubroutine special_calc_hydro
!***********************************************************************
    subroutine special_calc_magnetic(f,df,p)
!
!   Calculate a additional 'special' term on the right hand side of the 
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
      df(l1:l2,m,n,iax:iaz)=df(l1:l2,m,n,iax:iaz)-uuadvec_gaa
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_magnetic
!***********************************************************************
    subroutine special_calc_entropy(f,df,p)
!
!   Calculate a additional 'special' term on the right hand side of the 
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
      df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss)-uuadvec_gss
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_entropy
!***********************************************************************
    subroutine special_boundconds(f,bc)
!
!   Possibility of custom boundary condition.
!
!   06-oct-03/tony: coded
!
      use Cdata
!      
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      type (boundary_condition) :: bc
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(bc%bcname)
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
      integer :: nnghost
!
!  Just needs to the calculated at the first sub-timestep
!
      if (lfirst) then
!
!  Pre-calculate the average large scale speed of the flow
!
        fsum_tmp=0.
!
        do n=n1,n2;do m=m1,m2
          nnghost=n-nghost
          uphi=f(l1:l2,m,n,iuy)
          fsum_tmp(:,nnghost)=fsum_tmp(:,nnghost)+uphi*nygrid1
        enddo;enddo
!
! The sum has to be done processor-wise
! Sum over processors of same ipz, and different ipy
! --only relevant for 3D, but is here for generality
!
        call mpiallreduce_sum(fsum_tmp,uu_average,&
             (/nx,nz/),idir=2) !idir=2 is equal to old LSUMY=.true.
      endif
!
    endsubroutine special_before_boundary
!***********************************************************************
    subroutine special_after_timestep(f,df,dt_sub)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real :: dt_sub
!
      if (llast) then 
        if (lfargoadvection_as_shift) then
          call fourier_shift_fargo(f)
        else
          call advect_fargo(f)
        endif
      endif
!
!  Just for test purposes and comparison with the loop advection 
!  in Stone, J. et al., JCP 250, 509 (2005)
!
      if (lno_radial_advection) then 
        f(:,:,:,iux) = 0.
        df(:,:,:,iux) = 0.
      endif
!
      call keep_compiler_quiet(dt_sub)
!
    endsubroutine special_after_timestep
!***********************************************************************
    subroutine fourier_shift_fargo(f)
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
      use Fourier, only: fft_y_parallel
      use Cdata
      use Mpicomm
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,ny) :: a_re,a_im
      real, dimension (nx) :: phidot
      integer :: ivar,nnghost
!
!  Pencil uses linear velocity. Fargo will shift based on 
!  angular velocity. Get phidot from uphi. 
!
      do n=n1,n2
        nnghost=n-n1+1
        phidot=uu_average(:,nnghost)*rcyl_mn1
!
        do ivar=1,mvar
!
          a_re=f(l1:l2,m1:m2,n,ivar); a_im=0.
!
!  Forward transform. No need for computing the imaginary part. 
!  The transform is just a shift in y, so no need to compute 
!  the x-transform either. 
!
          call fft_y_parallel(a_re,a_im,SHIFT_Y=phidot*dt,lneed_im=.false.)
!
!  Inverse transform of the shifted array back into real space. 
!  No need again for either imaginary part of x-transform. 
!
          call fft_y_parallel(a_re,a_im,linv=.true.)
          f(l1:l2,m1:m2,n,ivar)=a_re
!
        enddo
      enddo
!
    endsubroutine fourier_shift_fargo
!********************************************************************
    subroutine advect_fargo(f)
!
!  Possibility to modify the f array after the evolution equations 
!  are solved.
!
!  In this case, do the fargo shift to the f and df-array, in 
!  real space. 
!
!  06-jul-06/tony: coded
!
      use Sub
      use Cdata
      use Mpicomm
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      real, dimension (nx,nygrid,mvar) :: faux_remap,faux_remap_shift
!
      real, dimension (nx,ny) :: faux,tmp2
      real, dimension (nx,nygrid) :: tmp
!
      integer :: ivar,ng,ig,mshift,cellshift,i,mserial
!
      integer, dimension (nx,nz) :: shift_intg
      real, dimension (nx,nz) :: shift_total,shift_frac
!
! For shift in real space, the shift is done in integer number of 
! cells in the azimuthal direction, so, take only the integer part 
! of the velocity for fargo advection. 
!
      do n=1,nz
        phidot_average(:,n) = uu_average(:,n)*rcyl_mn1
      enddo
!
! Define the integer circular shift
!
      shift_total = phidot_average*dt*dy_1(mpoint)
      shift_intg  = nint(shift_total)
!
! Do circular shift of cells
!
      do n=n1,n2
!
        do ivar=1,mvar
          faux=f(l1:l2,m1:m2,n,ivar)
          call remap_to_pencil_y(faux,tmp)
          faux_remap(:,:,ivar)=tmp
        enddo
!
        ng=n-n1+1
!
        do i=l1,l2
          ig=i-l1+1
          cellshift=shift_intg(ig,ng)
!
          do m=1,ny
            mserial=m+ipy*ny
            mshift=mserial-cellshift
            if (mshift .lt. 1 )     mshift = mshift + nygrid
            if (mshift .gt. nygrid) mshift = mshift - nygrid
!          
            do ivar=1,mvar
              faux_remap_shift(ig,mserial,ivar) = faux_remap(ig,mshift,ivar)
            enddo
          enddo
        enddo
!
        do ivar=1,mvar
          tmp=faux_remap_shift(:,:,ivar)
          call unmap_from_pencil_y(tmp, tmp2)
          f(l1:l2,m1:m2,n,ivar)=tmp2
        enddo
      enddo
!
! Fractional step
!
      shift_frac  = shift_total-shift_intg
      call fractional_shift(f,shift_frac)
!
    endsubroutine advect_fargo
!********************************************************************
    subroutine fractional_shift(f,shift_frac)
!
      use Deriv, only:der
      use Mpicomm
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df_frac
      real, dimension (nx,nz) :: shift_frac,uu_frac
      real, dimension(3) :: dt_sub
      real, dimension(nx) :: ushift,dfdy,facx
      integer :: j,itsubstep
!
! Set boundaries      
!
      if (nprocy==1) then 
        f(:,1:m1-1,:,:)  = f(:,m2i:m2,:,:)
        f(:,m2+1:my,:,:) = f(:,m1:m1i,:,:)
      else
        call initiate_isendrcv_bdry(f)
        call finalize_isendrcv_bdry(f)
      endif
!
      dt_sub=dt*beta_ts
!
      facx=1./(dt*dy_1(mpoint)*rcyl_mn1)
!
      do n=1,nz
        uu_frac(:,n)=shift_frac(:,n)*facx
      enddo
!
      do itsubstep=1,itorder
        if (itsubstep==1) then
          df_frac=0.
        else
          df_frac=alpha_ts(itsubstep)*df_frac
        endif
!
        do n=n1,n2;do m=m1,m2
!
          ushift=uu_frac(:,n-n1+1)
!
          do j=1,mvar
            call der(f,j,dfdy,2)
            df_frac(l1:l2,m,n,j)=df_frac(l1:l2,m,n,j)-ushift*dfdy
            f(l1:l2,m,n,j) = &
                 f(l1:l2,m,n,j)+ dt_sub(itsubstep)*df_frac(l1:l2,m,n,j)
          enddo
!
        enddo;enddo
!
      enddo
!
    endsubroutine fractional_shift
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

