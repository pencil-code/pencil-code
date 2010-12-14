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
  real :: nu_hyper3=1e-8
!
  namelist /special_init_pars/ dummy
!   
  namelist /special_run_pars/ lno_radial_advection, lfargoadvection_as_shift, nu_hyper3
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
!  Not implemented for the energy equation either
!
      if (pretend_lnTT.or.ltemperature) call stop_it("fargo advection not "//&
          "implemented for the temperature equation")
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
      real :: fac
      integer :: nnghost
!
!  Pre-calculate the average large scale speed of the flow
!
      fac=1.0/nygrid
      fsum_tmp=0.
!
      do m=m1,m2
      do n=n1,n2
        nnghost=n-nghost
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
          (/nx,nz/),idir=2) !idir=2 is equal to old LSUMY=.true.
!
    endsubroutine special_before_boundary
!***********************************************************************
    subroutine special_after_timestep(f,df,dt_sub)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real :: dt_sub
!
      if (lfargoadvection_as_shift) then
        call fourier_shift_fargo(f,df,dt_sub)
      else
        call advect_fargo(f,df,dt_sub)
      endif

    endsubroutine special_after_timestep
!***********************************************************************
    subroutine fourier_shift_fargo(f,df,dt_sub)
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
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,ny) :: a_re,a_im
      real :: dt_sub
      real, dimension (nxgrid) :: phidot_serial
      integer :: ivar,nnghost
!
      integer :: ixdo,ixup,izdo,izup,iproc_recv,jx,jz,iz_serial
      real, dimension (nx,nz) :: uu_average_recv
      real, dimension (nxgrid,nzgrid) :: uu_average_allprocs
!
      uu_average_recv=1.
      if (iproc/=root) then
        if (ipy==0) call mpisend_real(uu_average,(/nx,nz/),root,222)
      else
        do jx=0,nprocx-1 
          do jz=0,nprocz-1
            !iproc=ipx+nprocx*ipy+nprocx*nprocy*ipz
            iproc_recv=jx+nprocx*nprocy*jz
            if (iproc_recv/=root) then
              call mpirecv_real(uu_average_recv,(/nx,nz/),iproc_recv,222)
            else
              uu_average_recv=uu_average
            endif
            ixdo= jx   *nx + 1 ; izdo= jz   *nz + 1 
            ixup=(jx+1)*nx     ; izup=(jz+1)*nz
            uu_average_allprocs(ixdo:ixup,izdo:izup)=uu_average_recv
          enddo
        enddo
      endif
      call mpibcast_real(uu_average_allprocs,(/nxgrid,nzgrid/))
!
!  Pencil uses linear velocity. Fargo will shift based on 
!  angular velocity. Get phidot from uphi. 
!
      do n=n1,n2
        nnghost=n-n1+1
        iz_serial=ipz*nz + nnghost
        phidot_serial=uu_average_allprocs(:,iz_serial)/xgrid
!
        do ivar=1,mvar
!
          a_re=f(l1:l2,m1:m2,n,ivar); a_im=0.
!
!  Forward transform. No need for computing the imaginary part. 
!  The transform is just a shift in y, so no need to compute 
!  the x-transform either. 
!
          call fft_y_parallel(a_re,a_im,SHIFT_Y=phidot_serial*dt_sub,lneed_im=.false.)
!
!  Inverse transform of the shifted array back into real space. 
!  No need again for either imaginary part of x-transform. 
!
          call fft_y_parallel(a_re,a_im,linv=.true.)
          f(l1:l2,m1:m2,n,ivar)=a_re
!
!  Also shift df, unless we are at the last subtimestep
!
          if (lwrite_dvar.or..not.llast) then
            a_re=df(l1:l2,m1:m2,n,ivar); a_im=0.
            call fft_y_parallel(a_re,a_im,SHIFT_Y=phidot_serial*dt_sub,lneed_im=.false.)
            call fft_y_parallel(a_re,a_im,linv=.true.)
            df(l1:l2,m1:m2,n,ivar)=a_re
          endif
!
        enddo
      enddo
!
!  Just for test purposes and comparison with the loop advection 
!  in Stone, J. et al., JCP 250, 509 (2005)
!
      if (lno_radial_advection) then 
        f(:,:,:,iux) = 0.
        df(:,:,:,iux) = 0.
      endif
!
    endsubroutine fourier_shift_fargo
!********************************************************************
    subroutine advect_fargo(f,df,dt_sub)
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
      real, dimension (mx,my,mz,mvar) :: df
!
      real, dimension (ny,mvar) :: faux
      real, dimension (ny,mvar) :: dfaux
!
      real, dimension (ny,mvar) :: a_re,da_re
      integer :: ivar,ng,mg,ig,mshift,cellshift,i
      real :: dt_sub
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
      shift_total = phidot_average*dt_sub*dy_1(mpoint)
      shift_intg  = nint(shift_total)
      shift_frac  = shift_total-shift_intg
!
      do m=1,10
        f(l1:l2,m,n1:n2,iuz)=shift_total
      enddo
      do m=11,20
        f(l1:l2,m,n1:n2,iuz)=shift_intg
      enddo
      do m=21,30
        f(l1:l2,m,n1:n2,iuz)=shift_frac
      enddo
!
! Do circular shift of cells
!
      do n=n1,n2
        ng=n-n1+1
!
        do i=l1,l2
          ig=i-l1+1
          cellshift=shift_intg(ig,ng)
!
           faux =f(i,m1:m2,n,1:mvar)
          dfaux=df(i,m1:m2,n,1:mvar)
!
          do m=1,ny
            mshift=m-cellshift
            if (mshift .lt. 1 ) mshift = mshift + ny
            if (mshift .gt. ny) mshift = mshift - ny
!          
            do ivar=1,mvar
               a_re(m,ivar) = faux(mshift,ivar)
              !if (lwrite_dvar.or..not.llast) &
              da_re(m,ivar) = dfaux(mshift,ivar)
            enddo
          enddo
!         
          do ivar=1,mvar
            if (ivar.ne.iuz) then
              if (.not.((i.eq.l1).or.(i.eq.l2))) then
                f(i,m1:m2,n,ivar)=a_re(:,ivar)
                !if (lwrite_dvar.or..not.llast) &
                df(i,m1:m2,n,ivar)=da_re(:,ivar)
              endif
            endif
          enddo
        enddo
      enddo
!
! Fractional step
!
      !call fractional_step(f,df,shift_frac)
      call fractional_step_highorder(f,df,shift_frac,dt_sub)
!
    endsubroutine advect_fargo
!********************************************************************
    subroutine fractional_step_highorder(f,df,shift_frac,dt_full_step)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df      
      real, dimension (mx,my,mz,mvar) :: df_frac,d2f_frac      
      real, dimension (nx,nz) :: shift_frac,phidot_frac,uu_frac
      real :: dt_full_step,shift_cell,ushift
      real, dimension(3) :: dt_substep
      integer :: ivar,i,itsubstep
      real, dimension (ny) :: gradf,graddf
!
      dt_substep=dt_full_step*beta_ts
!
      !shift_total = phidot_average*dt_sub*dy_1(mpoint)
      !shift_intg  = nint(shift_total)
      !shift_frac  = shift_total-shift_intg
!
      !phidot_frac = shift_frac/(dt_full_step*dy_1(mpoint))
      !do n=1,nz
      !  uu_frac(:,n) = phidot_frac(:,n)*rcyl_mn
      !enddo
      uu_frac=shift_frac/dt_full_step

!
      do itsubstep=1,itorder
        if (itsubstep==1) then
          df_frac=0.
          d2f_frac=0.
        else
          df_frac=alpha_ts(itsubstep)*df_frac
          d2f_frac=alpha_ts(itsubstep)*d2f_frac          
        endif
!
        do n=n1,n2;do i=l1,l2
          !shift_cell=shift_frac(i-l1+1,n-n1+1)
          !shift_cell=uu_frac(i-l1+1,n-n1+1)
          ushift=uu_frac(i-l1+1,n-n1+1)
!
          do ivar=1,mvar
!
            call linterp( f(i,:,n,ivar),1.,gradf)
            call linterp(df(i,:,n,ivar),1.,graddf)
!
            df_frac(i,m1:m2,n,ivar) =  df_frac(i,m1:m2,n,ivar)  - ushift*gradf  
            d2f_frac(i,m1:m2,n,ivar) = d2f_frac(i,m1:m2,n,ivar) - ushift*graddf
!          
            if ((i.eq.l1).or.(i.eq.l2)) then
              df_frac(i,m1:m2,n,ivar)=0.
              df_frac(i,m1:m2,n,ivar)=0.
            endif

            f(i,m1:m2,n,ivar) = &
                 f(i,m1:m2,n,ivar)+ dt_substep(itsubstep)*df_frac(i,m1:m2,n,ivar)
!
            df(i,m1:m2,n,ivar) =&
                 df(i,m1:m2,n,ivar)+dt_substep(itsubstep)*d2f_frac(i,m1:m2,n,ivar)
!
          enddo
!
        enddo;enddo
!
      enddo
!
    endsubroutine fractional_step_highorder
!********************************************************************
    subroutine fractional_step(f,df,shift_frac)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
!
      real, dimension (nx,nz) :: shift_frac
      real :: shift_cell
      real, dimension (ny) :: g,dg
      integer :: i,ivar
!
! Calculate fractional shift in cells and interpolate using van leer
!
      do n=n1,n2;do i=l1,l2
        shift_cell=shift_frac(i-l1+1,n-n1+1)
!
        do ivar=1,mvar

          !call vanleer( f(i,:,n,ivar),-shift,g)
          !call vanleer(df(i,:,n,ivar),-shift,dg)

          call linterp( f(i,:,n,ivar),-shift_cell,g)
          call linterp(df(i,:,n,ivar),-shift_cell,dg)
!
          if (ivar .ne. iuz) then
!
          f(i,m1:m2,n,ivar)  =  f(i,m1:m2,n,ivar) + g
          df(i,m1:m2,n,ivar) = df(i,m1:m2,n,ivar) + dg
!          
          endif
        enddo
      enddo;enddo
!
    endsubroutine fractional_step
!********************************************************************
    subroutine linterp(f,shift_cell,g)
!
      real, dimension (my) :: f
      real, dimension (ny) :: g,der1f,der6f
      real :: shift_cell,fac
!
!      do m=m1,m2
!        g(m-m1+1) = f(m) + (f(m+1)-f(m))*shift
!      enddo
!
! Set boundaries      
!
      f(1:m1-1) = f(m2i:m2)
      f(m2+1:my) = f(m1:m1i)
!
      do m=m1,m2
        der1f(m-m1+1) = 1.0/60*       &
             (+ 45.0*(f(m+1)-f(m-1))  &
              -  9.0*(f(m+2)-f(m-2))  &
              +      (f(m+3)-f(m-3)))
!
        fac=dy_1(m)**6
        der6f(m-m1+1)=fac*(- 20.0* f(m) &
             + 15.0*(f(m+1)+f(m-1))     &
             -  6.0*(f(m+2)+f(m-2))     &
             +      (f(m+3)+f(m-3)))
      enddo
!
!      nu_hyper3=1e-10
      !shiftf = shift_cell*der1f + nu_hyper3*der6f
!      
      g = der1f !+ nu_hyper3*der6f
!
    endsubroutine linterp
!********************************************************************
    subroutine vanleer(f,shift_cell,g)
!
      real, dimension (my) :: f
      real, dimension (ny) :: g
      real :: shift_cell,gl,gm,gr,limiter,tmp,hh
!
! Set boundaries      
!
      f(1:m1-1) = f(m2i:m2)
      f(m2+1:my) = f(m1:m1i)
!
      do m=m1,m2
        gl=f(m-1); gm=f(m); gr=f(m+1)
!
        hh=gr-gl
!
        if (hh == 0) then 
          g(m-m1+1)=gm
        else
          tmp = (gm-gl)*(gr-gm)
          if (tmp > 0.) then 
            limiter=tmp
          else
            limiter=0.
          endif
          g(m-m1+1) = gm + 2*shift_cell/(gr-gl)*limiter
        endif
!
      enddo
!
    endsubroutine vanleer
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

