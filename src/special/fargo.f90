! $Id$
!
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
! Global arrays
  real, dimension (nx,nz) :: uu_average
  real, dimension (nx,nz) :: phidot_average
! "pencils"
  real, dimension (nx,3) :: uuadvec_guu,uuadvec_gaa
  real, dimension (nx) :: uuadvec_grho,uuadvec_glnrho,uuadvec_gss
  real, dimension (nx) :: uu_residual
!
  real, dimension (nxgrid) :: xgrid1
  real :: nygrid1,dummy
!
  namelist /special_run_pars/ dummy
!
  contains
!
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
    subroutine initialize_special(f)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  06-oct-03/tony: coded
!
      use Mpicomm, only: stop_it
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!
!  Make it possible to switch the algorithm off while still
!  having this file compiled, for debug purposes.
!
      if (lrun .and. .not. lfargo_advection) then
        if (lroot) then
          print*,""
          print*,"Switch"
          print*," lfargo_advection=T"
          print*,"in init_pars of start.in if you"
          print*,"want to use the fargo algorithm"
          print*,""
        endif
        call warning("initialize_special","")
      endif
!
      xgrid1=1./xgrid
      nygrid1=1./nygrid
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine pencil_criteria_special()
!
!  All pencils that this special module depends on are specified here.
!
!  18-07-06/tony: coded
!
      if (lfargo_advection) lpenc_requested(i_uu)=.true.
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
      if (lfargo_advection) then
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
        call keep_compiler_quiet(f)
!
      endif
!
    endsubroutine calc_pencils_special
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
!  Just needs to calculate at the first sub-timestep
!
      if (lfargo_advection.and.lfirst) then
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
