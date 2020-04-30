! $Id$
!
!  This module tries to solve for gradient matrix of particle velocities
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MPVAR CONTRIBUTION 9
! MPAUX CONTRIBUTION 9
! CPARAM logical, parameter :: lparticles_grad=.true.
!
!***************************************************************
module Particles_grad
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
  use Particles_sub
!
  implicit none
!
  include 'particles_grad.h'
!
  !namelist /particles_grad_init_pars/ 
!
  real :: fluid_mu,rhodust

  namelist /particles_grad_run_pars/ &
          fluid_mu,rhodust
!
  integer :: idiag_sigmap11max=0,idiag_sigmap12max=0,idiag_sigmap13max=0, &
             idiag_sigmap21max=0,idiag_sigmap22max=0,idiag_sigmap23max=0, &
             idiag_sigmap31max=0,idiag_sigmap32max=0,idiag_sigmap33max=0
!
  contains
!***********************************************************************
    subroutine register_particles_grad
!
!  Set up indices for access to the fp and dfp arrays.
!
!  17-sep-15/dhruba: coded
!
      use FArrayManager, only: farray_register_auxiliary
!
      if (lroot) call svn_id( &
          "$Id$")
!
!  Indices for particle gradients.
!
      call append_npvar('isigmap11',isigmap11)
      call append_npvar('isigmap12',isigmap12)
      call append_npvar('isigmap13',isigmap13)
      call append_npvar('isigmap21',isigmap21)
      call append_npvar('isigmap22',isigmap22)
      call append_npvar('isigmap23',isigmap23)
      call append_npvar('isigmap31',isigmap31)
      call append_npvar('isigmap32',isigmap32)
      call append_npvar('isigmap33',isigmap33)
!
! This module always demands gradu stored as an auxiliary array
!
      call farray_register_auxiliary('guij',iguij,vector=9)
      igradu11=iguij+0; igradu12=iguij+1; igradu13=iguij+2
      igradu21=iguij+3; igradu22=iguij+4; igradu23=iguij+5
      igradu31=iguij+6; igradu32=iguij+7; igradu33=iguij+8
!
    endsubroutine register_particles_grad
!***********************************************************************
    subroutine initialize_particles_grad(f)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  15-sep-15/dhruba: coded
!
      use SharedVariables, only: put_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f

      call keep_compiler_quiet(f)
!
    endsubroutine initialize_particles_grad
!***********************************************************************
    subroutine pencil_criteria_par_grad
!
!  All pencils that the Particles_grad module depends on are specified here.
!
!  17-sep-15/dhruba: coded
!
!
    endsubroutine pencil_criteria_par_grad
!***********************************************************************
    subroutine set_particle_grad(f,fp,npar_low,npar_high,ineargrid,init)
!
!  Set radius of new particles.
!
!  18-sep-09/nils: adapted from init_particles_radius
!
      use General, only: random_number_wrapper
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
      integer :: npar_low,npar_high
      integer, dimension(mpar_loc,3) :: ineargrid
      logical, optional :: init
!
      integer :: i, j, k, ij, kend, ind, ibin
      logical :: initial
      integer :: imn
      real, dimension(nx,3,3) :: gradu
      real, dimension(3,3) :: gradup
      real, dimension(9) :: gradup_lin
!
      initial=.false.
      if (present(init)) then
        if (init) initial=.true.
      endif
!
! The dust particles are initialized with the velocities of the flow. Hence
! the initial value of the gradient of dust particle velocities is 
! same as the interpolated value of gradient of gas velocity. 
! At this stage the gradient of gas velocity has not been calculated, although
! the gas velocity is known. So we first calculate the gradient of gas velocity
!
      do imn=1,ny*nz
        n=nn(imn)
        m=mm(imn)
        call gij(f,iuu,gradu,1)
        ij=iguij-1
        do i=1,3; do j=1,3
            ij=ij+1
            f(l1:l2,m,n,ij) = gradu(:,i,j) 
        enddo;enddo
      enddo
!
! Next impose periodic boundary condition on the gradu 
!
      call set_periodic_boundcond_on_aux(f,igradu11)
      call set_periodic_boundcond_on_aux(f,igradu12)
      call set_periodic_boundcond_on_aux(f,igradu13)
      call set_periodic_boundcond_on_aux(f,igradu21)
      call set_periodic_boundcond_on_aux(f,igradu22)
      call set_periodic_boundcond_on_aux(f,igradu23)
      call set_periodic_boundcond_on_aux(f,igradu31)
      call set_periodic_boundcond_on_aux(f,igradu32)
      call set_periodic_boundcond_on_aux(f,igradu33)
!
! Now interpolate gradu to the location of the particles. 
!
      do k=1,npar_loc
        call interpolate_linear(f,igradu11,igradu33, &
          fp(k,ixp:izp),gradup_lin,ineargrid(k,:),0,0)
        fp(k,isigmap11:isigmap33) =  gradup_lin
      enddo
!
    endsubroutine set_particle_grad
!***********************************************************************
    subroutine dsigmap_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of the gradient of particle velocities.
!
!  17-sep-15/dhruba: coded
!
      use Particles_number
      use Sub, only: matrix2linarray, linarray2matrix
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      type (pencil_case) :: p
      integer, dimension (mpar_loc,3) :: ineargrid
      logical :: lfirstcall=.true., lheader
      integer :: k,i,k1,k2
      real :: radius,one_over_tau
      real,dimension(3,3) :: sigmap,gradup,dsigma,dsigmap
      real, dimension(9) :: gradup_lin
      real,dimension(9) :: dsigmaplin
!
      intent (in) :: f
      intent (out) :: dfp
      intent (inout) :: fp
!
      real :: one_by_tau
!
!  Print out header information in first time step.
!
      lheader=lfirstcall.and.lroot.and.(.not.lpencil_check_at_work)
!
!  Identify module.
!
      if (lheader) print*,'dsigmap_dt_pencil: Calculate dsigmap/dt'
      lfirstcall=.false.
!
      if (npar_imn(imn) /= 0) then

        k1 = k1_imn(imn)
        k2 = k2_imn(imn)
!
        do k = k1,k2
          call interpolate_linear(f,igradu11,igradu33,&
               fp(k,ixp:izp),gradup_lin,ineargrid(k,:),0,ipar(k))
          call linarray2matrix(gradup_lin,gradup)
          radius=fp(k,iap)
          one_by_tau=(9./2.)*fluid_mu/(radius*radius*rhodust)
          call linarray2matrix(fp(k,isigmap11:isigmap33),sigmap)
          dsigma=one_by_tau*(gradup-sigmap)-matmul(sigmap,sigmap)
          call matrix2linarray(dsigmap, dsigmaplin)
          dfp(k,isigmap11:isigmap33) =  dfp(k,isigmap11:isigmap33) + dsigmaplin 
        enddo
      endif
!
    endsubroutine dsigmap_dt_pencil
!***********************************************************************
    subroutine dsigmap_dt(f,df,fp,dfp,ineargrid)
!
!  Evolution of particle gradient (actually merely the calculation of diagnostic 
! variables on particles because the actual evolution is calculated in a pencilized manner )
!
!  17-sep-15/dhruba: coded
!
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
!  Diagnostic output.
!
      if (ldiagnos) then
        if (idiag_sigmap11max/=0) &
            call max_par_name(fp(1:npar_loc,isigmap11),idiag_sigmap11max)
        if (idiag_sigmap12max/=0) &
            call max_par_name(fp(1:npar_loc,isigmap12),idiag_sigmap12max)
        if (idiag_sigmap13max/=0) &
            call max_par_name(fp(1:npar_loc,isigmap13),idiag_sigmap13max)
        if (idiag_sigmap21max/=0) &
            call max_par_name(fp(1:npar_loc,isigmap21),idiag_sigmap21max)
        if (idiag_sigmap22max/=0) &
            call max_par_name(fp(1:npar_loc,isigmap22),idiag_sigmap22max)
        if (idiag_sigmap23max/=0) &
            call max_par_name(fp(1:npar_loc,isigmap23),idiag_sigmap23max)
        if (idiag_sigmap31max/=0) &
            call max_par_name(fp(1:npar_loc,isigmap31),idiag_sigmap31max)
        if (idiag_sigmap32max/=0) &
            call max_par_name(fp(1:npar_loc,isigmap32),idiag_sigmap32max)
        if (idiag_sigmap33max/=0) &
            call max_par_name(fp(1:npar_loc,isigmap33),idiag_sigmap33max)
      endif
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dsigmap_dt
!***********************************************************************
    subroutine read_particles_grad_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
      integer :: pos
!
      !read(parallel_unit, NML=particles_grad_init_pars, IOSTAT=iostat)
!
    endsubroutine read_particles_grad_init_pars
!***********************************************************************
    subroutine write_particles_grad_init_pars(unit)
!
      integer, intent(in) :: unit
!
      !write(unit, NML=particles_grad_init_pars)
!
    endsubroutine write_particles_grad_init_pars
!***********************************************************************
    subroutine read_particles_grad_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=particles_grad_run_pars, IOSTAT=iostat)
!
    endsubroutine read_particles_grad_run_pars
!***********************************************************************
    subroutine write_particles_grad_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=particles_grad_run_pars)
!
    endsubroutine write_particles_grad_run_pars
!***********************************************************************
    subroutine rprint_particles_grad(lreset,lwrite)
!
!  Read and register print parameters relevant for particles grad.
!
!  22-aug-05/anders: coded
!
      use Diagnostics, only: parse_name
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname
!
!  Reset everything in case of reset.
!
      if (lreset) then
        idiag_sigmap11max=0; idiag_sigmap12max=0; idiag_sigmap13max=0
        idiag_sigmap21max=0; idiag_sigmap22max=0; idiag_sigmap23max=0
        idiag_sigmap31max=0; idiag_sigmap32max=0; idiag_sigmap33max=0
      endif
!
!  Run through all possible names that may be listed in print.in.
!
      if (lroot.and.ip<14) &
          print*, 'rprint_particles_grad: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'sigmap11max',idiag_sigmap11max)
        call parse_name(iname,cname(iname),cform(iname),'sigmap12max',idiag_sigmap12max)
        call parse_name(iname,cname(iname),cform(iname),'sigmap13max',idiag_sigmap13max)
        call parse_name(iname,cname(iname),cform(iname),'sigmap21max',idiag_sigmap21max)
        call parse_name(iname,cname(iname),cform(iname),'sigmap22max',idiag_sigmap22max)
        call parse_name(iname,cname(iname),cform(iname),'sigmap23max',idiag_sigmap23max)
        call parse_name(iname,cname(iname),cform(iname),'sigmap31max',idiag_sigmap31max)
        call parse_name(iname,cname(iname),cform(iname),'sigmap32max',idiag_sigmap32max)
        call parse_name(iname,cname(iname),cform(iname),'sigmap33max',idiag_sigmap33max)
      enddo
!
    endsubroutine rprint_particles_grad
!***********************************************************************
endmodule Particles_grad
