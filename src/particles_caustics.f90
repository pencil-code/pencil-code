! $Id$
!
!  This module writes information about the local state of the gas at
!  the positions of a selected number of particles.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_caustics=.true.
!
! MAUX CONTRIBUTION 9
! COMMUNICATED AUXILIARIES 9
! MPVAR CONTRIBUTION 6 
!
!***************************************************************
module Particles_caustics
!
  use Cdata
  use Messages
  use Particles_cdata
  use Particles_map
  use Particles_mpicomm
!
  implicit none
!
  include 'particles_caustics.h'
!
  real :: fake_eta=0., epsilondX=1e-4
!
  namelist /particles_caustics_init_pars/ &
    epsilondX
!
  namelist /particles_caustics_run_pars/ &
  fake_eta
!
  contains
!***********************************************************************
    subroutine register_particles_caustics()
!
!  Set up indices for access to the fp and dfp arrays
!
!  May-16/dhruba: coded
!
      use FArrayManager, only: farray_register_auxiliary
!
      if (lroot) call svn_id("particles_caustics")
!
!  Indices for \deltaX and \deltaV at particle positions
!
      idXp1=npvar+1
      pvarname(npvar+1)='idXp1'
      idXp2=npvar+2
      pvarname(npvar+2)='idXp2'
      idXp3=npvar+3
      pvarname(npvar+3)='idXp3'
      idVp1=npvar+4
      pvarname(npvar+4)='idVp1'
      idVp2=npvar+5
      pvarname(npvar+5)='idVp2'
      idVp3=npvar+6
      pvarname(npvar+6)='idVp3'
!
!  Increase npvar accordingly.
!
      npvar=npvar+6
!
!  Set indices for velocity gradient matrix at grid points
!
      call farray_register_auxiliary('guij',iguij,communicated=.true.,vector=9)
      igu11=iguij; igu12=iguij+1; igu13=iguij+2
      igu21=iguij+3; igu22=iguij+4; igu23=iguij+5
      igu31=iguij+6; igu32=iguij+7; igu33=iguij+8
!
!  Check that the fp and dfp arrays are big enough.
!
      if (npvar > mpvar) then
        if (lroot) write(0,*) 'npvar = ', npvar, ', mpvar = ', mpvar
        call fatal_error('register_particles','npvar > mpvar')
      endif
!
    endsubroutine register_particles_caustics
!***********************************************************************
    subroutine initialize_particles_caustics(f)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!
      use General, only: keep_compiler_quiet
      use FArrayManager
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!  Stop if there is no velocity to calculate derivatives
!
      if (.not. (lhydro.or.lhydro_kinematic)) &
        call fatal_error('initialize_particles_caustics','you must select either hydro or hydro_kinematic')
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_particles_caustics
!***********************************************************************
    subroutine init_particles_caustics(f,fp,ineargrid)
!
      use Sub, only: kronecker_delta,linarray2matrix,matrix2linarray
      use General, only: keep_compiler_quiet,random_number_wrapper
      use Mpicomm, only:  mpiallreduce_sum_int
      use Hydro, only: calc_gradu
      real, dimension (mx,my,mz,mfarray), intent (in) :: f
      real, dimension (mpar_loc,mparray), intent (out) :: fp
      integer, dimension (mpar_loc,3), intent (in) :: ineargrid
      integer, dimension(3) :: iXpdX
      real, dimension(3) :: uup1, uup2, Xp,dXp, XpdX
      real, dimension(nx,3:3) :: uij 
      real, dimension(9) :: Sij_lin
      real, dimension(3,3) :: Sijp
      real :: rno01 
      integer :: ipzero,ik,ii,jj,ij
!
!
     if (lroot) print*, 'init_particles_caustics: setting init. cond.'
!
! We set the gradient of the particle velocity field as that of the fluid velocity
! field. To do that we need to calculate the calculate the gradient of velocity 
! on a grid. This is stored as an auxiliary array within f (:,:,:,igu11:igu33) 
!
      call calc_gradu(f) 
      do ip=1,npar_loc
!
! Initialize the deltaX by a small random vector. 
! It is a factor epsilondX of grid spacing (we ignore non-uniform grids)
! Remember; dx_1 is 1/dx available from cdata.f90
!
        Xp = fp(ip,ixp:izp)
        do ii=1,3 
          call random_number_wrapper(rno01)
          fp(ip,idXp1+ii-1)=epsilondX*(2*rno01-1)*dx
        enddo
        dXp = fp(ip,idXp1:idXp3)
!
!  interpolate the gradu matrix to particle positions
!
        call interpolate_linear(f,igu11,igu33,fp(ip,ixp:izp),Sij_lin,ineargrid(ip,:), &
              0,ipar(ip))
        call linarray2matrix(Sij_lin,Sijp)
        do ii=1,3
          fp(ip,idVp1+ii-1) =  dot_product(Sijp(ii,:),dXp) 
        enddo
! loop over ip ends below
     enddo
!
    endsubroutine init_particles_caustics
!***********************************************************************
    subroutine dcaustics_dt(f,df,fp,dfp,ineargrid)
!
!
      use Diagnostics
      use Particles_sub, only: sum_par_name
!
      real, dimension (mx,my,mz,mfarray), intent (in) :: f
      real, dimension (mx,my,mz,mvar), intent (inout) :: df
      real, dimension (mpar_loc,mparray), intent (in) :: fp
      real, dimension (mpar_loc,mpvar), intent (inout) :: dfp
      integer, dimension (mpar_loc,3), intent (in) :: ineargrid
      logical :: lheader, lfirstcall=.true.
!
!  Print out header information in first time step.
!
      lheader=lfirstcall .and. lroot
      if (lheader) then
        print*,'dcaustics_dt: Calculate dcaustics_dt'
      endif
!
      if (ldiagnos) then
!
! No diagnostic is coded yet, but would be.
!
      endif
!
      if (lfirstcall) lfirstcall=.false.
!
    endsubroutine dcaustics_dt
!***********************************************************************
    subroutine dcaustics_dt_pencil(f,df,fp,dfp,p,ineargrid,k,taup1)
!
      use Sub, only: linarray2matrix,matrix2linarray
!
      real, dimension (mx,my,mz,mfarray),intent(in) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (mpar_loc,mparray), intent(in) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      integer, dimension (mpar_loc,3) :: ineargrid
      integer :: k
      real :: taup1
      intent (inout) :: df, dfp,ineargrid
      intent (in) :: k,taup1
      real, dimension(3) :: dXp,dVp 
      real, dimension(9) :: Sij_lin
      real, dimension(3,3) :: Sijp
      integer :: ii 
!
!  Identify module.
!
      if (headtt) then
        if (lroot) print*,'dcaustics_dt_pencil: calculate dcaustics_dt_pencil'
        if (lroot) print*,'called from dvvp_dt_pencil in the particles_dust module'
      endif
!
!
!  Loop over all particles in current pencil.
!
!
      dXp = fp(k,idXp1:idXp3)
      dVp = fp(k,idVp1:idVp3)
!
! solve (d/dt) \deltaX_k = \deltaV_k 
!
      dfp(k,idXp1:idXp3)= dVp 
!
!  interpolate the gradu matrix to particle positions
!
      call interpolate_linear(f,igu11,igu33,fp(k,ixp:izp),Sij_lin,ineargrid(k,:), &
            0,ipar(k))
      call linarray2matrix(Sij_lin,Sijp)
!
! solve (d/dt) \deltaV_k = (1/\taup)S_kj \deltaV_j 
!
      do ii=1,3
        dfp(k,idVp1+ii-1)= taup1*(dot_product(Sijp(ii,:),dXp) - dVp(ii)) 
      enddo
!
    endsubroutine dcaustics_dt_pencil
!***********************************************************************
    subroutine read_pcaustics_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=particles_caustics_init_pars, IOSTAT=iostat)
!
    endsubroutine read_pcaustics_init_pars
!***********************************************************************
    subroutine write_pcaustics_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=particles_caustics_init_pars)
!
    endsubroutine write_pcaustics_init_pars
!***********************************************************************
    subroutine read_pcaustics_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=particles_caustics_run_pars, IOSTAT=iostat)
!
    endsubroutine read_pcaustics_run_pars
!***********************************************************************
    subroutine write_pcaustics_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=particles_caustics_run_pars)
!
    endsubroutine write_pcaustics_run_pars
!***********************************************************************
    subroutine rprint_particles_caustics(lreset,lwrite)
!
!  Read and register print parameters relevant for particles.
!
!  may-2016/dhruba+akshay: coded
!
      use Diagnostics
      use General,   only: itoa
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname,inamez,inamey,inamex,inamexy,inamexz,inamer,inamerz
      integer :: k
      logical :: lwr
      character (len=intlen) :: srad
!
!  Write information to index.pro.
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
      if (lwr) then
        write(3,*) 'idXp1=',idXp1
        write(3,*) 'idXp2=',idXp2
        write(3,*) 'idXp3=',idXp3
        write(3,*) 'idVp1=',idVp1
        write(3,*) 'idVp2=',idVp2
        write(3,*) 'idVp3=',idVp3
      endif
!
!  Reset everything in case of reset.
!
      if (lreset) then
      endif
!
!  Run through all possible names that may be listed in print.in.
!
      if (lroot .and. ip<14) print*,'rprint_particles: run through parse list'
      do iname=1,nname
      enddo
!
    endsubroutine rprint_particles_caustics
!***********************************************************************
endmodule Particles_caustics
