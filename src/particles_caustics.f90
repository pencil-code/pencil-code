! $Id$
!
! This modules solves for the gradient matrix of flow velocities. 
! The gradient matrix Sigma obeys the equation:
!  (d/dt) Sigma = (1/(taup))*(S - Sigma) - Sigma^2
! where A is the flow gradient matrix at the position of the 
! particle; and taup is the Stokes time. 
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_caustics=.true.
!
! MAUX CONTRIBUTION  9
! COMMUNICATED AUXILIARIES 9
! MPVAR CONTRIBUTION 9
! MPAUX CONTRIBUTION 4
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
  integer :: dummy
  real :: TrSigma_Cutoff=-1e10
!
  namelist /particles_caustics_init_pars/ &
  dummy
!
  namelist /particles_caustics_run_pars/ &
       TrSigma_Cutoff
!
! Diagnostic variables
!
  integer :: idiag_TrSigmapm=0      ! DIAG_DOC: $\langle{\rm Tr}\left[\sigma\right]\rangle$
  integer :: idiag_blowupm=0        ! DIAG_DOC: Mean no. of times $\sigma$ falls below cutoff
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
      iSigmap11=npvar+1
      pvarname(npvar+1)='iSigmap11'
      isigmap12=npvar+2
      pvarname(npvar+2)='isigmap12'
      isigmap13=npvar+3
      pvarname(npvar+3)='isigmap13'
      isigmap21=npvar+4
      pvarname(npvar+4)='isigmap21'
      isigmap22=npvar+5
      pvarname(npvar+5)='isigmap22'
      isigmap23=npvar+6
      pvarname(npvar+6)='isigmap23'
      isigmap31=npvar+7
      pvarname(npvar+7)='isigmap31'
      isigmap32=npvar+8
      pvarname(npvar+8)='isigmap32'
      isigmap33=npvar+9
      pvarname(npvar+9)='isigmap33'
!
!  Increase npvar accordingly.
!
      npvar=npvar+9
!
! Now register the particle auxiliaries
!
      iPPp=mpvar+npaux+1
      pvarname(iPPp)='iPPp'
      iQQp=iPPp+1
      pvarname(iQQp)='iQQp'
      iRRp=iPPp+2
      pvarname(iRRp)='iRRp' 
      pvarname(iblowup) = 'iblowup'
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
        call fatal_error('register_particles_caustics','npvar > mpvar')
      endif
!
! Check the same for particle auxiliaries too. 
!
      if (npaux > mpaux) then
        if (lroot) write (0,*) 'npaux = ', npaux, ', mpaux = ', mpaux
        call fatal_error('register_particles_caustics: npaux > mpaux','')
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
! We set the gradient of the particle velocity field as zero. 
!
      do ip=1,npar_loc
        fp(ip,isigmap11:isigmap33) = 0.
      enddo
!
    endsubroutine init_particles_caustics
!***********************************************************************
    subroutine dcaustics_dt(f,df,fp,dfp,ineargrid)
!
      use Diagnostics
      use Particles_sub, only: sum_par_name
      use Sub, only : linarray2matrix,Inv2_3X3mat,det3X3mat
!
      real, dimension (mx,my,mz,mfarray), intent (in) :: f
      real, dimension (mx,my,mz,mvar), intent (inout) :: df
      real, dimension (mpar_loc,mparray), intent (inout) :: fp
      real, dimension (mpar_loc,mpvar), intent (inout) :: dfp
      integer, dimension (mpar_loc,3), intent (in) :: ineargrid
      logical :: lheader, lfirstcall=.true.
      real, dimension(3,3) :: Sigmap
      real, dimension(9) :: Sigma_lin
      real :: TrSigma,QSigma,detSigma
      integer :: ip
!
!  Print out header information in first time step.
!
      lheader=lfirstcall .and. lroot
      if (lheader) then
        print*,'dcaustics_dt: Calculate dcaustics_dt'
      endif
!
! Calculates the three invariants of the matrix sigma and stores them as auxiliary
! variable. 
!
      do ip=1,npar_loc
         Sigma_lin=fp(ip,isigmap11:isigmap33)
         call linarray2matrix(Sigma_lin,Sigmap)
         TrSigma=Sigmap(1,1)+Sigmap(2,2)+Sigmap(3,3)
         fp(ip,iPPp) = TrSigma
         call Inv2_3X3mat(Sigmap,QSigma)
         fp(ip,iQQp) = QSigma
         call det3X3mat(Sigmap,detSigma)
         fp(ip,iRRp) = detSigma
      enddo
      
        if (ldiagnos) then
!
! No diagnostic is coded yet, but would be.
           if (idiag_TrSigmapm/=0) &
                call sum_par_name(fp(1:npar_loc,iPPp),idiag_TrSigmapm)
           if (idiag_blowupm/=0) &
                call sum_par_name(fp(1:npar_loc,iblowup),idiag_blowupm)
           
           !
!          do k=1,ndustrad
!            if (idiag_npvzmz(k)/=0) call xysum_mn_name_z(p%npvz(:,k),idiag_npvzmz(k))
!            if (idiag_npcaustics(k)/=0) call xysum_mn_name_z(ncaustics(k),idiag_npcaustics(k))
!          enddo
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
      real, dimension(9) :: Sij_lin,Sigma_lin,dSigmap_lin
      real, dimension(3,3) :: Sijp,Sigmap,dSigmap
!
!  Identify module.
!
      if (headtt.and.lfirstpoint) then
        if (lroot) print*,'dcaustics_dt_pencil: calculate dcaustics_dt_pencil'
        if (lroot) print*,'called from dvvp_dt_pencil in the particles_dust module'
      endif
!
!
!  Loop over all particles in current pencil.
!
!
      Sigma_lin = fp(k,isigmap11:isigmap33)
      call linarray2matrix(Sigma_lin,Sigmap)
!
!  interpolate the gradu matrix to particle positions
!
      call interpolate_linear(f,igu11,igu33,fp(k,ixp:izp),Sij_lin,ineargrid(k,:), &
            0,ipar(k))
      call linarray2matrix(Sij_lin,Sijp)
!
! solve dynamical equation
!
      dSigmap = taup1*(Sijp-Sigmap) - matmul(Sigmap,Sigmap)
      call matrix2linarray(dSigmap,dSigmap_lin)
      dfp(k,isigmap11:isigmap33) = dSigmap_lin
!
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
        write(3,*) 'isigmap11=', isigmap11
        write(3,*) 'isigmap12=', isigmap13
        write(3,*) 'isigmap13=', isigmap13
        write(3,*) 'isigmap21=', isigmap21
        write(3,*) 'isigmap22=', isigmap22
        write(3,*) 'isigmap23=', isigmap23
        write(3,*) 'isigmap31=', isigmap31
        write(3,*) 'isigmap32=', isigmap32
        write(3,*) 'isigmap33=', isigmap33
        write(3,*) 'iPPp=', iPPp
        write(3,*) 'iQQp=', iQQp
        write(3,*) 'iRRp=', iRRp
        write(3,*) 'iblowup=', iblowup
      endif
!
!  Reset everything in case of reset.
!
      if (lreset) then
         idiag_TrSigmapm=0
      endif
!
!  Run through all possible names that may be listed in print.in.
!
      if (lroot .and. ip<14) print*,'rprint_particles: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'trsigmapm',idiag_TrSigmapm)
     enddo
!
    endsubroutine rprint_particles_caustics
!***********************************************************************
    subroutine reset_caustics(fp)
      real, dimension (mpar_loc,mparray), intent (out) :: fp
      integer :: ip
      real :: TrSigma, QSigma, RSigma
      real, dimension(9) :: Sigma_lin
      real, dimension(3,3) :: Sigmap

      do ip=1,npar_loc
         if (fp(ip,iPPp).lt.TrSigma_Cutoff) then
            fp(ip,iblowup)=fp(ip,iblowup)+1.
            fp(ip,isigmap11:isigmap33) = 0.
         endif
      enddo
    endsubroutine reset_caustics
!***********************************************************************
endmodule Particles_caustics
