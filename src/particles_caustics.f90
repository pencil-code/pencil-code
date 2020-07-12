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
! MPVAR CONTRIBUTION 10
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
  use Particles_sub
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
  integer :: idiag_lnVpm=0          ! DIAG_DOC: Mean of (logarithm of) Volume around an inertial particle
  integer, dimension(ndustrad)  :: idiag_ncaustics=0
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
!  Indices for 9 elements of the sigmap matrix
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
! One extra variable that calculate the evolution
! of (logarithm of) volume around an inertial particle. 
!
      call append_npvar('lnVp',ilnVp)
!
! Now register the particle auxiliaries
!
      call append_npaux('iPPp',iPPp)
      call append_npaux('iQQp',iQQp)
      call append_npaux('iRRp',iRRp)
      call append_npaux('iblowup',iblowup)
!
!  Set indices for velocity gradient matrix at grid points
!
      call farray_register_auxiliary('guij',iguij,communicated=.true.,array=9)
      igu11=iguij; igu12=iguij+1; igu13=iguij+2
      igu21=iguij+3; igu22=iguij+4; igu23=iguij+5
      igu31=iguij+6; igu32=iguij+7; igu33=iguij+8
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
      use General, only: keep_compiler_quiet,random_number_wrapper
      use Mpicomm, only:  mpiallreduce_sum_int
      use Hydro, only: calc_gradu
      real, dimension (mx,my,mz,mfarray), intent (in) :: f
      real, dimension (mpar_loc,mparray), intent (out) :: fp
      integer, dimension (mpar_loc,3), intent (in) :: ineargrid
!
!
      if (lroot) then
         print*, 'init_particles_caustics: setting init. cond.'
      endif
      call reinitialize_caustics(fp)
!
    endsubroutine init_particles_caustics
!***********************************************************************
    subroutine dcaustics_dt(f,df,fp,dfp,ineargrid)
!
      use Particles_radius, only: get_stbin
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
      integer :: ip,iStbin,k
      real :: prad
      real,dimension(ndustrad) :: blowup_Stdep,np_Stdep
!
!  Print out header information in first time step.
!
      lheader=lfirstcall .and. lroot
      if (lheader) then
        print*,'dcaustics_dt: Calculate dcaustics_dt'
        lfirstcall=.false.
      endif
!
! Calculates the three invariants of the matrix sigma and stores them as auxiliary
! variable. 
!
      blowup_Stdep=0.
      np_Stdep=0.
      do ip=1,npar_loc
         Sigma_lin=fp(ip,isigmap11:isigmap33)
         call linarray2matrix(Sigma_lin,Sigmap)
         TrSigma=Sigmap(1,1)+Sigmap(2,2)+Sigmap(3,3)
         fp(ip,iPPp) = TrSigma
         call Inv2_3X3mat(Sigmap,QSigma)
         fp(ip,iQQp) = QSigma
         call det3X3mat(Sigmap,detSigma)
         fp(ip,iRRp) = detSigma
         prad=fp(ip,iap)
         call get_stbin(prad,iStbin)
         blowup_Stdep(iStbin) = blowup_Stdep(iStbin)+fp(ip,iblowup)
         np_Stdep(iStbin) = np_Stdep(iStbin)+1
      enddo
!      
        if (ldiagnos) then
           if (idiag_TrSigmapm/=0) &
                call sum_par_name(fp(1:npar_loc,iPPp),idiag_TrSigmapm)
           if (idiag_blowupm/=0) &
                call sum_par_name(fp(1:npar_loc,iblowup),idiag_blowupm)
           if (idiag_lnVpm/=0) &
                call sum_par_name(fp(1:npar_loc,ilnVp),idiag_lnVpm)
           do k=1,ndustrad
              if (idiag_ncaustics(k)/=0) &
                   fname(idiag_ncaustics(k))=blowup_Stdep(k)/np_Stdep(k)
           enddo
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
      logical :: lheader,lfirstcall=.true.
      intent (inout) :: df, dfp,ineargrid
      intent (in) :: k,taup1
      real, dimension(9) :: Sij_lin,Sigma_lin,dSigmap_lin
      real, dimension(3,3) :: Sijp,Sigmap,dSigmap
!
!  Identify module.
!
      lheader=lfirstcall .and. lroot
      if (lheader) then
        print*,'dcaustics_dt_pencil: Calculate dcaustics_dt_pencil'
        print*,'called from dvvp_dt_pencil in the particles_dust module'
        lfirstcall=.false.
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
! update the (logarithm of) volume
!
      dfp(k,ilnVp) = Sigmap(1,1)+Sigmap(2,2)+Sigmap(3,3)
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
      use FArrayManager, only: farray_index_append
      use General,   only: itoa
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname
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
        call farray_index_append('isigmap11', isigmap11)
        call farray_index_append('isigmap12', isigmap13)
        call farray_index_append('isigmap13', isigmap13)
        call farray_index_append('isigmap21', isigmap21)
        call farray_index_append('isigmap22', isigmap22)
        call farray_index_append('isigmap23', isigmap23)
        call farray_index_append('isigmap31', isigmap31)
        call farray_index_append('isigmap32', isigmap32)
        call farray_index_append('isigmap33', isigmap33)
        call farray_index_append('ilnVp', ilnVp)
        call farray_index_append('iPPp', iPPp)
        call farray_index_append('iQQp', iQQp)
        call farray_index_append('iRRp', iRRp)
        call farray_index_append('iblowup', iblowup)
      endif
!
!  Reset everything in case of reset.
!
      if (lreset) then
         idiag_TrSigmapm=0
         idiag_blowupm=0
         idiag_ncaustics=0
         idiag_lnVpm=0
      endif
!
!  Run through all possible names that may be listed in print.in.
!
      if (lroot .and. ip<14) print*,'rprint_particles: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'trsigmapm',idiag_TrSigmapm)
        call parse_name(iname,cname(iname),cform(iname),'blowupm',idiag_blowupm)
        call parse_name(iname,cname(iname),cform(iname),'lnVpm',idiag_lnVpm)
        do k=1,ndustrad
           srad=itoa(k)
           call parse_name(iname,cname(iname),cform(iname), &
                'ncaustics'//trim(srad),idiag_ncaustics(k))
        enddo
     enddo
!
    endsubroutine rprint_particles_caustics
!***********************************************************************
    subroutine reinitialize_caustics(fp)
      real, dimension (mpar_loc,mparray), intent (out) :: fp
      integer :: ip
      if (lroot) then
         print*, 'The sigmap matrix always starts from zero,'
         print*, 'even when we restart from earlier runs.'
      endif
!
! We set the gradient of the particle velocity field as zero. 
!
      do ip=1,npar_loc
         fp(ip,isigmap11:isigmap33) = 0.
         fp(ip,ilnVp) = 0.
         fp(ip,iPPp) = 0.
         fp(ip,iQQp) = 0.
         fp(ip,iRRp) = 0.
         fp(ip,iblowup) = 0.
     enddo

   endsubroutine reinitialize_caustics
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
