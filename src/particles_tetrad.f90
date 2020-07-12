! $Id$
!
! Follow a particle. For this particle solve for three other particles
! that that are in the neighbourhood of the particle we were following.
! These four particles together forms a tetrad. 
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_tetrad=.true.
!
! MAUX CONTRIBUTION  9
! COMMUNICATED AUXILIARIES 9
! MPVAR CONTRIBUTION 18
! MPAUX CONTRIBUTION 2
!
!***************************************************************
module Particles_tetrad
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
  include 'particles_tetrad.h'
!
  integer :: dummy
  real :: rsmall = 1e-6
!
  namelist /particles_tetrad_init_pars/ &
  rsmall
!
  namelist /particles_tetrad_run_pars/ &
       rsmall
!
! Diagnostic variables
!
  integer :: idiag_TVolm            ! DIAG_DOC: Mean absolute volume of the tetrads
  integer :: idiag_TVolpm            ! DIAG_DOC: Mean of positive volume of the tetrads
  integer :: idiag_TVolnm           ! DIAG_DOC: Mean of negative volume of the tetrads
  integer :: idiag_VelVolm            ! DIAG_DOC: Mean absolute volume of the tetrads in velocity space
  integer :: idiag_VelVolpm           ! DIAG_DOC: Mean of positive volume of the tetrads in velocity space
  integer :: idiag_VelVolnm           ! DIAG_DOC: Mean of negative volume of the tetrads in velocity space.
  integer :: idiag_pspaceVolm         !DIAG_DOC : mean of the phase-space volume always positive.
contains
!***********************************************************************
    subroutine register_particles_tetrad()
!
!  Set up indices for access to the fp and dfp arrays
!
!  May-16/dhruba: coded
!
      use FArrayManager, only: farray_register_auxiliary
!
      if (lroot) call svn_id("particles_tetrad")
!
!  Indices for first particle
!
      
      call append_npvar('idR11',idR11)
      call append_npvar('idR12',idR12)
      call append_npvar('idR13',idR13)
      call append_npvar('idV11',idV11)
      call append_npvar('idV12',idV12)
      call append_npvar('idV13',idV13)
!
!  Indices for second particle
!
      
      call append_npvar('idR21',idR21)
      call append_npvar('idR22',idR22)
      call append_npvar('idR23',idR23)
      call append_npvar('idV21',idV21)
      call append_npvar('idV22',idV22)
      call append_npvar('idV23',idV23)
!
!  Indices for third particle
!
      
      call append_npvar('idR31',idR31)
      call append_npvar('idR32',idR32)
      call append_npvar('idR33',idR33)
      call append_npvar('idV31',idV31)
      call append_npvar('idV32',idV32)
      call append_npvar('idV33',idV33)
!
! One particle auxiliary that tracks the volume
! of the tetrad:
!
      call append_npaux('iVolp',iVolp)
!
! One MORE particle auxiliary that tracks the 
! volume  of the tetrad in velocity space:
!
      call append_npaux('iVelVolp',iVelVolp)
!
!  Set indices for velocity gradient matrix at grid points
!
      call farray_register_auxiliary('guij',iguij,communicated=.true.,array=9)
      igu11=iguij; igu12=iguij+1; igu13=iguij+2
      igu21=iguij+3; igu22=iguij+4; igu23=iguij+5
      igu31=iguij+6; igu32=iguij+7; igu33=iguij+8
!
    endsubroutine register_particles_tetrad
!***********************************************************************
    subroutine initialize_particles_tetrad(f)
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
        call fatal_error('initialize_particles_tetrad','you must select either hydro or hydro_kinematic')
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_particles_tetrad
!***********************************************************************
    subroutine init_particles_tetrad(f,fp,ineargrid)
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
         print*, 'init_particles_tetrad: setting init. cond.'
      endif
      call reinitialize_tetrad(fp)
!
    endsubroutine init_particles_tetrad
!***********************************************************************
    subroutine dtetrad_dt(f,df,fp,dfp,ineargrid)
!
      use Particles_radius, only: get_stbin
      use Sub, only: ScalarTripleProduct
      use Diagnostics
!
      real, dimension (mx,my,mz,mfarray), intent (in) :: f
      real, dimension (mx,my,mz,mvar), intent (inout) :: df
      real, dimension (mpar_loc,mparray), intent (inout) :: fp
      real, dimension (mpar_loc,mpvar), intent (inout) :: dfp
      integer, dimension (mpar_loc,3), intent (in) :: ineargrid
      logical :: lheader, lfirstcall=.true.
      integer :: ip,iStbin,k
      real :: Vol,VelVol
      real,dimension(3) :: dR1,dR2,dR3,dV1,dV2,dV3
!
!  Print out header information in first time step.
!
      lheader=lfirstcall .and. lroot
      if (lheader) then
        print*,'dtetrad_dt: Calculate dtetrad_dt'
        lfirstcall=.false.
      endif
!
! Calculate the volume and store it as auxiliary
! variable. 
!
      do ip=1,npar_loc
         dR1=fp(ip,idR11:idR13)
         dR2=fp(ip,idR21:idR23)
         dR3=fp(ip,idR31:idR33)
         call ScalarTripleProduct(dR1,dR2,dR3,Vol)
         fp(ip,iVolp) = Vol
!
! Calculate the phase-space volume here
!
         dV1=fp(ip,idV11:idV13)
         dV2=fp(ip,idV21:idV23)
         dV3=fp(ip,idV31:idV33)
         call ScalarTripleProduct(dV1,dV2,dV3,VelVol)
         fp(ip,iVelVolp) = VelVol
      enddo
!      
        if (ldiagnos) then
           if (idiag_TVolm/=0) &
                call sum_par_name(abs(fp(1:npar_loc,iVolp)),idiag_TVolm)
                call sum_par_name(abs(fp(1:npar_loc,iVelVolp)),idiag_VelVolm)
                call sum_par_name(&
                  abs(fp(1:npar_loc,iVelVolp)*fp(1:npar_loc,iVolp)),idiag_pspaceVolm)
        endif
!
      if (lfirstcall) lfirstcall=.false.
!
    endsubroutine dtetrad_dt
!***********************************************************************
    subroutine dtetrad_dt_pencil(f,df,fp,dfp,p,ineargrid,k,taup1)
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
      real, dimension(9) :: Sij_lin
      real, dimension(3,3) :: Sijp
      real,dimension(3) :: dR1,dR2,dR3,dV1,dV2,dV3
!
!  Identify module.
!
      lheader=lfirstcall .and. lroot
      if (lheader) then
        print*,'dtetrad_dt_pencil: Calculate dtetrad_dt_pencil'
        print*,'called from dvvp_dt_pencil in the particles_dust module'
        lfirstcall=.false.
      endif
!
      dR1=fp(k,idR11:idR13)
      dV1=fp(k,idV11:idV13)
      dR2=fp(k,idR21:idR23)
      dV2=fp(k,idV21:idV23)
      dR3=fp(k,idR31:idR33)
      dV3=fp(k,idV31:idV33)
!
!  interpolate the gradu matrix to particle positions
!
      call interpolate_linear(f,igu11,igu33,fp(k,ixp:izp),Sij_lin,ineargrid(k,:), &
            0,ipar(k))
      call linarray2matrix(Sij_lin,Sijp)
!
! solve dynamical equation
!
      dfp(k,idR11:idR13) = dV1
      dfp(k,idR21:idR23) = dV2
      dfp(k,idR31:idR33) = dV3
      
      dfp(k,idV11:idV13) = taup1*(matmul(Sijp,dR1)-dV1)
      dfp(k,idV21:idV23) = taup1*(matmul(Sijp,dR2)-dV2)
      dfp(k,idV31:idV33) = taup1*(matmul(Sijp,dR3)-dV3)
!
! 
    endsubroutine dtetrad_dt_pencil
!***********************************************************************
    subroutine read_ptetrad_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=particles_tetrad_init_pars, IOSTAT=iostat)
!
    endsubroutine read_ptetrad_init_pars
!***********************************************************************
    subroutine write_ptetrad_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=particles_tetrad_init_pars)
!
    endsubroutine write_ptetrad_init_pars
!***********************************************************************
    subroutine read_ptetrad_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=particles_tetrad_run_pars, IOSTAT=iostat)
!
    endsubroutine read_ptetrad_run_pars
!***********************************************************************
    subroutine write_ptetrad_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=particles_tetrad_run_pars)
!
    endsubroutine write_ptetrad_run_pars
!***********************************************************************
    subroutine rprint_particles_tetrad(lreset,lwrite)
!
!  Read and register print parameters relevant for particles.
!
!  may-2016/dhruba+akshay: coded
!
      use Diagnostics
      use General, only: itoa
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname
      integer :: k
      character (len=intlen) :: srad
!
!  Reset everything in case of reset.
!
      if (lreset) then
         idiag_TVolm=0
         idiag_TVolpm=0
         idiag_TVolnm=0
         idiag_VelVolm=0
         idiag_pspaceVolm=0
      endif
!
!  Run through all possible names that may be listed in print.in.
!
      if (lroot .and. ip<14) print*,'rprint_particles: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'tvolm',idiag_TVolm)
        call parse_name(iname,cname(iname),cform(iname),'tvel_volm',idiag_VelVolm)
        call parse_name(iname,cname(iname),cform(iname),'tpvolm',idiag_pspaceVolm)
     enddo
!
    endsubroutine rprint_particles_tetrad
!***********************************************************************
    subroutine reinitialize_tetrad(fp)
      real, dimension (mpar_loc,mparray), intent (out) :: fp
      integer :: ip
      if (lroot) then
         print*, 'The tetrads are always reinitialized,'
         print*, 'even when we restart from earlier runs.'
      endif
!
! 
!
      do ip=1,npar_loc
         fp(ip,idR11:idR13) = 0.
         fp(ip,idR11) = rsmall
         fp(ip,idV11:idV13) = 0.
!
         fp(ip,idR21:idR23) = 0.
         fp(ip,idR22) = rsmall
         fp(ip,idV21:idV23) = 0.
!
         fp(ip,idR31:idR33) = 0.
         fp(ip,idR33) = rsmall
         fp(ip,idV31:idV33) = 0.
         
!
         fp(ip,iVolp) = rsmall**3
     enddo

   endsubroutine reinitialize_tetrad
!***********************************************************************
endmodule Particles_tetrad
