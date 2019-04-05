! $Id$
!
! This module calculates the probability distribution function
! the particles that have moved a certain distance away from their
! initial position. 
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_persistence=.true.
!
! MPAUX CONTRIBUTION 4
!
!***************************************************************
module Particles_persistence
!
  use Cdata
  use Messages
  use Particles_cdata
  use Particles_map
  use Particles_mpicomm
!
  implicit none
!
  include 'particles_persistence.h'
!
  real,dimension(NPARTDISP) :: RR=0.0
  integer(kind=4),dimension(NPARTDISP) :: pinside=0

!
  namelist /particles_persistence_init_pars/ &
   RR 
!
  namelist /particles_persistence_run_pars/ &
   RR
!
  integer :: idiag_ds2pm=0
!
contains
!***********************************************************************
    subroutine register_particles_persistence
!
!  Set up indices for access to the fp and dfp arrays
!
!  May-16/dhruba: coded
!
      use FArrayManager, only: farray_register_auxiliary
!
      if (lroot) call svn_id( &
          "$Id$")
!
!  Indices to keep track of persistence particle positions
!
     ixp0=mpvar+npaux+1
     iyp0=mpvar+npaux+2
     izp0=mpvar+npaux+3
     ippersist=mpvar+npaux+4
     npaux=npaux+4
!
!  Check that the fp and dfp arrays are big enough.
!
     if (npaux > mpaux) then
        if (lroot) write (0,*) 'npaux = ', npaux, ', mpaux = ', mpaux
        call fatal_error('register_particles_mass: npaux > mpaux','')
     endif
!
    endsubroutine register_particles_persistence
!***********************************************************************
    subroutine initialize_particles_persist(f)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  13-nov-07/anders: coded
!
      use General, only: keep_compiler_quiet
      use FArrayManager
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_particles_persist
!***********************************************************************
    subroutine init_particles_persistence(fp)
!
      use Sub, only: kronecker_delta
      use General, only: keep_compiler_quiet,random_number_wrapper
      use Mpicomm, only:  mpiallreduce_sum_int
      real, dimension (mpar_loc,mparray), intent (out) :: fp
      real, dimension(nx,3:3) :: uij 
      integer, dimension (ncpus) :: my_particles=0,all_particles=0
      real :: random_number
      integer :: ipzero,ik,ii,jj,ij
!
      do ip=1,npar_loc
!
! Record the initial position 
!
         fp(ip,ixp0)=fp(ip,ixp)
         fp(ip,iyp0)=fp(ip,iyp)
         fp(ip,izp0)=fp(ip,izp)
         fp(ip,ippersist) = 0.
      enddo
!      
    endsubroutine init_particles_persistence
!***********************************************************************
    subroutine dpersist_dt(f,df,fp,dfp,ineargrid)
!
!
      use Diagnostics
      use Particles_sub, only: sum_par_name
!
      real, dimension (mx,my,mz,mfarray), intent (in) :: f
      real, dimension (mx,my,mz,mvar), intent (inout) :: df
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar), intent (inout) :: dfp
      integer, dimension (mpar_loc,3), intent (in) :: ineargrid
      logical :: lheader, lfirstcall=.true.
      real,dimension(mpar_loc) :: Ds2
      real, dimension(3) :: Xp,Xp0
      real :: dR,dR2
      integer :: plabel,j
!
!  Print out header information in first time step.
!
      lheader=lfirstcall .and. lroot
      if (lheader) then
        print*,'dpersist_dt: Calculate dpersist_dt'
      endif
!
! We are calculating whether a particle has moved (for the first time) 
! beyond a certain distance from its initial position. The distances are
! given in the array RR. If it has then we increment the persistence label 
! by unity. If the particle has moved beyond the maximum distance in the array
! RR then we do not check that particle again.
!
      do ip=1,npar_loc
         Xp = fp(ip,ixp:izp)
         Xp0 = fp(ip,ixp0:izp0)
         dR2=0
         do j=1,3
            dR2 = dR2+(Xp(j)-Xp0(j))**2
         enddo
         Ds2(ip)=dR2
         dR = sqrt(dR2)
         plabel=int(fp(ip,ippersist))+1
         if ((plabel .le. NPARTDISP ).and.(dR .gt. RR(plabel))) & 
            fp(ip,ippersist) = fp(ip,ippersist)+1.
!
! The next plabel may be different from the previous one
!
         plabel=int(fp(ip,ippersist))+1
         pinside(plabel) = pinside(plabel)+1
         
      enddo
      if (ldiagnos) then
        if (idiag_ds2pm/=0) call sum_par_name(Ds2(1:npar_loc),idiag_ds2pm)
      endif
!
      if (lfirstcall) lfirstcall=.false.
!
    endsubroutine dpersist_dt
!***********************************************************************
    subroutine read_ppersist_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=particles_persistence_init_pars, IOSTAT=iostat)
!
    endsubroutine read_ppersist_init_pars
!***********************************************************************
    subroutine write_ppersist_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=particles_persistence_init_pars)
!
    endsubroutine write_ppersist_init_pars
!***********************************************************************
    subroutine read_ppersist_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=particles_persistence_run_pars, IOSTAT=iostat)
!
    endsubroutine read_ppersist_run_pars
!***********************************************************************
    subroutine write_ppersist_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=particles_persistence_run_pars)
!
    endsubroutine write_ppersist_run_pars
!***********************************************************************
    subroutine rprint_particles_persist(lreset,lwrite)
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
         call farray_index_append('ixp0',ixp0)
         call farray_index_append('iyp0',iyp0)
         call farray_index_append('izp0',izp0)
         call farray_index_append('ippersist', ippersist)
      endif
!
!  Reset everything in case of reset
!
      if (lreset) then
         idiag_ds2pm=0
      endif
!
!  Run through all possible names that may be listed in print.in
!
      if (lroot .and. ip<14) print*,'rprint_particles_persist: run through parse list'
!
      do iname=1,nname
         call parse_name(iname,cname(iname),cform(iname),'ds2pm',idiag_ds2pm)
      enddo
!
!  Reset everything in case of reset.
!
      if (lreset) then
        idiag_ds2pm=0
      endif
!
    endsubroutine rprint_particles_persist
!***********************************************************************
endmodule Particles_persistence 
