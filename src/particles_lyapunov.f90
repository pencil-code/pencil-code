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
! CPARAM logical, parameter :: lparticles_lyapunov=.true.
!
! MAUX CONTRIBUTION 12
! COMMUNICATED AUXILIARIES 12
! MPVAR CONTRIBUTION 12
! PENCILS PROVIDED bbf(3)
!
!***************************************************************
module Particles_lyapunov
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
  include 'particles_lyapunov.h'
!
  logical :: lnoise2pvector=.false.,linit_largeb=.false.
  real :: fake_eta=0.,bamp=1e-2,kmode_forb=3
  integer :: idiag_bx2pm=0,idiag_by2pm=0,idiag_bz2pm=0
!
  namelist /particles_lyapunov_init_pars/ &
    bamp,linit_largeb,kmode_forb
!
  namelist /particles_lyapunov_run_pars/ &
  lnoise2pvector,fake_eta
!
  contains
!***********************************************************************
    subroutine register_particles_lyapunov()
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
!  Indices for velocity gradient matrix (Wij) at particle positions
!
      call append_npvar('iup11',iup11)
      call append_npvar('iup12',iup12)
      call append_npvar('iup13',iup13)
      call append_npvar('iup21',iup21)
      call append_npvar('iup22',iup22)
      call append_npvar('iup23',iup23)
      call append_npvar('iup31',iup31)
      call append_npvar('iup32',iup32)
      call append_npvar('iup33',iup33)
!
!  Indices for a passive vector at particle positions
!
      call append_npvar('ibpx',ibpx)
      call append_npvar('ibpy',ibpy)
      call append_npvar('ibpz',ibpz)
!
!  Set indices for velocity gradient matrix at grid points
!
      call farray_register_auxiliary('guij',iguij,communicated=.true.,vector=9)
      igu11=iguij; igu12=iguij+1; igu13=iguij+2
      igu21=iguij+3; igu22=iguij+4; igu23=iguij+5
      igu31=iguij+6; igu32=iguij+7; igu33=iguij+8
!
!  Set indices for bp data mapped back on grid; this is auxilliary variable 
!
      call farray_register_auxiliary('bbf',ibbf,communicated=.true.,vector=3)
      ibxf=ibbf; ibyf=ibbf+1; ibzf=ibbf+2
!
!  Check that the fp and dfp arrays are big enough.
!
      if (npvar > mpvar) then
        if (lroot) write(0,*) 'npvar = ', npvar, ', mpvar = ', mpvar
        call fatal_error('register_particles','npvar > mpvar')
      endif
!
    endsubroutine register_particles_lyapunov
!***********************************************************************
    subroutine initialize_particles_lyapunov(f)
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
!  Stop if there is no velocity to calculate derivatives
!
      if (.not. (lhydro.or.lhydro_kinematic)) &
        call fatal_error('initialize_particles_lyapunov','you must select either hydro or hydro_kinematic')
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_particles_lyapunov
!***********************************************************************
    subroutine init_particles_lyapunov(fp)
!
      use Sub, only: kronecker_delta
      use General, only: keep_compiler_quiet,random_number_wrapper
      use Mpicomm, only:  mpiallreduce_sum_int
      real, dimension (mpar_loc,mparray), intent (out) :: fp
      real, dimension(nx,3:3) :: uij 
      integer, dimension (ncpus) :: my_particles=0,all_particles=0
      real :: random_number,xcoord
      integer :: ipzero,ik,ii,jj,ij
!
      if (lroot) then 
         write(*,*) 'init_particles_lyapunov:'
         write(*,*) 'linit_largeb=',linit_largeb
      endif
      do ip=1,npar_loc
!
! Initialize the gradU matrix at particle position
! The initial condition is kroner delta.
!
        ij=0
        do ii=1,3; do jj=1,3 
          fp(ip,iup11+ij)=kronecker_delta(ii,jj)
          ij=ij+1
        enddo;enddo
!
! Assign random initial value for the passive vector 
!
        if (linit_largeb) then
          xcoord=fp(ip,ixp)
          fp(ip,ibpx+ik-1) = bamp*sin(kmode_forb*xcoord)
        else
          do ik=1,3
            call random_number_wrapper(random_number)
            fp(ip,ibpx+ik-1)= bamp*random_number
          enddo
        endif
      enddo
!
    endsubroutine init_particles_lyapunov
!***********************************************************************
    subroutine dlyapunov_dt(f,df,fp,dfp,ineargrid)
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
        print*,'dlyapunov_dt: Calculate dlyapunov_dt'
      endif
!
      if (ldiagnos) then
        if (idiag_bx2pm/=0) then 
          call sum_par_name(fp(1:npar_loc,ibpx)*fp(1:npar_loc,ibpx),idiag_bx2pm)
        endif
        if (idiag_by2pm/=0) then 
          call sum_par_name(fp(1:npar_loc,ibpy)*fp(1:npar_loc,ibpy),idiag_by2pm)
        endif
        if (idiag_bz2pm/=0) then 
          call sum_par_name(fp(1:npar_loc,ibpz)*fp(1:npar_loc,ibpz),idiag_bz2pm)
        endif
      endif
!
      if (lfirstcall) lfirstcall=.false.
!
    endsubroutine dlyapunov_dt
!***********************************************************************
    subroutine dlyapunov_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
      use Sub, only: linarray2matrix,matrix2linarray
!
      real, dimension (mx,my,mz,mfarray),intent(in) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (mpar_loc,mparray), intent(in) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      integer, dimension (mpar_loc,3) :: ineargrid
      intent (inout) :: df, dfp,ineargrid
      real, dimension(3) :: Bp,dBp
      real, dimension(9) :: Sij_lin,Wij_lin,dWij_lin
      real, dimension(3,3) :: Sijp,Wijp,dWijp
      integer :: ix0,iy0,iz0,i,j,ij,k
!
!  Identify module.
!
      if (headtt) then
        if (lroot) print*,'dlyapunov_dt_pencil: calculate dlyapunov_dt'
      endif
!
      if (npar_imn(imn)/=0) then
!
!  Loop over all particles in current pencil.
!
        do k=k1_imn(imn),k2_imn(imn)
!
          ix0=ineargrid(k,1)
          iy0=ineargrid(k,2)
          iz0=ineargrid(k,3)
!
!  interpolate the gradu matrix to particle positions
!
          call interpolate_linear(f,igu11,igu33,fp(k,ixp:izp),Sij_lin,ineargrid(k,:), &
            0,ipar(k))
          call linarray2matrix(Sij_lin,Sijp)
!
! solve (d/dt)W_ij = S_ik W_kj
!
          Wij_lin = fp(k,iup11:iup33)
          call linarray2matrix(Wij_lin,Wijp)
          dWijp=matmul(Sijp,Wijp)
          call matrix2linarray(dWijp,dWij_lin)
          dfp(k,iup11:iup33)= dfp(k,iup11:iup33)+dWij_lin
!
! solve (d/dt) B_i = S_ik B_k 
!
          Bp=fp(k,ibpx:ibpz)
          dBp=matmul(Sijp,Bp)
          dfp(k,ibpx:ibpz)=dfp(k,ibpx:ibpz)+dBp
        enddo
      endif
!
    endsubroutine dlyapunov_dt_pencil
!***********************************************************************
    subroutine read_plyapunov_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=particles_lyapunov_init_pars, IOSTAT=iostat)
!
    endsubroutine read_plyapunov_init_pars
!***********************************************************************
    subroutine write_plyapunov_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=particles_lyapunov_init_pars)
!
    endsubroutine write_plyapunov_init_pars
!***********************************************************************
    subroutine read_plyapunov_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=particles_lyapunov_run_pars, IOSTAT=iostat)
!
    endsubroutine read_plyapunov_run_pars
!***********************************************************************
    subroutine write_plyapunov_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=particles_lyapunov_run_pars)
!
    endsubroutine write_plyapunov_run_pars
!***********************************************************************
    subroutine particles_stochastic_lyapunov(fp)
      use General,only : gaunoise_number
      real, dimension (mpar_loc,mparray), intent(inout) :: fp
      integer :: ip,ik
      real,dimension(2) :: grandom
!
      do ip=1,npar_loc
        do ik=1,3
          call gaunoise_number(grandom)
          fp(ip,ibpx+ik-1) = fp(ip,ibpx+ik-1)+sqrt(fake_eta)*grandom(2)*sqrt(dt)
        enddo
      enddo
!    
    endsubroutine particles_stochastic_lyapunov
!***********************************************************************
    subroutine calc_pencils_par_lyapunov(f,p)
!
      use Sub, only: grad
!
!  This calculates the bbf data to a pencil.
!  Most basic pencils should come first, as others may depend on them.
!
!  16-feb-06/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      if (lpencil(i_bbf)) then
        if (ibbf/=0) then
          p%bbf=f(l1:l2,m,n,ibxf:ibzf)
        else
          p%bbf=0.0
        endif
      endif
!
    endsubroutine calc_pencils_par_lyapunov
!***********************************************************************
    subroutine rprint_particles_lyapunov(lreset,lwrite)
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
      integer :: iname,inamez,inamey,inamex,inamexy,inamexz,inamer,inamerz
      integer :: k
      character (len=intlen) :: srad
!
!  Reset everything in case of reset.
!
      if (lreset) then
        idiag_bx2pm=0;idiag_by2pm=0;idiag_bz2pm=0
      endif
!
!  Run through all possible names that may be listed in print.in.
!
      if (lroot .and. ip<14) print*,'rprint_particles: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'bx2pm',idiag_bx2pm)
        call parse_name(iname,cname(iname),cform(iname),'by2pm',idiag_by2pm)
        call parse_name(iname,cname(iname),cform(iname),'bz2pm',idiag_bz2pm)
      enddo
!
    endsubroutine rprint_particles_lyapunov
!***********************************************************************
endmodule Particles_lyapunov
