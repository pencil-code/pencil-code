! $Id: particles_number.f90,v 1.1 2005-11-24 15:35:40 ajohan Exp $
!
!  This module takes care of everything related to particle number.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MPVAR CONTRIBUTION 1
! CPARAM logical, parameter :: lparticles_number=.true.
!
!***************************************************************
module Particles_number

  use Cdata
  use Particles_cdata
  use Particles_sub

  implicit none

  include 'particles_number.h'

  real :: nptilde0
  integer, dimension (mpar_loc) :: ineighbour
  integer, dimension (mx,my,mz) :: ishepherd
  character (len=labellen), dimension(ninit) :: initnptilde='nothing'

  namelist /particles_number_init_pars/ &
      initnptilde, nptilde0

  namelist /particles_number_run_pars/ &
      nptilde0

  contains

!***********************************************************************
    subroutine register_particles_number()
!
!  Set up indices for access to the fp and dfp arrays.
!
!  24-nov-05/anders: adapted
!
      use Messages, only: fatal_error, cvs_id
!
      logical, save :: first=.true.
!
      if (.not. first) call fatal_error('register_particles_number: called twice','')
      first = .false.
!
      if (lroot) call cvs_id( &
           "$Id: particles_number.f90,v 1.1 2005-11-24 15:35:40 ajohan Exp $")
!
!  Index for particle internal number.
!
      inptilde=npvar+1
!
!  Increase npvar accordingly.
!
      npvar=npvar+1
!
!  Check that the fp and dfp arrays are big enough.
!
      if (npvar > mpvar) then
        if (lroot) write(0,*) 'npvar = ', npvar, ', mpvar = ', mpvar
        call fatal_error('register_particles: npvar > mpvar','')
      endif
!
    endsubroutine register_particles_number
!***********************************************************************
    subroutine initialize_particles_number(lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  24-nov-05/anders: adapted
!
      logical :: lstarting
!
    endsubroutine initialize_particles_number
!***********************************************************************
    subroutine init_particles_number(f,fp)
!
!  Initial internal particle number.
!
!  24-nov-05/anders: adapted
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mpar_loc,mpvar) :: fp
!
      integer :: j
!
      do j=1,ninit
        
        select case(initnptilde(j))

        case('nothing')
          if (lroot.and.j==1) print*, 'init_particles_number: nothing'

        case('constant')
          if (lroot) print*, 'init_particles_number: constant internal number'
          fp(1:npar_loc,inptilde)=nptilde0

        endselect

      enddo
!
    endsubroutine init_particles_number
!***********************************************************************
    subroutine dnptilde_dt(f,df,fp,dfp)
!
!  Evolution of internal particle number.
!
!  24-oct-05/anders: adapted
!
      use Messages, only: fatal_error
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
!
      real :: deltavp
      integer :: j, k, l, m, n
      logical :: lheader, lfirstcall=.true.
!
!      intent (in) :: f, fp
!      intent (out) :: dfp
!
!  Print out header information in first time step.
!
      lheader=lfirstcall .and. lroot
!
!  Identify module and boundary conditions.
!
      if (lheader) print*,'dnptilde_dt: Calculate dnptilde_dt'
!
!
!
      call nearest_gridpoint_map(f,fp)

      do l=l1,l2; do m=m1,m2; do n=n1,n2
        print*, 'AAA', l, m, n, f(l,m,n,inp)
        k=ishepherd(l,m,n)
        do while (ineighbour(k)/=0)
          j=k
          do while (ineighbour(j)/=0)
            j=ineighbour(j)
            print*, k, j
            deltavp=sqrt( &
                (fp(k,ivpx)-fp(j,ivpx))**2 + &
                (fp(k,ivpy)-fp(j,ivpy))**2 + &
                (fp(k,ivpz)-fp(j,ivpz))**2 )
            dfp(k,inptilde)=dfp(k,inptilde)+pi*(fp(j,iap)+fp(k,iap))**2*fp(j,inptilde)*fp(k,inptilde)*deltavp
          enddo
          k=ineighbour(k)
        enddo
      enddo; enddo; enddo
      stop

!
!  Diagnostic output
!
      if (ldiagnos) then
!        if (idiag_apm/=0) call sum_par_name(fp(1:npar_loc,iap),idiag_apm)
!        if (idiag_ap2m/=0) call sum_par_name(fp(1:npar_loc,iap)**2,idiag_ap2m)
      endif
!
      lfirstcall=.false.
!
    endsubroutine dnptilde_dt
!***********************************************************************
    subroutine nearest_gridpoint_map(f,fp)
!
!  Attach information about present particles to each grid point.
!
!  24-oct-05/anders: coded
!
      use Messages, only: fatal_error
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mpar_loc,mpvar) :: fp
!
      integer :: k, ix0, iy0, iz0
!
      intent (in) :: fp
!
      f(l1:l2,m1:m2,n1:n2,inp)=0.0
      ineighbour=0
      ishepherd=0
!
      do k=1,npar_loc
        call find_closest_gridpoint(fp(k,ixp:izp),ix0,iy0,iz0)
        f(ix0,iy0,iz0,inp)=f(ix0,iy0,iz0,inp)+1
        ineighbour(k)=ishepherd(ix0,iy0,iz0)
        ishepherd(ix0,iy0,iz0)=k
      enddo
!
    endsubroutine nearest_gridpoint_map
!***********************************************************************
    subroutine read_particles_num_init_pars(unit,iostat)
!    
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_number_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_number_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_num_init_pars
!***********************************************************************
    subroutine write_particles_num_init_pars(unit)
!    
      integer, intent (in) :: unit
!
      write(unit,NML=particles_number_init_pars)
!
    endsubroutine write_particles_num_init_pars
!***********************************************************************
    subroutine read_particles_num_run_pars(unit,iostat)
!    
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_number_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_number_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_num_run_pars
!***********************************************************************
    subroutine write_particles_num_run_pars(unit)
!    
      integer, intent (in) :: unit
!
      write(unit,NML=particles_number_run_pars)
!
    endsubroutine write_particles_num_run_pars
!***********************************************************************
    subroutine rprint_particles_number(lreset,lwrite)
!   
!  Read and register print parameters relevant for internal particle number.
!
!  24-aug-05/anders: adapted
!
      use Cdata
      use Sub, only: parse_name
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname
      logical :: lwr
!
!  Write information to index.pro
! 
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
      if (lwr) write(3,*) 'iap=', inptilde
!
!  Reset everything in case of reset
!
!      if (lreset) then
!        idiag_apm=0; idiag_ap2m=0
!      endif
!
!  Run through all possible names that may be listed in print.in
!
!      if (lroot.and.ip<14) &
!          print*, 'rprint_particles_radius: run through parse list'
!      do iname=1,nname
!        call parse_name(iname,cname(iname),cform(iname),'apm',idiag_apm)
!        call parse_name(iname,cname(iname),cform(iname),'ap2m',idiag_ap2m)
!      enddo
!
    endsubroutine rprint_particles_number
!***********************************************************************

endmodule Particles_number
