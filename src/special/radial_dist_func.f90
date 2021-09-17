! $Id$
!
!  Solve the lucky droplet model for many realizations.
!  The different realizations correspond to "meshpoints".
!  To add the contributions for each step, we use the usual
!  time step in the Pencil Code, so t is just the step, and
!  the accumulated collision times (after 125 steps or so)
!  for all realizations at the same time are the values in
!  the f-array.
!
!  16-apr-20/axel: adapted from nospecial.f90
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
!
! Declare index of variables
!
  integer, parameter :: npgrid=256, nbin=npgrid/2
  integer :: ispecial=0, rdf_stride_outer=1, rdf_stride_inner=1
  logical :: lloffset_search=.false.
!
  ! input parameters
  real :: gam_lucky=fourthird, efficiency_exponent=0., efficiency_prefactor=1.
  real :: rcrit=0., lambda_star1=1, runit=1.
  character(len=50) :: init_qq='read512cubed'
  character(LEN=labellen) :: file_in='../np_ap1.dat'
  namelist /special_init_pars/ &
    rdf_stride_outer, rdf_stride_inner, file_in, lloffset_search
!
  ! run parameters
  logical :: lMFT=.false., lrA=.false., lrB=.false.
  namelist /special_run_pars/ &
    gam_lucky, lMFT, lrA, lrB, efficiency_exponent, efficiency_prefactor, &
    rcrit, lambda_star1, runit
!
! other variables (needs to be consistent with reset list below)
!
  integer :: idiag_rad =0 ! DIAG_DOC: $r/r_\ast$
  integer :: idiag_tauk=0 ! DIAG_DOC: $\tau_k$
  integer :: idiag_tt1m=0 ! DIAG_DOC: $\langle T \rangle$
  integer :: idiag_qq1m=0 ! DIAG_DOC: $\langle \ln T \rangle$
  integer :: idiag_qq2m=0 ! DIAG_DOC: $\langle \ln T^2 \rangle$
  integer :: idiag_qq3m=0 ! DIAG_DOC: $\langle \ln T^3 \rangle$
  integer :: idiag_qq4m=0 ! DIAG_DOC: $\langle \ln T^4 \rangle$
!
  contains
!
!***********************************************************************
    subroutine register_special()
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!  19-feb-2019/axel: coded
!
      use FArrayManager
      use SharedVariables, only: put_shared_variable
!
!  Identify CVS/SVN version information.
!
      if (lroot) call svn_id( &
           "$Id$")
!
!  Set ichemistry to consecutive numbers nvar+1, nvar+2, ..., nvar+nchemspec.
!
      call farray_register_auxiliary('special',ispecial,array=nbin+1)
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  19-feb-2019/axel: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!  Initialize any module variables which are parameter dependent
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f)
!
!  initialise special condition; called from start.f90
!  16-apr-2020/axel: coded
!
      use Mpicomm
      use Sub
!
      real, dimension (npgrid,npgrid,npgrid) :: np
      real, dimension (mx,my,mz,mfarray) :: f
      real, parameter :: length=4*pi
      real :: distance, delx2, dely2, delz2
      integer :: ll,mm,nn, i, l, noffset, moffset, loffset, idist, ibin, nnp, mnp, lnp
      integer :: loffset_search, l1_search
      real, dimension(npgrid) :: xx = (/ (real(i*length/npgrid), i=0, npgrid-1) /)
      real, dimension(npgrid) :: yy = (/ (real(i*length/npgrid), i=0, npgrid-1) /)
      real, dimension(npgrid) :: zz = (/ (real(i*length/npgrid), i=0, npgrid-1) /)
      real, dimension(nbin+1) :: rdf, rdf_sum
!
      intent(inout) :: f
!
!  Initial condition; same for every population.
!
      select case (init_qq)
        case ('nothing'); if (lroot) print*,'init_qq: nothing'
        case ('zero'); f(:,:,:,ispecial)=0.
        case ('read512cubed')
          open(1,FILE=trim(file_in),FORM='unformatted')
          read(1) np
          close(1)
!
!  Go through all starting points. On each processor, the indices l,m,n
!  access only points within the f-array, but not the full array np.
!  Therefore we define a processor-dependent offset for each direction.
!
          f=0.
          noffset=nz*ipz
          moffset=ny*ipy
          loffset=nx*ipx
!
!  return to this point from below unless np(lnp,mnp,nnp)/=0.
!
          loffset_search=0
2000      continue
          l1_search=l1+loffset_search
!
!  By default, lloffset_search=F, so loffset_search will remain 0,
!  so the l loop with start with l1.
!
          do n=n1,n2,rdf_stride_outer
          do m=m1,m2,rdf_stride_outer
          do l=l1_search,l2,rdf_stride_outer
            nnp=n+noffset
            mnp=m+moffset
            lnp=l+loffset
            if (np(lnp,mnp,nnp)/=0.) then
!
!  Go through all points of the np array.
!  The indices l,m,n are the local indices per processor, but
!  xx,yy,zz are global, so they need to be accessed by indices
!  lnp, mnp, nnp. The indices ll,mm,nn are always global.
!
              do nn=1,npgrid,rdf_stride_inner
              do mm=1,npgrid,rdf_stride_inner
              do ll=1,npgrid,rdf_stride_inner
                if (np(ll,mm,nn)/=0.) then
                  delx2=min((xx(lnp)-xx(ll))**2, (xx(lnp)-xx(ll)-length)**2, (xx(lnp)-xx(ll)+length)**2)
                  dely2=min((yy(mnp)-yy(mm))**2, (yy(mnp)-yy(mm)-length)**2, (yy(mnp)-yy(mm)+length)**2)
                  delz2=min((zz(nnp)-zz(nn))**2, (zz(nnp)-zz(nn)-length)**2, (zz(nnp)-zz(nn)+length)**2)
                  distance=sqrt(delx2+dely2+delz2)
                  idist=int(nbin*distance/(.5*length))+1
!
!  Select all the points in a shell, or rather, find the shell into which
!  the number of pairs, np(lnp,mnp,nnp)*np(ll,mm,nn), is added.
!
                  if (idist>=1.and.idist<=nbin) then
                    f(l,m,n,idist) =f(l,m,n,idist) +np(lnp,mnp,nnp)*np(ll,mm,nn)
                    f(l,m,n,nbin+1)=f(l,m,n,nbin+1)+np(lnp,mnp,nnp)*np(ll,mm,nn)
                  endif
                endif
              enddo
              enddo
              enddo
            else
              if (lloffset_search) then
                loffset_search=loffset_search+1
                goto 2000
              endif
            endif
          enddo
          enddo
          print*,'iproc,n=',iproc,n
          enddo
!
        case default
          !
          !  Catch unknown values
          !
          if (lroot) print*,'init_qq: No such value for init_qq: ', trim(init_qq)
          call stop_it("")
      endselect
!
!  compute total sum
!
      do ibin=1,nbin+1
        rdf(ibin)=sum(f(:,:,:,ibin))
      enddo
      open (1, file=trim(directory_snap)//'/rdf.txt', form='formatted')
      write(1,'(1p,8e10.3)') rdf
      close(1)
  !
  !  Summing up the results from the different processors.
  !  The result is available only on root.
  !
  call mpireduce_sum(rdf,rdf_sum,nbin+1)
  !
  if (lroot) then
    open(1,file=trim(datadir)//'/rdf.txt')
    write(1,'(1p,8e10.3)') rdf_sum
    close(1)
  endif
!
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special()
!
!  All pencils that this special module depends on are specified here.
!
!  18-07-06/tony: coded
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine pencil_interdep_special(lpencil_in)
!
!  Interdependency among pencils provided by this module are specified here.
!
!  18-07-06/tony: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
!  Calculate Hydro pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!   24-nov-04/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine dspecial_dt(f,df,p)
!
!  calculate right hand side of ONE OR MORE extra coupled PDEs
!  along the 'current' Pencil, i.e. f(l1:l2,m,n) where
!  m,n are global variables looped over in equ.f90
!
!  Due to the multi-step Runge Kutta timestepping used one MUST always
!  add to the present contents of the df array.  NEVER reset it to zero.
!
!  several precalculated Pencils of information are passed if for
!  efficiency.
!
!  16-apr-2020/axel: coded
!
      use General, only: random_number_wrapper
      use Diagnostics
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: rrr, tauk, lamk
      real, dimension (nx) :: tt, qq, rad, ttauk
      real :: r, r2, rrunit
      type (pencil_case) :: p
!
      intent(in) :: f,p
      intent(inout) :: df
!
      call keep_compiler_quiet(p)
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dSPECIAL_dt'
!
!  Define state vector
!
      if (t==tstart) then
        tt=.0
        qq=.0
      else
        tt=f(l1:l2,m,n,ispecial)
        qq=alog(tt)
      endif
!
!  particle radius, and radius squared
!
      r=t**onethird
      r2=r**2
!
!  Selection of different combinations of rA and rB
!  Black line in paper corresponds to: lrA=T, lrB=T
!
      if (lrA.and.lrB) then
        lamk=(r+1.)**2*(r2-1.)
      elseif (lrA.and..not.lrB) then
        lamk=(r+1.)**2* r2
      elseif (.not.lrA.and.lrB) then
        lamk= r2      *(r2-1.)
      else
        lamk=t**gam_lucky
      endif
!
!  Efficiency prescription: include a powerlaw prescription
!  for t>tstar with exponent efficiency_exponent.
!  Allow also efficiency_prefactor, just to demonstrate that
!  it has no effect on the final P(T) curve, as expected.
!  Entering the inverse, lambda_star1, is also useful and more
!  intuitive (lambda_star1=1123.18s for standard parameters)
!
      rrunit=r*runit
      if (efficiency_exponent/=0. .and. rrunit>rcrit) then
        lamk=lamk*(rrunit/rcrit)**efficiency_exponent
      endif
      lamk=lamk*efficiency_prefactor/lambda_star1
!
!  Produce exponentially distributed random numbers,
!  but can also do mean-field theory as a test.
!
      if (lMFT) then
        tauk=1./lamk
      else
        call random_number_wrapper(rrr)
        tauk=-alog(rrr)/lamk
      endif
      df(l1:l2,m,n,ispecial)=df(l1:l2,m,n,ispecial)+tauk
!
!  diagnostics
!
      if (ldiagnos) then
        rad=rrunit
        ttauk=tauk
        if (idiag_rad /=0) call sum_mn_name(rad  ,idiag_rad)
        if (idiag_tauk/=0) call sum_mn_name(tauk ,idiag_tauk)
        if (idiag_tt1m/=0) call sum_mn_name(tt   ,idiag_tt1m)
        if (idiag_qq1m/=0) call sum_mn_name(qq   ,idiag_qq1m)
        if (idiag_qq2m/=0) call sum_mn_name(qq**2,idiag_qq2m)
        if (idiag_qq3m/=0) call sum_mn_name(qq**3,idiag_qq3m)
        if (idiag_qq4m/=0) call sum_mn_name(qq**4,idiag_qq4m)
      endif
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine read_special_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=special_init_pars, IOSTAT=iostat)
!
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=special_init_pars)
!
    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=special_run_pars, IOSTAT=iostat)
!
    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=special_run_pars)
!
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
!  reads and registers print parameters relevant to special
!
!  19-feb-2019/axel: coded
!
      use Diagnostics
      use FArrayManager, only: farray_index_append
      use Sub
!
!   SAMPLE IMPLEMENTATION
!
      integer :: iname
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_rad=0; idiag_tauk=0; idiag_tt1m=0
        idiag_qq1m=0; idiag_qq2m=0; idiag_qq3m=0; idiag_qq4m=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'rad' ,idiag_rad)
        call parse_name(iname,cname(iname),cform(iname),'tauk',idiag_tauk)
        call parse_name(iname,cname(iname),cform(iname),'tt1m',idiag_tt1m)
        call parse_name(iname,cname(iname),cform(iname),'qq1m',idiag_qq1m)
        call parse_name(iname,cname(iname),cform(iname),'qq2m',idiag_qq2m)
        call parse_name(iname,cname(iname),cform(iname),'qq3m',idiag_qq3m)
        call parse_name(iname,cname(iname),cform(iname),'qq4m',idiag_qq4m)
      enddo
!
    endsubroutine rprint_special
!***********************************************************************
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
