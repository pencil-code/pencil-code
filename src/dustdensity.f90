! $Id: dustdensity.f90,v 1.21 2003-12-29 17:09:26 ajohan Exp $

!  This module is used both for the initial condition and during run time.
!  It contains dnd_dt and init_nd, among other auxiliary routines.

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Dustdensity

  use Cparam
  use Cdata
  use Dustvelocity

  implicit none
  
  real, dimension(nx,ndustspec,ndustspec) :: dkern
  real, dimension(ndustspec) :: cdiffnd=0
  real :: nd_const=1.,dkern_cst=1.,dust_to_gas_ratio=0.,rhod0=1.,nd00=0.
  real :: cdiffnd_all
  character (len=labellen) :: initnd='zero'
  logical :: lcalcdkern=.true.

  namelist /dustdensity_init_pars/ &
       rhod0, initnd, dust_to_gas_ratio, nd_const, dkern_cst, nd00

  namelist /dustdensity_run_pars/ &
       rhod0, cdiffnd, cdiffnd_all, lcalcdkern

  ! diagnostic variables (needs to be consistent with reset list below)
  integer :: i_rhodm
  integer, dimension(ndustspec) :: i_ndm=0

  contains

!***********************************************************************
    subroutine register_dustdensity()
!
!  Initialise variables which should know that we solve the
!  compressible hydro equations: ind; increase nvar accordingly.
!
!   4-jun-02/axel: adapted from hydro
!
      use Mpicomm, only: stop_it
      use Sub
      use General, only: chn
!
      logical, save :: first=.true.
      integer :: idust
      character (len=4) :: sdust
!
      if (.not. first) call stop_it('register_dustdensity: called twice')
      first = .false.
!
      ldustdensity = .true.
!
      do idust=1,ndustspec
        if (idust .eq. 1) then
          ind(1) = iuud(1)+3         ! indix to access lam
        else
          ind(idust) = ind(idust-1)+4
        endif  
        nvar = nvar+1                 ! add 1 variable pr. dust layer
!
        if ((ip<=8) .and. lroot) then
          print*, 'register_dustdensity: idust = ', idust
          print*, 'register_dustdensity: nvar = ', nvar
          print*, 'register_dustdensity: ind = ', ind(idust)
        endif
      enddo
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id: dustdensity.f90,v 1.21 2003-12-29 17:09:26 ajohan Exp $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('register_dustdensity: nvar > mvar')
      endif
!
!  Writing files for use with IDL
!
      do idust=1,ndustspec
        call chn(idust,sdust)
        if (ndustspec .eq. 1) sdust = ''
        if (lroot) then
          if (maux == 0) then
            if (nvar < mvar) write(4,*) ',nd'//trim(sdust)//' $'
            if (nvar == mvar) write(4,*) ',nd'//trim(sdust)
          else
            write(4,*) ',nd'//trim(sdust)//' $'
          endif
          write(15,*) 'nd'//trim(sdust)//' = fltarr(mx,my,mz,1)*one'
        endif
      enddo
!
    endsubroutine register_dustdensity
!***********************************************************************
    subroutine initialize_dustdensity()
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  24-nov-02/tony: coded 
!  08-dec-03/anders: Copy *_all parameters to whole array
      integer :: idust
!
      select case (initnd)     
      
      case('kernel_cst')
        dkern(:,:,:) = dkern_cst
        lcalcdkern=.false.

      endselect
!
!  If *_all set, make all empty *(:) = *_all
!
      if (cdiffnd_all .ne. 0.) then
        if (lroot .and. ip<6) &
            print*, 'initialize_dustdensity: cdiffnd_all=',cdiffnd_all
        do idust=1,ndustspec
          if (cdiffnd(idust) .eq. 0.) cdiffnd(idust) = cdiffnd_all
        enddo
      endif
!
    endsubroutine initialize_dustdensity
!***********************************************************************
    subroutine init_nd(f,xx,yy,zz)
!
!  initialise nd; called from start.f90
!
!  7-nov-01/wolf: coded
! 28-jun-02/axel: added isothermal
!
      use Mpicomm
      use Sub
      use IO
      use Global
      use Gravity
      use Initcond
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
      integer :: i,j,idust
!
!  different initializations of nd (called from start).
!
      select case(initnd)
 
      case('zero'); if(lroot) print*,'init_nd: zero nd'
      case('first')
        f(:,:,:,ind(1):ind(ndustspec)) = 0.
        f(:,:,:,ind(1)) = nd00
      case('MRN77')
      case('const_nd'); f(:,:,:,ind) = nd_const
      case('frac_of_gas_loc')
        if (dust_to_gas_ratio .lt. 0.) &
            call stop_it("init_nd: Negative dust_to_gas_ratio!")
        do idust=1,ndustspec
          f(:,:,:,ind(idust)) = dust_to_gas_ratio/mg(idust)*exp(f(:,:,:,ilnrho))
        enddo
      case('frac_of_gas_glo')
        if (dust_to_gas_ratio .lt. 0.) &
            call stop_it("init_nd: Negative dust_to_gas_ratio!")
        do i=1,mx
          do j=1,my
            do idust=1,ndustspec
              f(i,j,:,ind(idust)) = &
                  dust_to_gas_ratio/mg(idust)*exp(f(4,4,:,ilnrho))
            enddo
          enddo
        enddo
      case('kernel_cst')
        print*, 'init_nd: Test case with constant kernel'
        f(:,:,:,ind(1):ind(ndustspec)) = 0.
        f(:,:,:,ind(1)) = nd00
      case default
!
!  Catch unknown values
!
        if (lroot) print*, 'init_nd: No such value for initnd: ', &
            trim(initnd)
        call stop_it("")

      endselect
!
!  sanity check
!
      if ( notanumber(f(:,:,:,ind)) ) then
        STOP "init_nd: Imaginary dustdensity values"
      endif
!
      if(ip==0) print*,xx,yy,zz ! keep compiler quiet
!
    endsubroutine init_nd
!***********************************************************************
    subroutine dnd_dt(f,df,uu,uud,divud,gnd)
!
!  continuity equation
!  calculate dnd/dt = - u.gradnd - nd*divud
!
!   7-jun-02/axel: incoporated from subroutine pde
!
      use Sub
      use Density, only: cs0
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3,3) :: udij
      real, dimension (nx,3,ndustspec) :: gnd
      real, dimension (nx,ndustspec) :: nd,uu,uud,divud
      real, dimension (nx) :: ugnd,gnd2,del2nd,udiudj,dndfac,deltaud
      real :: diffnd
      integer :: i,j,k,idust
      logical :: lfirstpoint2
!
      intent(in)  :: f,uu,uud,divud
      intent(out) :: df,gnd
!
!  Abbreviations
!
      do idust=1,ndustspec
        nd(:,idust) = f(l1:l2,m,n,ind(idust))
      enddo
!
!  Calculate kernel of coagulation equation
!
      if (lcalcdkern) then
        do i=1,ndustspec
          do j=i,ndustspec
            call dot_mn (f(l1:l2,m,n,iudx(i):iudz(i)), &
                f(l1:l2,m,n,iudx(j):iudz(j)),udiudj)
            deltaud = sqrt(udiudj)
            dkern(:,i,j) = scolld(i,j)*deltaud
            dkern(:,j,i) = dkern(:,i,j)
          enddo
        enddo
      endif
!
!  Grain growth
!
      do i=1,ndustspec
        do j=i,ndustspec
          dndfac = -dkern(:,i,j)*nd(:,i)*nd(:,j)
          df(l1:l2,m,n,ind(i)) = df(l1:l2,m,n,ind(i)) + dndfac
          df(l1:l2,m,n,ind(j)) = df(l1:l2,m,n,ind(j)) + dndfac
          forall (k=1:ndustspec, &
              mg(i) + mg(j) .ge. mglow(k) .and. mg(i) + mg(j) .lt. mghig(k)) 
            df(l1:l2,m,n,ind(k)) = &
                df(l1:l2,m,n,ind(k)) - dndfac*(mg(i)+mg(j))/mg(k)
          endforall
        enddo
      enddo
!
!  Loop over dust layers
!
      do idust=1,ndustspec
!
!  identify module and boundary conditions
!
        if (headtt.or.ldebug) print*,'dnd_dt: SOLVE dnd_dt'
        if (headtt) call identify_bcs('nd',ind(idust))
!
!  calculate dustdensity gradient and advection term
!
        call grad(f,ind(idust),gnd(:,:,idust))
        call dot_mn(uud(:,idust),gnd(:,:,idust),ugnd)
!
!  continuity equation
!
        df(l1:l2,m,n,ind(idust)) = df(l1:l2,m,n,ind(idust)) - &
            ugnd - f(l1:l2,m,n,ind(idust))*divud(:,idust)
!
!  mass diffusion, in units of dxmin*cs0
!
        if (cdiffnd(idust) /= 0.) then
          diffnd=cdiffnd(idust)*dxmin*cs0
          call del2(f,ind(idust),del2nd)
          call dot2_mn(gnd(:,:,idust),gnd2)
          df(l1:l2,m,n,ind(idust)) = &
              df(l1:l2,m,n,ind(idust)) + diffnd*(del2nd+nd(:,idust)*gnd2)
          maxdiffus=amax1(maxdiffus,diffnd)
        endif
!
        if (ldiagnos) then
          if (i_ndm(idust)/=0) call sum_mn_name(nd(:,idust),i_ndm(idust))
          if (i_rhodm/=0) then
            if (lfirstpoint .and. idust .ne. 1) then
              lfirstpoint2 = .true.
              lfirstpoint = .false.
            endif
            call sum_mn_name(nd(:,idust)*mg(idust),i_rhodm)
            lfirstpoint = lfirstpoint2
          endif
        endif
!
!  End loop over dust layers
!
      enddo
!
    endsubroutine dnd_dt
!***********************************************************************
    subroutine rprint_dustdensity(lreset,lwrite)
!
!  reads and registers print parameters relevant for compressible part
!
!   3-may-02/axel: coded
!  27-may-02/axel: added possibility to reset list
!
      use Sub
      use General, only: chn
!
      integer :: iname,idust
      logical :: lreset,lwr
      logical, optional :: lwrite
      character (len=4) :: sdust,sdustspec,snd1
!
!  Write information to index.pro that should not be repeated for idust
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite

      if (lwr) then
        write(3,*) 'ndustspec=',ndustspec
        write(3,*) 'nname=',nname
      endif
!
!  reset everything in case of reset
!
      if (lreset) then
        i_ndm = 0
        i_rhodm = 0
      endif
!
!  Define arrays for multiple dust species
!
      if (lwr) then
        call chn(ndustspec,sdustspec)
        write(3,*) 'ind=intarr('//trim(sdustspec)//')'
        write(3,*) 'i_ndm=intarr('//trim(sdustspec)//')'
      endif
!
!  Loop over dust species (for species-dependent diagnostics)
!
      do idust=1,ndustspec
        call chn(idust-1,sdust)
        if (ndustspec .eq. 1) sdust=''
!
!  iname runs through all possible names that may be listed in print.in
!
        if(lroot.and.ip<14) print*,'rprint_dustdensity: run through parse list'
        do iname=1,nname
          call parse_name(iname,cname(iname),cform(iname), &
              'ndm'//trim(sdust),i_ndm(idust))
        enddo
!
!  write column where which variable is stored
!
        if (lwr) then
          if (i_ndm(idust) .ne. 0) &
              write(3,*) 'i_ndm('//trim(sdust)//')=',i_ndm(idust)
        endif
!
!  End loop over dust layers
!
      enddo
!
!  Non-species-dependent diagnostics
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'rhodm',i_rhodm)
      enddo
      if (lwr) then
        if (i_rhodm .ne. 0) write(3,*) 'i_rhodm=',i_rhodm
      endif
!
!  Write dust index in short notation
!      
      call chn(ind(1),snd1)
      if (lwr) then
        write(3,*) 'ind=indgen('//trim(sdustspec)//')*4 + '//trim(snd1)
      endif
!
    endsubroutine rprint_dustdensity
!***********************************************************************

endmodule Dustdensity
