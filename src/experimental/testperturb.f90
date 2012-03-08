! $Id$
!
!  test perturbation method
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: ltestperturb = .true.
!
!***************************************************************
module TestPerturb
!
  use Sub
  use Cdata
  use Messages

  implicit none
!
!  cosine and sine function for setting test fields and analysis
!
  real, dimension(mz) :: cz,sz,Bk1cz,Bk1sz,B1k1cz,B1k1sz,B1cz,B1sz
  integer :: njtest=0
!
!  input parameters
!
  character (len=labellen) :: itestfield='B11-B21'
  integer :: nt_testperturb=0, it_testperturb_finalize=0
  real :: ktestfield=1., ktestfield1=1., Btest0=1.
  real :: dummy
!
  include 'testperturb.h'
!
  namelist /testperturb_init_pars/ &
      dummy
!
  namelist /testperturb_run_pars/ &
      itestfield,ktestfield,Btest0,nt_testperturb
!
  integer :: idiag_alp11=0      ! DIAG_DOC: $\alpha_{11}$
  integer :: idiag_alp21=0      ! DIAG_DOC: $\alpha_{21}$
  integer :: idiag_alp31=0      ! DIAG_DOC: $\alpha_{31}$
  integer :: idiag_alp12=0      ! DIAG_DOC: $\alpha_{12}$
  integer :: idiag_alp22=0      ! DIAG_DOC: $\alpha_{22}$
  integer :: idiag_alp32=0      ! DIAG_DOC: $\alpha_{32}$
  integer :: idiag_eta11=0      ! DIAG_DOC: $\eta_{113}k$
  integer :: idiag_eta21=0      ! DIAG_DOC: $\eta_{213}k$
  integer :: idiag_eta31=0      ! DIAG_DOC: $\eta_{313}k$
  integer :: idiag_eta12=0      ! DIAG_DOC: $\eta_{123}k$
  integer :: idiag_eta22=0      ! DIAG_DOC: $\eta_{223}k$
  integer :: idiag_eta32=0      ! DIAG_DOC: $\eta_{323}k$
!
  contains
!***********************************************************************
    subroutine register_testperturb()
!
!  Initialise variables
!
!  2-july-02/nils: coded
!
      use Mpicomm, only: stop_it
!
!  Identify version number.
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_testperturb
!***********************************************************************
    subroutine initialize_testperturb()
!
!  22-mar-08/axel: adapted from testfield_z
!
!  pre-calculate cos(kz) and sin(kz), as well as 1/ktestfield
!
      if (lroot.and.ip<=12) &
          print*,'initialize_testperturb: ktestfield=',ktestfield
!
!  set cosine and sine function for setting test fields and analysis
!
      cz=cos(ktestfield*z)
      sz=sin(ktestfield*z)
!
!  Also calculate its inverse, but only if different from zero
!
      if (ktestfield==0) then
        ktestfield1=1.
      else
        ktestfield1=1./ktestfield
      endif
!
!  other abbreviations
!
      Bk1cz=Btest0*ktestfield1*cz
      Bk1sz=Btest0*ktestfield1*sz
!
!  scale with inverse Btest0
!
      B1cz=cz/Btest0
      B1sz=sz/Btest0
!
!  scale with inverse Btest0 and with inverse ktestfield
!
      B1k1cz=ktestfield1*B1cz
      B1k1sz=ktestfield1*B1sz
!
!  determine number of testfields, njest
!
        select case (itestfield)
          case ('B11-B21+B=0'); njtest=3
          case ('B11-B21'); njtest=2
          case ('B11-B22'); njtest=4
          case ('B=0'); njtest=1
        case default
          call fatal_error('testperturb','undefined itestfield value')
        endselect
!
!  write testfield information to a file (for convenient post-processing)
!
      if (lroot) then
        open(1,file=trim(datadir)//'/testperturb_info.dat',STATUS='unknown')
        write(1,'(3a)') "itestfield='",trim(itestfield)//"'"
        write(1,'(a,f5.2)') 'ktestfield=',ktestfield
        write(1,'(a,i5)') 'nt_testperturb=',nt_testperturb
        close(1)
      endif
!
!  print a warning if nt_testperturb exceeds it1
!
      if (nt_testperturb>it1) then
        if (lroot) print*,'WARNING: nt_testperturb>it1 may be problematic'
      endif
!
    endsubroutine initialize_testperturb
!***********************************************************************
    subroutine read_testperturb_init_pars(unit,iostat)
!
!  read initial testperturb parameters
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=testperturb_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=testperturb_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_testperturb_init_pars
!***********************************************************************
    subroutine write_testperturb_init_pars(unit)
!
!  write initial testperturb parameters
!
      integer, intent(in) :: unit
!
      write(unit,NML=testperturb_init_pars)
!
    endsubroutine write_testperturb_init_pars
!***********************************************************************
    subroutine read_testperturb_run_pars(unit,iostat)
!
!  read run testperturb parameters
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=testperturb_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=testperturb_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_testperturb_run_pars
!***********************************************************************
    subroutine write_testperturb_run_pars(unit)
!
!  write run testperturb parameters
!
      integer, intent(in) :: unit
!
      write(unit,NML=testperturb_run_pars)
!
    endsubroutine write_testperturb_run_pars
!***********************************************************************
    subroutine testperturb_begin(f,df)
!
!  Calculate the response to perturbations after one timestep.
!
!  19-mar-08/axel: coded
!
      use Diagnostics
      use Hydro, only: calc_pencils_hydro
      use Timestep, only: time_step

      real, dimension (mx,my,mz,mfarray) :: fsave,f
      real, dimension (mx,my,mz,mvar) :: df
!
      real, dimension (nx,3) :: uu,bb,uxb
      logical :: headtt_save
      integer :: jtest,it_testperturb
      real :: tsave
      type (pencil_case) :: p
!
!  print identifier
!
      if (headtt.or.ldebug) print*, 'testperturb_begin'
!
!  continue only if lout=.true.
!
      if (ip<13) print*,'testperturb_begin: it,t,lout=',it,t,lout
      if (lout) then
!
!  Save current f-array and current time
!
        fsave=f
        tsave=t
!
!  set timestep for final analysis step in testperturb_finalize(f)
!
        it_testperturb_finalize=it+nt_testperturb-1
        if (ip<14) print*,'testperturb_begin: it,it_testperturb_finalize=', &
                                              it,it_testperturb_finalize
!
!  do each of the test fields one after another
!
        do jtest=1,njtest
          select case (itestfield)
            case ('B11-B22'); call add_A0test_B11_B22(f,jtest)
            case ('B=0') !(dont do anything)
          case default
            call fatal_error('testperturb_begin','undefined itestfield value')
          endselect
!
!  Time advance for nt_test timesteps
!
          do it_testperturb=1,nt_testperturb
            call time_step(f,df,p)
            if (ip<13) print*,'testperturb_begin: stay in the loop; t=',t
          enddo
          if (ip<14) print*,'testperturb_begin: do analysis; t,jtest=',t,jtest
!
!  calculate emf for calculation of alpha_ij and eta_ij tensor.
!  Note that this result applies to the end of this timestep, and does not
!  agree with that calcuted in pdf for the beginning of the timestep.
!
          lfirstpoint=.true.
          headtt_save=headtt
          do m=m1,m2
          do n=n1,n2
            call calc_pencils_hydro(f,p)
            call curl(f,iaa,bb)
            call cross_mn(p%uu,bb,uxb)
!
!  Note that we have not subtracted the contributions from the EMF
!  of the mean magnetic and velocity fields. Need to check.
!
            select case (jtest)
            case (1)
              call sum_mn_name(B1cz(n)*uxb(:,1),idiag_alp11)
              call sum_mn_name(B1cz(n)*uxb(:,2),idiag_alp21)
              call sum_mn_name(B1cz(n)*uxb(:,3),idiag_alp31)
              call sum_mn_name(-B1k1sz(n)*uxb(:,1),idiag_eta11)
              call sum_mn_name(-B1k1sz(n)*uxb(:,2),idiag_eta21)
              call sum_mn_name(-B1k1sz(n)*uxb(:,3),idiag_eta31)
            case (2)
              call sum_mn_name(B1sz(n)*uxb(:,1),idiag_alp11,ipart=2)
              call sum_mn_name(B1sz(n)*uxb(:,2),idiag_alp21,ipart=2)
              call sum_mn_name(B1sz(n)*uxb(:,3),idiag_alp31,ipart=2)
              call sum_mn_name(B1k1cz(n)*uxb(:,1),idiag_eta11,ipart=2)
              call sum_mn_name(B1k1cz(n)*uxb(:,2),idiag_eta21,ipart=2)
              call sum_mn_name(B1k1cz(n)*uxb(:,3),idiag_eta31,ipart=2)
            case (3)
              call sum_mn_name(B1cz(n)*uxb(:,1),idiag_alp12)
              call sum_mn_name(B1cz(n)*uxb(:,2),idiag_alp22)
              call sum_mn_name(B1cz(n)*uxb(:,3),idiag_alp32)
              call sum_mn_name(-B1k1sz(n)*uxb(:,1),idiag_eta12)
              call sum_mn_name(-B1k1sz(n)*uxb(:,2),idiag_eta22)
              call sum_mn_name(-B1k1sz(n)*uxb(:,3),idiag_eta32)
            case (4)
              call sum_mn_name(B1sz(n)*uxb(:,1),idiag_alp12,ipart=2)
              call sum_mn_name(B1sz(n)*uxb(:,2),idiag_alp22,ipart=2)
              call sum_mn_name(B1sz(n)*uxb(:,3),idiag_alp32,ipart=2)
              call sum_mn_name(B1k1cz(n)*uxb(:,1),idiag_eta12,ipart=2)
              call sum_mn_name(B1k1cz(n)*uxb(:,2),idiag_eta22,ipart=2)
              call sum_mn_name(B1k1cz(n)*uxb(:,3),idiag_eta32,ipart=2)
            endselect
            headtt=.false.
            lfirstpoint=.false.
          enddo
          enddo
!
!  reset f to what it was before.
!
          f=fsave
          t=tsave
        enddo
!
!  Since this routine was called before any regular call to the
!  timestepping routine, we must reset headtt to what it was before.
!
        headtt=headtt_save
      else
        if (ip<13) print*,'testperturb_begin: skip; lout,it=',lout,it
      endif
!
    endsubroutine testperturb_begin
!***********************************************************************
    subroutine testperturb_finalize(f)
!
!  Finalize testperturb method by subtracting contribution
!  emf of the reference state.
!
!  22-mar-08/axel: adapted from testperturb_begin
!
      use Diagnostics
      use Equ, only: diagnostic
      use Hydro, only: calc_pencils_hydro

      real, dimension (mx,my,mz,mfarray) :: f
!
      real, dimension (nx,3) :: uu,bb,uxb
      logical :: headtt_save
      integer :: jtest
      real :: tmp
      type (pencil_case) :: p
!
!  print identifier
!
      if (headtt.or.ldebug) print*, 'testperturb_finalize'
!
!  continue only if it >= it_testperturb
!
      if (it<it_testperturb_finalize) then
        if (ip<14) print*,'testperturb_finalize: return; t=',t
      else
        if (ip<14) print*,'testperturb_finalize: final analysis; t=',t
!
!  calculate emf for calculation of alpha_ij and eta_ij tensor
!
        lfirstpoint=.true.
        do m=m1,m2
        do n=n1,n2
          call calc_pencils_hydro(f,p)
          call curl(f,iaa,bb)
          call cross_mn(p%uu,bb,uxb)
!
!  alpha terms
!
          tmp=-(B1cz(n)+B1sz(n))
          call sum_mn_name(tmp*uxb(:,1),idiag_alp11,ipart=3)
          call sum_mn_name(tmp*uxb(:,2),idiag_alp21,ipart=3)
          call sum_mn_name(tmp*uxb(:,3),idiag_alp31,ipart=3)
          call sum_mn_name(tmp*uxb(:,1),idiag_alp12,ipart=3)
          call sum_mn_name(tmp*uxb(:,2),idiag_alp22,ipart=3)
          call sum_mn_name(tmp*uxb(:,3),idiag_alp32,ipart=3)
!
!  eta terms
!
          tmp=B1k1sz(n)-B1k1cz(n)
          call sum_mn_name(tmp*uxb(:,1),idiag_eta11,ipart=3)
          call sum_mn_name(tmp*uxb(:,2),idiag_eta21,ipart=3)
          call sum_mn_name(tmp*uxb(:,3),idiag_eta31,ipart=3)
          call sum_mn_name(tmp*uxb(:,1),idiag_eta12,ipart=3)
          call sum_mn_name(tmp*uxb(:,2),idiag_eta22,ipart=3)
          call sum_mn_name(tmp*uxb(:,3),idiag_eta32,ipart=3)
!
          headtt=.false.
          lfirstpoint=.false.
        enddo
        enddo
      endif
!
    endsubroutine testperturb_finalize
!***********************************************************************
    subroutine add_A0test_B11_B22 (f,jtest)
!
!  add testfields into f array:
!
!         ( 0 )   (     0      )     ( B*coskz )
!  B^11 = ( 0 ) x ( -B/k*sinkz )  =  (    0    )
!         ( dz)   (     0      )     (    0    )
!
!         ( 0 )   (     0      )     ( B*sinkz )
!  B^21 = ( 0 ) x ( +B/k*coskz )  =  (    0    )
!         ( dz)   (     0      )     (    0    )
!
!         ( 0 )   ( +B/k*sinkz )     (    0    )
!  B^12 = ( 0 ) x (      0    )   =  ( B*coskz )
!         ( dz)   (      0    )      (    0    )
!
!         ( 0 )   ( -B/k*coskz )     (    0    )
!  B^22 = ( 0 ) x (      0    )   =  ( B*sinkz )
!         ( dz)   (      0    )      (    0    )
!
!  22-mar-08/axel: adapted from testfield_z
!
      use Cdata
      use Mpicomm, only: stop_it
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: jtest
!
      intent(inout) :: f
      intent(in)  :: jtest
!
!  add testfield contribution to the f array
!
      select case (jtest)
      case (1); f(:,:,:,iay)=f(:,:,:,iay)-spread(spread(Bk1sz,1,mx),2,my)
      case (2); f(:,:,:,iay)=f(:,:,:,iay)+spread(spread(Bk1cz,1,mx),2,my)
      case (3); f(:,:,:,iax)=f(:,:,:,iax)+spread(spread(Bk1sz,1,mx),2,my)
      case (4); f(:,:,:,iax)=f(:,:,:,iax)-spread(spread(Bk1cz,1,mx),2,my)
      case default; call stop_it('add_A0test_B11_B22: jtest incorrect')
      endselect
!
    endsubroutine add_A0test_B11_B22
!***********************************************************************
    subroutine rprint_testperturb(lreset,lwrite)
!
!  reads and registers print parameters relevant to testperturb
!
!  20-mar-08/axel: adapted from shear.f90
!
      use Cdata
      use Diagnostics
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
        idiag_alp11=0; idiag_alp21=0; idiag_alp31=0
        idiag_alp12=0; idiag_alp22=0; idiag_alp32=0
        idiag_eta11=0; idiag_eta21=0; idiag_eta31=0
        idiag_eta12=0; idiag_eta22=0; idiag_eta32=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'alp11',idiag_alp11)
        call parse_name(iname,cname(iname),cform(iname),'alp21',idiag_alp21)
        call parse_name(iname,cname(iname),cform(iname),'alp31',idiag_alp31)
        call parse_name(iname,cname(iname),cform(iname),'alp12',idiag_alp12)
        call parse_name(iname,cname(iname),cform(iname),'alp22',idiag_alp22)
        call parse_name(iname,cname(iname),cform(iname),'alp32',idiag_alp32)
        call parse_name(iname,cname(iname),cform(iname),'eta11',idiag_eta11)
        call parse_name(iname,cname(iname),cform(iname),'eta21',idiag_eta21)
        call parse_name(iname,cname(iname),cform(iname),'eta31',idiag_eta31)
        call parse_name(iname,cname(iname),cform(iname),'eta12',idiag_eta12)
        call parse_name(iname,cname(iname),cform(iname),'eta22',idiag_eta22)
        call parse_name(iname,cname(iname),cform(iname),'eta32',idiag_eta32)
      enddo
!
!  write column where which testperturb variable is stored
!
      if (lwr) then
        write(3,*) 'idiag_alp11=',idiag_alp11
        write(3,*) 'idiag_alp21=',idiag_alp21
        write(3,*) 'idiag_alp31=',idiag_alp31
        write(3,*) 'idiag_alp12=',idiag_alp12
        write(3,*) 'idiag_alp22=',idiag_alp22
        write(3,*) 'idiag_alp32=',idiag_alp32
        write(3,*) 'idiag_eta11=',idiag_eta11
        write(3,*) 'idiag_eta21=',idiag_eta21
        write(3,*) 'idiag_eta31=',idiag_eta31
        write(3,*) 'idiag_eta12=',idiag_eta12
        write(3,*) 'idiag_eta22=',idiag_eta22
        write(3,*) 'idiag_eta32=',idiag_eta32
      endif
!
    endsubroutine rprint_testperturb
!***********************************************************************
endmodule TestPerturb
