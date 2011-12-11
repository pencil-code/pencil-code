! $Id$
!
!  This modules implements viscous heating and diffusion terms
!  here smagorinsky viscosity
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lviscosity = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 1
!
!***************************************************************
module Viscosity
!
  use Cparam
  use Cdata
  use Density
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  character (len=labellen) :: ivisc='smagorinsky'
  integer :: idiag_dtnu=0
  real :: maxeffectivenu,nu_mol,C_smag=0.0
  logical :: lvisc_first=.false.
!
  integer :: dummy
  namelist /viscosity_init_pars/ dummy
!
  namelist /viscosity_run_pars/ nu, lvisc_first,ivisc,c_smag
!
  contains
!***********************************************************************
    subroutine register_viscosity()
!
!  19-nov-02/tony: coded
!  16-july-04/nils: adapted from visc_shock
!
      use Cdata
      use FArrayManager
      use Sub
!
      lvisc_LES  = .true.
      lvisc_smagorinsky = .true.
!
      call farray_register_auxiliary('smagorinsky',ismagorinsky,communicated=.true.)
!
!  identify version number
!
      if (lroot) call svn_id( &
           "$Id$")
!
!  Writing files for use with IDL
!
      if (naux < maux) aux_var(aux_count)=',smagorinsky $'
      if (naux == maux) aux_var(aux_count)=',smagorinsky'
      aux_count=aux_count+1
      if (lroot) write(15,*) 'smagorinsky = fltarr(mx,my,mz)*one'
!
    endsubroutine register_viscosity
!***********************************************************************
    subroutine initialize_viscosity()
!
!  20-nov-02/tony: coded
!  16-july-04/nils: adapted from visc_shock
!
       use Cdata
!
        if (headtt.and.lroot) print*,'viscosity: nu=',nu
        if (headtt.and.lroot) print*,'viscosity: c_smag=',c_smag
!
        lneed_sij=.true.
!
    endsubroutine initialize_viscosity
!*******************************************************************
    subroutine rprint_viscosity(lreset,lwrite)
!
!  Writes ismagorinsky to index.pro file
!
!  16-july-04/nils: adapted from visc_shock
!
      use Cdata
      use Diagnostics
      use Sub
!
      logical :: lreset
      logical, optional :: lwrite
      integer :: iname
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_nu_LES=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if (lroot.and.ip<14) print*,'rprint_ionization: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'nu_LES',idiag_nu_LES)
      enddo
!
!  write column where which viscosity variable is stored
!
      if (present(lwrite)) then
        if (lwrite) then
!
        endif
      endif
!
      call keep_compiler_quiet(lreset)
!
    endsubroutine rprint_viscosity
!***********************************************************************
    subroutine pencil_criteria_viscosity()
!
!  All pencils that the Viscosity module depends on are specified here.
!
!  21-11-04/anders: coded
!
    endsubroutine pencil_criteria_viscosity
!***********************************************************************
    subroutine pencil_interdep_viscosity(lpencil_in)
!
!  Interdependency among pencils from the Viscosity module is specified here.
!
!  21-11-04/anders: coded
!
      use Cdata
!
      logical, dimension (npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_viscosity
!***********************************************************************
    subroutine calc_pencils_viscosity(f,p)
!
!  Calculate Viscosity pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  21-11-04/anders: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f,p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_viscosity
!***********************************************************************
    subroutine calc_viscosity(f)
!
!  calculate nu_smag for full domain
!
!  16-july-04/nils: coded
!
      use Diagnostics
      use Mpicomm
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,2) :: tmp
      real, dimension (nx,ny,nz) :: sij2
      real, dimension (nx) :: nu_smag
      integer :: i,j,ncount,mcount
!
!
!
      sij2=0
      do i=1,2
         do j=i+1,3
!
! find uij and uji (stored in tmp in order to save space)
!
            call der_2nd_nof(f(:,:,:,iux+i-1),tmp(:,:,:,1),j) ! uij
            call der_2nd_nof(f(:,:,:,iux+j-1),tmp(:,:,:,2),i) ! uji
!
! find off diagonal part of sij2 for full array
!
            sij2=sij2+tmp(l1:l2,m1:m2,n1:n2,1)*tmp(l1:l2,m1:m2,n1:n2,2) &
                 +0.5*tmp(l1:l2,m1:m2,n1:n2,1)**2 &
                 +0.5*tmp(l1:l2,m1:m2,n1:n2,2)**2
         enddo
      enddo
!
!  add diagonal part of sij2 for full array
!
      tmp(:,:,:,1)=0
      call der_2nd_nof(f(:,:,:,iux-1+1),tmp(:,:,:,1),1)
      sij2=sij2+tmp(l1:l2,m1:m2,n1:n2,1)**2
      call der_2nd_nof(f(:,:,:,iux-1+2),tmp(:,:,:,1),2)
      sij2=sij2+tmp(l1:l2,m1:m2,n1:n2,1)**2
      call der_2nd_nof(f(:,:,:,iux-1+3),tmp(:,:,:,1),3)
      sij2=sij2+tmp(l1:l2,m1:m2,n1:n2,1)**2
!
      f(l1:l2,m1:m2,n1:n2,ismagorinsky)=(C_smag*dxmax)**2.*sqrt(2*sij2)
!
! Diagnostics
!
      if (ldiagnos) then
         if (idiag_nu_LES /= 0) then
            itype_name(idiag_nu_LES)=ilabel_sum
            do mcount=m1,m2
               do ncount=n1,n2
                  nu_smag=f(l1:l2,mcount,ncount,ismagorinsky)
                  call sum_mn_name(nu_smag,idiag_nu_LES)
               enddo
            enddo
         endif
      endif
!
!
!
!
    endsubroutine calc_viscosity
!***********************************************************************
    subroutine der_2nd_nof(var,tmp,j)
!
!  24-nov-03/nils: coded
!
      use Cdata
!
      real, dimension (mx,my,mz) :: var
      real, dimension (mx,my,mz) :: tmp
      integer :: j
!
      intent (in) :: var,j
      intent (out) :: tmp
!
      tmp=0.
!
      if (j==1 .and. nxgrid/=1) then
          tmp(     1,:,:) = (- 3.*var(1,:,:) &
                            + 4.*var(2,:,:) &
                            - 1.*var(3,:,:))/(2.*dx)
          tmp(2:mx-1,:,:) = (- 1.*var(1:mx-2,:,:) &
                            + 1.*var(3:mx  ,:,:))/(2.*dx)
          tmp(    mx,:,:) = (+ 1.*var(mx-2,:,:) &
                            - 4.*var(mx-1,:,:) &
                            + 3.*var(mx  ,:,:))/(2.*dx)
      endif
!
      if (j==2 .and. nygrid/=1) then
          tmp(:,     1,:) = (- 3.*var(:,1,:) &
                            + 4.*var(:,2,:) &
                            - 1.*var(:,3,:))/(2.*dy)
          tmp(:,2:my-1,:) = (- 1.*var(:,1:my-2,:) &
                            + 1.*var(:,3:my  ,:))/(2.*dy)
          tmp(:,    my,:) = (+ 1.*var(:,my-2,:) &
                            - 4.*var(:,my-1,:) &
                            + 3.*var(:,my  ,:))/(2.*dy)
      endif
!
      if (j==3 .and. nzgrid/=1) then
          tmp(:,:,     1) = (- 3.*var(:,:,1) &
                            + 4.*var(:,:,2) &
                            - 1.*var(:,:,3))/(2.*dz)
          tmp(:,:,2:mz-1) = (- 1.*var(:,:,1:mz-2) &
                            + 1.*var(:,:,3:mz  ))/(2.*dz)
          tmp(:,:,    mz) = (+ 1.*var(:,:,mz-2) &
                            - 4.*var(:,:,mz-1) &
                            + 3.*var(:,:,mz  ))/(2.*dz)
      endif
!
    endsubroutine der_2nd_nof
!
!!**********************************************************************
    subroutine calc_viscous_heat(f,df,glnrho,divu,rho1,cs2,TT1,shock,Hmax)
!
!  calculate viscous heating term for right hand side of entropy equation
!
!  20-nov-02/tony: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx)   :: rho1,TT1,cs2,Hmax
      real, dimension (nx)   :: sij2, divu,shock
      real, dimension (nx,3) :: glnrho
!
!  traceless strain matrix squared
!
      call multm2_mn(sij,sij2)
!
      select case (ivisc)
        case default
          if (headtt) print*,'no heating: ivisc=',ivisc
      endselect
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(cs2)
      call keep_compiler_quiet(divu)
      call keep_compiler_quiet(glnrho)
      call keep_compiler_quiet(shock)
      call keep_compiler_quiet(Hmax)
!
    endsubroutine calc_viscous_heat
!***********************************************************************
    subroutine calc_viscous_force(f,df,glnrho,divu,rho,rho1,shock,gshock,bij)
!
!  calculate viscous heating term for right hand side of entropy equation
!
!   16-jul-04/nils: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3,3) :: bij
      real, dimension (nx,3) :: glnrho,del2u,del6u,graddivu,fvisc,sglnrho
      real, dimension (nx,3) :: gshock,sgradnu_smag
      real, dimension (nx,3) :: nusglnrho,tmp1,tmp2,gradnu_smag
      real, dimension (nx) :: murho1,rho,rho1,divu,shock,SS12,sij2
      integer :: i
!
      intent (in) :: f, glnrho, rho1,rho
      intent (out) :: df
!
!  viscosity operator
!  rho1 is pre-calculated in equ
!
      select case (ivisc)
!
        case ('smagorinsky')
          !
          !  viscous force: nu_smag*(del2u+graddivu/3+2S.glnrho)+2S.gradnu_smag
          !  where nu_smag=(C_smag*dxmax)**2*sqrt(SS)
          !
          if (headtt) print*,'viscous force: Smagorinsky'
          if (ldensity) then
             call del2v_etc(f,iuu,del2u,GRADDIV=graddivu)
             call multmv_mn(sij,glnrho,sglnrho)
             call multsv_mn(f(l1:l2,m,n,ismagorinsky),sglnrho,nusglnrho)
             tmp1=del2u+1./3.*graddivu
             call multsv_mn(f(l1:l2,m,n,ismagorinsky),tmp1,tmp2)
             call grad(f,ismagorinsky,gradnu_smag)
             call multmv_mn(sij,gradnu_smag,sgradnu_smag)
             fvisc=2*nusglnrho+tmp2+2*sgradnu_smag
             diffus_nu=diffus_nu+f(l1:l2,m,n,ismagorinsky)*dxyz_2
          else
             if (lfirstpoint) &
                  print*,"ldensity better be .true. for ivisc='smagorinsky'"
          endif
       case ('smagorinsky_simplified')
          if (headtt) print*, 'for ivisc=smagorinsky_simplified use visc_const'
          !
          !  viscous force: nu_smag*(del2u+graddivu/3+2S.glnrho)
          !  where nu_smag=(C_smag*dxmax)**2*sqrt(SS)
          !
          if (headtt) print*,'viscous force: Smagorinsky_simplified'
          if (ldensity) then
            call del2v_etc(f,iuu,del2u,GRADDIV=graddivu)
            call multmv_mn(sij,glnrho,sglnrho)
            call multsv_mn(f(l1:l2,m,n,ismagorinsky),sglnrho,nusglnrho)
            tmp1=del2u+1./3.*graddivu
            call multsv_mn(f(l1:l2,m,n,ismagorinsky),tmp1,tmp2)
            call grad(f,ismagorinsky,gradnu_smag)
            call multmv_mn(2*sij,gradnu_smag,sgradnu_smag)
            fvisc=2*nusglnrho+tmp2+sgradnu_smag
            diffus_nu=diffus_nu+f(l1:l2,m,n,ismagorinsky)*dxyz_2
         else
            if (lfirstpoint) print*,&
                 "ldensity better be true for ivisc='smagorinsky_simplified'"
         endif
      case default
         !
         !  Catch unknown values
         !
         if (lroot) print*, 'No such value for ivisc: ', trim(ivisc)
         call stop_it('calc_viscous_forcing')
!
      endselect
      !
      ! Add ordinary viscosity if nu /= 0
      !
      if (nu /= 0.) then
         !
         !  viscous force: nu*(del2u+graddivu/3+2S.glnrho)
         !  -- the correct expression for nu=const
         !
         fvisc=fvisc+2*nu*sglnrho+nu*(del2u+1./3.*graddivu)
         diffus_nu=diffus_nu+nu*dxyz_2
      endif
!
      df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)+fvisc
!
!  set viscous time step
!
      if (ldiagnos.and.idiag_dtnu/=0) then
        call max_mn_name(diffus_nu/cdtv,idiag_dtnu,l_dt=.true.)
      endif
!
      call keep_compiler_quiet(divu)
      call keep_compiler_quiet(shock)
      call keep_compiler_quiet(gshock)
!
    endsubroutine calc_viscous_force
!***********************************************************************
endmodule Viscosity
