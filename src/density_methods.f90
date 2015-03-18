module DensityMethods
!
!  11-mar-15/MR:  Created to avoid circular dependencies with EquationOfState.
!
  use Cparam
  use Cdata

  implicit none

  include 'density_methods.h'  

  public :: putrho, putlnrho

  interface putrho
    module procedure putrho_s
    module procedure putrho_v
  endinterface
!
  interface putlnrho
    module procedure putlnrho_s
    module procedure putlnrho_v
  endinterface

  real, dimension(:,:), pointer :: reference_state

  contains
!***********************************************************************
    subroutine initialize_density_methods

      use SharedVariables, only: get_shared_variable

      call get_shared_variable('reference_state',reference_state, &
                               caller='initialize_density_methods')

    endsubroutine initialize_density_methods
!***********************************************************************
    function getrho_s(f,lf)

      real                :: getrho_s
      real,    intent(in) :: f
      integer, intent(in) :: lf 

      if (ldensity_nolog) then
        if (lreference_state) then
          getrho_s=f+reference_state(lf-l1+1,iref_rho)
        else
          getrho_s=f
        endif
      else
        getrho_s=exp(f)
      endif

    endfunction getrho_s
!***********************************************************************
    subroutine getrho_1d(f,rho)

      real, dimension(mx), intent(in) :: f
      real, dimension(nx), intent(out):: rho

      if (ldensity_nolog) then
        if (lreference_state) then
          rho=f(l1:l2)+reference_state(:,iref_rho)
        else
          rho=f(l1:l2)
        endif
      else
        rho=exp(f(l1:l2))
      endif

    endsubroutine getrho_1d
!***********************************************************************
    subroutine getlnrho_1d(f,lnrho)

      real, dimension(mx), intent(in) :: f
      real, dimension(nx), intent(out):: lnrho

      if (ldensity_nolog) then
        if (lreference_state) then
          lnrho=log(f(l1:l2)+reference_state(:,iref_rho))
        else
          lnrho=log(f(l1:l2))
        endif
      else
        lnrho=f(l1:l2)
      endif

    endsubroutine getlnrho_1d
!***********************************************************************
    subroutine getlnrho_2dxy(f,lnrho)

      real, dimension(mx,my), intent(in) :: f
      real, dimension(mx,my), intent(out):: lnrho

      if (ldensity_nolog) then
        if (lreference_state) then
          lnrho(l1:l2,:)=log(f(l1:l2,:)+transpose(spread(reference_state(:,iref_rho),1,my)))  !!!
        else
          lnrho=log(f)
        endif
      else
        lnrho=f
      endif

    endsubroutine getlnrho_2dxy
!***********************************************************************
    subroutine getlnrho_2dyz(f,lnrho,ix)

      real, dimension(my,mz), intent(in) :: f
      real, dimension(my,mz), intent(out):: lnrho
      integer,                intent(in) :: ix

      if (ldensity_nolog) then
        if (lreference_state) then
          lnrho=log(f+reference_state(ix,iref_rho))
        else
          lnrho=log(f)
        endif
      else
        lnrho=f
      endif

    endsubroutine getlnrho_2dyz
!***********************************************************************
    subroutine getrho_2dxy(f,rho)

      real, dimension(mx,my), intent(in) :: f
      real, dimension(mx,my), intent(out):: rho

      if (ldensity_nolog) then
        if (lreference_state) then
          rho(l1:l2,:)=f(l1:l2,:)+transpose(spread(reference_state(:,iref_rho),1,my))  !!!
        else
          rho=f
        endif
      else
        rho=exp(f)
      endif

    endsubroutine getrho_2dxy
!***********************************************************************
    subroutine getrho_2dyz(f,rho,ix)

      real, dimension(my,mz), intent(in) :: f
      real, dimension(my,mz), intent(out):: rho
      integer,                intent(in) :: ix

      if (ldensity_nolog) then
        if (lreference_state) then
          rho=f+reference_state(ix,iref_rho)
        else
          rho=f
        endif
      else
        rho=exp(f)
      endif

    endsubroutine getrho_2dyz 
!***********************************************************************
    subroutine putrho_v(f,rho)

      real, dimension(mx), intent(out):: f
      real, dimension(nx), intent(in) :: rho
!
      if (ldensity_nolog) then
        f(l1:l2)=rho
        if (lreference_state) &
          f(l1:l2)=f(l1:l2)-reference_state(:,iref_rho)
      else
        f(l1:l2)=log(rho)
      endif

    endsubroutine putrho_v
!***********************************************************************
    subroutine putrho_s(f,rho)

      real, dimension(mx,my), intent(out):: f
      real,                   intent(in) :: rho
!
      integer :: m

      if (ldensity_nolog) then
        f(l1:l2,:)=rho
        if (lreference_state) then
          do m=1,my
            f(l1:l2,m)=f(l1:l2,m)-reference_state(:,iref_rho)
          enddo
        endif
      else
        f(l1:l2,:)=log(rho)
      endif

    endsubroutine putrho_s
!***********************************************************************
    subroutine putlnrho_v(f,lnrho)

      real, dimension(mx), intent(out):: f
      real, dimension(nx), intent(in) :: lnrho
!
      if (ldensity_nolog) then
        f(l1:l2)=exp(lnrho)
        if (lreference_state) &
          f(l1:l2)=f(l1:l2)-reference_state(:,iref_rho)
      else
        f(l1:l2)=lnrho
      endif

    endsubroutine putlnrho_v
!***********************************************************************
    subroutine putlnrho_s(f,lnrho)

      real, dimension(mx,my), intent(out):: f
      real,                   intent(in) :: lnrho
!
      if (ldensity_nolog) then
        f(l1:l2,:)=exp(lnrho)
        if (lreference_state) &
          f(l1:l2,:)=f(l1:l2,:)-transpose(spread(reference_state(:,iref_rho),1,my))
      else
        f(l1:l2,:)=lnrho
      endif
    endsubroutine putlnrho_s
!***********************************************************************
    subroutine getdlnrho_z(f,in,dlnrho)

      integer,                       intent(in) :: in
      real, dimension(mx,my,-in:in), intent(in) :: f
      real, dimension(mx,my),        intent(out):: dlnrho

      if (ldensity_nolog) then
        dlnrho = f(:,:,in)-f(:,:,-in)
        if (lreference_state) then
          dlnrho(l1:l2,:) = dlnrho(l1:l2,:)/(f(l1:l2,:,0)+transpose(spread(reference_state(:,iref_rho),1,my)))   !!!
        else
          dlnrho = dlnrho/f(:,:,0)
        endif
      else
        dlnrho = f(:,:,in)-f(:,:,-in)
      endif
!
    endsubroutine getdlnrho_z
!***********************************************************************
    subroutine getdlnrho_x(f,il,dx2,dlnrho,ix)

      integer,                       intent(in) :: il,ix
      real, dimension(-il:il,my,mz), intent(in) :: f
      real,                          intent(in) :: dx2
      real, dimension(my,mz),        intent(out):: dlnrho

      if (ldensity_nolog) then
        dlnrho = f(il,:,:)-f(-il,:,:)
        if (lreference_state) then
          dlnrho = (dlnrho + dx2*reference_state(ix,iref_grho)) &
                   /(f(0,:,:)+reference_state(ix,iref_rho))
        else
          dlnrho = dlnrho/f(0,:,:)
        endif
      else
        dlnrho = f(il,:,:)-f(-il,:,:)
      endif

    endsubroutine getdlnrho_x
!***********************************************************************
endmodule DensityMethods
