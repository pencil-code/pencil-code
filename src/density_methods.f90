module DensityMethods
!
!  11-mar-15/MR:  Created to avoid circular dependencies with EquationOfState.
!
  use Cdata

  implicit none

  include 'density_methods.h'  

  interface getrho1
    module procedure getrho1_1d
  endinterface
!
  interface getrho
    module procedure getrho_1d
    module procedure getrho_2dyz
    module procedure getrho_2d
    module procedure getrho_3d
  endinterface
!
  interface getlnrho
    module procedure getlnrho_1d_x
    module procedure getlnrho_1d_y
    module procedure getlnrho_2dyz
    module procedure getlnrho_2d
  endinterface
!
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

      if (lreference_state) &
        call get_shared_variable('reference_state',reference_state, &
                                 caller='initialize_density_methods')

    endsubroutine initialize_density_methods
!***********************************************************************
    subroutine getrho1_1d(f,rho1)
!
!  Fetches inverse of density from x-dependent 1D array f.
!
!   4-oct.17/MR: derived from getrho_1d.
!
      real, dimension(mx), intent(in) :: f
      real, dimension(nx), intent(out):: rho1

      if (ldensity_nolog) then
        if (lreference_state) then
          rho1=1./(f(l1:l2)+reference_state(:,iref_rho))
        else
          rho1=1./f(l1:l2)
        endif
      else
        rho1=exp(-f(l1:l2))
      endif

    endsubroutine getrho1_1d
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
    subroutine getlnrho_1d_x(f,lnrho)

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

    endsubroutine getlnrho_1d_x
!***********************************************************************
    subroutine getlnrho_1d_y(f,ix,lnrho)

      real, dimension(my), intent(in) :: f
      real, dimension(ny), intent(out):: lnrho
      integer,             intent(in) :: ix

      if (ldensity_nolog) then
        if (lreference_state) then
          lnrho=log(f(m1:m2)+reference_state(ix,iref_rho))
        else
          lnrho=log(f(m1:m2))
        endif
      else
        lnrho=f(m1:m2)
      endif

    endsubroutine getlnrho_1d_y
!***********************************************************************
    subroutine getlnrho_2d(f,lnrho)

      real, dimension(:,:), intent(in) :: f
      real, dimension(size(f,1),size(f,2)), intent(out):: lnrho

      if (ldensity_nolog) then
        if (lreference_state) then
          lnrho(l1:l2,:)= log(f(l1:l2,:) &
                         +spread(reference_state(:,iref_rho),2,size(f,2)))  !!!
        else
          lnrho=log(f)
        endif
      else
        lnrho=f
      endif

    endsubroutine getlnrho_2d
!***********************************************************************
    subroutine getlnrho_2dyz(f,ix,lnrho)

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
    subroutine getrho_2d(f,rho)

      real, dimension(:,:), intent(in) :: f
      real, dimension(size(f,1),size(f,2)), intent(out):: rho

      if (ldensity_nolog) then
        if (lreference_state) then
          rho(l1:l2,:)= f(l1:l2,:) &
                       +spread(reference_state(:,iref_rho),2,size(f,2))  !!!
        else
          rho=f
        endif
      else
        rho=exp(f)
      endif

    endsubroutine getrho_2d
!***********************************************************************
    subroutine getrho_3d(f,rho)

      real, dimension(:,:,:), intent(in) :: f
      real, dimension(size(f,1),size(f,2),size(f,3)), intent(out):: rho

      if (ldensity_nolog) then
        if (lreference_state) then
          rho(l1:l2,:,:) = f(l1:l2,:,:) &
                          +spread(spread(reference_state(:,iref_rho),2,size(f,2)),3,size(f,3))  !!!
        else
          rho=f
        endif
      else
        rho=exp(f)
      endif

    endsubroutine getrho_3d
!***********************************************************************
    subroutine getrho_2dyz(f,ix,rho)

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
          f(l1:l2,:)=f(l1:l2,:)-spread(reference_state(:,iref_rho),2,my)
      else
        f(l1:l2,:)=lnrho
      endif
    endsubroutine putlnrho_s
!***********************************************************************
    subroutine getderlnrho_z(f,iz,derlnrho)
!
!  Evaluates derlnrho as d_z ln(rho) for all x,y at z-position iz.
!
!  30-sep-16/MR: coded
!
      use Deriv, only: der

      integer,                             intent(in) :: iz
      real, dimension(:,:,:,:),            intent(in) :: f
      real, dimension(size(f,1),size(f,2)),intent(out):: derlnrho

      integer :: n_save, m_save

      n_save=n; m_save=m         ! save global n,m as we might be inside an mn-loop
      n=iz
      do m=1,size(f,2)
        call der(f(:,:,:,ilnrho),derlnrho(:,m),3)        ! = d \[ln]rho / dz
      enddo
      n=n_save; m=m_save

      if (ldensity_nolog) then
        if (lreference_state) then
          derlnrho(l1:l2,:) = derlnrho(l1:l2,:)/(f(l1:l2,:,iz,ilnrho) &
                             +spread(reference_state(:,iref_rho),2,my))   !!!
        else
          derlnrho = derlnrho/f(:,:,iz,ilnrho)
        endif
      endif
!
    endsubroutine getderlnrho_z
!***********************************************************************
    subroutine getdlnrho_x(f,rl,il,rho,dlnrho)

      integer,                   intent(in) :: rl,il
      real, dimension(mx,my,mz), intent(in) :: f
      real, dimension(my,mz),    intent(in) :: rho
      real, dimension(my,mz),    intent(out):: dlnrho

      integer :: id

      dlnrho = f(rl+il,:,:)-f(rl-il,:,:)

      if (ldensity_nolog) then
        if (lreference_state) then
          if (rl <= (nx+1)/2) then
            id = -il
          else
            id = il
          endif
          dlnrho = dlnrho + dx2_bound(id)*reference_state(rl,iref_grho) !!!
        endif
        dlnrho = dlnrho/rho
      endif

    endsubroutine getdlnrho_x
!***********************************************************************
    subroutine getdlnrho_y(f,rm,im,dlnrho)

      integer,                   intent(in) :: rm,im
      real, dimension(mx,my,mz), intent(in) :: f
      real, dimension(mx,mz),    intent(out):: dlnrho

      dlnrho = f(:,rm+im,:)-f(:,rm-im,:)

      if (ldensity_nolog) then
        if (lreference_state) then
          dlnrho(l1:l2,:) = dlnrho(l1:l2,:)/(f(l1:l2,rm,:) &
                           +spread(reference_state(:,iref_rho),2,mz))   !!!
        else
          dlnrho = dlnrho/f(:,rm,:)
        endif
      endif
!
    endsubroutine getdlnrho_y
!***********************************************************************
    subroutine getdlnrho_z(f,rn,in,dlnrho)

      integer,                              intent(in) :: rn,in
      real, dimension(:,:,:),               intent(in) :: f
      real, dimension(size(f,1),size(f,2)), intent(out):: dlnrho

      dlnrho = f(:,:,rn+in)-f(:,:,rn-in)          ! = Delta \rho or Delta log(\rho)
      if (ldensity_nolog) then
        if (lreference_state) then
          dlnrho(l1:l2,:) = dlnrho(l1:l2,:)/(f(l1:l2,:,rn) &
                           +spread(reference_state(:,iref_rho),2,my))   !!!
        else
          dlnrho = dlnrho/f(:,:,rn)
        endif
      endif
!
    endsubroutine getdlnrho_z
!***********************************************************************
endmodule DensityMethods
