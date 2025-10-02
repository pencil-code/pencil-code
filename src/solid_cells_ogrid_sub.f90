module solid_cells_ogrid_sub
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use solid_cells_ogrid_cdata
!
  implicit none

private

public :: grad_ogrid, grad_other_ogrid, del2_ogrid, g2ij_ogrid, dot2_mn_ogrid
public :: gij_ogrid, div_mn_ogrid, traceless_strain_ogrid, multm2_sym_mn_ogrid
public :: curl_mn_ogrid, dot_mn_ogrid, dot2_0_ogrid, u_dot_grad_vec_ogrid
public :: u_dot_grad_scl_ogrid, doupwind_ogrid, multmv_mn_ogrid
public :: gij_etc_ogrid, der_ogrid, der_other_ogrid, der2_ogrid
public :: derij_ogrid, der6_ogrid, deri_3d_inds_ogrid
public :: bval_from_neumann_arr_ogrid, bval_from_neumann_SBP
public :: bval_from_neumann_bdry5, set_ghosts_onesided_ogrid
public :: der_ogrid_SBP, der_ogrid_bdry5, der2_ogrid_bdry5
public :: der2_ogrid_bdry5_alt, der2_ogrid_SBP, der_ijm_ogrid_SBP
public :: der_ijn_ogrid_SBP, der_ijm_ogrid_bdry5, der_ijn_ogrid_bdry5
public :: der_ogrid_SBP_experimental, der2_ogrid_SBP_experimental

  contains

!***********************************************************************
    subroutine grad_ogrid(f,k,g)
!
!  Calculate gradient of a scalar, get vector.
!  f is the ogrid, k the index of the field to be differentiated      
!  and g the pencil (x direction) of the gradient of the scalar field k
!
!  07-feb-17/Jorgen: Adapted from sub.f90
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid) :: f
      real, dimension (nx_ogrid,3) :: g
      integer :: k
!
      intent(in) :: k
      intent(out) :: g
      intent(inout) :: f
!     
      call der_ogrid(f,k,g(:,1),1)
      call der_ogrid(f,k,g(:,2),2)
      call der_ogrid(f,k,g(:,3),3)
!
    endsubroutine grad_ogrid
!***********************************************************************
    subroutine grad_other_ogrid(f,g)
!
!  For non 'mvar' variable calculate gradient of a scalar, get vector
!
!  26-nov-02/tony: coded
!
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid) :: f
      real, dimension (nx_ogrid,3) :: g
!
      intent(in) :: f
      intent(out) :: g
!
!  Uses overloaded der routine.
!
      call der_other_ogrid(f,g(:,1),1); 
      call der_other_ogrid(f,g(:,2),2);
      call der_other_ogrid(f,g(:,3),3);
!
    endsubroutine grad_other_ogrid
!***********************************************************************
    subroutine del2_ogrid(f,k,del2f)
!
!  Calculate del2 of a scalar, get scalar.
!
!  07-feb-17/Jorgen: Adapted from sub.f90
!
      intent(in) :: f,k
      intent(out) :: del2f
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid) :: f
      real, dimension (nx_ogrid) :: del2f,d2fdx,d2fdy,d2fdz,tmp
      integer :: k
!
      call der2_ogrid(f,k,d2fdx,1)
      call der2_ogrid(f,k,d2fdy,2)
      call der2_ogrid(f,k,d2fdz,3)
      del2f=d2fdx+d2fdy+d2fdz
!
!  Since we have cylindrical coordinates
      call der_ogrid(f,k,tmp,1)
      del2f=del2f+tmp*rcyl_mn1_ogrid
!
    endsubroutine del2_ogrid
!***********************************************************************
    subroutine g2ij_ogrid(f,k,g)
!
!  Calculates the Hessian, i.e. all second derivatives of a scalar.
!
!  07-feb-17/Jorgen: Adapted from sub.f90
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid) :: f
      real, dimension (nx_ogrid,3,3) :: g
      real, dimension (nx_ogrid) :: tmp
      integer :: i,j,k
!
      intent(in) :: f,k
      intent(out) :: g
!
!  Run though all 9 possibilities, treat diagonals separately.
!
      do j=1,3
        call der2_ogrid(f,k,tmp,j)
        g(:,j,j)=tmp
        do i=j+1,3
          call derij_ogrid(f,k,tmp,i,j)
          g(:,i,j)=tmp
          g(:,j,i)=tmp
        enddo
      enddo
      if(SBP.or.BDRY5) then
        call fatal_error('solid_cells: g2ij_ogrid',&
          'not implemented with summation by parts property')
      endif
!
    endsubroutine g2ij_ogrid
!***********************************************************************
    subroutine dot2_mn_ogrid(a,b)
!
!  Dot product with itself, to calculate max and rms values of a vector.
!
!  07-feb-17/Jorgen: Adapted from sub.f90
!
      real, dimension (nx_ogrid,3) :: a
      real, dimension (nx_ogrid) :: b
!
      intent(in) :: a
      intent(out) :: b

      b=a(:,1)**2+a(:,2)**2+a(:,3)**2
!
    endsubroutine dot2_mn_ogrid
!***********************************************************************
    subroutine gij_ogrid(f,k,g)!,nder)
!
!  Calculate gradient of a vector, return matrix.
!
!   07-feb-17/Jorgen: Adapted from sub.f90
!   16-okt-17/Jorgen: Possibility of gradient squadred is removed, as 
!                     the expression is not tested properly, and it should
!                     most likely contain a +(1/r)df/dr term
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid) :: f
      real, dimension (nx_ogrid,3,3) :: g
      real, dimension (nx_ogrid) :: tmp
      integer :: i,j,k,k1!,nder
!
      intent(in) :: f,k!,nder
      intent(out) :: g
!
      k1=k-1
      do i=1,3 
        do j=1,3
!          if (nder == 1) then
            call der_ogrid(f,k1+i,tmp,j)
!          elseif (nder == 2) then
!            call der2_ogrid(f,k1+i,tmp,j)
!          endif
          g(:,i,j)=tmp
        enddo 
      enddo
!
    endsubroutine gij_ogrid
!***********************************************************************
    subroutine div_mn_ogrid(aij,b,a)
!
!  Calculate divergence from derivative matrix.
!
!  07-feb-17/Jorgen: Adapted from sub.f90
!
      real, dimension (nx_ogrid,3,3) :: aij
      real, dimension (nx_ogrid,3) :: a
      real, dimension (nx_ogrid) :: b
!
      intent(in) :: aij,a
      intent(out) :: b
!
      b=aij(:,1,1)+aij(:,2,2)+aij(:,3,3)
!
!  Adjustments for other cylindrical coordinate system
!
      b=b+rcyl_mn1_ogrid*a(:,1)
!
    endsubroutine div_mn_ogrid
!***********************************************************************
    subroutine traceless_strain_ogrid(uij,divu,sij,uu)
!
!  Calculates traceless rate-of-strain tensor sij from derivative tensor uij
!  and divergence divu within each pencil;
!  curvilinear co-ordinates require velocity argument uu, so this is not optional here
!  In-place operation is possible, i.e. uij and sij may refer to the same array.
!
!  07-feb-17/Jorgen: Adapted from sub.f90
!
    real, dimension (nx_ogrid,3,3)   :: uij, sij
    real, dimension (nx_ogrid)       :: divu
    real, dimension (nx_ogrid,3)     :: uu
!
    integer :: i,j
!
    intent(in)  :: uij, divu
    intent(out) :: sij
!
    do j=1,3
      sij(:,j,j)=uij(:,j,j)-(1./3.)*divu
      do i=j+1,3
        sij(:,i,j)=.5*(uij(:,i,j)+uij(:,j,i))
        sij(:,j,i)=sij(:,i,j)
      enddo
    enddo
    sij(:,1,2)=sij(:,1,2)-.5*rcyl_mn1_ogrid*uu(:,2)
    sij(:,2,2)=sij(:,2,2)+rcyl_mn1_ogrid*uu(:,1)
    sij(:,2,1)=sij(:,1,2)
!
    endsubroutine traceless_strain_ogrid
!***********************************************************************
    subroutine multm2_sym_mn_ogrid(a,b)
!
!  Symmetric matrix squared, gives scalar.
!
!  07-feb-17/Jorgen: Adapted from sub.f90
!
      real, dimension (nx_ogrid,3,3), intent(in) :: a
      real, dimension (nx_ogrid), intent(out) :: b
!
      integer :: i, j
!
      b = a(:,1,1)**2
      do i = 2, 3
        b = b + a(:,i,i)**2
        do j = 1, i-1
          b = b + 2 * a(:,i,j)**2
        enddo
      enddo
!
    endsubroutine multm2_sym_mn_ogrid
!***********************************************************************
    subroutine curl_mn_ogrid(aij,b,a)
!
!  Calculate curl from derivative matrix.
!
!  07-feb-17/Jorgen: Adapted from sub.f90
!
      real, dimension (nx_ogrid,3,3), intent (in) :: aij
      real, dimension (nx_ogrid,3),   intent (in), optional :: a
      real, dimension (nx_ogrid,3),   intent (out) :: b
      integer :: i1=1,i2=2,i3=3,i4=4,i5=5,i6=6,i7=7
!
      b(:,1)=aij(:,3,2)-aij(:,2,3)
      b(:,2)=aij(:,1,3)-aij(:,3,1)
      b(:,3)=aij(:,2,1)-aij(:,1,2)
!
!  Adjustments for cylindrical coordinate system.
!  If we go all the way to the center, we need to put a regularity condition.
!  We do this here currently up to second order, and only for curl_mn.
!
      b(:,3)=b(:,3)+a(:,2)*rcyl_mn1_ogrid
      if (rcyl_mn_ogrid(1)==0.) b(i1,3)=(360.*b(i2,3)-450.*b(i3,3)+400.*b(i4,3) &
                                  -225.*b(i5,3)+72.*b(i6,3)-10.*b(i7,3))/147.
!
    endsubroutine curl_mn_ogrid
!***********************************************************************
    subroutine dot_mn_ogrid(a,b,c)
!
!  Dot product, c=a.b, on pencil arrays
!
!  07-feb-17/Jorgen: Adapted from sub.f90
!
      real, dimension (:,:) :: a,b
      real, dimension (:) :: c
!
      intent(in) :: a,b
      intent(inout) :: c
!
      integer :: i
!
      c=0.
      do i=1,size(a,2)
        c=c+a(:,i)*b(:,i)
      enddo
!
    endsubroutine dot_mn_ogrid
!***********************************************************************
    subroutine dot2_0_ogrid(a,b)
!
!  Dot product, c=a.b, of two simple 3-d arrays.
!
!  07-fab-15/Jorgen: Copied from sub.f90
!
      real, dimension (:) :: a
      real :: b
!
      intent(in) :: a
      intent(out) :: b
!
      b = dot_product(a,a)
!
    endsubroutine dot2_0_ogrid
!***********************************************************************
    subroutine u_dot_grad_vec_ogrid(f,k,gradf,uu,ugradf,upwind)
!
!  u.gradu
!
!  07-feb-17/Jorgen: Adapted from sub.f90
!
      use General, only: loptest
!
      intent(in) :: f,k,gradf,uu,upwind
      intent(out) :: ugradf
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid) :: f
      real, dimension (nx_ogrid,3,3) :: gradf
      real, dimension (nx_ogrid,3) :: uu,ff,ugradf,grad_f_tmp
      real, dimension (nx_ogrid) :: tmp
      integer :: j,k
      logical, optional :: upwind
!
      if (k<1 .or. k>mfarray_ogrid) then
        call fatal_error('u_dot_grad_vec_ogrid','variable index is out of bounds')
        return
      endif
!
      do j=1,3
        grad_f_tmp = gradf(:,j,:)
        call u_dot_grad_scl_ogrid(f,k+j-1,grad_f_tmp,uu,tmp,UPWIND=loptest(upwind))
        ugradf(:,j)=tmp
      enddo
!
!  Adjust to cylindrical coordinate system
!  The following now works for general u.gradA.
!
      ff=f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,k:k+2)
      ugradf(:,1)=ugradf(:,1)-rcyl_mn1_ogrid*(uu(:,2)*ff(:,2))
      ugradf(:,2)=ugradf(:,2)+rcyl_mn1_ogrid*(uu(:,1)*ff(:,2))
!
    endsubroutine u_dot_grad_vec_ogrid
!***********************************************************************
    subroutine u_dot_grad_scl_ogrid(f,k,gradf,uu,ugradf,upwind)
!
!  Do advection-type term u.grad f_k.
!  Assumes gradf to be known, but takes f and k as arguments to be able
!  to calculate upwind correction.
!
!  07-feb-17/Jorgen: Adapted from sub.f90
!
      use General, only: loptest
!
      intent(in) :: f,k,gradf,uu,upwind
      intent(out) :: ugradf
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid) :: f
      real, dimension (nx_ogrid,3) :: uu,gradf
      real, dimension (nx_ogrid) :: ugradf
      integer :: k
      logical, optional :: upwind
!
      if (k<1 .or. k>mfarray_ogrid) then
        call fatal_error('u_dot_grad_scl_ogrid','variable index is out of bounds')
        return
      endif
!
      call dot_mn_ogrid(uu,gradf,ugradf)
!
!  Upwind correction
!
      if (present(upwind)) then
        if (upwind) call doupwind_ogrid(f,k,uu,ugradf)
      endif
!
    endsubroutine u_dot_grad_scl_ogrid
!***********************************************************************
    subroutine doupwind_ogrid(f,k,uu,ugradf)
!
!  Calculates upwind correction, works incrementally on ugradf
!
!  27-feb-17/Jorgen: Adapted from sub.f90
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid), intent(IN)    :: f
      integer :: k
      real, dimension (nx_ogrid,3),                         intent(IN)    :: uu
      real, dimension (nx_ogrid),                           intent(INOUT) :: ugradf
!
      real, dimension (nx_ogrid,3) :: del6f
      integer                      :: ii
!
      do ii=1,3
        call der6_ogrid(f,k,del6f(:,ii),ii)
        del6f(:,ii) = abs(uu(:,ii))*del6f(:,ii)
      enddo
!
      del6f(:,2) = rcyl_mn1_ogrid*del6f(:,2)
!
!
      ugradf = ugradf-sum(del6f,2)
!
    endsubroutine doupwind_ogrid
!***********************************************************************
    subroutine multmv_mn_ogrid(a,b,c)
!
!  Matrix multiplied with vector, gives vector.
!
!  C_i = A_{i,j} B_j
!
!  07-feb-17/Jorgen: Adapted from sub.f90
!
      real, dimension (nx_ogrid,3,3) :: a
      real, dimension (nx_ogrid,3) :: b,c
      real, dimension (nx_ogrid) :: tmp
      integer :: i,j
!
      intent(in) :: a,b
      intent(out) :: c
!
      do i=1,3
        j=1
        tmp=a(:,i,j)*b(:,j)
        do j=2,3
          tmp=tmp+a(:,i,j)*b(:,j)
        enddo
        c(:,i)=tmp
      enddo
!
    endsubroutine multmv_mn_ogrid
!***********************************************************************
    subroutine gij_etc_ogrid(f,iref,aa,aij,bij,del2,graddiv,lcovariant_derivative)
!
!  Calculate B_i,j = eps_ikl A_l,jk and A_l,kk.
!
!  05-apr-17/Jorgen - Adapted from gij_etc in sub.f90
!
      use Deriv, only: der2,derij
      use General, only: loptest
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid), intent (in) :: f
      integer, intent (in) :: iref
      logical, intent (in), optional :: lcovariant_derivative
      real, dimension (nx_ogrid,3), intent (in)   :: aa
      real, dimension (nx_ogrid,3,3), intent (in) :: aij
      real, dimension (nx_ogrid,3,3), intent (out), optional :: bij
      real, dimension (nx_ogrid,3), intent (out), optional :: del2,graddiv
!
!  Locally used variables.
!
      real, dimension (nx_ogrid,3,3,3) :: d2A
      integer :: iref1,i,j
!
!  Reference point of argument.
!
      iref1=iref-1
!
!  Calculate all mixed and non-mixed second derivatives
!  of the vector potential (A_k,ij).
!
!  Do not calculate both d^2 A_k/(dx dy) and d^2 A_k/(dy dx).
!  (This wasn't spotted by me but by a guy from SGI...)
!  Note: for non-cartesian coordinates there are different correction terms,
!  see below.
!
      do i=1,3
        do j=1,3
          call der2_ogrid(f,iref1+i,d2A(:,j,j,i),j)
        enddo
        call derij_ogrid(f,iref1+i,d2A(:,2,3,i),2,3); d2A(:,3,2,i)=d2A(:,2,3,i)
        call derij_ogrid(f,iref1+i,d2A(:,3,1,i),3,1); d2A(:,1,3,i)=d2A(:,3,1,i)
        call derij_ogrid(f,iref1+i,d2A(:,1,2,i),1,2); d2A(:,2,1,i)=d2A(:,1,2,i)
      enddo
!      d2A(:,2,1,i)=d2A(:,2,1,i)-aij(:,i,2)*rcyl_mn1_ogrid
!
!  Calculate optionally b_i,j = eps_ikl A_l,kj,
!  del2_i = A_i,jj and graddiv_i = A_j,ji .
!
      if (present(bij)) then
        call fatal_error('gij_etc_ogrid','Vorticity calculation not yet implemented on ogrid')
! !
!         bij(:,1,:)=d2A(:,2,:,3)-d2A(:,3,:,2)
!         bij(:,2,:)=d2A(:,3,:,1)-d2A(:,1,:,3)
!         bij(:,3,:)=d2A(:,1,:,2)-d2A(:,2,:,1)
! !  Corrections for cylindrical coordinates.
!           bij(:,3,2)=bij(:,3,2)+ aij(:,2,2)*rcyl_mn1_ogrid
! !          !bij(:,3,1)=bij(:,3,1)+(aij(:,2,1)+aij(:,1,2))*rcyl_mn1-aa(:,2)*rcyl_mn2
! !  FAG:Partial correction to -d2A(:,2,1,1) already applied above +aij(:,i,2)*rcyl_mn1
! !
!           bij(:,3,1)=bij(:,3,1)+aij(:,2,1)*rcyl_mn1_ogrid-aa(:,2)*rcyl_mn2_ogrid
!           if (loptest(lcovariant_derivative)) then
!             !bij(:,1,1)=bij(:,1,1)
!             bij(:,1,2)=bij(:,1,2)+(aij(:,3,1)-aij(:,1,3))*rcyl_mn1_ogrid
!             !bij(:,1,3)=bij(:,1,3)
!             !bij(:,2,1)=bij(:,2,1)
!             bij(:,2,2)=bij(:,2,2)+(aij(:,3,2)-aij(:,2,3))*rcyl_mn1_ogrid
!             !bij(:,2,3)=bij(:,2,3)
!             !bij(:,3,1)=bij(:,3,1)
!             !bij(:,3,2)=bij(:,3,2)
!             bij(:,3,3)=bij(:,3,3)+aij(:,2,3)*rcyl_mn1_ogrid
!           endif
      endif
!
!  Calculate del2 and graddiv, if requested.
!
      if (present(graddiv)) then
        graddiv(:,:)=d2A(:,1,:,1)+d2A(:,2,:,2)+d2A(:,3,:,3)
!  Since we have cylindrical coordinates
        graddiv(:,1)=graddiv(:,1)+rcyl_mn1_ogrid*(aij(:,1,1)-aij(:,2,2)) - rcyl_mn2_ogrid*aa(:,1)
        graddiv(:,2)=graddiv(:,2)+rcyl_mn1_ogrid*aij(:,1,2)
        graddiv(:,3)=graddiv(:,3)+rcyl_mn1_ogrid*aij(:,1,3)
      endif
!
      if (present(del2)) then
        del2(:,:)=d2A(:,1,1,:)+d2A(:,2,2,:)+d2A(:,3,3,:)
!  Since we have cylindrical coordinates
        del2(:,1)= del2(:,1) + rcyl_mn1_ogrid*(aij(:,1,1)-2*aij(:,2,2)) - rcyl_mn2_ogrid*aa(:,1)
        del2(:,2)= del2(:,2) + rcyl_mn1_ogrid*(aij(:,2,1)+2*aij(:,1,2)) - rcyl_mn2_ogrid*aa(:,2)
      endif
!
    endsubroutine gij_etc_ogrid
!***********************************************************************
!***********************************************************************
!***********************************************************************
! FROM DERIV.F90
!***********************************************************************
! ROUTINES
!   der_ogrid
!   der2_ogrid
!   derij_ogrid
!   der6_ogrid
!   deri_3d_inds
!***********************************************************************
    subroutine der_ogrid(f, k, df, j)
!
!  calculate derivative df_k/dx_j 
!  accurate to 6th order, explicit, periodic
!  replace cshifts by explicit construction -> x6.5 faster!
!
!  07-feb-17/Jorgen: Adapted from deriv.f90
!
      real, dimension(mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid), intent(in) :: f
      real, dimension(nx_ogrid), intent(out) :: df
      integer, intent(in) :: j, k
!
      integer :: i
      real, parameter :: a = 1.0 / 60.0
      real, dimension(nx_ogrid) :: fac
!
      if (j==1) then
        if (nxgrid_ogrid/=1) then
          if(lfirst_proc_x) then
            if(SBP) then
              call der_ogrid_SBP(f(l1_ogrid:l1_ogrid+8,m_ogrid,n_ogrid,k),df(1:6))
              if (.not. SBP_reduced_to2) then
                i=6
              else
                i=4
              endif
            elseif(BDRY5) then
              call der_ogrid_bdry5(f(l1_ogrid:l1_ogrid+5,m_ogrid,n_ogrid,k),df(1:3))
              i=3
            else
              i=0
            endif
          else
            i=0
          endif
          fac = a * dx_1_ogrid(l1_ogrid:l2_ogrid)
          df(1+i:nx_ogrid)=fac(1+i:nx_ogrid) * &
                          (+ 45.0*(f(l1_ogrid+1+i:l2_ogrid+1,m_ogrid,n_ogrid,k)-f(l1_ogrid-1+i:l2_ogrid-1,m_ogrid,n_ogrid,k)) &
                           -  9.0*(f(l1_ogrid+2+i:l2_ogrid+2,m_ogrid,n_ogrid,k)-f(l1_ogrid-2+i:l2_ogrid-2,m_ogrid,n_ogrid,k)) &
                           +      (f(l1_ogrid+3+i:l2_ogrid+3,m_ogrid,n_ogrid,k)-f(l1_ogrid-3+i:l2_ogrid-3,m_ogrid,n_ogrid,k)))
        else
          df=0.
          if (ip<=5) print*, 'der_main: Degenerate case in x-direction'
        endif
      elseif (j==2) then
        if (nygrid_ogrid/=1) then
          fac = a * dy_1_ogrid(m_ogrid)
          df=fac*(+ 45.0*(f(l1_ogrid:l2_ogrid,m_ogrid+1,n_ogrid,k)-f(l1_ogrid:l2_ogrid,m_ogrid-1,n_ogrid,k)) &
                  -  9.0*(f(l1_ogrid:l2_ogrid,m_ogrid+2,n_ogrid,k)-f(l1_ogrid:l2_ogrid,m_ogrid-2,n_ogrid,k)) &
                  +      (f(l1_ogrid:l2_ogrid,m_ogrid+3,n_ogrid,k)-f(l1_ogrid:l2_ogrid,m_ogrid-3,n_ogrid,k)))
          ! Since we have cylindrical coordinates
          df = df * rcyl_mn1_ogrid
        else
          df=0.
          if (ip<=5) print*, 'der_main: Degenerate case in y-direction'
        endif
      elseif (j==3) then
        if (nzgrid_ogrid/=1) then
          fac = a * dz_1_ogrid(n_ogrid)
          df=fac*(+ 45.0*(f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid+1,k)-f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid-1,k)) &
                  -  9.0*(f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid+2,k)-f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid-2,k)) &
                  +      (f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid+3,k)-f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid-3,k)))
        else
          df=0.
          if (ip<=5) print*, 'der_main: Degenerate case in z-direction'
        endif
      endif
      !if(m_ogrid == 16 .and. j==1) then
      !  print*, ''
      !  print*, 'df,k'
      !  print*, df,k
      !endif
!
    endsubroutine der_ogrid
!***********************************************************************
    subroutine der_other_ogrid(f,df,j)
!
!  Along one pencil in NON f variable
!  calculate derivative of a scalar, get scalar
!  accurate to 6th order, explicit, periodic
!  replace cshifts by explicit construction -> x6.5 faster!
!
!  26-nov-02/tony: coded - duplicate der_main but without k subscript
!                          then overload the der interface.
!  25-jun-04/tobi+wolf: adapted for non-equidistant grids
!  21-feb-07/axel: added 1/r and 1/pomega factors for non-coord basis
!  30-sep-16/MR: allow results dimensions > nx
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid) :: f
      real, dimension (nx_ogrid) :: df
      integer :: j,i
!
      intent(in)  :: f,j
      intent(out) :: df

      real, dimension (nx_ogrid) :: fac
!
      if (j==1) then
        if (nxgrid_ogrid/=1) then
          if(lfirst_proc_x) then
            if(SBP) then
              call der_ogrid_SBP(f(l1_ogrid:l1_ogrid+8,m_ogrid,n_ogrid),df(1:6))
              if (.not. SBP_reduced_to2) then
                i=6
              else
                i=4
              endif
            elseif(BDRY5) then
              call der_ogrid_bdry5(f(l1_ogrid:l1_ogrid+5,m_ogrid,n_ogrid),df(1:3))
              i=3
            else
              i=0
            endif
          else
            i=0
          endif
          fac=(1./60)*dx_1_ogrid(l1_ogrid:l2_ogrid)
          df=fac*(+ 45.0*(f(l1_ogrid+1:l2_ogrid+1,m_ogrid,n_ogrid)-f(l1_ogrid-1:l2_ogrid-1,m_ogrid,n_ogrid)) &
                  -  9.0*(f(l1_ogrid+2:l2_ogrid+2,m_ogrid,n_ogrid)-f(l1_ogrid-2:l2_ogrid-2,m_ogrid,n_ogrid)) &
                  +      (f(l1_ogrid+3:l2_ogrid+3,m_ogrid,n_ogrid)-f(l1_ogrid-3:l2_ogrid-3,m_ogrid,n_ogrid)))
        else
          df=0.
          if (ip<=5) print*, 'der_other: Degenerate case in x-direction'
        endif
      else
         if (j==2) then
          if (nygrid_ogrid/=1) then
            fac=(1./60)*dy_1_ogrid(m_ogrid)
            df=fac*(+ 45.0*(f(l1_ogrid:l2_ogrid,m_ogrid+1,n_ogrid)-f(l1_ogrid:l2_ogrid,m_ogrid-1,n_ogrid)) &
                    -  9.0*(f(l1_ogrid:l2_ogrid,m_ogrid+2,n_ogrid)-f(l1_ogrid:l2_ogrid,m_ogrid-2,n_ogrid)) &
                    +      (f(l1_ogrid:l2_ogrid,m_ogrid+3,n_ogrid)-f(l1_ogrid:l2_ogrid,m_ogrid-3,n_ogrid)))
          ! Since we have cylindrical coordinates
            df=df*rcyl_mn1_ogrid
          else
            df=0.
            if (ip<=5) print*, 'der_other: Degenerate case in y-direction'
          endif
        elseif (j==3) then
          if (nzgrid_ogrid/=1) then
            fac=(1./60)*dz_1_ogrid(n_ogrid)
            df=fac*(+ 45.0*(f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid+1)-f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid-1)) &
                    -  9.0*(f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid+2)-f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid-2)) &
                    +      (f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid+3)-f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid-3)))
          else
            df=0.
            if (ip<=5) print*, 'der_other: Degenerate case in z-direction'
          endif
        endif
      endif
!
    endsubroutine der_other_ogrid
!***********************************************************************
    subroutine der2_ogrid(f,k,df2,j)
!
!  calculate 2nd derivative d^2f_k/dx_j^2
!  accurate to 6th order, explicit, periodic
!  replace cshifts by explicit construction -> x6.5 faster!
!
!  07-feb-17/Jorgen: Adapted from deriv.f90
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid) :: f
      real, dimension (nx_ogrid) :: df2,fac,df
      integer :: j,k,i
!
      real :: der2_coef0, der2_coef1, der2_coef2, der2_coef3

      intent(in)  :: f,k,j
      intent(out) :: df2
!
      der2_coef0=-490./180.; der2_coef1=270./180.
      der2_coef2=-27./180.; der2_coef3=2./180.

      if (j==1) then
        if (nxgrid_ogrid/=1) then
          if(lfirst_proc_x) then
            if(SBP) then
              call der2_ogrid_SBP(f(l1_ogrid:l1_ogrid+8,m_ogrid,n_ogrid,k),df2(1:6))
              if (.not. SBP_reduced_to2) then
                i=6
              else
                i=4
              endif
            elseif(BDRY5) then
              call der2_ogrid_bdry5(f(l1_ogrid:l1_ogrid+6,m_ogrid,n_ogrid,k),df2(1:3))
              !call der2_ogrid_bdry5_alt(f,df2(1:3),k)
              i=3
            else
              i=0
            endif
          else
            i=0
          endif
          fac=dx_1_ogrid(l1_ogrid:l2_ogrid)**2
          df2(1+i:nx_ogrid)=fac(1+i:nx_ogrid) * &
                  (der2_coef0* f(l1_ogrid  +i:l2_ogrid  ,m_ogrid,n_ogrid,k) &
                  +der2_coef1*(f(l1_ogrid+1+i:l2_ogrid+1,m_ogrid,n_ogrid,k)+f(l1_ogrid-1+i:l2_ogrid-1,m_ogrid,n_ogrid,k)) &
                  +der2_coef2*(f(l1_ogrid+2+i:l2_ogrid+2,m_ogrid,n_ogrid,k)+f(l1_ogrid-2+i:l2_ogrid-2,m_ogrid,n_ogrid,k)) &
                  +der2_coef3*(f(l1_ogrid+3+i:l2_ogrid+3,m_ogrid,n_ogrid,k)+f(l1_ogrid-3+i:l2_ogrid-3,m_ogrid,n_ogrid,k)))
          if (.not.lequidist_ogrid(j)) then
            call der_ogrid(f,k,df,j)
            df2=df2+dx_tilde_ogrid(l1_ogrid:l2_ogrid)*df
          endif
        else
          df2=0.
        endif
      elseif (j==2) then
        if (nygrid_ogrid/=1) then
          fac=dy_1_ogrid(m_ogrid)**2
          df2=fac*(der2_coef0* f(l1_ogrid:l2_ogrid,m_ogrid  ,n_ogrid,k) &
                  +der2_coef1*(f(l1_ogrid:l2_ogrid,m_ogrid+1,n_ogrid,k)+f(l1_ogrid:l2_ogrid,m_ogrid-1,n_ogrid,k)) &
                  +der2_coef2*(f(l1_ogrid:l2_ogrid,m_ogrid+2,n_ogrid,k)+f(l1_ogrid:l2_ogrid,m_ogrid-2,n_ogrid,k)) &
                  +der2_coef3*(f(l1_ogrid:l2_ogrid,m_ogrid+3,n_ogrid,k)+f(l1_ogrid:l2_ogrid,m_ogrid-3,n_ogrid,k)))
          df2=df2*rcyl_mn2_ogrid
          if (.not.lequidist_ogrid(j)) then
            call der_ogrid(f,k,df,j)
            df2=df2+dy_tilde_ogrid(m_ogrid)*df
          endif
        else
          df2=0.
        endif
      elseif (j==3) then
        if (nzgrid_ogrid/=1) then
          fac=dz_1_ogrid(n_ogrid)**2
          df2=fac*( der2_coef0* f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid    ,k) &
                   +der2_coef1*(f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid+1,k)+f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid-1,k)) &
                   +der2_coef2*(f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid+2,k)+f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid-2,k)) &
                   +der2_coef3*(f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid+3,k)+f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid-3,k)))
          if (.not.lequidist_ogrid(j)) then
            call der_ogrid(f,k,df,j)
            df2=df2+dz_tilde_ogrid(n_ogrid)*df
          endif
        else
          df2=0.
        endif
      endif
!
    endsubroutine der2_ogrid
!***********************************************************************
    subroutine derij_ogrid(f,k,df,i,j)
!
!  calculate 2nd derivative with respect to two different directions
!  input: scalar, output: scalar
!  accurate to 6th order, explicit, periodic
!
!  17-feb-17/Jorgen: Adapted from deriv.f90
!  05-apr-17/Jorgen: Added summation by parts operator near cylinder surface
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid) :: f
      real, dimension (nx_ogrid) :: df,fac
      integer :: i,j,k,ii
!
      intent(in) :: f,k,i,j
      intent(out) :: df
!
      if (lbidiagonal_derij_ogrid) then
        !
        ! Use bidiagonal mixed-derivative operator, i.e.
        ! employ only the three neighbouring points on each of the four
        ! half-diagonals. This gives 6th-order mixed derivatives as the
        ! version below, but involves just 12 points instead of 36.
        !
        if ((i==1.and.j==2).or.(i==2.and.j==1)) then
          if (nxgrid_ogrid/=1.and.nygrid_ogrid/=1) then
            fac=(1./720.)*dx_1_ogrid(l1_ogrid:l2_ogrid)*dy_1_ogrid(m_ogrid)
            df=fac*( &
                        270.*( f(l1_ogrid+1:l2_ogrid+1,m_ogrid+1,n_ogrid,k)-f(l1_ogrid-1:l2_ogrid-1,m_ogrid+1,n_ogrid,k)  &
                              +f(l1_ogrid-1:l2_ogrid-1,m_ogrid-1,n_ogrid,k)-f(l1_ogrid+1:l2_ogrid+1,m_ogrid-1,n_ogrid,k)) &
                       - 27.*( f(l1_ogrid+2:l2_ogrid+2,m_ogrid+2,n_ogrid,k)-f(l1_ogrid-2:l2_ogrid-2,m_ogrid+2,n_ogrid,k)  &
                              +f(l1_ogrid-2:l2_ogrid-2,m_ogrid-2,n_ogrid,k)-f(l1_ogrid+2:l2_ogrid+2,m_ogrid-2,n_ogrid,k)) &
                       +  2.*( f(l1_ogrid+3:l2_ogrid+3,m_ogrid+3,n_ogrid,k)-f(l1_ogrid-3:l2_ogrid-3,m_ogrid+3,n_ogrid,k)  &
                              +f(l1_ogrid-3:l2_ogrid-3,m_ogrid-3,n_ogrid,k)-f(l1_ogrid+3:l2_ogrid+3,m_ogrid-3,n_ogrid,k)) &
                   )
          else
            df=0.
            if (ip<=5) print*, 'derij: Degenerate case in x- or y-direction'
          endif
        elseif ((i==2.and.j==3).or.(i==3.and.j==2)) then
          if (nygrid_ogrid/=1.and.nzgrid_ogrid/=1) then
            fac=(1./720.)*dy_1_ogrid(m_ogrid)*dz_1_ogrid(n_ogrid)
            df=fac*( &
                        270.*( f(l1_ogrid:l2_ogrid,m_ogrid+1,n_ogrid+1,k)-f(l1_ogrid:l2_ogrid,m_ogrid+1,n_ogrid-1,k)  &
                              +f(l1_ogrid:l2_ogrid,m_ogrid-1,n_ogrid-1,k)-f(l1_ogrid:l2_ogrid,m_ogrid-1,n_ogrid+1,k)) &
                       - 27.*( f(l1_ogrid:l2_ogrid,m_ogrid+2,n_ogrid+2,k)-f(l1_ogrid:l2_ogrid,m_ogrid+2,n_ogrid-2,k)  &
                              +f(l1_ogrid:l2_ogrid,m_ogrid-2,n_ogrid-2,k)-f(l1_ogrid:l2_ogrid,m_ogrid-2,n_ogrid+2,k)) &
                       +  2.*( f(l1_ogrid:l2_ogrid,m_ogrid+3,n_ogrid+3,k)-f(l1_ogrid:l2_ogrid,m_ogrid+3,n_ogrid-3,k)  &
                              +f(l1_ogrid:l2_ogrid,m_ogrid-3,n_ogrid-3,k)-f(l1_ogrid:l2_ogrid,m_ogrid-3,n_ogrid+3,k)) &
                   )
          else
            df=0.
            if (ip<=5) print*, 'derij: Degenerate case in y- or z-direction'
          endif
        elseif ((i==3.and.j==1).or.(i==1.and.j==3)) then
          if (nzgrid_ogrid/=1.and.nxgrid_ogrid/=1) then
            fac=(1./720.)*dz_1_ogrid(n_ogrid)*dx_1_ogrid(l1_ogrid:l2_ogrid)
            df=fac*( &
                        270.*( f(l1_ogrid+1:l2_ogrid+1,m_ogrid,n_ogrid+1,k)-f(l1_ogrid-1:l2_ogrid-1,m_ogrid,n_ogrid+1,k)  &
                              +f(l1_ogrid-1:l2_ogrid-1,m_ogrid,n_ogrid-1,k)-f(l1_ogrid+1:l2_ogrid+1,m_ogrid,n_ogrid-1,k)) &
                       - 27.*( f(l1_ogrid+2:l2_ogrid+2,m_ogrid,n_ogrid+2,k)-f(l1_ogrid-2:l2_ogrid-2,m_ogrid,n_ogrid+2,k)  &
                              +f(l1_ogrid-2:l2_ogrid-2,m_ogrid,n_ogrid-2,k)-f(l1_ogrid+2:l2_ogrid+2,m_ogrid,n_ogrid-2,k)) &
                       +  2.*( f(l1_ogrid+3:l2_ogrid+3,m_ogrid,n_ogrid+3,k)-f(l1_ogrid-3:l2_ogrid-3,m_ogrid,n_ogrid+3,k)  &
                              +f(l1_ogrid-3:l2_ogrid-3,m_ogrid,n_ogrid-3,k)-f(l1_ogrid+3:l2_ogrid+3,m_ogrid,n_ogrid-3,k)) &
                   )
          else
            df=0.
            if (ip<=5) print*, 'derij: Degenerate case in x- or z-direction'
          endif
        endif
!
      else                      ! not using bidiagonal mixed derivatives
        !
        ! This is the old, straight-forward scheme
        ! Note that the summation by parts operators are only implemented for this scheme
        !
        if ((i==1.and.j==2).or.(i==2.and.j==1)) then
          if (nxgrid_ogrid/=1.and.nygrid_ogrid/=1) then
            if(lfirst_proc_x) then
              if(SBP) then
                if (.not. SBP_reduced_to2) then
                  ii=6
                else
                  ii=4
                endif
                call der_ijm_ogrid_SBP(f,df(1:6),k)
              elseif(BDRY5) then
                ii=3
                call der_ijm_ogrid_bdry5(f,df(1:3),k)
              else
                ii=0
              endif
            else
              ii=0
            endif
            fac=(1./60.**2)*dx_1_ogrid(l1_ogrid:l2_ogrid)*dy_1_ogrid(m_ogrid)
            df(1+ii:nx_ogrid)=fac(1+ii:nx_ogrid)*( &
              45.*((45.*(f(l1_ogrid+1+ii:l2_ogrid+1,m_ogrid+1,n_ogrid,k)-f(l1_ogrid-1+ii:l2_ogrid-1,m_ogrid+1,n_ogrid,k))  &
                    -9.*(f(l1_ogrid+2+ii:l2_ogrid+2,m_ogrid+1,n_ogrid,k)-f(l1_ogrid-2+ii:l2_ogrid-2,m_ogrid+1,n_ogrid,k))  &
                       +(f(l1_ogrid+3+ii:l2_ogrid+3,m_ogrid+1,n_ogrid,k)-f(l1_ogrid-3+ii:l2_ogrid-3,m_ogrid+1,n_ogrid,k))) &
                  -(45.*(f(l1_ogrid+1+ii:l2_ogrid+1,m_ogrid-1,n_ogrid,k)-f(l1_ogrid-1+ii:l2_ogrid-1,m_ogrid-1,n_ogrid,k))  &
                    -9.*(f(l1_ogrid+2+ii:l2_ogrid+2,m_ogrid-1,n_ogrid,k)-f(l1_ogrid-2+ii:l2_ogrid-2,m_ogrid-1,n_ogrid,k))  &
                       +(f(l1_ogrid+3+ii:l2_ogrid+3,m_ogrid-1,n_ogrid,k)-f(l1_ogrid-3+ii:l2_ogrid-3,m_ogrid-1,n_ogrid,k))))&
              -9.*((45.*(f(l1_ogrid+1+ii:l2_ogrid+1,m_ogrid+2,n_ogrid,k)-f(l1_ogrid-1+ii:l2_ogrid-1,m_ogrid+2,n_ogrid,k))  &
                    -9.*(f(l1_ogrid+2+ii:l2_ogrid+2,m_ogrid+2,n_ogrid,k)-f(l1_ogrid-2+ii:l2_ogrid-2,m_ogrid+2,n_ogrid,k))  &
                       +(f(l1_ogrid+3+ii:l2_ogrid+3,m_ogrid+2,n_ogrid,k)-f(l1_ogrid-3+ii:l2_ogrid-3,m_ogrid+2,n_ogrid,k))) &
                  -(45.*(f(l1_ogrid+1+ii:l2_ogrid+1,m_ogrid-2,n_ogrid,k)-f(l1_ogrid-1+ii:l2_ogrid-1,m_ogrid-2,n_ogrid,k))  &
                    -9.*(f(l1_ogrid+2+ii:l2_ogrid+2,m_ogrid-2,n_ogrid,k)-f(l1_ogrid-2+ii:l2_ogrid-2,m_ogrid-2,n_ogrid,k))  &
                       +(f(l1_ogrid+3+ii:l2_ogrid+3,m_ogrid-2,n_ogrid,k)-f(l1_ogrid-3+ii:l2_ogrid-3,m_ogrid-2,n_ogrid,k))))&
                 +((45.*(f(l1_ogrid+1+ii:l2_ogrid+1,m_ogrid+3,n_ogrid,k)-f(l1_ogrid-1+ii:l2_ogrid-1,m_ogrid+3,n_ogrid,k))  &
                    -9.*(f(l1_ogrid+2+ii:l2_ogrid+2,m_ogrid+3,n_ogrid,k)-f(l1_ogrid-2+ii:l2_ogrid-2,m_ogrid+3,n_ogrid,k))  &
                       +(f(l1_ogrid+3+ii:l2_ogrid+3,m_ogrid+3,n_ogrid,k)-f(l1_ogrid-3+ii:l2_ogrid-3,m_ogrid+3,n_ogrid,k))) &
                  -(45.*(f(l1_ogrid+1+ii:l2_ogrid+1,m_ogrid-3,n_ogrid,k)-f(l1_ogrid-1+ii:l2_ogrid-1,m_ogrid-3,n_ogrid,k))  &
                    -9.*(f(l1_ogrid+2+ii:l2_ogrid+2,m_ogrid-3,n_ogrid,k)-f(l1_ogrid-2+ii:l2_ogrid-2,m_ogrid-3,n_ogrid,k))  &
                       +(f(l1_ogrid+3+ii:l2_ogrid+3,m_ogrid-3,n_ogrid,k)-f(l1_ogrid-3+ii:l2_ogrid-3,m_ogrid-3,n_ogrid,k))))&
                   )
          else
            df=0.
            if (ip<=5) print*, 'derij: Degenerate case in x- or y-direction'
          endif
        elseif ((i==2.and.j==3).or.(i==3.and.j==2)) then
          if (nygrid_ogrid/=1.and.nzgrid_ogrid/=1) then
            fac=(1./60.**2)*dy_1_ogrid(m_ogrid)*dz_1_ogrid(n_ogrid)
            df=fac*( &
              45.*((45.*(f(l1_ogrid:l2_ogrid,m_ogrid+1,n_ogrid+1,k)-f(l1_ogrid:l2_ogrid,m_ogrid-1,n_ogrid+1,k))  &
                    -9.*(f(l1_ogrid:l2_ogrid,m_ogrid+2,n_ogrid+1,k)-f(l1_ogrid:l2_ogrid,m_ogrid-2,n_ogrid+1,k))  &
                       +(f(l1_ogrid:l2_ogrid,m_ogrid+3,n_ogrid+1,k)-f(l1_ogrid:l2_ogrid,m_ogrid-3,n_ogrid+1,k))) &
                  -(45.*(f(l1_ogrid:l2_ogrid,m_ogrid+1,n_ogrid-1,k)-f(l1_ogrid:l2_ogrid,m_ogrid-1,n_ogrid-1,k))  &
                    -9.*(f(l1_ogrid:l2_ogrid,m_ogrid+2,n_ogrid-1,k)-f(l1_ogrid:l2_ogrid,m_ogrid-2,n_ogrid-1,k))  &
                       +(f(l1_ogrid:l2_ogrid,m_ogrid+3,n_ogrid-1,k)-f(l1_ogrid:l2_ogrid,m_ogrid-3,n_ogrid-1,k))))&
              -9.*((45.*(f(l1_ogrid:l2_ogrid,m_ogrid+1,n_ogrid+2,k)-f(l1_ogrid:l2_ogrid,m_ogrid-1,n_ogrid+2,k))  &
                    -9.*(f(l1_ogrid:l2_ogrid,m_ogrid+2,n_ogrid+2,k)-f(l1_ogrid:l2_ogrid,m_ogrid-2,n_ogrid+2,k))  &
                       +(f(l1_ogrid:l2_ogrid,m_ogrid+3,n_ogrid+2,k)-f(l1_ogrid:l2_ogrid,m_ogrid-3,n_ogrid+2,k))) &
                  -(45.*(f(l1_ogrid:l2_ogrid,m_ogrid+1,n_ogrid-2,k)-f(l1_ogrid:l2_ogrid,m_ogrid-1,n_ogrid-2,k))  &
                    -9.*(f(l1_ogrid:l2_ogrid,m_ogrid+2,n_ogrid-2,k)-f(l1_ogrid:l2_ogrid,m_ogrid-2,n_ogrid-2,k))  &
                       +(f(l1_ogrid:l2_ogrid,m_ogrid+3,n_ogrid-2,k)-f(l1_ogrid:l2_ogrid,m_ogrid-3,n_ogrid-2,k))))&
                 +((45.*(f(l1_ogrid:l2_ogrid,m_ogrid+1,n_ogrid+3,k)-f(l1_ogrid:l2_ogrid,m_ogrid-1,n_ogrid+3,k))  &
                    -9.*(f(l1_ogrid:l2_ogrid,m_ogrid+2,n_ogrid+3,k)-f(l1_ogrid:l2_ogrid,m_ogrid-2,n_ogrid+3,k))  &
                       +(f(l1_ogrid:l2_ogrid,m_ogrid+3,n_ogrid+3,k)-f(l1_ogrid:l2_ogrid,m_ogrid-3,n_ogrid+3,k))) &
                  -(45.*(f(l1_ogrid:l2_ogrid,m_ogrid+1,n_ogrid-3,k)-f(l1_ogrid:l2_ogrid,m_ogrid-1,n_ogrid-3,k))  &
                    -9.*(f(l1_ogrid:l2_ogrid,m_ogrid+2,n_ogrid-3,k)-f(l1_ogrid:l2_ogrid,m_ogrid-2,n_ogrid-3,k))  &
                       +(f(l1_ogrid:l2_ogrid,m_ogrid+3,n_ogrid-3,k)-f(l1_ogrid:l2_ogrid,m_ogrid-3,n_ogrid-3,k))))&
                   )
          else
            df=0.
            if (ip<=5) print*, 'derij: Degenerate case in y- or z-direction'
          endif
        elseif ((i==3.and.j==1).or.(i==1.and.j==3)) then
          if (nzgrid_ogrid/=1.and.nxgrid_ogrid/=1) then
            if(lfirst_proc_x) then
              if(SBP) then
                if (.not. SBP_reduced_to2) then
                  ii=6
                else
                  ii=4
                endif
                call der_ijn_ogrid_SBP(f,df(1:6),k)
              elseif(BDRY5) then
                ii=3
                call der_ijn_ogrid_bdry5(f,df(1:3),k)
              else
                ii=0
              endif
            else
              ii=0
            endif
            fac=(1./60.**2)*dz_1_ogrid(n_ogrid)*dx_1_ogrid(l1_ogrid:l2_ogrid)
            df(1+ii:nx_ogrid)=fac(1+ii:nx_ogrid)*( &
              45.*((45.*(f(l1_ogrid+1+ii:l2_ogrid+1,m_ogrid,n_ogrid+1,k)-f(l1_ogrid-1+ii:l2_ogrid-1,m_ogrid,n_ogrid+1,k))  &
                    -9.*(f(l1_ogrid+2+ii:l2_ogrid+2,m_ogrid,n_ogrid+1,k)-f(l1_ogrid-2+ii:l2_ogrid-2,m_ogrid,n_ogrid+1,k))  &
                       +(f(l1_ogrid+3+ii:l2_ogrid+3,m_ogrid,n_ogrid+1,k)-f(l1_ogrid-3+ii:l2_ogrid-3,m_ogrid,n_ogrid+1,k))) &
                  -(45.*(f(l1_ogrid+1+ii:l2_ogrid+1,m_ogrid,n_ogrid-1,k)-f(l1_ogrid-1+ii:l2_ogrid-1,m_ogrid,n_ogrid-1,k))  &
                    -9.*(f(l1_ogrid+2+ii:l2_ogrid+2,m_ogrid,n_ogrid-1,k)-f(l1_ogrid-2+ii:l2_ogrid-2,m_ogrid,n_ogrid-1,k))  &
                       +(f(l1_ogrid+3+ii:l2_ogrid+3,m_ogrid,n_ogrid-1,k)-f(l1_ogrid-3+ii:l2_ogrid-3,m_ogrid,n_ogrid-1,k))))&
              -9.*((45.*(f(l1_ogrid+1+ii:l2_ogrid+1,m_ogrid,n_ogrid+2,k)-f(l1_ogrid-1+ii:l2_ogrid-1,m_ogrid,n_ogrid+2,k))  &
                    -9.*(f(l1_ogrid+2+ii:l2_ogrid+2,m_ogrid,n_ogrid+2,k)-f(l1_ogrid-2+ii:l2_ogrid-2,m_ogrid,n_ogrid+2,k))  &
                       +(f(l1_ogrid+3+ii:l2_ogrid+3,m_ogrid,n_ogrid+2,k)-f(l1_ogrid-3+ii:l2_ogrid-3,m_ogrid,n_ogrid+2,k))) &
                  -(45.*(f(l1_ogrid+1+ii:l2_ogrid+1,m_ogrid,n_ogrid-2,k)-f(l1_ogrid-1+ii:l2_ogrid-1,m_ogrid,n_ogrid-2,k))  &
                    -9.*(f(l1_ogrid+2+ii:l2_ogrid+2,m_ogrid,n_ogrid-2,k)-f(l1_ogrid-2+ii:l2_ogrid-2,m_ogrid,n_ogrid-2,k))  &
                       +(f(l1_ogrid+3+ii:l2_ogrid+3,m_ogrid,n_ogrid-2,k)-f(l1_ogrid-3+ii:l2_ogrid-3,m_ogrid,n_ogrid-2,k))))&
                 +((45.*(f(l1_ogrid+1+ii:l2_ogrid+1,m_ogrid,n_ogrid+3,k)-f(l1_ogrid-1+ii:l2_ogrid-1,m_ogrid,n_ogrid+3,k))  &
                    -9.*(f(l1_ogrid+2+ii:l2_ogrid+2,m_ogrid,n_ogrid+3,k)-f(l1_ogrid-2+ii:l2_ogrid-2,m_ogrid,n_ogrid+3,k))  &
                       +(f(l1_ogrid+3+ii:l2_ogrid+3,m_ogrid,n_ogrid+3,k)-f(l1_ogrid-3+ii:l2_ogrid-3,m_ogrid,n_ogrid+3,k))) &
                  -(45.*(f(l1_ogrid+1+ii:l2_ogrid+1,m_ogrid,n_ogrid-3,k)-f(l1_ogrid-1+ii:l2_ogrid-1,m_ogrid,n_ogrid-3,k))  &
                    -9.*(f(l1_ogrid+2+ii:l2_ogrid+2,m_ogrid,n_ogrid-3,k)-f(l1_ogrid-2+ii:l2_ogrid-2,m_ogrid,n_ogrid-3,k))  &
                       +(f(l1_ogrid+3+ii:l2_ogrid+3,m_ogrid,n_ogrid-3,k)-f(l1_ogrid-3+ii:l2_ogrid-3,m_ogrid,n_ogrid-3,k))))&
                   )
            !fac=(1./60.**2)*dz_1_ogrid(n_ogrid)*dx_1_ogrid(l1_ogrid:l2_ogrid)
            !df=fac*( &
            !  45.*((45.*(f(l1_ogrid+1:l2_ogrid+1,m_ogrid,n_ogrid+1,k)-f(l1_ogrid-1:l2_ogrid-1,m_ogrid,n_ogrid+1,k))  &
            !        -9.*(f(l1_ogrid+2:l2_ogrid+2,m_ogrid,n_ogrid+1,k)-f(l1_ogrid-2:l2_ogrid-2,m_ogrid,n_ogrid+1,k))  &
            !           +(f(l1_ogrid+3:l2_ogrid+3,m_ogrid,n_ogrid+1,k)-f(l1_ogrid-3:l2_ogrid-3,m_ogrid,n_ogrid+1,k))) &
            !      -(45.*(f(l1_ogrid+1:l2_ogrid+1,m_ogrid,n_ogrid-1,k)-f(l1_ogrid-1:l2_ogrid-1,m_ogrid,n_ogrid-1,k))  &
            !        -9.*(f(l1_ogrid+2:l2_ogrid+2,m_ogrid,n_ogrid-1,k)-f(l1_ogrid-2:l2_ogrid-2,m_ogrid,n_ogrid-1,k))  &
            !           +(f(l1_ogrid+3:l2_ogrid+3,m_ogrid,n_ogrid-1,k)-f(l1_ogrid-3:l2_ogrid-3,m_ogrid,n_ogrid-1,k))))&
            !  -9.*((45.*(f(l1_ogrid+1:l2_ogrid+1,m_ogrid,n_ogrid+2,k)-f(l1_ogrid-1:l2_ogrid-1,m_ogrid,n_ogrid+2,k))  &
            !        -9.*(f(l1_ogrid+2:l2_ogrid+2,m_ogrid,n_ogrid+2,k)-f(l1_ogrid-2:l2_ogrid-2,m_ogrid,n_ogrid+2,k))  &
            !           +(f(l1_ogrid+3:l2_ogrid+3,m_ogrid,n_ogrid+2,k)-f(l1_ogrid-3:l2_ogrid-3,m_ogrid,n_ogrid+2,k))) &
            !      -(45.*(f(l1_ogrid+1:l2_ogrid+1,m_ogrid,n_ogrid-2,k)-f(l1_ogrid-1:l2_ogrid-1,m_ogrid,n_ogrid-2,k))  &
            !        -9.*(f(l1_ogrid+2:l2_ogrid+2,m_ogrid,n_ogrid-2,k)-f(l1_ogrid-2:l2_ogrid-2,m_ogrid,n_ogrid-2,k))  &
            !           +(f(l1_ogrid+3:l2_ogrid+3,m_ogrid,n_ogrid-2,k)-f(l1_ogrid-3:l2_ogrid-3,m_ogrid,n_ogrid-2,k))))&
            !     +((45.*(f(l1_ogrid+1:l2_ogrid+1,m_ogrid,n_ogrid+3,k)-f(l1_ogrid-1:l2_ogrid-1,m_ogrid,n_ogrid+3,k))  &
            !        -9.*(f(l1_ogrid+2:l2_ogrid+2,m_ogrid,n_ogrid+3,k)-f(l1_ogrid-2:l2_ogrid-2,m_ogrid,n_ogrid+3,k))  &
            !           +(f(l1_ogrid+3:l2_ogrid+3,m_ogrid,n_ogrid+3,k)-f(l1_ogrid-3:l2_ogrid-3,m_ogrid,n_ogrid+3,k))) &
            !      -(45.*(f(l1_ogrid+1:l2_ogrid+1,m_ogrid,n_ogrid-3,k)-f(l1_ogrid-1:l2_ogrid-1,m_ogrid,n_ogrid-3,k))  &
            !        -9.*(f(l1_ogrid+2:l2_ogrid+2,m_ogrid,n_ogrid-3,k)-f(l1_ogrid-2:l2_ogrid-2,m_ogrid,n_ogrid-3,k))  &
            !           +(f(l1_ogrid+3:l2_ogrid+3,m_ogrid,n_ogrid-3,k)-f(l1_ogrid-3:l2_ogrid-3,m_ogrid,n_ogrid-3,k))))&
            !       )
          else
            df=0.
            if (ip<=5) print*, 'derij: Degenerate case in x- or z-direction'
          endif
        endif
!
      endif                     ! bidiagonal derij

!  Since we have cylindrical coordinates
      if ( i+j==3 .or. i+j==5 ) df=df*rcyl_mn1_ogrid
!
    endsubroutine derij_ogrid
!***********************************************************************
    subroutine der6_ogrid(f, k, df, j)
!
!  Calculats D^(6)*dx^5/60, which is the upwind correction of centered derivatives.
!
!   27-feb-17/Jorgen: Adapted from deriv.f90
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid) :: f
      real, dimension (nx_ogrid) :: df,fac
      integer :: j,k
!
      intent(in)  :: f,k,j
      intent(out) :: df
!
      if (j==1) then
        if (nxgrid_ogrid/=1) then
          fac=(1.0/60)*dx_1_ogrid(l1_ogrid:l2_ogrid)
          !df(1:i)=0
          df(1:nx_ogrid)=fac(1:nx_ogrid)* &
                               (- 20.0* f(l1_ogrid  :l2_ogrid  ,m_ogrid,n_ogrid,k) &
                                + 15.0*(f(l1_ogrid+1:l2_ogrid+1,m_ogrid,n_ogrid,k)+f(l1_ogrid-1:l2_ogrid-1,m_ogrid,n_ogrid,k)) &
                                -  6.0*(f(l1_ogrid+2:l2_ogrid+2,m_ogrid,n_ogrid,k)+f(l1_ogrid-2:l2_ogrid-2,m_ogrid,n_ogrid,k)) &
                                +      (f(l1_ogrid+3:l2_ogrid+3,m_ogrid,n_ogrid,k)+f(l1_ogrid-3:l2_ogrid-3,m_ogrid,n_ogrid,k)))
        else
          df=0.
        endif
!
!  Settin df(1:i) = 0 means setting the upwind correction to zero for the points 
!  closest to the surface and at the surface. This is necessary to be consistent
!  with the boundary closures, and to not use values of f_ogrid inside the cylinder
!
        if(lfirst_proc_x) then
          if(SBP) then
            if (.not. SBP_reduced_to2) then
              df(1:6)=0.
            else
              df(1:4)=0.
            endif
          elseif(BDRY5) then
            df(1:3)=0.
          endif
        endif
      elseif (j==2) then
        if (ny_ogrid/=1) then
          fac=(1.0/60)*dy_1_ogrid(m_ogrid)
          df=fac*(- 20.0* f(l1_ogrid:l2_ogrid,m_ogrid  ,n_ogrid,k) &
                  + 15.0*(f(l1_ogrid:l2_ogrid,m_ogrid+1,n_ogrid,k)+f(l1_ogrid:l2_ogrid,m_ogrid-1,n_ogrid,k)) &
                  -  6.0*(f(l1_ogrid:l2_ogrid,m_ogrid+2,n_ogrid,k)+f(l1_ogrid:l2_ogrid,m_ogrid-2,n_ogrid,k)) &
                  +      (f(l1_ogrid:l2_ogrid,m_ogrid+3,n_ogrid,k)+f(l1_ogrid:l2_ogrid,m_ogrid-3,n_ogrid,k)))
        else
          df=0.
        endif
      elseif (j==3) then
        if (nz_ogrid/=1) then
          fac=(1.0/60)*dz_1_ogrid(n_ogrid)
          df=fac*(- 20.0* f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid  ,k) &
                  + 15.0*(f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid+1,k)+f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid-1,k)) &
                  -  6.0*(f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid+2,k)+f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid-2,k)) &
                  +      (f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid+3,k)+f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid-3,k)))
        else
          df=0.
        endif
      endif
!
    endsubroutine der6_ogrid
!! !***********************************************************************
!!     subroutine der_upwnd(f, k, df, j)
!! !
!! !  Calculats d(f)/dx_j using 5th order upwind meathod, with 3rd order closure near surface
!! !  Only to be used for irho at the moment
!! !
!! !   25-sep-17/Jorgen: Coded
!! !
!!       real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid) :: f
!!       real, dimension (nx_ogrid) :: df,fac
!!       integer :: j,k,i
!! !
!!       intent(in)  :: f,k,j
!!       intent(out) :: df
!! !
!!         
!! !TODO: FIX THIS
!! !TODO: FIX THIS
!!       if (j==1) then
!!         if (nxgrid_ogrid/=1) then
!!           fac=(1.0/60)*dx_1_ogrid(l1_ogrid:l2_ogrid)
!!           fac(2) = fac(2)*10
!!           if(k==irho) then
!!             if(lfirst_proc_x) then
!!               df(1)=0.
!!               df(2)=fac(2)*(3.0*f(l1_ogrid+1  ,m_ogrid,n_ogrid,k) &
!!                            -6.0*f(l1_ogrid+1+1,m_ogrid,n_ogrid,k) &
!!                            +2.0*f(l1_ogrid+1-1,m_ogrid,n_ogrid,k) &
!!                            +    f(l1_ogrid+1+2,m_ogrid,n_ogrid,k))
!!               i=2
!!             else
!!               i=0
!!           else
!!             call fatal_error('der_upwnd','Upwinding only implemented for density')
!!           i3=l1_ogrid+2
!!           df(1+i:nx_ogrid)=fac(1+i:nx_ogrid)* &
!!                                (  20.0*f(i3  :l2_ogrid  ,m_ogrid,n_ogrid,k) &
!!                                 - 60.0*f(i3+1:l2_ogrid+1,m_ogrid,n_ogrid,k) &
!!                                 + 30.0*f(i3-1:l2_ogrid-1,m_ogrid,n_ogrid,k) &
!!                                 + 15.0*f(i3+2:l2_ogrid+2,m_ogrid,n_ogrid,k) &
!!                                 -  3.0*f(i3-2:l2_ogrid-2,m_ogrid,n_ogrid,k) &
!!                                 -  2.0*f(i3+3:l2_ogrid+3,m_ogrid,n_ogrid,k) )
!!         else
!!           df=0.
!!         endif
!!       elseif (j==2) then
!!         if (nygrid_ogrid/=1) then
!!           fac=(1.0/60)*dy_1_ogrid(m_ogrid)
!!           df=fac*(- 20.0* f(l1_ogrid:l2_ogrid,m_ogrid  ,n_ogrid,k) &
!!                   + 15.0*(f(l1_ogrid:l2_ogrid,m_ogrid+1,n_ogrid,k)+f(l1_ogrid:l2_ogrid,m_ogrid-1,n_ogrid,k)) &
!!                   -  6.0*(f(l1_ogrid:l2_ogrid,m_ogrid+2,n_ogrid,k)+f(l1_ogrid:l2_ogrid,m_ogrid-2,n_ogrid,k)) &
!!                   +      (f(l1_ogrid:l2_ogrid,m_ogrid+3,n_ogrid,k)+f(l1_ogrid:l2_ogrid,m_ogrid-3,n_ogrid,k)))
!!         else
!!           df=0.
!!         endif
!!       elseif (j==3) then
!!         if (nzgrid_ogrid/=1) then
!!           fac=(1.0/60)*dz_1_ogrid(n_ogrid)
!!           df=fac*(- 20.0* f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid  ,k) &
!!                   + 15.0*(f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid+1,k)+f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid-1,k)) &
!!                   -  6.0*(f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid+2,k)+f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid-2,k)) &
!!                   +      (f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid+3,k)+f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid-3,k)))
!!         else
!!           df=0.
!!         endif
!!       endif
!! !
!!     endsubroutine der6_ogrid
!! !***********************************************************************
    subroutine deri_3d_inds_ogrid(f,df,inds,j,lignored,lnometric)
!
!  dummy routine for compatibility
!
!  27-feb-17/Jorgen: Adapted from deriv.f90
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid) :: f
      real, dimension (nx_ogrid)                           :: df
      integer                             :: j
      logical,                   optional :: lignored, lnometric
      integer, dimension(nx_ogrid)              :: inds
!
      intent(in)  :: f,j,inds,lignored,lnometric
      intent(out) :: df
!
      call fatal_error('deri_3d_inds_ogrid','Upwinding not implemented for nonuniform grids')
!
! dummy computation to avoid compiler warnings of unused variables
      if (present(lignored).and.present(lnometric)) &
          df  = inds + f(l1_ogrid:l2_ogrid,1,1,1) + j
!
    endsubroutine deri_3d_inds_ogrid
!************************************************************************
    subroutine set_ghosts_onesided_ogrid(ivar)
!
!   Set ghost points for onesided boundary conditions with Dirichlet BC
!   on the cylidner surface.
!   Only works for the radial direction.
!
!   16-feb-17/Jorgen: Adapted from deriv.f90.
!
      integer :: k,ivar,i
!
      do i=1,nghost
        k=l1_ogrid-i
        f_ogrid(k,:,:,ivar)=7*f_ogrid(k+1,:,:,ivar) &
                          -21*f_ogrid(k+2,:,:,ivar) &
                          +35*f_ogrid(k+3,:,:,ivar) &
                          -35*f_ogrid(k+4,:,:,ivar) &
                          +21*f_ogrid(k+5,:,:,ivar) &
                           -7*f_ogrid(k+6,:,:,ivar) &
                             +f_ogrid(k+7,:,:,ivar)
      enddo

    endsubroutine set_ghosts_onesided_ogrid
!***********************************************************************
    subroutine bval_from_neumann_arr_ogrid
!
!  Calculates the boundary value from the Neumann BC d f/d x_i = val employing
!  one-sided difference formulae. val depends on x,y.
!
!  16-feb-17/Jorgen: Adapted from deriv.f90
!  jan-19/Eva: added BC for chemistry
!
!
      real :: val=0.
      real, dimension (my_ogrid, mz_ogrid) :: mdot_c,diff_coeff
      integer :: k,j,i,count_dos,i_CO2
      real :: MC_MO2 = 12.0107/(2.*15.9994), MC_MCO2 = 12.0107/(12.0107 + 2.*15.9994)
      real :: MCO_MO2 = (12.0107+15.9994)/(2.*15.9994), MCO_MCO2 = (12.0107+15.9994)/(12.0107 + 2.*15.9994)
      real :: A_CO2, A_O2, B_CO2, B_O2, C_CO2, C_O2, rho_Rate_CO2, rho_Rate_O2
      real :: YO2_chosen_l=0., YO2_chosen_r=0., YO2_computed_l, YO2_computed_r, YO2_computed=0, Y_CO2=0
      real :: Y_min, Y_max, q, residual_l, residual_r

      k=l1_ogrid
      if (lexpl_rho) then
        f_ogrid(k,:,:,irho) = (-val*60.*(1./dx_1_ogrid(l1_ogrid))+ 360.*f_ogrid(k+1,:,:,irho) &
                                                 - 450.*f_ogrid(k+2,:,:,irho) &
                                                 + 400.*f_ogrid(k+3,:,:,irho) &
                                                 - 225.*f_ogrid(k+4,:,:,irho) &
                                                 +  72.*f_ogrid(k+5,:,:,irho) &
                                                 -  10.*f_ogrid(k+6,:,:,irho) )/147.
      endif
      if (lchemistry) then
        if (.not. lreac_heter) then
          do j = 1,nchemspec
              f_ogrid(k,:,:,ichemspec(j)) = (-val*60.*1./dx_1_ogrid(l1_ogrid) + 360.*f_ogrid(k+1,:,:,ichemspec(j)) &
                                             - 450.*f_ogrid(k+2,:,:,ichemspec(j)) &
                                             + 400.*f_ogrid(k+3,:,:,ichemspec(j)) &
                                             - 225.*f_ogrid(k+4,:,:,ichemspec(j)) &
                                             +  72.*f_ogrid(k+5,:,:,ichemspec(j)) &
                                             -  10.*f_ogrid(k+6,:,:,ichemspec(j)) )/147.
! Set-up ghost points based on non-symmetric stencils
              f_ogrid(k-1,:,:,ichemspec(j)) = (-val*60.*1./dx_1_ogrid(l1_ogrid) -  77.*f_ogrid(k  ,:,:,ichemspec(j)) &
                                               + 150.*f_ogrid(k+1,:,:,ichemspec(j)) &
                                               - 100.*f_ogrid(k+2,:,:,ichemspec(j)) &
                                               +  50.*f_ogrid(k+3,:,:,ichemspec(j)) &
                                               -  15.*f_ogrid(k+4,:,:,ichemspec(j)) &
                                               +   2.*f_ogrid(k+5,:,:,ichemspec(j)) )/10.

              f_ogrid(k-2,:,:,ichemspec(j)) = (-val*60.*1./dx_1_ogrid(l1_ogrid) -  24.*f_ogrid(k-1,:,:,ichemspec(j)) &
                                               -  35.*f_ogrid(k  ,:,:,ichemspec(j)) &
                                               +  80.*f_ogrid(k+1,:,:,ichemspec(j)) &
                                               -  30.*f_ogrid(k+2,:,:,ichemspec(j)) &
                                               +   8.*f_ogrid(k+3,:,:,ichemspec(j)) &
                                               -   1.*f_ogrid(k+4,:,:,ichemspec(j)) )/(-2.)

              f_ogrid(k-3,:,:,ichemspec(j)) = (-val*60.*1./dx_1_ogrid(l1_ogrid) +   9.*f_ogrid(k-2,:,:,ichemspec(j)) &
                                               -  45.*f_ogrid(k-1,:,:,ichemspec(j)) &
                                               +  45.*f_ogrid(k+1,:,:,ichemspec(j)) &
                                               -   9.*f_ogrid(k+2,:,:,ichemspec(j)) &
                                               +   1.*f_ogrid(k+3,:,:,ichemspec(j)) )
          enddo
        else

! if heterogeneous reactions grad(Y_k) = (-mdot_c*Y_k-m_k)/(rho*D_k)
!
          if (lnonegative_Yk) then
            i_CO2 = (mvar-nchemspec) + ichem_CO2
            q=(3.-sqrt(5.))/2.
            diff_coeff = f_ogrid(k,:,:,iviscosity)*Pr_number1
            do i = m1_ogrid,m2_ogrid
                count_dos = 0
                Y_min = 0.0
                Y_max = 0.2
                if (f_ogrid(k,i,4,i_O2) > 0) then
                  rho_Rate_O2 = -heter_reaction_rate(i,4,ichem_O2)/f_ogrid(k,i,4,i_O2)
                else
                  rho_Rate_O2 = 0.
                endif
                if (f_ogrid(k,i,4,i_CO2) > 0) then
                  rho_Rate_CO2 = -heter_reaction_rate(i,4,ichem_CO2)/f_ogrid(k,i,4,i_CO2)
                else
                  rho_Rate_CO2 = 0.
                endif
                A_CO2 = -MC_MCO2*rho_Rate_CO2
                A_O2 = -2.*MC_MO2*rho_Rate_O2

                if ((A_CO2 .ne. 0.) .and. (A_O2 .ne. 0.)) then

                  B_CO2 = f_ogrid(k,i,4,irho)*diff_coeff(i,4)*Lewis_coef1(ichem_CO2)*(-147./60.*dx_1_ogrid(l1_ogrid)) - rho_Rate_CO2
                  B_O2 = f_ogrid(k,i,4,irho)*diff_coeff(i,4)*Lewis_coef1(ichem_O2)*(-147./60.*dx_1_ogrid(l1_ogrid)) - rho_Rate_O2

                  C_CO2 = f_ogrid(k,i,4,irho)*diff_coeff(i,4)*Lewis_coef1(ichem_CO2)*dx_1_ogrid(l1_ogrid)/60.*( &
                        + 360.*f_ogrid(k+1,i,4,i_CO2) &
                        - 450.*f_ogrid(k+2,i,4,i_CO2) &
                        + 400.*f_ogrid(k+3,i,4,i_CO2) &
                        - 225.*f_ogrid(k+4,i,4,i_CO2) &
                        +  72.*f_ogrid(k+5,i,4,i_CO2) &
                        -  10.*f_ogrid(k+6,i,4,i_CO2))
                  C_O2 = f_ogrid(k,i,4,irho)*diff_coeff(i,4)*Lewis_coef1(ichem_O2)*dx_1_ogrid(l1_ogrid)/60.*( &
                       + 360.*f_ogrid(k+1,i,4,i_O2) &
                       - 450.*f_ogrid(k+2,i,4,i_O2) &
                       + 400.*f_ogrid(k+3,i,4,i_O2) &
                       - 225.*f_ogrid(k+4,i,4,i_O2) &
                       +  72.*f_ogrid(k+5,i,4,i_O2) &
                       -  10.*f_ogrid(k+6,i,4,i_O2))

                  YO2_chosen_l = Y_min + q*(Y_max-Y_min)
                  YO2_chosen_r = Y_max - q*(Y_max-Y_min)
                  call compute_Yresidual(YO2_chosen_l,YO2_computed,Y_CO2,residual_l,A_O2,B_O2,C_O2,A_CO2,B_CO2,C_CO2,&
                                         MC_MO2,MC_MCO2,rho_Rate_O2,rho_Rate_CO2,i)
                  call compute_Yresidual(YO2_chosen_r,YO2_computed,Y_CO2,residual_r,A_O2,B_O2,C_O2,A_CO2,B_CO2,C_CO2,&
                                         MC_MO2,MC_MCO2,rho_Rate_O2,rho_Rate_CO2,i)
                  do while ((residual_l > toler) .and. (residual_r > toler))
                      if (residual_l .le. residual_r) then
                          Y_max = YO2_chosen_r
                          residual_r = residual_l
                          YO2_chosen_r = YO2_chosen_l
                          YO2_chosen_l = Y_min + q*(Y_max-Y_min)
                          call compute_Yresidual(YO2_chosen_l,YO2_computed,Y_CO2,residual_l,A_O2,B_O2,C_O2,A_CO2,B_CO2,C_CO2,&
                                                 MC_MO2,MC_MCO2,rho_Rate_O2,rho_Rate_CO2,i)
                      elseif (residual_l > residual_r) then
                          Y_min = YO2_chosen_l
                          residual_l = residual_r
                          YO2_chosen_l = YO2_chosen_r
                          YO2_chosen_r = Y_max - q*(Y_max-Y_min)
                          call compute_Yresidual(YO2_chosen_r,YO2_computed,Y_CO2,residual_r,A_O2,B_O2,C_O2,A_CO2,B_CO2,C_CO2,&
                                                 MC_MO2,MC_MCO2,rho_Rate_O2,rho_Rate_CO2,i)
                      endif
                      count_dos = count_dos + 1
                      if (count_dos > 10000) then
                          print*,'residual_l,residual_r,YO2_l,YO2_r',residual_l,residual_r,YO2_chosen_l,YO2_chosen_r
                          call fatal_error('bval_from_neumann_arr_ogrid','no solutions for YO2 in the assumed range')
                      endif
                  end do
                  f_ogrid(k,i,:,i_CO2) = Y_CO2
                  f_ogrid(k,i,:,i_O2 ) = YO2_computed
                  heter_reaction_rate(i,:,nchemspec+1) = -(2.*MC_MO2*rho_Rate_O2*YO2_computed + MC_MCO2*rho_Rate_CO2*Y_CO2)
                  heter_reaction_rate(i,:,ichem_CO2) = -rho_Rate_CO2*Y_CO2
                  heter_reaction_rate(i,:,ichem_O2) = -rho_Rate_O2*YO2_computed
                  heter_reaction_rate(i,:,ichem_CO) = 2.*(MCO_MO2*rho_Rate_O2*YO2_computed + MCO_MCO2*rho_Rate_CO2*Y_CO2)
                endif
            enddo
            mdot_c = heter_reaction_rate(:,:,nchemspec+1)

            do j = 1,nchemspec
                if (((j .ne. 4) .and. (j .ne. 1)) .or. &
                ((j .eq. 4) .and. (A_O2 .eq. 0.)) .or. ((j .eq. 1) .and. (A_CO2 .eq. 0.))) then
                  diff_coeff = f_ogrid(k,:,:,iviscosity)*Pr_number1*Lewis_coef1(j)
                  f_ogrid(k,:,:,ichemspec(j)) = (heter_reaction_rate(:,:,j)*60.*1./dx_1_ogrid(l1_ogrid)&
                                              /(f_ogrid(k,:,:,irho)*diff_coeff) &
                                              + 360.*f_ogrid(k+1,:,:,ichemspec(j)) &
                                              - 450.*f_ogrid(k+2,:,:,ichemspec(j)) &
                                              + 400.*f_ogrid(k+3,:,:,ichemspec(j)) &
                                              - 225.*f_ogrid(k+4,:,:,ichemspec(j)) &
                                              +  72.*f_ogrid(k+5,:,:,ichemspec(j)) &
                                              -  10.*f_ogrid(k+6,:,:,ichemspec(j))) &
                                              /(147.-mdot_c*60.*1./dx_1_ogrid(l1_ogrid)/(f_ogrid(k,:,:,irho)*diff_coeff))

                endif
            enddo
          else
            mdot_c = heter_reaction_rate(:,:,nchemspec+1)
            do j = 1,nchemspec
              diff_coeff = f_ogrid(k,:,:,iviscosity)*Pr_number1*Lewis_coef1(j)
              f_ogrid(k,:,:,ichemspec(j)) = (heter_reaction_rate(:,:,j)*60.*1./dx_1_ogrid(l1_ogrid) &
                                          / (f_ogrid(k,:,:,irho)*diff_coeff)   &
                                          + 360.*f_ogrid(k+1,:,:,ichemspec(j)) &
                                          - 450.*f_ogrid(k+2,:,:,ichemspec(j)) &
                                          + 400.*f_ogrid(k+3,:,:,ichemspec(j)) &
                                          - 225.*f_ogrid(k+4,:,:,ichemspec(j)) &
                                          +  72.*f_ogrid(k+5,:,:,ichemspec(j)) &
                                          -  10.*f_ogrid(k+6,:,:,ichemspec(j))) &
                                          /(147.-mdot_c*60.*1./dx_1_ogrid(l1_ogrid)/(f_ogrid(k,:,:,irho)*diff_coeff))
            enddo
          endif
        endif
      endif
!
    endsubroutine bval_from_neumann_arr_ogrid
!***********************************************************************
    subroutine compute_Yresidual(YO2,YO2_computed,Y_CO2,residual,A_O2,&
               B_O2,C_O2,A_CO2,B_CO2,C_CO2,MC_MO2,MC_MCO2,rho_Rate_O2,rho_Rate_CO2,i)
!
      real, intent(in) :: A_CO2, A_O2, B_CO2, B_O2, C_CO2, C_O2
      real, intent(out) :: YO2_computed, Y_CO2, residual
      real, intent(in) :: rho_Rate_O2, rho_Rate_CO2, MC_MO2, MC_MCO2, YO2
      real  :: B1, B2, YO2_YCO2, B_CO2_rest, F_O2
      integer :: i
!
        B_CO2_rest = - 2.*MC_MO2*rho_Rate_O2
        B1 = B_CO2 + B_CO2_rest*YO2
        Y_CO2 = (- B1 - sqrt(B1**2 - 4.*A_CO2*C_CO2))/(2.*A_CO2)
        if (Y_CO2 > 1. .or. Y_CO2 < 0.) Y_CO2 = Y_CO2 + sqrt(B1**2 - 4.*A_CO2*C_CO2)/A_CO2
!
        B2 = - MC_MCO2*rho_Rate_CO2
        YO2_YCO2 = (-A_CO2*Y_CO2**2 - B_CO2*Y_CO2 - C_CO2)/B_CO2_rest
        F_O2 = A_O2*YO2**2 + B_O2*YO2 + B2*YO2_YCO2 + C_O2

        residual = abs(F_O2)
        YO2_computed = YO2
!
    endsubroutine compute_Yresidual
!***********************************************************************
    subroutine bval_from_neumann_SBP(f_og)
!
!  Calculates the boundary value from the Neumann BC d f/d x_i = val employing
!  one-sided difference formulae. val depends on x,y.
!  Only implemented for df/dx_i = 0 at boundary.
!
!  16-feb-17/Jorgen: Adapted from deriv.f90
!  07-nov-18/Jonas: Added treatment for cases with temperature gradient
!  jan-19/Eva: added BC for chemistry
!
      use SharedVariables, only: get_shared_variable
!
      integer :: k, j, i, l
      real, dimension (mx_ogrid, my_ogrid, mz_ogrid,mfarray_ogrid), intent(inout) ::  f_og
      real, dimension (my_ogrid, mz_ogrid) :: df_surf, dR, diff_coeff, mdot_c
      real, dimension (my_ogrid, mz_ogrid) :: first_term
!
      k=l1_ogrid
!
      if (lexpl_rho) then
        if (iTT==0) then
          f_og   (k,:,:,irho) = -(D1_SBP(1,2)*f_og   (k+1,:,:,irho) + &
              D1_SBP(1,3)*f_og   (k+2,:,:,irho) + &
              D1_SBP(1,4)*f_og   (k+3,:,:,irho) + &
              D1_SBP(1,5)*f_og   (k+4,:,:,irho) + &
              D1_SBP(1,6)*f_og   (k+5,:,:,irho) + &
              D1_SBP(1,7)*f_og   (k+6,:,:,irho) + &
              D1_SBP(1,8)*f_og   (k+7,:,:,irho) + &
              D1_SBP(1,9)*f_og   (k+8,:,:,irho) )/D1_SBP(1,1)
        else if (iTT>0) then
          call der_ogrid_SBP_surf(f_og,df_surf,iTT)
          f_og   (k,:,:,irho) = -(D1_SBP(1,2)*f_og   (k+1,:,:,irho) + &
              D1_SBP(1,3)*f_og   (k+2,:,:,irho) + &
              D1_SBP(1,4)*f_og   (k+3,:,:,irho) + &
              D1_SBP(1,5)*f_og   (k+4,:,:,irho) + &
              D1_SBP(1,6)*f_og   (k+5,:,:,irho) + &
              D1_SBP(1,7)*f_og   (k+6,:,:,irho) + &
              D1_SBP(1,8)*f_og   (k+7,:,:,irho) + &
              D1_SBP(1,9)*f_og   (k+8,:,:,irho) )/ &
              (D1_SBP(1,1)+ df_surf / f_og(k,:,:,iTT))
        else
          call fatal_error('solid_cells_ogrid','temperature SBP index not found')
        endif
      endif
!
! grad(Y_k) = 0 when homogeneous reactions only
!
      if (lchemistry .and. (.not. lreac_heter)) then
        do j = 1,nchemspec
          f_og   (k,:,:,ichemspec(j)) = -(D1_SBP(1,2)*f_og   (k+1,:,:,ichemspec(j)) + &
              D1_SBP(1,3)*f_og   (k+2,:,:,ichemspec(j)) + &
              D1_SBP(1,4)*f_og   (k+3,:,:,ichemspec(j)) + &
              D1_SBP(1,5)*f_og   (k+4,:,:,ichemspec(j)) + &
              D1_SBP(1,6)*f_og   (k+5,:,:,ichemspec(j)) + &
              D1_SBP(1,7)*f_og   (k+6,:,:,ichemspec(j)) + &
              D1_SBP(1,8)*f_og   (k+7,:,:,ichemspec(j)) + &
              D1_SBP(1,9)*f_og   (k+8,:,:,ichemspec(j)) )/D1_SBP(1,1)
        enddo
!
! if heterogeneous reactions grad(Y_k) = (-mdot_c*Y_k-m_k)/(rho*D_k) 
!
      elseif (lchemistry) then
        mdot_c = heter_reaction_rate(:,:,nchemspec+1)
        do j = 1,nchemspec
          diff_coeff = f_og(l1_ogrid,:,:,iviscosity)*Pr_number1*Lewis_coef1(j)
          first_term = -(D1_SBP(1,2)*f_og(k+1,:,:,ichemspec(j)) + &
              D1_SBP(1,3)*f_og   (k+2,:,:,ichemspec(j)) + &
              D1_SBP(1,4)*f_og   (k+3,:,:,ichemspec(j)) + &
              D1_SBP(1,5)*f_og   (k+4,:,:,ichemspec(j)) + &
              D1_SBP(1,6)*f_og   (k+5,:,:,ichemspec(j)) + &
              D1_SBP(1,7)*f_og   (k+6,:,:,ichemspec(j)) + &
              D1_SBP(1,8)*f_og   (k+7,:,:,ichemspec(j)) + &
              D1_SBP(1,9)*f_og   (k+8,:,:,ichemspec(j)))* &
              dx_1_ogrid(l1_ogrid)*diff_coeff*f_og(k,:,:,irho)
!       if (nchemspec == 2 .or. nchemspec == 3 .or. nchemspec == 5) then
          if (.not. llin_BC) then
            f_og(k,:,:,ichemspec(j)) = (first_term - heter_reaction_rate(:,:,j))/&
                  (D1_SBP(1,1)*diff_coeff*f_og(k,:,:,irho)*dx_1_ogrid(l1_ogrid) + mdot_c)
          else
            f_og(k,:,:,ichemspec(j)) = -(heter_reaction_rate(:,:,j) + f_og(k+1,:,:,ichemspec(j))*&
                                        diff_coeff*f_og(k,:,:,irho)*dx_1_ogrid(l1_ogrid))/&
                                       (mdot_c - diff_coeff*f_og(k,:,:,irho)*dx_1_ogrid(l1_ogrid))
          endif
!
          if (lnonegative_Yk) then
            do i = m1_ogrid,m2_ogrid
              if (f_og(k,i,4,ichemspec(j)) .lt. 0) then
                heter_reaction_rate(i,:,j) = first_term(i,:)
                heter_reaction_rate(i,:,ichem_CO2) = -2.*((12.0107+15.9994)&
                  /(2.*15.9994)*heter_reaction_rate(i,:,ichem_O2)+(12.0107+15.9994)&
                  /(12.0107+2.*15.9994)*heter_reaction_rate(i,:,ichem_CO2))
                mdot_c(i,:) = -(heter_reaction_rate(i,:,ichem_CO) + &
                               heter_reaction_rate(i,:,ichem_CO2) + heter_reaction_rate(i,:,ichem_O2))
                do l = 1,nchemspec
                  diff_coeff(i,:) = f_og(l1_ogrid,i,:,iviscosity)*Pr_number1*Lewis_coef1(l)
                  first_term(i,:) = -(D1_SBP(1,2)*f_og(k+1,i,:,ichemspec(l)) + &
                      D1_SBP(1,3)*f_og   (k+2,i,:,ichemspec(l)) + &
                      D1_SBP(1,4)*f_og   (k+3,i,:,ichemspec(l)) + &
                      D1_SBP(1,5)*f_og   (k+4,i,:,ichemspec(l)) + &
                      D1_SBP(1,6)*f_og   (k+5,i,:,ichemspec(l)) + &
                      D1_SBP(1,7)*f_og   (k+6,i,:,ichemspec(l)) + &
                      D1_SBP(1,8)*f_og   (k+7,i,:,ichemspec(l)) + &
                      D1_SBP(1,9)*f_og   (k+8,i,:,ichemspec(l)))* &
                      dx_1_ogrid(l1_ogrid)*diff_coeff(i,:)*f_og(k,i,:,irho)
                  f_og   (k,i,:,ichemspec(l)) = (first_term(i,:) - heter_reaction_rate(i,:,l))/ &
                      (D1_SBP(1,1)*diff_coeff(i,:)*f_og(k,i,:,irho)*dx_1_ogrid(l1_ogrid) + mdot_c(i,:))
                enddo
              endif
            enddo
          endif
        ! The below is the proper expression for the line above if gradM is accounted for in equations
        ! (D1_SBP(1,1)*diff_coeff*f_og(k,:,:,irho)*dx_1_ogrid(l1_ogrid) + mdot_c - diff_coeff*f_og(k,:,:,irho)*grad_lnR)
        ! If this is used remember to compute grad_lnR :
        !        call der_ogrid_SBP_surf(f_og,dR,iRR)
        !        grad_lnR(m1_ogrid:m2_ogrid,n1_ogrid:n2_ogrid) &
        !            = dR(m1_ogrid:m2_ogrid,n1_ogrid:n2_ogrid)&
        !            /f_og(k,m1_ogrid:m2_ogrid,n1_ogrid:n2_ogrid,iRR)
        enddo
      endif
!
    endsubroutine bval_from_neumann_SBP
!***********************************************************************
    subroutine bval_from_neumann_bdry5(f_og)
!
!  Calculates the boundary value from the Neumann BC d f/d x_i = val employing
!  one-sided difference formulae. val depends on x,y.
!  Only implemented for df/dx_i = 0 at boundary.
!
!  14-okt-17/Jorgen: Adapted from deriv.f90
!
     real, dimension (mx_ogrid, my_ogrid, mz_ogrid,mfarray_ogrid), intent(inout) ::  f_og
     integer :: j
!
     if(lexpl_rho) then
       f_og   (l1_ogrid,:,:,irho)  = ( 5.0000000000000000 *f_og   (l1_ogrid+1,:,:,irho) &
           -5.0000000000000000 *f_og   (l1_ogrid+2,:,:,irho) &
           +3.3333333333333333 *f_og   (l1_ogrid+3,:,:,irho) &
           -1.2500000000000000 *f_og   (l1_ogrid+4,:,:,irho) &
           +0.20000000000000000*f_og   (l1_ogrid+5,:,:,irho))&
           /( 2.2833333333333333)   
     endif
!
     if(lchemistry) then
       do j = 1,nchemspec
         f_og   (l1_ogrid,:,:,ichemspec(j))  = &
             (5.0000000000000000 *f_og   (l1_ogrid+1,:,:,ichemspec(j)) &
             -5.0000000000000000 *f_og   (l1_ogrid+2,:,:,ichemspec(j)) &
             +3.3333333333333333 *f_og   (l1_ogrid+3,:,:,ichemspec(j)) &
             -1.2500000000000000 *f_og   (l1_ogrid+4,:,:,ichemspec(j)) &
             +0.20000000000000000*f_og   (l1_ogrid+5,:,:,ichemspec(j)))&
             /( 2.2833333333333333)  
       enddo
     endif
     
    endsubroutine bval_from_neumann_bdry5
!***********************************************************************
    subroutine der_ogrid_SBP(f,df)
! 
!  Summation by parts boundary condition for first derivative.
!  Only implemented in radial direction.
!
!  21-mar-17/Jorgen: Coded
      real, dimension(l1_ogrid:l1_ogrid+8), intent(in) :: f
      real, dimension(6), intent(out) :: df
      integer :: i

      do i=1,6
        df(i)=dx_1_ogrid(l1_ogrid+i-1)*(D1_SBP(i,1)*f(l1_ogrid  ) + &
                                        D1_SBP(i,2)*f(l1_ogrid+1) + &
                                        D1_SBP(i,3)*f(l1_ogrid+2) + &
                                        D1_SBP(i,4)*f(l1_ogrid+3) + &
                                        D1_SBP(i,5)*f(l1_ogrid+4) + &
                                        D1_SBP(i,6)*f(l1_ogrid+5) + &
                                        D1_SBP(i,7)*f(l1_ogrid+6) + &
                                        D1_SBP(i,8)*f(l1_ogrid+7) + &
                                        D1_SBP(i,9)*f(l1_ogrid+8) )
      enddo

    endsubroutine der_ogrid_SBP
!***********************************************************************
    subroutine der_ogrid_SBP_surf(f,df_surf,ii)
! 
!  Summation by parts boundary condition for first derivative.
!  Only implemented in radial direction.
!
!  21-mar-17/Jorgen: Coded
!  06-nov-18/jonas: adapted
!
      real, dimension(mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid), intent(in) :: f
      real, dimension(my_ogrid,mz_ogrid), intent(out) :: df_surf
      integer :: k, ii

      k = l1_ogrid
      df_surf(:,:) = D1_SBP(1,1)*f(k,:,:,ii) + &
           D1_SBP(1,2)*f(k+1,:,:,ii) + &
           D1_SBP(1,3)*f(k+2,:,:,ii) + &
           D1_SBP(1,4)*f(k+3,:,:,ii) + &
           D1_SBP(1,5)*f(k+4,:,:,ii) + &
           D1_SBP(1,6)*f(k+5,:,:,ii) + &
           D1_SBP(1,7)*f(k+6,:,:,ii) + &
           D1_SBP(1,8)*f(k+7,:,:,ii) + &
           D1_SBP(1,9)*f(k+8,:,:,ii)
      
    endsubroutine der_ogrid_SBP_surf
!***********************************************************************
    subroutine der_ogrid_bdry5(f,df)
! 
!  Fifth order boundary closures for first derivatives 
!  Only implemented in radial direction.
!
!  13-okt-17/Jorgen: Coded
      real, dimension(l1_ogrid:l1_ogrid+5), intent(in) :: f
      real, dimension(3), intent(out) :: df

      df(1) = dx_1_ogrid(l1_ogrid)*(&
              -2.2833333333333333   *f(l1_ogrid  ) &
              +5.0000000000000000   *f(l1_ogrid+1) &
              -5.0000000000000000   *f(l1_ogrid+2) &
              +3.3333333333333333   *f(l1_ogrid+3) &
              -1.2500000000000000   *f(l1_ogrid+4) &
              +0.20000000000000000  *f(l1_ogrid+5) )
      
      df(2) = dx_1_ogrid(l1_ogrid+1)*(&
              -0.20000000000000000  *f(l1_ogrid  ) &
              -1.0833333333333333   *f(l1_ogrid+1) &
              +2.0000000000000000   *f(l1_ogrid+2) &
              -1.0000000000000000   *f(l1_ogrid+3) &
              +0.33333333333333333  *f(l1_ogrid+4) &
              -0.050000000000000000 *f(l1_ogrid+5) )
      
      df(3) = dx_1_ogrid(l1_ogrid+2)*(&
              +0.050000000000000000 *f(l1_ogrid  ) &
              -0.50000000000000000  *f(l1_ogrid+1) &
              -0.33333333333333333  *f(l1_ogrid+2) &
              +1.0000000000000000   *f(l1_ogrid+3) &
              -0.25000000000000000  *f(l1_ogrid+4) &
              +0.033333333333333333 *f(l1_ogrid+5) )

    endsubroutine der_ogrid_bdry5
!***********************************************************************
    subroutine der2_ogrid_bdry5(f,df)
! 
!  Fifth order boundary closures for second derivatives 
!  Only implemented in radial direction.
!
!  13-okt-17/Jorgen: Coded
      real, dimension(l1_ogrid:l1_ogrid+6), intent(in) :: f
      real, dimension(3), intent(out) :: df

      df(1) = dx_1_ogrid(l1_ogrid)**2*(&
              +4.5111111111111111 *f(l1_ogrid  ) &
              -17.39999999999999  *f(l1_ogrid+1) &
              +29.25000000000000  *f(l1_ogrid+2) &
              -28.22222222222222  *f(l1_ogrid+3) &
              +16.50000000000000  *f(l1_ogrid+4) &
              -5.4000000000000000 *f(l1_ogrid+5) &
              +0.76111111111111111*f(l1_ogrid+6) )

      df(2) = dx_1_ogrid(l1_ogrid+1)**2*(&
              +0.76111111111111111  *f(l1_ogrid  ) &
              -0.81666666666666667  *f(l1_ogrid+1) &
              -1.4166666666666667   *f(l1_ogrid+2) &
              +2.6111111111111111   *f(l1_ogrid+3) &
              -1.5833333333333333   *f(l1_ogrid+4) &
              +0.51666666666666667  *f(l1_ogrid+5) &
              -0.072222222222222222 *f(l1_ogrid+6) )
                
      df(3) = dx_1_ogrid(l1_ogrid+2)**2*(&
              -0.072222222222222222*f(l1_ogrid  ) &
              +1.26666666666666667 *f(l1_ogrid+1) &
              -2.33333333333333333 *f(l1_ogrid+2) &
              +1.11111111111111111 *f(l1_ogrid+3) &
              +0.083333333333333333*f(l1_ogrid+4) &
              -0.066666666666666666*f(l1_ogrid+5) &
              +0.011111111111111111*f(l1_ogrid+6) )
    endsubroutine der2_ogrid_bdry5
!***********************************************************************
    subroutine der2_ogrid_bdry5_alt(f,df2,k)
! 
!  Fifth order boundary closures for second derivatives 
!  Only implemented in radial direction.
!  Alternative edition, that uses first derivatives twice, instead
!  of explicit use of second derivatives
!
!  13-okt-17/Jorgen: Coded
      real, dimension(mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid), intent(in) :: f
      real, dimension(3), intent(out) :: df2
      integer, intent(in) :: k
      real, dimension(6) :: df

      call der_ogrid_bdry5(f(l1_ogrid:l1_ogrid+5,m_ogrid,n_ogrid,k),df(1:3))
      df(4:6)=(1/60.)*dx_1_ogrid(l1_ogrid+3:l1_ogrid+5) *( &
                      + 45.0*(f(l1_ogrid+4:l1_ogrid+6,m_ogrid,n_ogrid,k)-f(l1_ogrid+2:l1_ogrid+4,m_ogrid,n_ogrid,k)) &
                      -  9.0*(f(l1_ogrid+5:l1_ogrid+7,m_ogrid,n_ogrid,k)-f(l1_ogrid+1:l1_ogrid+3,m_ogrid,n_ogrid,k)) &
                      +      (f(l1_ogrid+6:l1_ogrid+8,m_ogrid,n_ogrid,k)-f(l1_ogrid  :l1_ogrid+2,m_ogrid,n_ogrid,k)))
      df2(1) = dx_1_ogrid(l1_ogrid)*(&
              -2.2833333333333333   *df(1) &
              +5.0000000000000000   *df(2) &
              -5.0000000000000000   *df(3) &
              +3.3333333333333333   *df(4) &
              -1.2500000000000000   *df(5) &
              +0.20000000000000000  *df(6) )
      
      df2(2) = dx_1_ogrid(l1_ogrid+1)*(&
              -0.20000000000000000  *df(1) &
              -1.0833333333333333   *df(2) &
              +2.0000000000000000   *df(3) &
              -1.0000000000000000   *df(4) &
              +0.33333333333333333  *df(5) &
              -0.050000000000000000 *df(6) )
      
      df2(3) = dx_1_ogrid(l1_ogrid+2)*(&
              +0.050000000000000000 *df(1) &
              -0.50000000000000000  *df(2) &
              -0.33333333333333333  *df(3) &
              +1.0000000000000000   *df(4) &
              -0.25000000000000000  *df(5) &
              +0.033333333333333333 *df(6) )

    endsubroutine der2_ogrid_bdry5_alt
!***********************************************************************
    subroutine der2_ogrid_SBP(f,df2)
! 
!  Summation by parts boundary condition for second derivative.
!  Only implemented in radial direction.
!
!  21-mar-17/Jorgen: Coded
      real, dimension(l1_ogrid:l1_ogrid+8), intent(in) :: f
      real, dimension(6), intent(out) :: df2
      integer :: i

      do i=1,6
        df2(i)=(dx_1_ogrid(l1_ogrid+i-1)**2)*(D2_SBP(i,1)*f(l1_ogrid  ) + &
                                              D2_SBP(i,2)*f(l1_ogrid+1) + &
                                              D2_SBP(i,3)*f(l1_ogrid+2) + &
                                              D2_SBP(i,4)*f(l1_ogrid+3) + &
                                              D2_SBP(i,5)*f(l1_ogrid+4) + &
                                              D2_SBP(i,6)*f(l1_ogrid+5) + &
                                              D2_SBP(i,7)*f(l1_ogrid+6) + &
                                              D2_SBP(i,8)*f(l1_ogrid+7) + &
                                              D2_SBP(i,9)*f(l1_ogrid+8) )
      enddo

    endsubroutine der2_ogrid_SBP
!***********************************************************************
    subroutine der_ijm_ogrid_SBP(f,df,k)
! 
!  Cross derivatives for SBP boundary closures 
!  Only implemented in radial direction.
!
!  10-okt-17/Jorgen: Coded
      real, dimension(mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid), intent(in) :: f
      real, dimension(6), intent(out) :: df
      real, dimension(6)  :: df_mn1,df_mn_1,df_mn2,df_mn_2,df_mn3,df_mn_3
      integer, intent(in) :: k

      call der_ogrid_SBP(f(l1_ogrid:l1_ogrid+8,m_ogrid+1,n_ogrid,k),df_mn1)
      call der_ogrid_SBP(f(l1_ogrid:l1_ogrid+8,m_ogrid-1,n_ogrid,k),df_mn_1)
      call der_ogrid_SBP(f(l1_ogrid:l1_ogrid+8,m_ogrid+2,n_ogrid,k),df_mn2)
      call der_ogrid_SBP(f(l1_ogrid:l1_ogrid+8,m_ogrid-2,n_ogrid,k),df_mn_2)
      call der_ogrid_SBP(f(l1_ogrid:l1_ogrid+8,m_ogrid+3,n_ogrid,k),df_mn3)
      call der_ogrid_SBP(f(l1_ogrid:l1_ogrid+8,m_ogrid-3,n_ogrid,k),df_mn_3)

      

     df(1:6)=(1./60.)*dy_1_ogrid(m_ogrid)*( &
                        45.*(df_mn1-df_mn_1) &
                        -9.*(df_mn2-df_mn_2) &
                           +(df_mn3-df_mn_3))

     endsubroutine der_ijm_ogrid_SBP
!*********************************************************************** 
    subroutine der_ijn_ogrid_SBP(f,df,k)
! 
!  Cross derivatives for SBP boundary closures 
!  Only implemented in radial direction.
!
!  10-okt-17/Jorgen: Coded
      real, dimension(mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid), intent(in) :: f
      integer, intent(in) :: k
      real, dimension(6), intent(out) :: df
      real, dimension(6)  :: df_mn1,df_mn_1,df_mn2,df_mn_2,df_mn3,df_mn_3

      call der_ogrid_SBP(f(l1_ogrid:l1_ogrid+8,m_ogrid,n_ogrid+1,k),df_mn1)
      call der_ogrid_SBP(f(l1_ogrid:l1_ogrid+8,m_ogrid,n_ogrid-1,k),df_mn_1)
      call der_ogrid_SBP(f(l1_ogrid:l1_ogrid+8,m_ogrid,n_ogrid+2,k),df_mn2)
      call der_ogrid_SBP(f(l1_ogrid:l1_ogrid+8,m_ogrid,n_ogrid-2,k),df_mn_2)
      call der_ogrid_SBP(f(l1_ogrid:l1_ogrid+8,m_ogrid,n_ogrid+3,k),df_mn3)
      call der_ogrid_SBP(f(l1_ogrid:l1_ogrid+8,m_ogrid,n_ogrid-3,k),df_mn_3)

      
      df(1:6)=(1./60.)*dz_1_ogrid(n_ogrid)*( &
                       45.*(df_mn1-df_mn_1) &
                       -9.*(df_mn2-df_mn_2) &
                          +(df_mn3-df_mn_3))

     endsubroutine der_ijn_ogrid_SBP
!*********************************************************************** 
    subroutine der_ijm_ogrid_bdry5(f,df,k)
! 
!  Cross derivatives for fifth order boundary closures 
!  Only implemented in radial direction.
!
!  13-okt-17/Jorgen: Coded
      real, dimension(mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid), intent(in) :: f
      integer, intent(in) :: k
      real, dimension(3), intent(out) :: df
      real, dimension(3)  :: df_mn1,df_mn_1,df_mn2,df_mn_2,df_mn3,df_mn_3

      call der_ogrid_bdry5(f(l1_ogrid:l1_ogrid+5,m_ogrid+1,n_ogrid,k),df_mn1)
      call der_ogrid_bdry5(f(l1_ogrid:l1_ogrid+5,m_ogrid-1,n_ogrid,k),df_mn_1)
      call der_ogrid_bdry5(f(l1_ogrid:l1_ogrid+5,m_ogrid+2,n_ogrid,k),df_mn2)
      call der_ogrid_bdry5(f(l1_ogrid:l1_ogrid+5,m_ogrid-2,n_ogrid,k),df_mn_2)
      call der_ogrid_bdry5(f(l1_ogrid:l1_ogrid+5,m_ogrid+3,n_ogrid,k),df_mn3)
      call der_ogrid_bdry5(f(l1_ogrid:l1_ogrid+5,m_ogrid-3,n_ogrid,k),df_mn_3)

      
      df(1:3)=(1./60.)*dy_1_ogrid(m_ogrid)*( &
                      45.*(df_mn1-df_mn_1) &
                      -9.*(df_mn2-df_mn_2) &
                         +(df_mn3-df_mn_3))

     endsubroutine der_ijm_ogrid_bdry5
!*********************************************************************** 
    subroutine der_ijn_ogrid_bdry5(f,df,k)
! 
!  Cross derivatives for fifth order boundary closures 
!  Only implemented in radial direction.
!
!  13-okt-17/Jorgen: Coded
      real, dimension(mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid), intent(in) :: f
      integer, intent(in) :: k
      real, dimension(3), intent(out) :: df
      real, dimension(3)  :: df_mn1,df_mn_1,df_mn2,df_mn_2,df_mn3,df_mn_3

      call der_ogrid_bdry5(f(l1_ogrid:l1_ogrid+5,m_ogrid,n_ogrid+1,k),df_mn1)
      call der_ogrid_bdry5(f(l1_ogrid:l1_ogrid+5,m_ogrid,n_ogrid-1,k),df_mn_1)
      call der_ogrid_bdry5(f(l1_ogrid:l1_ogrid+5,m_ogrid,n_ogrid+2,k),df_mn2)
      call der_ogrid_bdry5(f(l1_ogrid:l1_ogrid+5,m_ogrid,n_ogrid-2,k),df_mn_2)
      call der_ogrid_bdry5(f(l1_ogrid:l1_ogrid+5,m_ogrid,n_ogrid+3,k),df_mn3)
      call der_ogrid_bdry5(f(l1_ogrid:l1_ogrid+5,m_ogrid,n_ogrid-3,k),df_mn_3)

      
      df(1:3)=(1./60.)*dz_1_ogrid(n_ogrid)*( &
                       45.*(df_mn1-df_mn_1) &
                       -9.*(df_mn2-df_mn_2) &
                          +(df_mn3-df_mn_3))

     endsubroutine der_ijn_ogrid_bdry5
!*********************************************************************** 
    subroutine der_ogrid_SBP_experimental(f,k,df)
! 
!  Summation by parts boundary condition for first derivative.
!  Only implemented in radial direction.
!  This experimental routine is for the outer boundary
!
!  12-may-17/Jorgen: Coded
!
      real, dimension(9,my_ogrid,mz_ogrid,mfarray_ogrid), intent(in) :: f
      real, dimension(6), intent(out) :: df
      integer, intent(in) :: k
      real :: fac
      integer :: i

      do i=1,6
        fac = dx_1_ogrid(mx_ogrid-i+1)
        df(7-i)=fac*        (D1_SBP(i,1)*f(9,m_ogrid,n_ogrid,k) + &
                             D1_SBP(i,2)*f(8,m_ogrid,n_ogrid,k) + &
                             D1_SBP(i,3)*f(7,m_ogrid,n_ogrid,k) + &
                             D1_SBP(i,4)*f(6,m_ogrid,n_ogrid,k) + &
                             D1_SBP(i,5)*f(5,m_ogrid,n_ogrid,k) + &
                             D1_SBP(i,6)*f(4,m_ogrid,n_ogrid,k) + &
                             D1_SBP(i,7)*f(3,m_ogrid,n_ogrid,k) + &
                             D1_SBP(i,8)*f(2,m_ogrid,n_ogrid,k) + &
                             D1_SBP(i,9)*f(1,m_ogrid,n_ogrid,k) )
      enddo

    endsubroutine der_ogrid_SBP_experimental
!***********************************************************************
    subroutine der2_ogrid_SBP_experimental(f,k,df2)
! 
!  Summation by parts boundary condition for second derivative.
!  Only implemented in radial direction.
!  This experimental routine is for the outer boundary
!
!  12-may-17/Jorgen: Coded
      real, dimension(9,my_ogrid,mz_ogrid,mfarray_ogrid), intent(in) :: f
      real, dimension(6), intent(out) :: df2
      integer, intent(in) :: k
      real :: fac
      integer :: i

      do i=1,6
        fac = dx_1_ogrid(mx_ogrid-i+1)
        df2(7-i)=(fac**2)       *(D2_SBP(i,1)*f(9,m_ogrid,n_ogrid,k) + &
                                  D2_SBP(i,2)*f(8,m_ogrid,n_ogrid,k) + &
                                  D2_SBP(i,3)*f(7,m_ogrid,n_ogrid,k) + &
                                  D2_SBP(i,4)*f(6,m_ogrid,n_ogrid,k) + &
                                  D2_SBP(i,5)*f(5,m_ogrid,n_ogrid,k) + &
                                  D2_SBP(i,6)*f(4,m_ogrid,n_ogrid,k) + &
                                  D2_SBP(i,7)*f(3,m_ogrid,n_ogrid,k) + &
                                  D2_SBP(i,8)*f(2,m_ogrid,n_ogrid,k) + &
                                  D2_SBP(i,9)*f(1,m_ogrid,n_ogrid,k) )
      enddo

    endsubroutine der2_ogrid_SBP_experimental
!*********************************************************************** 
  end module solid_cells_ogrid_sub
