! $Id: timestep.f90 9840 2008-09-05 07:29:37Z ajohan $
!
!  This module add solid (as in no-fluid) cells in the domain.
!  This can be used e.g. in order to simulate a cylinder in a cross flow.
!
module Solid_Cells

  use Cparam
  use Cdata
  
  implicit none
  
  include 'solid_cells.h'

  integer, parameter           :: max_items=10
  integer                      :: ncylinders,nrectangles
  real, dimension(max_items,4) :: cylinder
  real, dimension(max_items,6) :: rectangle
  integer, dimension(mx,my,mz,3)  :: ba,ba_shift
  
  contains
!***********************************************************************
    subroutine initialize_solid_cells
!
!  DOCUMENT ME!
!
!  19-nov-2008/nils: coded
!
      lsolid_cells=.true.
!
!  Define the geometry of the solid object.
!  This shold probably be included as a geometry.local file such that
!  one can define complex geometries on a case to case basis.
!  Alternatively one will here end up with a terribly long series
!  of case checks.
!
      cylinder(1,:)=(/16.85e-3,0,0,0/)
      ncylinders=1
!
      call find_solid_cell_boundaries
      call calculate_shift_matrix
!
    endsubroutine initialize_solid_cells
!***********************************************************************  
    subroutine update_solid_cells(f)
!
!  Set the boundary values of the solid area such that we get a 
!  correct fluid-solid interface.
!
!  19-nov-2008/nils: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i,j,k,idir,xind,yind,zind
!
      do i=l1,l2
        do j=m1,m2
          do k=n1,n2
            do idir=1,3
              if (ba_shift(i,j,k,idir).ne.0) then
                xind=i
                yind=j
                zind=k
!
                if (idir==1) then
                  xind=i-ba_shift(i,j,k,idir)
                elseif (idir==2) then
                  yind=j-ba_shift(i,j,k,idir)
                elseif (idir==3) then
                  zind=k-ba_shift(i,j,k,idir)
                else
                  print*,'No such idir!...exiting!'
                  stop
                endif
                f(i,j,k,iux:iuz)=-f(xind,yind,zind,iux:iuz)
                if (ilnrho>0) f(i,j,k,ilnrho) = f(xind,yind,zind,ilnrho)
              endif
            enddo
          enddo
        enddo
      enddo
!
    endsubroutine update_solid_cells
!***********************************************************************  
    subroutine freeze_solid_cells(df)
!
!  If we are in a solid cell set df=0 for all variables
!
!  19-nov-2008/nils: coded
!
      real, dimension (mx,my,mz,mvar) :: df
      integer :: i,j,k
! 
      do i=l1,l2
        if (&
             (ba(i,m,n,1).ne.0).or.&
             (ba(i,m,n,2).ne.0).or.&
             (ba(i,m,n,3).ne.0)) then
          df(i,m,n,:)=0
        endif
      enddo
!
    endsubroutine freeze_solid_cells
!***********************************************************************
    subroutine find_solid_cell_boundaries
!
!  Find the boundaries of the geometries such that we can set the
!  ghost points inside the solid geometry in order to achieve the
!  correct no-slip boundaries.
!
!  19-nov-2008/nils: coded
!
      integer :: i,j,k,icyl
      real :: x2,y2,xval_p,xval_m,yval_p,yval_m
!
!  Initialize ba
!
      ba=0
!
!  Loop over all cylinders (this should actually be a loop over all
!  geometries!)
!
      do icyl=1,ncylinders
!
!  First we look in x-direction
!
        do j=m1,m2
          x2=cylinder(icyl,1)**2-(y(j)-cylinder(icyl,3))**2
          if (x2>0) then
            xval_p=cylinder(icyl,2)+sqrt(x2)
            xval_m=cylinder(icyl,2)-sqrt(x2)            
            do i=l1,l2
              if (x(i)<xval_p .and. x(i)>xval_m) then
                if (x(i+1)>xval_p) then
                  ba(i,j,:,1)=-1
                elseif (x(i+2)>xval_p) then
                  ba(i,j,:,1)=-2
                elseif (x(i+3)>xval_p) then
                  ba(i,j,:,1)=-3
                elseif (x(i-1)<xval_m) then
                  ba(i,j,:,1)=1
                elseif (x(i-2)<xval_m) then
                  ba(i,j,:,1)=2
                elseif (x(i-3)<xval_m) then
                  ba(i,j,:,1)=3
                else
                  ba(i,j,:,1)=9
                endif
              endif
            enddo
          endif
        enddo
!
!  Then we look in y-direction
!
        do i=l1,l2
          y2=cylinder(icyl,1)**2-(x(i)-cylinder(icyl,2))**2
          if (y2>0) then
            yval_p=cylinder(icyl,3)+sqrt(y2)
            yval_m=cylinder(icyl,3)-sqrt(y2)            
            do j=m1,m2
              if (y(j)<yval_p .and. y(j)>yval_m) then
                if (y(j+1)>yval_p) then
                  ba(i,j,:,2)=-1
                elseif (y(j+2)>yval_p) then
                  ba(i,j,:,2)=-2
                elseif (y(j+3)>yval_p) then
                  ba(i,j,:,2)=-3
                elseif (y(j-1)<yval_m) then
                  ba(i,j,:,2)=1
                elseif (y(j-2)<yval_m) then
                  ba(i,j,:,2)=2
                elseif (y(j-3)<yval_m) then
                  ba(i,j,:,2)=3
                else
                  ba(i,j,:,2)=9
                endif
              endif
            enddo
          endif
        enddo
      enddo
!
    endsubroutine find_solid_cell_boundaries
!***********************************************************************
    subroutine calculate_shift_matrix
!
!  Set up the shift matrix
!
!  19-nov-2008/nils: coded
!
      integer :: i,j,k,idir
      integer :: sgn
!
      ba_shift=0
!
      do i=l1,l2
        do j=m1,m2
          do k=n1,n2
            do idir=1,3
!
!  If ba is non-zero find the shift matrix
!
              if (ba(i,j,k,idir).ne.0 .and. ba(i,j,k,idir).ne.9) then
                sgn=-ba(i,j,k,idir)/abs(ba(i,j,k,idir))
                ba_shift(i,j,k,idir)=2*ba(i,j,k,idir)+sgn
              endif
            enddo
          enddo
        enddo
      enddo
!
    endsubroutine calculate_shift_matrix
!***********************************************************************  
  endmodule Solid_Cells
