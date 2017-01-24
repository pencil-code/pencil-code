!
!
!   New solid cells module
!   O-grid around cylinder or sphere, coupled to cartesian
!   grid outside the o-grid by interpolation.
!
!   Very fine resolution of boundary layer
!
!
module Solid_Cells_O_Grid

  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  contains

    type :: ogrid
      real, dimension(3) :: x
      real :: r_comp, r_inter
      real :: dr, dtheta
    endtype ogrid 

      
!***********************************************************************
subroutine initialize_grid()
end subroutine initialize_grid
!***********************************************************************
subroutine initialize_solid()
end subroutine initialize_solid()
!***********************************************************************
subroutine recive_flow_info()
end subroutine recive_flow_info
!***********************************************************************
subroutine interpolate_cartesian_to_ogrid()
end subroutine interpolate_cartesian_to_ogrid
!***********************************************************************
subroutine timestep_solid()
end subroutine timestep_solid
!***********************************************************************
subroutine resolve_boundary()
end subroutine resolve_boundary
!***********************************************************************
subroutine resolve_flow()
end subroutine resolve_flow
!***********************************************************************
subroutine interpolate_ogrid_to_cartesian()
end subroutine interpolate_ogrid_to_cartesian
!***********************************************************************
subroutine send_flow_info()
end subroutine send_flow_info()
!***********************************************************************
subroutine compute_draglift()
end subroutine compute_draglift
!***********************************************************************
subroutine print_solid()
end subroutine print_solid
!***********************************************************************

! Particles could check first if thery are within the o-grids computational domain (closer than interpolation domain)
! when the particle positions are interpolated. If this is the case they should use the o-grid for interpolation,
! but if they are outside, they should use the cartesian grid
! Need a variable R_grid for this, that is, the radius of the o-grid. Mabye a different name than R_grid.

logical function within_ogrid_comp(xxp)
!
! Check if current possition is within the o-grids computational domain
!
  use Cdata ! or something, need to get solid and ogrid from somewhere
  real, dimension(3), intent(in) :: xxp
  
  within_ogrid_comp = .false.

  if(sum((xxp(1:2)-solid%x(1:2))**2) < ogrid%r_comp) then
    within_ogrid_comp = .true.
  endif
  
endfunction within_ogrid_comp

