  module Geometrical_types
!
! $Id$
!
! MODULE_DOC: Collection of geometrical object types.
! MODULE_DOC: (Presently only rectangular toroid)
!
! 16-May-20/MR: coded
!
  implicit none
!
  private
!
  public :: torus_rect, torus_init
  public :: torus_precess, torus_wobble
  public :: torus_extend_r, torus_extend_z
  public :: torus_rect_unfmt_read, torus_rect_unfmt_write
!
  type torus_rect
    real, dimension(3) :: center
    real :: th=0.,ph=0.
    real :: r_in=0., thick=0., height=0.
    real :: Omega_prec, extr_rate, extz_rate
    real, dimension(3) :: wob_amp, wob_om, wob_phase

    real, dimension(3) :: center0
    real :: th0,ph0,r_in0,height0
    contains
    !  procedure, pass :: precess => my_precess
    procedure :: torus_rect_unfmt_write
    procedure :: torus_rect_unfmt_read
    generic :: write (unformatted) => torus_rect_unfmt_write
    generic :: read (unformatted) => torus_rect_unfmt_read

  endtype torus_rect

  contains
!
!***********************************************************************
    subroutine torus_init(torus)

    type(torus_rect) :: torus

    torus%center0=torus%center
    torus%th0    =torus%th
    torus%ph0    =torus%ph
    torus%r_in0  =torus%r_in
    torus%height0=torus%height

    endsubroutine torus_init
!***********************************************************************
    subroutine torus_precess(torus,t)

    type(torus_rect) :: torus
    double precision :: t

    if (torus%Omega_prec==0.) return
    torus%ph=torus%ph0+torus%Omega_prec*t

    endsubroutine torus_precess
!***********************************************************************
    subroutine torus_wobble(torus,t)

    type(torus_rect) :: torus
    double precision :: t

    if (all(torus%wob_om.eq.0.)) return

    torus%center=torus%center0+torus%wob_amp*cos(torus%wob_om*t+torus%wob_phase)

    endsubroutine torus_wobble
!***********************************************************************
    subroutine torus_extend_r(torus,t)

    type(torus_rect) :: torus
    double precision :: t

    if (torus%extr_rate.eq.0.) return

    torus%r_in=max(torus%r_in0*(1.+torus%extr_rate*t),0.d0)

    endsubroutine torus_extend_r
!***********************************************************************
    subroutine torus_extend_z(torus,t)

    type(torus_rect) :: torus
    double precision :: t

    if (torus%extz_rate.eq.0.) return

    torus%height=max(torus%height0*(1.+torus%extz_rate*t),0.d0)

    endsubroutine torus_extend_z
!***********************************************************************
    subroutine torus_rect_unfmt_write(torus,unit,iostat,iomsg)

    class(torus_rect),intent(IN)   :: torus
    integer,          intent(IN)   :: unit
    integer,          intent(OUT)  :: iostat
    character(LEN=*), intent(INOUT):: iomsg

    write(unit=unit,iostat=iostat,iomsg=iomsg) torus%center, &
                                               torus%th,torus%ph, &
                                               torus%r_in, torus%thick, &
                                               torus%height
    endsubroutine torus_rect_unfmt_write
!***********************************************************************
    subroutine torus_rect_unfmt_read(torus,unit,iostat,iomsg)

    class(torus_rect),intent(INOUT):: torus
    integer,          intent(IN)   :: unit
    integer,          intent(OUT)  :: iostat
    character(LEN=*), intent(INOUT):: iomsg

    read(unit=unit,iostat=iostat,iomsg=iomsg) torus%center, &
                                              torus%th,torus%ph, &
                                              torus%r_in, torus%thick, &
                                              torus%height
    endsubroutine torus_rect_unfmt_read
!***********************************************************************
  endmodule Geometrical_types

