  private
!
  public :: initialize_density_methods, getrho, getlnrho, getdlnrho, getrho_s, getdlnrho_y
  public :: putrho, putlnrho
!
  interface getdlnrho
    module procedure getdlnrho_x
!    module procedure getdlnrho_y
    module procedure getdlnrho_z
  endinterface
