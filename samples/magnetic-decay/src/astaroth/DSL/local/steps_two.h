input real AC_dt
input int AC_step_num
ComputeSteps AC_rhs(boundconds)
{
        twopass_solve_intermediate(AC_step_num,AC_dt)
        twopass_solve_final(AC_step_num)
}
BoundConds boundconds{
  #include "boundconds.h"
}
//TP: periodic in XY sym in Z
//BoundConds boundconds{
//  periodic(BOUNDARY_XY)
//  bc_sym_z(BOUNDARY_Z_TOP, 1.0,AC_top,UUX,false)
//  bc_sym_z(BOUNDARY_Z_TOP, 1.0,AC_top,UUY,false)
//  bc_sym_z(BOUNDARY_Z_TOP, 1.0,AC_top,UUZ,false)
//
//
//  bc_sym_z(BOUNDARY_Z_TOP, 1.0,AC_top,AAX,false)
//  bc_sym_z(BOUNDARY_Z_TOP, 1.0,AC_top,AAY,false)
//  bc_sym_z(BOUNDARY_Z_TOP, 1.0,AC_top,AAZ,false)
//
//  bc_sym_z(BOUNDARY_Z_TOP, 1.0,AC_top,SS,false)
//  bc_sym_z(BOUNDARY_Z_TOP, 1.0,AC_top,RHO,false)
//
//  bc_sym_z(BOUNDARY_Z_BOT, 1.0,AC_bot,UUX,false)
//  bc_sym_z(BOUNDARY_Z_BOT, 1.0,AC_bot,UUY,false)
//  bc_sym_z(BOUNDARY_Z_BOT, 1.0,AC_bot,UUZ,false)
//
//
//  bc_sym_z(BOUNDARY_Z_BOT, 1.0,AC_bot,AAX,false)
//  bc_sym_z(BOUNDARY_Z_BOT, 1.0,AC_bot,AAY,false)
//  bc_sym_z(BOUNDARY_Z_BOT, 1.0,AC_bot,AAZ,false)
//
//  bc_sym_z(BOUNDARY_Z_BOT, 1.0,AC_bot,SS,false)
//  bc_sym_z(BOUNDARY_Z_BOT, 1.0,AC_bot,RHO,false)
//}
////TP example 2
//BoundConds boundconds{
//  periodic(BOUNDARY_XY)
//  bc_sym_z(BOUNDARY_Z_TOP, 1.0,AC_top,UUY,false)
//  bc_sym_z(BOUNDARY_Z_TOP, 1.0,AC_top,UUZ,false)
//  bc_sym_z(BOUNDARY_Z_BOT, 1.0,AC_bot,UUY,false)
//  bc_sym_z(BOUNDARY_Z_BOT, 1.0,AC_bot,UUZ,false)
//
//  bc_sym_z(BOUNDARY_Z_TOP, 1.0,AC_top,AAX,false)
//  bc_sym_z(BOUNDARY_Z_TOP, 1.0,AC_top,AAY,false)
//  bc_sym_z(BOUNDARY_Z_TOP, 1.0,AC_top,AAZ,false)
//
//  bc_ss_flux(BOUNDARY_Z_TOP, AC_top)
//
//
//
//  bc_sym_z(BOUNDARY_Z_BOT, 1.0,AC_bot,AAX,false)
//  bc_sym_z(BOUNDARY_Z_BOT, 1.0,AC_bot,AAY,false)
//  bc_sym_z(BOUNDARY_Z_BOT, 1.0,AC_bot,AAZ,false)
//
//  bc_sym_z(BOUNDARY_Z_BOT, 1.0,AC_bot,SS,false)
//
//  bc_steady_z(BOUNDARY_Z_TOP, AC_top,UUX)
//  bc_steady_z(BOUNDARY_Z_BOT, AC_bot,UUX)
//
//  bc_ism(BOUNDARY_Z_BOT, AC_bot,RHO)
//  bc_ism(BOUNDARY_Z_TOP, AC_top,RHO)
//
//}

