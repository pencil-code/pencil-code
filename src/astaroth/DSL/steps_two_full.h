#include "../shock/kernels.ac"
#include "../density/mass_conservation.h"
#include "../selfgravity.h"
#include "../magnetic/before_boundary.h"
#include "../alphadisk/after_timestep.h"
#include "../gravitational_waves.h"
#include "../hydro/hydro_after_boundary_conservative.h"
#include "../ioncalc.h"
#include "../newton_cooling.h"
#include "../polymer.h"

input real AC_dt
input PC_SUB_STEP_NUMBER AC_step_num
input bool AC_lrmv
input real AC_t

ComputeSteps AC_rhs(boundconds)
{
	shock_1_divu(AC_step_num)
	shock_2_max(AC_step_num)
	shock_3_smooth(AC_step_num)
        twopass_solve_intermediate(AC_step_num,AC_dt,AC_t,AC_lrmv)
        twopass_solve_final(AC_step_num)
}
ComputeSteps AC_calculate_timestep(boundconds)
{
	shock_1_divu(AC_step_num)
	shock_2_max(AC_step_num)
	shock_3_smooth(AC_step_num)
	twopass_solve_intermediate(PC_FIRST_SUB_STEP,AC_dt,AC_t,AC_lrmv)
}

ComputeSteps AC_calc_selfgravity_rhs(boundconds)
{
	selfgravity_calc_rhs()
}
ComputeSteps AC_calc_final_potential(boundconds)
{
	calc_final_potential(AC_t)
}

ComputeSteps AC_sor_step(boundconds)
{
	selfgravity_sor_step(0)
	selfgravity_sor_step(1)
}

ComputeSteps AC_before_boundary_steps(boundconds)
{
	get_current_total_mass(AC_lrmv)
	fix_mass_drift(AC_lrmv)
	magnetic_before_boundary_reductions()
	hydro_before_boundary(AC_step_num)
	hydro_after_boundary_conservative(AC_t)
	calc_axion_integral(AC_t)
	ioncalc()
	prep_ode_right()
	calc_poly_fr()
}
ComputeSteps AC_after_timestep(boundconds)
{
	after_timestep_alphadisk()
}
ComputeSteps AC_gravitational_waves_solve_and_stress(boundconds)
{
	gravitational_waves_solve_and_stress(AC_t,AC_dt)
}

ComputeSteps AC_integrate_tau(boundconds)
{
  calc_kappar_and_dtau()
  integrate_tau_up()
  integrate_tau_down()
  calc_tau()
}
BoundConds boundconds{
  #include "boundconds.h"
#if LNEWTON_COOLING
  ac_const_bc(BOUNDARY_Y_BOT,TAU_BELOW,0.0)
  ac_const_bc(BOUNDARY_Y_TOP,TAU_ABOVE,0.0)
#endif
}
//TP: periodic in XY sym in Z
//BoundConds boundconds{
//  periodic(BOUNDARY_XY)
//  bc_sym_z(BOUNDARY_Z_TOP, 1.0,AC_top,UUX,false)
//  bc_sym_z(BOUNDARY_Z_TOP, 1.0,AC_top,UUY,false)
//  bc_sym_z(BOUNDARY_Z_TOP, 1.0,AC_top,UUZ,false)
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
//}

