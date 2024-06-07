#include "../shock/kernels.ac"
ComputeSteps AC_rhs_0(boundconds)
{
shock_1_divu(),
  shock_2_smooth(),
        twopass_solve_intermediate(0,AC_dt),
        twopass_solve_final(0)
}
ComputeSteps AC_rhs_1(boundconds)
{
shock_1_divu(),
  shock_2_smooth(),
        twopass_solve_intermediate(1,AC_dt),
        twopass_solve_final(1)
}
ComputeSteps AC_rhs_2(boundconds)
{
shock_1_divu(),
  shock_2_smooth(),
        twopass_solve_intermediate(2,AC_dt),
        twopass_solve_final(2)
}
BoundConds boundconds{
  periodic(BOUNDARY_XYZ)
}
