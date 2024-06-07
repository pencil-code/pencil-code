#include "../shock/kernels.ac"
ComputeSteps AC_rhs_0(boundconds) {
shock_1_divu(),
   shock_2_smooth(),
  singlepass_solve(0,AC_dt),
}
ComputeSteps AC_rhs_1(boundconds) {
shock_1_divu(),
   shock_2_smooth(),
  singlepass_solve(1,AC_dt),
}
ComputeSteps AC_rhs_2(boundconds) {
shock_1_divu(),
   shock_2_smooth(),
  singlepass_solve(2,AC_dt),
}
BoundConds boundconds{
  periodic(BOUNDARY_XYZ)
}
