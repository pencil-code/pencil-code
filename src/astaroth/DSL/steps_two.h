ComputeSteps AC_rhs_0(boundconds)
{
        twopass_solve_intermediate(0,AC_dt),
        twopass_solve_final(0)
}
ComputeSteps AC_rhs_1(boundconds)
{
        twopass_solve_intermediate(1,AC_dt),
        twopass_solve_final(1)
}
ComputeSteps AC_rhs_2(boundconds)
{
        twopass_solve_intermediate(2,AC_dt),
        twopass_solve_final(2)
}
BoundConds boundconds{
  periodic(BOUNDARY_XYZ)
}
