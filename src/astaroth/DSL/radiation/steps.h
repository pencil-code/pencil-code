input int AC_frequency_bin
ComputeSteps get_source_function_and_opacity(boundconds)
{
        source_function(AC_frequency_bin)
        opacity(AC_frequency_bin)
}

ComputeSteps Qintrinsic_steps(boundconds)
{
        integrate_Q_up()
        integrate_Q_down()
}
ComputeSteps Qextrinsic_steps(boundconds)
{
        revise_Q_up()
        revise_Q_down()
        sum_up_rays()
}

source_function_as_intensity(AcBoundary boundary, Field Q)
{
        const int3 normal = get_normal(boundary)
        const int3 boundary_point = get_boundary(normal)
        int3 domain = boundary_point
        int3 ghost  = boundary_point
        for i in 0:1
        {
                domain = domain - normal
                ghost  = ghost  + normal
                Q[ghost.x][ghost.y][ghost.z] = -SRAD[ghost.x][ghost.y][ghost.z];
        }
}
BoundConds Qintrinsic_bcs
{
  ac_const_bc(BOUNDARY_Z,Q_UP,0.0)
  ac_const_bc(BOUNDARY_Z,Q_DOWN,0.0)
  ac_const_bc(BOUNDARY_Z,TAU_UP,0.0)
  ac_const_bc(BOUNDARY_Z,TAU_DOWN,0.0)
  ac_fixed_bc(BOUNDARY_XYZ,SRAD)
  ac_fixed_bc(BOUNDARY_XYZ,F_KAPPARHO)
}
