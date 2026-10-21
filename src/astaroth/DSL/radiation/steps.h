input int AC_frequency_bin
input int AC_ray_idir
ComputeSteps get_source_function_and_opacity(boundconds)
{
        source_function(AC_frequency_bin)
        opacity(AC_frequency_bin)
}

ComputeSteps Qrad_reset(boundconds)
{
        reset_Qrad()
}

//Intrinsic step: all rays start with zero intensity and optical depth (Qintrinsic)
ComputeSteps Q_up_intrinsic(Qintrinsic_bcs)
{
        integrate_Q_up(AC_ray_idir)
}
ComputeSteps Q_down_intrinsic(Qintrinsic_bcs)
{
        integrate_Q_down(AC_ray_idir)
}
ComputeSteps Q_xp_intrinsic(Qhor_intrinsic_bcs)
{
        integrate_Q_xp(AC_ray_idir)
}
ComputeSteps Q_xm_intrinsic(Qhor_intrinsic_bcs)
{
        integrate_Q_xm(AC_ray_idir)
}
ComputeSteps Q_yp_intrinsic(Qhor_intrinsic_bcs)
{
        integrate_Q_yp(AC_ray_idir)
}
ComputeSteps Q_ym_intrinsic(Qhor_intrinsic_bcs)
{
        integrate_Q_ym(AC_ray_idir)
}

//Revision with the boundary intensity and accumulation into Qrad
ComputeSteps Q_up_extrinsic(boundconds)
{
        revise_and_accumulate_Q_up(AC_ray_idir)
}
ComputeSteps Q_down_extrinsic(boundconds)
{
        revise_and_accumulate_Q_down(AC_ray_idir)
}
ComputeSteps Q_xp_extrinsic(boundconds)
{
        revise_and_accumulate_Q_xp(AC_ray_idir)
}
ComputeSteps Q_xm_extrinsic(boundconds)
{
        revise_and_accumulate_Q_xm(AC_ray_idir)
}
ComputeSteps Q_yp_extrinsic(boundconds)
{
        revise_and_accumulate_Q_yp(AC_ray_idir)
}
ComputeSteps Q_ym_extrinsic(boundconds)
{
        revise_and_accumulate_Q_ym(AC_ray_idir)
}

//Boundary heating rate Q=I-S of the vertical rays in the upstream ghost zone (radboundary_xy_set):
//'0': I=0, 'S': I=S, 'c': I=Irad_z (layer of constant temperature); (1) bottom, (2) top.
rad_intensity_bc_z(AcBoundary boundary, Field Q)
{
        const int3 normal = get_normal(boundary)
        const int3 boundary_point = get_boundary(normal)
        const int ib = (normal.z < 0) ? 0 : 1
        const int ibc = AC_ibc_rad_z__mod__radiation[ib]
        const int3 ghost = boundary_point + normal
        real Q0 = 0.0
        if (ibc == 0) Q0 = -SRAD[ghost.x][ghost.y][ghost.z]
        if (ibc == 2) Q0 = AC_irad_z__mod__radiation[ib]-SRAD[ghost.x][ghost.y][ghost.z]
        Q[ghost.x][ghost.y][ghost.z] = Q0
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

//Zero also in the periodic directions; the periodicity enters in revise_periodic_ray
BoundConds Qhor_intrinsic_bcs
{
  ac_const_bc(BOUNDARY_XY,Q_HOR,0.0)
  ac_const_bc(BOUNDARY_XY,TAU_HOR,0.0)
  ac_fixed_bc(BOUNDARY_XYZ,SRAD)
  ac_fixed_bc(BOUNDARY_XYZ,F_KAPPARHO)
}
