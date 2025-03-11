#!/bin/csh
set differ=cmp
echo "====================="
echo "xy.h5"
$differ data/averages/xy.h5    reference_averages_xy.h5
echo "====================="
echo "lnTT_xy2.h5"
$differ   data/slices/lnTT_xy2.h5  reference_slices_lnTT_xy2.h5
echo "====================="
echo "lnTT_xy.h5"
$differ   data/slices/lnTT_xy.h5   reference_slices_lnTT_xy.h5
echo "====================="
echo "lnTT_xz.h5"
$differ   data/slices/lnTT_xz.h5   reference_slices_lnTT_xz.h5
echo "====================="
echo "lnTT_yz.h5"
$differ   data/slices/lnTT_yz.h5   reference_slices_lnTT_yz.h5
echo "====================="
echo "rho_xy2.h5"
$differ   data/slices/rho_xy2.h5   reference_slices_rho_xy2.h5
echo "====================="
echo "rho_xy.h5"
$differ   data/slices/rho_xy.h5    reference_slices_rho_xy.h5
echo "====================="
echo "rho_yz.h5"
$differ   data/slices/rho_xz.h5    reference_slices_rho_xz.h5
echo "====================="
echo "rho_yz.h5"
$differ   data/slices/rho_yz.h5    reference_slices_rho_yz.h5
echo "====================="
