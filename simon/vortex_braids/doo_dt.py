# doo_dt.py
"""
Compare the terms from domega_dt.
"""

import os as os
import numpy as np
import pencilnew as pn
import time
xder = pn.math.derivatives.xder
yder = pn.math.derivatives.yder
zder = pn.math.derivatives.zder

run_dir = 'n256_tw2_dz16_oo01_nu4e-5_cl_peri_peri_highTimeRes'

Omega = 0.1
#Omega = 0.01

# Prepare the integrated quantities.
curl_uu_x_oo = []
curl_uu_x_oo_0 = []
nu_Lap_oo = []
t = []
oo_list = []

file_idx = 0
while True:
    print('compute VAR{0}'.format(file_idx))

    try:
        var = pn.read.var(datadir=os.path.join(run_dir, 'data'), var_file='VAR{0}'.format(file_idx))
        file_idx += 1
    except:
        file_idx += 1
        if file_idx < 21:
            continue
        else:
            print('Stopping the calculation.')
            break

    params = pn.read.param(param2=True, datadir=os.path.join(run_dir, 'data'))
    nu = params.nu

    oo = pn.math.derivatives.curl(var.uu, var.dx, var.dy, var.dz, x=var.x, coordinate_system='cylindrical')

    oo_0 = np.zeros_like(oo)
    oo_0[2] = Omega
    uu_x_oo = pn.math.cross(var.uu, oo)
    curl_uu_x_oo.append(pn.math.derivatives.curl(uu_x_oo, var.dx, var.dy, var.dz, x=var.x, coordinate_system='cylindrical')[:, 3:-3, 3:-3, 3:-3])
    uu_x_oo_0 = pn.math.cross(var.uu, oo_0)
    curl_uu_x_oo_0.append(pn.math.derivatives.curl(uu_x_oo_0, var.dx, var.dy, var.dz, x=var.x, coordinate_system='cylindrical')[:, 3:-3, 3:-3, 3:-3])
    nu_Lap_oo.append(-nu*pn.math.derivatives.curl3(var.uu, var.dx, var.dy, var.dz, x=var.x, coordinate_system='cylindrical')[:, 3:-3, 3:-3, 3:-3])

    oo_list.append(oo[:, 3:-3, 3:-3, 3:-3])

    t.append(var.t)

curl_uu_x_oo = np.array(curl_uu_x_oo)
curl_uu_x_oo_0 = np.array(curl_uu_x_oo_0)
nu_Lap_oo = np.array(nu_Lap_oo)
oo = np.array(oo_list)
t = np.array(t)

doo_dt = (oo[2:] - oo[:-2])/(t[2:, np.newaxis, np.newaxis, np.newaxis, np.newaxis] - t[:-2, np.newaxis, np.newaxis, np.newaxis, np.newaxis])
doo_dt_rhs = curl_uu_x_oo + 2*curl_uu_x_oo_0 + nu_Lap_oo

