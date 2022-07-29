# der.py
#
"""
This is a wrapper for the actual derivatives, currently 6th order with
ghost zones included (pencil-code style).
"""


from pencil.math.derivatives.der_6th_order_w_ghosts import (
    xder5_6th,
    yder5_6th,
    zder5_6th,
)
from pencil.math.derivatives.der_4th_order_w_ghosts import (
    xder3_4th,
    yder3_4th,
    zder3_4th,
)
from pencil.math.derivatives.der_nonequi import (
    xder_6th,
    yder_6th,
    zder_6th,
    xder2_6th,
    yder2_6th,
    zder2_6th,
    xder6_2nd,
    yder6_2nd,
    zder6_2nd,
)

xder = xder_6th
yder = yder_6th
zder = zder_6th

xder2 = xder2_6th
yder2 = yder2_6th
zder2 = zder2_6th

xder3 = xder3_4th
yder3 = yder3_4th
zder3 = zder3_4th

xder5 = xder5_6th
yder5 = yder5_6th
zder5 = zder5_6th

xder6 = xder6_2nd
yder6 = yder6_2nd
zder6 = zder6_2nd
