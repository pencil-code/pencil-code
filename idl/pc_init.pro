;  Compile the derivative routines for data that have ghost zones
;  For analysis purposes, one may want to use other routines (called
;  xder_6th, yder_6th, ..., zder2_6th in Axel's idl/pro repo).
;
@xder_6th_ghost
@yder_6th_ghost
@zder_6th_ghost
@xder2_6th_ghost
@yder2_6th_ghost
@zder2_6th_ghost
@xderyder_6th_ghost
@xderzder_6th_ghost
@yderzder_6th_ghost
@xder6_6th_ghost
@yder6_6th_ghost
@zder6_6th_ghost
;
;  The following avoids a mysterious bug when using esrg_legend later
;  (box size was wrong, because lenstr(['a','b']) would be wrong,
;  because xyout would write all letters onto one posision) ..
;
@lenstr
;
pro pc_init
  print, "Compiled default derivative routines - 6th order finite differences with ghost zones"
  print, "Starting Pencil-Code Session"
end
