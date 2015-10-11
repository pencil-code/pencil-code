; $Id$

;;;;;;;;;;;;;;;
;;;  r.pro  ;;;
;;;;;;;;;;;;;;;

;;; Read the data produced on one processor
;;; You should have run `start.pro' once before.
;;; $Id$

if ((n_elements(started) le 0) or (n_elements(read_all) gt 0)) then begin
  message, "You need to run start.pro first: use `.rnew start'"
endif
undefine, read_all

default, varfile, 'var.dat'

; Read data
pc_read_var, obj=vars, varfile=varfile, dim=dim, grid=grid, param=param, datadir=datadir

; Shortcuts
tags = tag_names (vars)
num_tags = n_elements (tags)
for pos = 0, num_tags-1 do dummy = execute (tags[pos]+' = vars.'+tags[pos])
nx = dim.nx
ny = dim.ny
nz = dim.nz
mx = dim.mx
my = dim.my
mz = dim.mz

xx = spread(x, [1,2], [my,mz])
yy = spread(y, [0,2], [mx,mz])
zz = spread(z, [0,1], [mx,my])
rr = sqrt(xx^2+yy^2+zz^2)

END

