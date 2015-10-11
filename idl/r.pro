; Reads the data from a varfile.
; You should have run 'start.pro' once before.
; $Id$

if ((n_elements(started) le 0) or (n_elements(read_all) gt 0)) then begin
  message, "You need to run start.pro first: use '.rnew start'"
endif
undefine, read_all

default, varfile, 'var.dat'

; Read data
pc_read_var_raw, obj=var, tags=tags, time=t, varfile=varfile, dim=dim, grid=grid, start_param=par, run_param=par2, datadir=datadir

; Shortcuts
var_names = tag_names (tags)
num_tags = n_elements (var_names)
offset = 0
while (offset le num_tags-1) do begin
  indices = tags.(offset)
  if (size (indices, /type) eq 3) then begin
    dummy = execute (var_names[offset]+' = var[*,*,*,indices]')
    components = n_elements (indices)
    if (components gt 1) then offset += components
  end
  offset += 1
end
var = 0
x = grid.x
y = grid.y
z = grid.z
dx = grid.dx
dy = grid.dy
dz = grid.dz
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

if (par.lhydro) then begin
  nu=par2.nu
endif
if (par.ldensity) then begin
  cs0=par2.cs0
endif
if (par.lentropy) then begin
  hcond0=par2.hcond0
  hcond1=par2.hcond1
  hcond2=par2.hcond2
  luminosity=par2.luminosity
  wheat=par2.wheat
  cool=par2.cool
  wcool=par2.wcool
  Fbot=par2.Fbot
endif
if (par.lmagnetic) then begin
  eta=par2.eta
end

END

