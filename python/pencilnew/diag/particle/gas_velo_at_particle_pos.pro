;###############################################################
;#
;# calculate the gas velocity at particle position
;# also stored are: paticle position, velocities, gas velocities and gird cell indexes
;#
;###############################################################

pro gas_velo_at_particle_pos, datadir=datadir, destination=destination, $
          filename_prefix=filename_prefix, doforthelastNvar=doforthelastNvar, overwrite=overwrite

;#
;# setup and default
if (not keyword_set(doforthelastNvar)) then doforthelastNvar=4
if (not keyword_set(destination)) then destination='gas_velo_at_particle_pos'
if (not keyword_set(filename_prefix)) then filename_prefix=destination
if (not keyword_set(overwrite)) then overwrite=0

;#
;# Preparation
pen_dir = FILEPATH('pc', ROOT_DIR=datadir)
destination = FILEPATH(destination, ROOT_DIR=pen_dir)
FILE_MKDIR, destination

;# get varlist, eq. to pvarlist and count number of vars
varlist = file_basename(file_Search(FILEPATH('proc0', ROOT_DIR=datadir),'VAR*'))
varlist = strsplit(varlist,'VAR',/EXTRACT)
varlist = uint(sort(varlist.toarray()))   
numvar = size(varlist,/N_ELEMENTS)

print, "~~  gas velocities at particles position for "+str(numvar)+" snapshots."

;# get particle number
print, "~~ reading pdim"
pc_read_pdim, obj=pdim, datadir=datadir, /quiet
numPar = pdim.npar

;# read grid
print, "~~ reading grid"
pc_read_grid, obj=grid, datadir=datadir, /quiet

;# get fast access on grid quantities
dx = grid.DX
dy = grid.DY
dz = grid.DZ
Lx = grid.LX
Ly = grid.LY
Lz = grid.LZ

;# read dim
print, "~~ reading dim"
pc_read_dim, obj=dim, /quiet, datadir=datadir

;# Prepare storage object
gas_velo_at_particle_pos = CREATE_STRUCT('time', 	0.0D, $
					 'par_pos', 	dblarr(numPar,3), $
					 'par_velo', 	dblarr(numPar,3), $
					 'par_idx', 	dblarr(numPar,3), $
					 'npar', 	dblarr(numPar), $
					 'gas_velo', 	dblarr(numPar,3) ) 

;#
;# MAIN CALCULATIONS
FOR i=numvar-doforthelastNvar-1, numvar-1 DO BEGIN
  print, ' '
  print, "####### "+str(i)+"/"+str(numvar-1)+" ########"
  
  filename = FILEPATH(filename_prefix+'_'+str(i)+'.sav', ROOT_DIR=destination) 
  
  ;# check if this file already exists
  IF file_test(filename) EQ 1 && overwrite EQ 0 THEN BEGIN
    print, "~~ "+filename+".sav already existing!"
  ENDIF ELSE BEGIN

    ;# read needed PVAR and VAR file
    print, "~~ reading PVAR file PVAR"+strtrim(string(varlist[i]),1)
    pc_read_pvar, obj=pp, datadir=datadir, varfile="PVAR"+strtrim(string(varlist[i]),1), /quiet
    print, "~~ reading VAR file VAR"+strtrim(string(varlist[i]),1)
    pc_read_var, obj=ff, datadir=datadir, varfile="VAR"+strtrim(string(varlist[i]),1);#, variables=["uu"], /quiet, /magic
    print, ff.uu[0:2, 0:2, 3, *]
    
    ;# store known values to gas_velo_at_particle_pos
    print, "~~ Storing known information into gas_velo_at_particle_pos: snapshot time, particle positions, particle velocities, particle index"
    gas_velo_at_particle_pos.time = pp.t
    gas_velo_at_particle_pos.par_pos = [pp.xx[*,0], pp.xx[*,1], pp.xx[*,2]]
    gas_velo_at_particle_pos.par_velo = [pp.vv[*,0], pp.vv[*,1], pp.vv[*,2]]
    gas_velo_at_particle_pos.par_idx = [floor((pp.xx[*,0]+Lx/2)/dx), $
					floor((pp.xx[*,1]+Ly/2)/dy), $
					floor((pp.xx[*,2]+Lz/2)/dz)]
    np = pc_noghost(ff.np, dim=dim)
    gas_velo_at_particle_pos.npar = np[gas_velo_at_particle_pos.par_idx[*,0], $
				      gas_velo_at_particle_pos.par_idx[*,1], $
				      gas_velo_at_particle_pos.par_idx[*,2]]
    
    ;# prepare grid fix for 2d runs!
    IF dim.NX EQ 1 THEN BEGIN
      print, "~~ detected 2d run, need to prepare fake X grid for pc_gas_velocity_at_particle"
      gridx = indgen(ceil(n_elements(grid.x)/2.))
      gridx = [reverse(-gridx[1:-1]),[0],gridx[1:-1]]
    ENDIF ELSE BEGIN
      gridx = grid.x
    ENDELSE
    
    IF dim.NY EQ 1 THEN BEGIN
      print, "~~ detected 2d run, need to prepare fake Y grid for pc_gas_velocity_at_particle"
      gridy = indgen(ceil(n_elements(grid.y)/2.))
      gridy = [reverse(-gridy[1:-1]),[0],gridy[1:-1]]
    ENDIF ELSE BEGIN
      gridy = grid.y
    ENDELSE
    
    IF dim.NZ EQ 1 THEN BEGIN
      print, "~~ detected 2d run, need to prepare fake Z grid for pc_gas_velocity_at_particle"
      gridz = indgen(ceil(n_elements(grid.z)/2.))
      gridz = [reverse(-gridz[1:-1]),[0],gridz[1:-1]]
    ENDIF ELSE BEGIN
      gridz = grid.z
    ENDELSE
    
;#     print, 'X:'
;#     print, gridx
;#     print, 'Y:'
;#     print, gridy
;#     print, 'Z:'
;#     print, gridz
    
    ;# calculating gas velocity at particle position via AJ routine
    print, "~~ calculating gas velocity at particle position" 
    gas_velo_vec = pc_gas_velocity_at_particle(pp.xx, ff.uu, gridx, gridy, gridz, $
						tsc=1, datadir=datadir, quiet=1, dim=dim)
    gas_velo_at_particle_pos.gas_velo = [gas_velo_vec[*,0], gas_velo_vec[*,1], gas_velo_vec[*,2]]
    
    ;# save gas_velo_at_particle_pos structure
    print, '~~ Saving gas_velo_at_particle_pos object under '+filename
    save, gas_velo_at_particle_pos, filename=filename
    
  ENDELSE
ENDFOR
    
print, ' '
print, '~~ Done! Done! Done! :)'
print, ' '

END