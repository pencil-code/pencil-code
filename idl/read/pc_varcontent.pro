;
;  $Id$
;
; VARCONTENT STRUCTURE DESCRIPTION
;
; variable (string)
;   Human readable name for the variable
;
; idlvar (string)
;   Name of the variable (usually in the IDL global namespace)
;   in which the variables data will be stored
;
; idlinit (string)
;   IDL command to initialise the storage variable ready for
;   reading in from a file
;
; idlvarloc (string)
;   As idlvar but used when two grid sizes are used eg. global mesh
;   and processor mesh (local -> loc). And used in processes such
;   as in rall.pro.  Eg. uses mesh sizes of (mxloc,myloc,mzloc)
;
; idlinitloc (string)
;   Again as idlinit but used when two mesh sizes are required at once.
;   see idlvarloc
;
function pc_varcontent, datadir=datadir, dim=dim, ivar=ivar, param=param, par2=param2, $
    run2D=run2D, scalar=scalar, noaux=noaux, quiet=quiet
COMPILE_OPT IDL2,HIDDEN
;
;  Read grid dimensions, input parameters and location of datadir.
;
if (not keyword_set(datadir)) then datadir = pc_get_datadir()
if (n_elements(dim) eq 0) then pc_read_dim, obj=dim, datadir=datadir, quiet=quiet
if (n_elements(ivar) eq 0) then ivar=-1
if (n_elements(param) eq 0) then pc_read_param, obj=param, datadir=datadir, dim=dim, quiet=quiet
if (n_elements(par2) eq 0) then pc_read_param, param2=param2, datadir=datadir, dim=dim, quiet=quiet
default, noaux, 0

; 
;  Read the positions of variables in the f-array from the index file.
;
indices_file = datadir+'/index.pro'
num_lines = file_lines (indices_file)
index_pro = strarr (num_lines)
openr, lun, indices_file, /get_lun
readf, lun, index_pro
close, lun
free_lun, lun

; Read in and make accessible all 'nXXX' variables.
default, ntestfield, 0
default, ntestflow, 0
default, ntestscalar, 0
for line = 1, num_lines do begin
  str = stregex (index_pro[line-1], '^ *n[^= ]+[= ]+[0-9]+ *$', /extract)
  if (not execute (str)) then $
      message, 'pc_varcontent: there was a problem with "'+indices_file+'" at line '+str (line)+'.', /info
end

;
;  For EVERY POSSIBLE variable in a snapshot file, store a
;  description of the variable in an indexed array of structures
;  where the indexes line up with those in the saved f array.
;
;  Note: auxiliary variables should go to the table below the folling one.
;
indices = [ $
  { name:'iuu', label:'Velocity', dims:3 }, $
  { name:'ipp', label:'Pressure', dims:1 }, $
  { name:'ippp', label:'Pressure as auxiliary variable', dims:1 }, $
  { name:'iss', label:'Entropy', dims:1 }, $
  { name:'isss', label:'Entropy as auxiliary variable', dims:1 }, $
  { name:'icp', label:'Specific heat as auxiliary variable', dims:1 }, $
  { name:'icv', label:'Specific heat as auxiliary variable', dims:1 }, $
  { name:'igamma', label:'Ratio of specific heat as auxiliary variable', dims:1 }, $
  { name:'inabad', label:'adiabatic logarithmic temperature gradient as auxiliary variable', dims:1 }, $
  { name:'ics', label:'Sound speed as auxiliary variable', dims:1 }, $
  { name:'ilnrho', label:'Log density', dims:1 }, $
  { name:'irho', label:'Density', dims:1 }, $
  { name:'irho_b', label:'Base density', dims:1 }, $
  { name:'irhs', label:'RHS', dims:3 }, $
  { name:'iss_b', label:'Base Entropy', dims:1 }, $
  { name:'iaa', label:'Magnetic vector potential', dims:3 }, $
  { name:'iaphi', label:'A_phi', dims:1 }, $
  { name:'ibphi', label:'B_phi', dims:1 }, $
  { name:'ibb', label:'Magnetic field', dims:3 }, $
  { name:'ijj', label:'Current density', dims:3 }, $
  { name:'iemf', label:'Current density', dims:3 }, $
  { name:'ip11', label:'Polymer Tensor 11', dims:1 }, $
  { name:'ip12', label:'Polymer Tensor 12', dims:1 }, $
  { name:'ip13', label:'Polymer Tensor 13', dims:1 }, $
  { name:'ip22', label:'Polymer Tensor 22', dims:1 }, $
  { name:'ip23', label:'Polymer Tensor 23', dims:1 }, $
  { name:'ip33', label:'Polymer Tensor 33', dims:1 }, $
  { name:'igij', label:'Gravitational Metric', dims:3 }, $
  { name:'iuut', label:'Integrated velocity', dims:3 }, $
  { name:'iaatest', label:'Testfield vector potential', dims:ntestfield }, $
  { name:'iuutest', label:'Testflow', dims:ntestflow }, $
  { name:'icctest', label:'Testflow', dims:ntestscalar }, $
;  { name:'iuxb', label:'Testfield vector potential', dims:ntestfield }, $  ; is this art or can it be removed?
  { name:'iuun', label:'Velocity of neutrals', dims:3 }, $
  { name:'ispitzer', label:'Heat flux vector according to Spitzer', dims:3 }, $
  { name:'ilnrhon', label:'Log density of neutrals', dims:1 }, $
  { name:'ifx', label:'Radiation vector', dims:3 }, $
  { name:'ie', label:'Radiation scalar', dims:1 }, $
  { name:'icc', label:'Passive scalar', dims:1 }, $
  { name:'ilncc', label:'Log passive scalar', dims:1 }, $
  { name:'iXX_chiral', label:'XX chiral', dims:1 }, $
  { name:'iYY_chiral', label:'YY chiral', dims:1 }, $
  { name:'ispecial', label:'Special', dims:1 }, $
  { name:'ispec_3vec', label:'Special vector', dims:3 }, $
  { name:'iphi', label:'Electric potential', dims:1 }, $
  { name:'iLam', label:'Gauge potential', dims:1 }, $
  { name:'iecr', label:'Cosmic ray energy density', dims:1 }, $
  { name:'ifcr', label:'Cosmic ray energy flux', dims:3 }, $
  { name:'igtheta5', label:'Chemical potential gradient', dims:3 }, $
  { name:'itheta5', label:'Chemical potential', dims:1 }, $
  { name:'imu5', label:'Cosmic ray energy density', dims:1 }, $
  { name:'iam', label:'Meanfield dynamo', dims:3 }, $
  { name:'ipsi_real', label:'Wave function real part', dims:1 }, $
  { name:'ipsi_imag', label:'Wave function imaginary part', dims:1 }, $
  { name:'ialpm', label:'alpm', dims:1 }, $
  { name:'ietat', label:'etat', dims:1 }, $
  { name:'ieta', label:'Dust resistivity', dims:1 }, $
  { name:'izeta', label:'Ionization rate', dims:1 }, $
  { name:'ichemspec', label:'Chemical species mass fraction', dims:1 }, $
  { name:'iuud', label:'Dust velocity', dims:3 }, $
  { name:'ind', label:'Dust number density', dims:1 }, $
  { name:'imd', label:'Dust density', dims:1 }, $
  { name:'imi', label:'Dust mi ?something?', dims:1 }, $
  { name:'ilnTT', label:'Log temperature', dims:1 }, $
  { name:'iTT', label:'Temperature', dims:1 }, $
  { name:'ieth', label:'Thermal energy', dims:1 } $
  ; don't forget to add a comma above when extending
]

; Auxiliary variables:
indices_aux = [ $
  { name:'iee', label:'Electric field', dims:3 }, $
  { name:'iQrad', label:'Radiative heating rate', dims:1 }, $
  { name:'ikapparho', label:'Opacity', dims:1 }, $
  { name:'iKR_Frad', label:'Radiative flux scaled with kappa*rho', dims:3 }, $
  { name:'iyH', label:'Hydrogen ionization fraction', dims:1 }, $
  { name:'ishock', label:'Shock profile', dims:1 }, $
  { name:'ishock_perp', label:'B-perpendicular shock profile', dims:1 }, $
  { name:'icooling', label:'Cooling term', dims:1 }, $
  { name:'icooling2', label:'Applied cooling term', dims:1 }, $
  { name:'idetonate', label:'Detonation energy', dims:1 }, $
  { name:'inp', label:'Particle number', dims:1 }, $
  { name:'irhop', label:'Particle mass density', dims:1 }, $
  { name:'iuup', label:'Particle velocity field', dims:3 }, $
  { name:'ifgx', label:'Gas terms for stiff drag forces', dims:3 }, $
  { name:'ipviscx', label:'Particle viscosity field', dims:3 }, $
  { name:'ipotself', label:'Self gravity potential', dims:1 }, $
  { name:'igpotselfx', label:'Gradient of self gravity potential', dims:3 }, $
  { name:'ivisc_heat', label:'Viscous dissipation', dims:1 }, $
  { name:'ihypvis', label:'Hyperviscosity', dims:3 }, $
  { name:'ihypres', label:'Hyperresistivity', dims:3 }, $
  { name:'ippaux', label:'Auxiliary pressure', dims:1 }, $
  { name:'ispecaux', label:'Special auxiliary variable', dims:1 }, $
  { name:'ipsi', label:'Streamfunction', dims:1 }, $
  { name:'isigma', label:'Column density', dims:1 }, $
  { name:'imdot', label:'Mass accretion rate', dims:1 }, $
  { name:'itmid', label:'Midplane temperature', dims:1 }, $
  { name:'ipotturb', label:'Turbulent potential', dims:1 } $
  ; don't forget to add a comma above when extending
]
naux=n_elements(indices_aux)

; Inconsistent names (IDL-name is inconsistent with name in the main code):
inconsistent = [ $
  { name:'ifx', inconsistent_name:'ff' }, $
  { name:'ichemspec', inconsistent_name:'YY' }, $
  { name:'idetonate', inconsistent_name:'det' }, $
  { name:'ifgx', inconsistent_name:'ffg' }, $
  { name:'ipviscx', inconsistent_name:'pvisc' }, $
  { name:'igpotselfx', inconsistent_name:'gpotself' }, $
  { name:'ihypvis', inconsistent_name:'hyv' }, $
  { name:'ihypres', inconsistent_name:'hyr' } $
  ; don't forget to add a comma above when extending
]

; Inconsistent names in special modules (IDL-name is inconsistent with name in the main code):
inconsistent_special = [ $
  { name:'ikappar', inconsistent_name:'kappar' }, $
  { name:'ilambda', inconsistent_name:'lambda' }  $
  ; don't forget to add a comma above when extending
]

; Special variables:
file_special = datadir+'/index_special.pro'
if (file_test (file_special)) then begin
  openr, lun, file_special, /get_lun
  line = ''
  line_pos = 0
  num_inconsistent = n_elements (inconsistent_special)
  while (not eof (lun)) do begin
    readf, lun, line
    line_pos += 1
    ; Backwards-compatibility for old runs with alphadisk, flux_limdiff, streamfunction, or turbpotential.
    for pos = 0, num_inconsistent-1 do begin
      search = inconsistent_special[pos].inconsistent_name
      replace = inconsistent_special[pos].name
      str = stregex (line, '^ *'+search+' *(=.*)$', /extract, /sub)
      line = replace+str[1]
    endfor
    ; Parse line with number of components.
    str = stregex (line, '^ *n[^= ]+[= ]+[0-9]+ *$', /extract)
    if (not execute (str)) then $
        message, 'pc_varcontent: there was a problem with "'+file_special+'" at line '+str (line_pos)+'.', /info
    ; Parse line with "ispecial = ..." or similar.
    str = stregex (line, '^ *(i[^= ]+)[= ]+.*$', /extract, /sub)
    if (str[1] ne '') then begin
      indices = [ indices, { name:str[1], label:'Special', dims:1 } ]
      index_pro = [ index_pro, line ]
    endif
  endwhile
  close, lun
  free_lun, lun
endif

;
;  The number of variables in the snapshot file depends on whether we
;  are writing auxiliary data or not. Auxiliary variables can be turned
;  off by hand by setting noaux=1, e.g. for reading derivative snapshots.
;
if (not keyword_set (noaux)) then begin

  if (keyword_set(param2)) then $
    lpar2aux=keyword_set(param2.lwrite_aux) $
  else $
    lpar2aux=0

  if ( (keyword_set(param.lwrite_aux) and lpar2aux) or $
       (keyword_set(param.lwrite_aux) and ivar eq 0) or $
       (lpar2aux and ivar gt 0) ) then $
    indices = [ indices, indices_aux ]
endif
;
;  Predefine some variable types used regularly.
;
INIT_DATA = [ 'make_array (mx,my,mz,', 'type=type_idl)' ]
INIT_DATA_LOC = [ 'make_array (mxloc,myloc,mzloc,', 'type=type_idl)' ]

;
;  For 2-D runs with lwrite_2d=T. Data has been written by the code without
;  ghost zones in the missing direction. We add ghost zones here anyway so
;  that the array can be treated exactly like 3-D data.
;
if (keyword_set(run2D)) then begin
  if (dim.nx eq 1) then begin
    ; 2-D run in (y,z) plane.
    INIT_DATA_LOC = [ 'make_array (myloc,mzloc,', 'type=type_idl)' ]
  endif else if (dim.ny eq 1) then begin
    ; 2-D run in (x,z) plane.
    INIT_DATA_LOC = [ 'make_array (mxloc,mzloc,', 'type=type_idl)' ]
  endif else begin
    ; 2-D run in (x,y) plane.
    INIT_DATA_LOC = [ 'make_array (mxloc,myloc,', 'type=type_idl)' ]
  endelse
endif

;
;  Parse variables and count total number of variables.
;
totalvars = 0L
num_tags = n_elements (indices)
num_vars = 0
for tag = 1, num_tags do begin
  search = indices[tag-1].name
  dims = indices[tag-1].dims
  add_vars = dims
  ; Backwards-compatibility for old runs with dustdensity or dustvelocity.
  matches = stregex (index_pro, '^ *'+search+' *= *intarr *\( *([0-9]+) *\) *$', /extract, /sub)
  line = max (where (matches[0,*] ne ''))
  if (line ge 0) then begin
    offset_matches = stregex (index_pro, '^ *'+search+' *\[ *0 *\] *= *([0-9]+) *$', /extract, /sub)
    offset_line = max (where (offset_matches[0,*] ne ''))
    if (offset_line ge 0) then begin
      offset = long (offset_matches[1,offset_line])
      index_pro[where (matches[0,*] ne '')] = ''
      index_pro[line] = search+' = indgen ('+str (matches[1,line])+')*'+str (dims)+' + '+str (offset)
      dummy = execute ('n'+strmid (search, 1)+' = '+str (matches[1,line]))
    endif
  endif
  ; Identify f-array variables with multiple components.
  matches = stregex (index_pro, '^ *'+search+' *= *(indgen *\( *[0-9]+ *\).*)$', /extract, /sub)
  line = max (where (matches[0,*] ne ''))
  if (line ge 0) then begin
    if (not execute (index_pro[line])) then $
        message, 'pc_varcontent: there was a problem with "'+indices_file+'" at line '+str (line)+'.', /info
    if (not execute ('num_subtags = n'+strmid (search, 1))) then $
        message, 'pc_varcontent: there was a problem with reading n"'+strmid (search, 1)+'" at line '+str (line)+'.', /info
    if (search eq 'ichemspec') then begin
      matches = [ index_pro[line], '[ '+strjoin (str (ichemspec), ',')+' ]' ]
      add_vars *= num_subtags
    endif else begin
      matches = [ index_pro[line], matches[1,line] ]
      add_vars *= num_subtags
    endelse
  endif else begin
    ; Regular f-array variables.
    matches = stregex (index_pro, '^ *'+search+' *= *([0-9]+|\[[0-9][0-9, ]+\]) *$', /extract, /sub)
  endelse
  line = max (where (matches[0,*] ne ''))
  if (line lt 0) then continue
  exec_str = 'pos = '+matches[1,line]
  if (not execute (exec_str)) then $
      message, 'pc_varcontent: there was a problem with "'+indices_file+'" at line '+str (line)+'.', /info
  if (pos[0] le 0) then continue
  ; Append f-array variable to valid varcontent.
  num_vars += 1
  if (size (selected, /type) eq 0) then begin
    selected = [ tag-1 ]
    executes = [ exec_str ]
    position = [ pos[0] ]
  end else begin
    selected = [ selected, tag-1 ]
    executes = [ executes, exec_str ]
    position = [ position, pos[0] ]
  end
  totalvars += add_vars
endfor

;
;  Make an array of structures in which to store their descriptions.
;
varcontent = replicate ({ varcontent_all, variable:'UNKNOWN', idlvar:'dummy', idlinit:'0', $
    idlvarloc:'dummy_loc', idlinitloc:'0', skip:0 }, totalvars)

;
;  Fill varcontent array.
;
selected = selected[sort (position)]
executes = executes[sort (position)]
for var = 0, num_vars-1 do begin
  tag = selected[var]
  dims = indices[tag].dims
  if (dims eq 1) then joint = '' else joint = str (dims)+','
  replace = where (inconsistent[*].name eq indices[tag].name)
  name = strmid (indices[tag].name, 1)
  dummy = execute (executes[var])
  num_components = n_elements (pos)
  if (strpos (executes[var], 'indgen') ge 0) then begin
    joint += str (num_components)+','
    skip = num_components * dims
    num_components = 1
  endif else begin
    skip = dims
  endelse
  for component = 1, num_components do begin
    if (pos[component-1] gt 0) then begin
      idl_var = name
      if (replace[0] ge 0) then idl_var = inconsistent[replace[0]].inconsistent_name
      if (num_components gt 1) then idl_var += str (component)
      varcontent[pos[component-1]-1].variable = indices[tag].label + ' ('+idl_var+')'
      varcontent[pos[component-1]-1].idlvar = idl_var
      varcontent[pos[component-1]-1].idlinit = strjoin (INIT_DATA, joint)
      varcontent[pos[component-1]-1].idlvarloc = idl_var+'_loc'
      varcontent[pos[component-1]-1].idlinitloc = strjoin (INIT_DATA_LOC, joint)
      varcontent[pos[component-1]-1].skip = skip - 1
    endif
  endfor
endfor

;
;  Turn vector quantities into scalars if requested.
;
if (keyword_set(scalar)) then begin
  for i = 0L, totalvars-1L do begin
    if (varcontent[i].skip eq 2) then begin
      varcontent[i+2].variable  = varcontent[i].variable + ' 3rd component' 
      varcontent[i+1].variable  = varcontent[i].variable + ' 2nd component' 
      varcontent[i  ].variable  = varcontent[i].variable + ' 1st component' 
      varcontent[i+2].idlvar    = varcontent[i].idlvar + '3' 
      varcontent[i+1].idlvar    = varcontent[i].idlvar + '2' 
      varcontent[i  ].idlvar    = varcontent[i].idlvar + '1' 
      varcontent[i+2].idlvarloc = varcontent[i].idlvarloc + '3' 
      varcontent[i+1].idlvarloc = varcontent[i].idlvarloc + '2' 
      varcontent[i  ].idlvarloc = varcontent[i].idlvarloc + '1' 
      varcontent[i:i+2].idlinit    = strjoin (INIT_DATA)
      varcontent[i:i+2].idlinitloc = strjoin (INIT_DATA_LOC)
      varcontent[i:i+2].skip       = 0
      i=i+2
    endif   
  endfor
endif
;
return, varcontent
;
end
