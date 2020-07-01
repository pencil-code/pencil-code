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
function pc_varcontent, datadir=datadir, dim=dim, param=param, par2=run_param, $
                        run2D=run2D, scalar=scalar, noaux=noaux, quiet=quiet, down=down, single=single
;
;    /single: enforces single precision of returned data.
;      /down: data read from downsampled snapshot.
;
COMPILE_OPT IDL2,HIDDEN
;
;  Read grid dimensions, input parameters and location of datadir.
;
datadir = pc_get_datadir(datadir)
default, down, 0
if (n_elements(dim) eq 0) then pc_read_dim, obj=dim, datadir=datadir, quiet=quiet, down=down
if (n_elements(param) eq 0) then pc_read_param, obj=param, datadir=datadir, dim=dim, quiet=quiet
if (n_elements(run_param) eq 0) then pc_read_param, obj=run_param, /param2, datadir=datadir, dim=dim, quiet=quiet
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
default, ntestlnrho, 0
default, n_np_ap, 0

for line = 1, num_lines do begin
  str = stregex (index_pro[line-1], '^ *n[^= ]+[= ]+[0-9]+ *$', /extract)
  if (not execute (str)) then $
      message, 'pc_varcontent: there was a problem with "'+indices_file+'" at line '+str (line)+'.', /info
endfor

mvar=dim.mvar & maux=dim.maux
;
;  For EVERY POSSIBLE variable in a snapshot file, store a
;  description of the variable in an indexed array of structures
;  where the indexes line up with those in the saved f array.
;
;  Note: Integrated variables and variables which can be both integrated and auxiliary *must* be included here.
;        Auxiliary variables should go to the table below the following one, but work also here.
;
indices = [ $
  { name:'iuu', label:'Velocity', dims:3 }, $
  { name:'iadv_der_uu', label:'Advective acceleration as auxiliary variable', dims:3}, $
  { name:'ipp', label:'Pressure', dims:1 }, $
  { name:'ippp', label:'Pressure as auxiliary variable', dims:1 }, $
  { name:'iss', label:'Entropy', dims:1 }, $
  { name:'icp', label:'Specific heat as auxiliary variable', dims:1 }, $
  { name:'icv', label:'Specific heat as auxiliary variable', dims:1 }, $
  { name:'igamma', label:'Ratio of specific heat as auxiliary variable', dims:1 }, $
  { name:'inabad', label:'nabla adiabatic as auxiliary variable', dims:1 }, $
  { name:'idelta', label:'delta as auxiliary variable', dims:1 }, $
  { name:'iviscosity', label:'viscosity as auxiliary variable', dims:1 }, $
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
  { name:'ikappar', label:'kappar', dims:1 }, $
  { name:'itau', label:'tau', dims:1 }, $
  { name:'iggT', label:'ggT', dims:1 }, $
  { name:'iggX', label:'ggX', dims:1 }, $
  { name:'ihhT', label:'hhT', dims:1 }, $
  { name:'ihhX', label:'hhX', dims:1 }, $
  { name:'iggTim', label:'ggTim', dims:1 }, $
  { name:'iggXim', label:'ggXim', dims:1 }, $
  { name:'ihhTim', label:'hhTim', dims:1 }, $
  { name:'ihhXim', label:'hhXim', dims:1 }, $
  { name:'ihij', label:'hij', dims:6 }, $
  { name:'igij', label:'gij', dims:6 }, $
  { name:'ip11', label:'Polymer Tensor 11', dims:1 }, $
  { name:'ip12', label:'Polymer Tensor 12', dims:1 }, $
  { name:'ip13', label:'Polymer Tensor 13', dims:1 }, $
  { name:'ip22', label:'Polymer Tensor 22', dims:1 }, $
  { name:'ip23', label:'Polymer Tensor 23', dims:1 }, $
  { name:'ip33', label:'Polymer Tensor 33', dims:1 }, $
  { name:'iuut', label:'Integrated velocity', dims:3 }, $
  { name:'iaatest', label:'Testmethod vector potential', dims:ntestfield }, $
  { name:'iuutest', label:'Testmethod velocity', dims:ntestflow }, $
  { name:'icctest', label:'Testmethod scalar', dims:ntestscalar }, $
  { name:'ilnrhotest', label:'Testmethod log(rho)', dims:ntestlnrho }, $
  { name:'iuun', label:'Velocity of neutrals', dims:3 }, $
  { name:'ispitzer', label:'Heat flux vector according to Spitzer', dims:3 }, $
  { name:'iqq', label:'heatflux vector', dims:3 }, $
  { name:'ilnrhon', label:'Log density of neutrals', dims:1 }, $
  { name:'ifx', label:'Radiation vector', dims:3 }, $
  { name:'ie', label:'Radiation scalar', dims:1 }, $
  { name:'icc', label:'Passive scalar', dims:1 }, $
  { name:'ilncc', label:'Log passive scalar', dims:1 }, $
  { name:'iacc', label:'Active Scalar', dims:1 }, $
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
  { name:'imuS', label:'Chemical potential', dims:1 }, $
  { name:'imu5', label:'Chiral chemical potential', dims:1 }, $
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
  { name:'ieth', label:'Thermal energy', dims:1 }, $
  { name:'igpx', label:'Pressure gradient x', dims:1 }, $
  { name:'igpy', label:'Pressure gradient y', dims:1 }, $
  { name:'iRR', label:'Specific gas constant', dims:1 }, $
  { name:'iss_run_aver', label:'Running mean of entropy', dims:1 } $
  ; don't forget to add a comma above when extending
]

indices_vector = [ $
  { name:'iuu', components:['iux','iuy','iuz'] }, $
  { name:'iaa', components:['iax','iay','iaz'] }, $
  { name:'ibb', components:['ibx','iby','ibz'] }, $
  { name:'ijj', components:['ijx','ijy','ijz'] } $
  ; don't forget to add a comma above when extending
]

; Auxiliary variables: (see also explanation above)
indices_aux = [ $
  { name:'iee', label:'Electric field', dims:3 }, $
  { name:'iQrad', label:'Radiative heating rate', dims:1 }, $
  { name:'ikapparho', label:'Opacity', dims:1 }, $
  { name:'isss', label:'Entropy as auxiliary variable', dims:1 }, $
  { name:'iKR_Frad', label:'Radiative flux scaled with kappa*rho', dims:3 }, $
  { name:'iyH', label:'Hydrogen ionization fraction', dims:1 }, $
  { name:'ishock', label:'Shock profile', dims:1 }, $
  { name:'ishock_perp', label:'B-perpendicular shock profile', dims:1 }, $
  { name:'icooling', label:'ISM cooling term', dims:1 }, $
  { name:'inetheat', label:'Net applied ISM heating term', dims:1 }, $
  { name:'idetonate', label:'Detonation energy', dims:1 }, $
  { name:'inp', label:'Particle number', dims:1 }, $
  { name:'inp_ap', label:'Particle number, size dependent', dims:n_np_ap }, $
  { name:'iphiuu', label:'Potential of curl-free part of velocity field', dims:1 }, $
  { name:'irhop', label:'Particle mass density', dims:1 }, $
  { name:'iuup', label:'Particle velocity field', dims:3 }, $
  { name:'ifgx', label:'Gas terms for stiff drag forces', dims:3 }, $
  { name:'ipviscx', label:'Particle viscosity field', dims:3 }, $
  { name:'ipotself', label:'Self gravity potential', dims:1 }, $
  { name:'igpotselfx', label:'Gradient of self gravity potential', dims:3 }, $
  { name:'ivisc_heat', label:'Viscous dissipation', dims:1 }, $
  { name:'ihypvis', label:'Hyperviscosity', dims:3 }, $
  { name:'ihypres', label:'Hyperresistivity', dims:3 }, $
  { name:'ihcond', label:'Thermal conductivity', dims:1 }, $
  { name:'iglhc', label:'Gradient of thermal conductivity', dims:3 }, $
  { name:'ippaux', label:'Auxiliary pressure', dims:1 }, $
  { name:'ispecaux', label:'Special auxiliary variable', dims:1 }, $
  { name:'iStr', label:'Str', dims:6 }, $
  { name:'iStT', label:'StT', dims:1 }, $
  { name:'iStX', label:'StX', dims:1 }, $
  { name:'iStTim', label:'StTim', dims:1 }, $
  { name:'iStXim', label:'StXim', dims:1 }, $
  { name:'isld_char', label:'SLD characteristic speed', dims:1 }, $
  { name:'ipsi', label:'Streamfunction', dims:1 }, $
  { name:'isigma', label:'Column density', dims:1 }, $
  { name:'imdot', label:'Mass accretion rate', dims:1 }, $
  { name:'itmid', label:'Midplane temperature', dims:1 }, $
  { name:'ipotturb', label:'Turbulent potential', dims:1 }, $
  { name:'iff', label:'Forcing function', dims:3 }, $
  { name:'itauascalar', label:'Relaxation time', dims:1 }, $
  { name:'issat', label:'Supersaturation', dims:1 }, $
  { name:'icondensationRate', label:'Condensation rate', dims:1 }, $
  { name:'iwaterMixingRatio', label:'Water mixing ratio', dims:1 }, $
  { name:'inusmag', label:'Smagorinsky viscosity', dims:1 }, $
  { name:'ietasmag', label:'Smagorinsky diffusivity', dims:1 }, $
  { name:'iuxbtest', label:'Testfield EMF', dims:ntestfield }, $
  { name:'ijxbtest', label:'Testfield Lorentz force', dims:ntestfield }, $
  { name:'iugutest', label:'Testflow advective acc.', dims:ntestfield }, $
  { name:'iughtest', label:'Testflow enthalpy advection', dims:ntestfield } $
  ; don't forget to add a comma above when extending
]
;
; Inconsistent names (IDL-name is inconsistent with name in the main code):
; E.g., in Fortran we use "ifx", but some IDL scrips expect "ff" in varcontent.
; Note: the initial "i" is automatically removed and hence *not* inconsistent.
;
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

; Inconsistent names in special modules (see also explanation above):
inconsistent_special = [ $
  { name:'ikappar', inconsistent_name:'kappar' }, $   ; seems not inconsistent
  { name:'ilambda', inconsistent_name:'lambda' }  $   ; seems not inconsistent
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
  if maux gt 0 then $
    if keyword_set(param.lwrite_aux) or down then $
      indices = [ indices, indices_aux ] $
    else if is_defined(run_param) then $
      if  keyword_set(run_param.lwrite_aux) then indices = [ indices, indices_aux ] else maux=0
endif
;
;  Predefine some variable types used regularly.
;
if keyword_set(single) then type='4' else type='type_idl'
INIT_DATA = [ 'make_array (mx,my,mz,', 'type='+type+')' ]
;
;  For 2-D runs with lwrite_2d=T. Data has been written by the code without
;  ghost zones in the missing direction. We add ghost zones here anyway so
;  that the array can be treated exactly like 3-D data.
;
if (keyword_set(run2D)) then begin
  INIT_DATA_LOC = [ $
    'reform(make_array (dim.nx eq 1 ? 1 : mxloc,dim.ny eq 1 ? 1 : myloc,dim.nz eq 1 ? 1 : mzloc,', $
                   'type=type_idl))' ]
endif else $
  INIT_DATA_LOC = [ 'make_array (mxloc,myloc,mzloc,', 'type=type_idl)' ]
;
;  Parse variables and count total number of variables.
;
num_tags = n_elements(indices)
num_vars = 0

offsetv = down and (mvar eq 0) ? '-pos[0]+1' : ''    ; corrects index for downsampled varfile if no MVAR variables are contained
						     ; as indices in index.pro refer to the varfile not to the downsampled varfile
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
    found = where(search eq indices_vector[*].name, num_found)
    if (num_found ge 1) then search = indices_vector[found].components[0]
    matches = stregex (index_pro, '^ *'+search+' *= *([0-9]+|\[[0-9][0-9, ]+\]) *$', /extract, /sub)
  endelse
  line = max (where (matches[0,*] ne ''))
  if (line lt 0) then continue

  exec_str = 'pos = ' + matches[1,line]
  if (not execute (exec_str)) then $
      message, 'pc_varcontent: there was a problem with "'+indices_file+'" at line '+str (line)+'.', /info
  if (pos[0] le 0) then continue
  ; Append f-array variable to valid varcontent.
  num_vars += 1
  ncomps = n_elements(pos)

  if (size (selected, /type) eq 0) then begin
    selected = [ tag-1 ]
    executes = [ exec_str + offsetv ]
    position = [ pos[0] ]
    components = [ ncomps ]
  end else begin
    selected = [ selected, tag-1 ]
    executes = [ executes, exec_str + offsetv ]
    position = [ position, pos[0] ]
    components = [ components, ncomps ]
  end
endfor
;
; Reorder to be ascending w.r.t. position.
;
inds=sort(position)
selected = selected[inds]
executes = executes[inds]
components = components[inds]
;
; in the *ordered* list of hits
; only the first mvar+maux entries matter
;
totalvars = 0L
for var=0,num_vars-1 do begin
  tag = selected[var]
  totalvars += indices[tag].dims*components[var]
  if totalvars eq mvar+maux then begin
    selected = selected[0:var]
    executes = executes[0:var]
    num_vars=var+1
    break
  endif
endfor
;
;  Make an array of structures in which to store their descriptions.
;
varcontent = replicate ({ varcontent_all, variable:'UNKNOWN', idlvar:'dummy', idlinit:'0', $
                          idlvarloc:'dummy_loc', idlinitloc:'0', skip:0 }, totalvars)
;
;  Fill varcontent array.
;
vc_pos = 0
for var = 0, num_vars-1 do begin

  tag = selected[var]
  dims = indices[tag].dims
  if (dims eq 1) then joint = '' else joint = str (dims)+','
  replace = where (inconsistent[*].name eq indices[tag].name)
  name = strmid (indices[tag].name, 1)
  dummy = execute (executes[var])
  num_components = components[var]

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
      varcontent[vc_pos+component-1].variable = indices[tag].label + ' ('+idl_var+')'
      varcontent[vc_pos+component-1].idlvar = idl_var
      varcontent[vc_pos+component-1].idlinit = strjoin (INIT_DATA, joint)
      varcontent[vc_pos+component-1].idlvarloc = idl_var+'_loc'
      varcontent[vc_pos+component-1].idlinitloc = strjoin (INIT_DATA_LOC, joint)
      varcontent[vc_pos+component-1].skip = skip - 1
    endif
  endfor

  vc_pos += skip
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
