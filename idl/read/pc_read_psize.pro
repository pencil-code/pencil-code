;+
;  Read particle size distribution from file, as output by the Pencil Code.
;-
;  Author: Anders Johansen (anders@astro.lu.se)
;  Date  : 06-feb-11
;
pro pc_read_psize, object=object, dim=dim, datadir=datadir, filename=filename, $
    quiet=quiet, single=single, help=help
COMPILE_OPT IDL2,HIDDEN

common pc_precision, zero, one, precision, data_type, data_bytes, type_idl

  if (keyword_set(help)) then begin
    doc_library, 'pc_read_psize'
    return      
  endif
;
;  Defaul settings.
;
default, quiet, 0
default, single, 0
;
;  Default datadir, dim, and filename.
;
datadir = pc_get_datadir(datadir)
default, filename, 'particle_size_dist.dat'
if (size (dim, /type) ne 8) then pc_read_dim, obj=dim, datadir=datadir, /quiet
;
;  Set precision.
;
pc_set_precision,  datadir=datadir, dim=dim, /quiet
;
;  Count the number of output entries in the file.
;
nlines = file_lines(datadir+'/'+filename)
;
;  Open the file and read the data.
;
openr, lun, datadir+'/'+filename, /get_lun
;
;  Read information on size of particle bins.
;
log_ap_min_dist=one
log_ap_max_dist=one
nbin_ap_dist=0
;
readf, lun, log_ap_min_dist, log_ap_max_dist, nbin_ap_dist
;
;  Read information on particle bins.
;
log_ap=make_array(nbin_ap_dist, type=single ? 4 : type_idl)
log_ap_low=make_array(nbin_ap_dist, type=single ? 4 : type_idl)
log_ap_high=make_array(nbin_ap_dist, type=single ? 4 : type_idl)
;
readf, lun, log_ap
readf, lun, log_ap_low
readf, lun, log_ap_high
;
;  Read output times and measured particle size distribution.
;
ap_dist_tmp=make_array(nbin_ap_dist, type=type_idl)
t_tmp=one
ap_dist=make_array(nbin_ap_dist,nlines-4, type=single ? 4 : type_idl)
t=make_array(nlines-4, type=single ? 4 : type_idl)
;
for i=0,nlines-5 do begin
  readf, lun, t_tmp, ap_dist_tmp
  t[i]=t_tmp
  ap_dist[*,i]=ap_dist_tmp
endfor
;
close, lun
free_lun, lun
;
;  Inform about the size of the arrays.
;
if (not quiet) then print, 'Read particle size distribution at ' + $
    strtrim(n_elements(t),2)+' times'
;
;  Put data and parameters in object.
;
makeobject="object = create_struct(name=objectname," + $
    "['t','ap_dist','log_ap_min_dist','log_ap_max_dist','nbin_ap_dist', " + $
    " 'log_ap','log_ap_low','log_ap_high'], " + $
    "t, ap_dist, log_ap_min_dist, log_ap_max_dist, nbin_ap_dist, "+ $
    "log_ap, log_ap_low, log_ap_high)"
;
if (execute(makeobject) ne 1) then begin
  message, 'ERROR creating object: ' + makeobject, /info
  undefine, object
endif
;
end
