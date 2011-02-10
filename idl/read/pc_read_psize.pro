;
;  $Id$
;
;  Read particle size distribution from file, as output by the Pencil Code.
;
;  Author: Anders Johansen (anders@astro.lu.se)
;  Date  : 06-feb-11
;
pro pc_read_psize, object=object, datadir=datadir, filename=filename, $
    quiet=quiet
COMPILE_OPT IDL2,HIDDEN
COMMON pc_precision, zero, one
;
;  Defaul settings.
;
default, quiet, 0
;
;  Default datadir and filename.
;
if (not keyword_set(datadir)) then datadir=pc_get_datadir()
default, filename, 'particle_size_dist.dat'
;
;  Set the precision.
;
pc_set_precision, dim=dim, /quiet
;
;  Count the number of output entries in the file.
;
nlines = file_lines(datadir+'/'+filename)
;
;  Open the file and read the data.
;
openr, 1, datadir+'/'+filename
;
;  Read information on size of particle bins.
;
log_ap_min_dist=one
log_ap_max_dist=one
nbin_ap_dist=0
;
readf, 1, log_ap_min_dist, log_ap_max_dist, nbin_ap_dist
;
;  Read information on particle bins.
;
log_ap=fltarr(nbin_ap_dist)*one
log_ap_low=fltarr(nbin_ap_dist)*one
log_ap_high=fltarr(nbin_ap_dist)*one
;
readf, 1, log_ap
readf, 1, log_ap_low
readf, 1, log_ap_high
;
;  Read output times and measured particle size distribution.
;
ap_dist_tmp=fltarr(nbin_ap_dist)*one
t_tmp=one
ap_dist=fltarr(nbin_ap_dist,nlines-4)*one
t=fltarr(nlines-4)*one
;
for i=0,nlines-5 do begin
  readf, 1, t_tmp, ap_dist_tmp
  t[i]=t_tmp
  ap_dist[*,i]=ap_dist_tmp
endfor
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
