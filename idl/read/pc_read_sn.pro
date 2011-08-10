;
; $Id$
;
;  Read time_series.dat and sort data into structure or variables
;
;  Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
;  $Date: 2007-08-03 09:53:26 $
;  $Revision: 1.4 $
;
;  14-nov-02/wolf: coded
;  27-nov-02/tony: ported to routine of standard structure
;
;  Requires: input_table.pro
;
;  Usage:
;    pc_read_ts, OBJECT=ts1           ; Read all into a structure (recommended)
;    plot, ts1.t, ts1.urms
;    help, ts1, /STRUCTURE
;
;    pc_read_ts, T=time, DT=timestep  ; Name variables (limited to those coded)
;
;    pc_read_ts, /HELP
;
;
function parse_tsheader, hline
COMPILE_OPT IDL2,HIDDEN
;
;  Split header line into the individual variable names.
;
  line = strmid(hline,strpos(hline,'#')+1)
  labels = ['']
  ;
  ; strsplit() is not available in IDL prior to 5.3 (and str_sep() was
  ; obsoleted after 5.2..), so we do this manually:
  ;
  while (line ne '') do begin
    repeat begin
      line = strmid(line,1)     ; remove first character
    endrep until (strmid(line,0,1) ne '-')
    endlb = strpos(line,'-')
    labels = [labels, strmid(line,0,endlb)]
    line = strmid(line,endlb)
  endwhile
  ;
  ;  eliminate empty labels
  ;
  good = where(labels ne '')
  return, labels[good]
end
; ---------------------------------------------------------------------- ;
function list_idx, label, list
COMPILE_OPT IDL2,HIDDEN
;
;  Return index if label is contained in list, else 0
;
  return, (where(list eq label))[0]
end
; ---------------------------------------------------------------------- ;
function in_list, label, list
COMPILE_OPT IDL2,HIDDEN
;
;  Return 1 if label is contained in list, else 0
;
  return, (list_idx(label,list)+1) ne 0
end
; ---------------------------------------------------------------------- ;
pro pc_read_sn,n=n, object=object, filename=filename, $ 
                 datadir=datadir,PRINT=PRINT,QUIET=QUIET,HELP=HELP
COMPILE_OPT IDL2,HIDDEN

; If no meaningful parameters are given show some help!
  IF ( keyword_set(HELP) ) THEN BEGIN
    print, "Usage: "
    print, ""
    print, "pc_read_sn,  t=t,                                                                           "
    print, "             object=object,                                       " 
    print, "             datadir=datadir, /PRINT, /DOUBLE, /QUIET, /HELP      "
    print, "                                                                  "
    print, "Read time series data from time_series.dat into separate "
    print, "variables or a structure.          "
    print, "                                              "
    print, "  datadir: specify the root data directory.               [string]"
    print, "           Default is './data'           "
    print, " filename: specify the filename of the time series data   [string]"
    print, "           Default is 'time_series.dat'                           "
    print, "                                                                  "
    print, "        n: number of entries (valid - not commented out) [integer]"
    print, ""
    print, "   object: optional structure in which to return all   [structure]"
    print, "           the above as tags   "
    print, ""
    print, "  /DOUBLE: instruction to read all values as double precision    "
    print, "   /PRINT: instruction to print all variables to standard output "
    print, "   /QUIET: instruction not to print any 'helpful' information    "
    print, "    /HELP: display this usage information, and exit              "
    return
  ENDIF

; Default data directory

if (not keyword_set(datadir)) then datadir=pc_get_datadir()
default, filename, 'sn_series.dat'
;
; Initialize / set default returns for ALL variables
;
n=0

; Get a unit number
GET_LUN, file

; Build the full path and filename
fullfilename=datadir+'/'+filename         

;
;  read header
;
; Check for existance and read the data
if (file_test(fullfilename)) then begin
    if ( not keyword_set(QUIET) ) THEN print, 'Reading ' + fullfilename + '...'
    openr, file, fullfilename
    line = ''
    repeat begin
        readf, file, line
        line = strtrim(line)
        hashpos = strpos(line,'#')
        ;; identify header line as / *#--/
    endrep until ((hashpos ge 0) and (strmid(line,hashpos+1,2) eq '--'))
    point_lun,-file,fileposition
    close,file
    FREE_LUN,file
end else begin
    FREE_LUN,file
    message, 'ERROR: cannot find file ' + fullfilename
end

;
;  read table
;
newheader=line
full_data=0.
full_labels=''

while (fileposition ne -1) do begin
  labels = parse_tsheader(newheader)
  ncols = n_elements(labels)
  newheader='^#--'
 
  data = input_table(fullfilename,DOUBLE=DOUBLE,  $
                     STOP_AT=newheader,FILEPOSITION=fileposition)
  if ((size(data))[1] ne ncols) then begin
    message, /INFO, 'Inconsistency: label number different from column number'
  endif

  ; Merge read data into full data set
  if ((size(full_labels))[0] eq 0) then begin
    ;If it's the first time just chunk the data into place
    full_labels=labels
    full_data=data
    if ((size(full_data))[0] eq 1) then $
        full_data=reform(full_data,n_elements(full_data),1)
    if ((size(full_data))[0] eq 0) then $
        full_data=reform([full_data],1,1)
    end else begin
    col_index=intarr(ncols)
    for i=0,ncols-1 do begin
      if (not in_list(labels[i],full_labels)) then begin
         old_full_labels=full_labels
         old_ncols=n_elements(old_full_labels)
         full_labels=strarr(old_ncols+1)
         full_labels[0:old_ncols-1]=old_full_labels[*]
         old_full_labels=0. 

         full_labels[old_ncols]=labels[i]
	 col_index[i]=old_ncols
      endif else begin
        col_index[i]=list_idx(labels[i],full_labels)
      endelse
    endfor
         
    old_ncols=(size(full_data))[1]
    
    old_nrows=(size(full_data))[2]

    new_ncols=(size(full_labels))[1]
    new_nrows=old_nrows + (size(data))[2]
    
    old_full_data=full_data
    if (keyword_set(double)) then begin
      full_data = dblarr(new_ncols, new_nrows)
      full_data[*,*] = !VALUES.D_NAN
    endif else begin
      full_data = fltarr(new_ncols, new_nrows)
      full_data[*,*] = !VALUES.F_NAN
    endelse
    
    full_data[0:old_ncols-1,0:old_nrows-1] = old_full_data[0:old_ncols-1,0:old_nrows-1]
    old_full_data=0.  ; Clear the allocated memory 
    
    for i=0,ncols-1 do begin
      full_data[col_index[i],old_nrows:new_nrows-1]=data[i,*]
    endfor   
  endelse
endwhile


;
;  assemble the data
;
ncols = n_elements(full_labels)
cmd = 'object = {'           ; build a command to execute
for i=0,ncols-1 do begin
  if (i eq 0) then pref=' ' else pref=', '
  cmd = cmd + pref + full_labels[i] $
                   + ': reform(full_data[' $
                   + strtrim(i,2) + ',*])'
endfor
cmd = cmd + ' }'

if (execute(cmd) ne 1) then $
  message, 'There was a problem executing <' + cmd + '>', /INFO

; Unwrap and quantities that may have been separately requested from object
n = (size(data))[1]
;if (in_list('t',full_labels))    then t = object.t

; If requested print a summary
if keyword_set(PRINT) then begin
    print, 'For GLOBAL calculation domain:'
    print, '    NO SUMMARY INFORMATION CONFIGURED - edit pc_read_ts.pro'
endif

end
; ---------------------------------------------------------------------- ;
