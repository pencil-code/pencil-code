; $Id: pc_read_ts.pro,v 1.1 2003-08-02 15:38:25 mee Exp $
;
;  Read time_series.dat
;
;  Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
;  $Date: 2003-08-02 15:38:25 $
;  $Revision: 1.1 $
;
;  14-nov-02/wolf: coded
;  27-nov-02/tony: ported to routine of standard structure
;
;  REQUIRES: input_table.pro (WD)
;  
;
pro pc_read_ts,n=n,it=it,t=t,dt=dt,dtc=dtc,urms=urms,ekin=ekin,eth=eth,rhom=rhom,ssm=ssm, $
                 object=object, $ 
                 datadir=datadir,PRINT=PRINT,QUIET=QUIET,HELP=HELP

; If no meaningful parameters are given show some help!
  IF ( keyword_set(HELP) ) THEN BEGIN
    print, "Usage: "
    print, ""
    print, "pc_read_ts,  t=t,                                                                           "
    print, "             object=object,                                                                 "
    print, "             /DOUBLE,                                                                       "
    print, "             datadir=datadir, /PRINT, /QUIET, /HELP                                         "
    print, "                                                                                            "
    print, "Read time series data from time_series.dat into separate variables or a structure.          "
    print, "                                                                                            "
    print, "  datadir: specify the root data directory. Default is './data'                     [string]"
    print, "                                                                                            "
    print, "        n: number of entries (valid - not commented out)                           [integer]"
    print, "       it: array of time step numbers                                            [float(it)]"
    print, "        t: array containing time in code units                                   [float(it)]"
    print, "       dt: array of time step sizes                                              [float(it)]"
    print, "      dtc: time step limit by CFL condition                                      [float(it)]"
    print, "     urms: RMS velocity                                                          [float(it)]"
    print, "     ekin: total kinetic energy                                                  [float(it)]"
    print, "      eth: total thermal energy                                                  [float(it)]"
    print, "     rhom:                                                                       [float(it)]"
    print, "      ssm:                                                                       [float(it)]"
    print, ""
    print, "   object: optional structure in which to return all the above as tags           [structure]"
    print, ""
    print, "  /DOUBLE: instruction to read all values as double precision                               "
    print, "   /PRINT: instruction to print all variables to standard output                            "
    print, "   /QUIET: instruction not to print any 'helpful' information                               "
    print, "    /HELP: display this usage information, and exit                                         "
    return
  ENDIF

; Default data directory

default, datadir, 'data'

;
; Initialize / set default returns for ALL variables
;
n=0
it=0.
t=0.
dt=0.
dtc=0.
urms=0.
ekin=0.
eth=0.
rhom=0.
ssm=0.

; Get a unit number
GET_LUN, file

; Build the full path and filename
filename=datadir+'/time_series.dat'         

;
;  read header
;
; Check for existance and read the data
dummy=findfile(filename, COUNT=found)
if (found gt 0) then begin
    if ( not keyword_set(QUIET) ) THEN print, 'Reading ' + filename + '...'
    openr, file, filename
    line = ''
    repeat begin
        readf, file, line
        line = strtrim(line)
        hashpos = strpos(line,'#')
        ;; identify header line as / *#--/
    endrep until ((hashpos ge 0) and (strmid(line,hashpos+1,2) eq '--'))
    close,file
end else begin
    message, 'ERROR: cannot find file ' + filename
end

labels = parse_tsheader(line)
ncols = n_elements(labels)
;
;  read table
;
data = input_table(filename)
if ((size(data))[1] ne ncols) then begin
  message, /INFO, 'Inconsistency: label numer different from column number'
endif
;
;  assemble the data
;
cmd = 'object = {'           ; build a command to execute
for i=0,ncols-1 do begin
  if (i eq 0) then pref=' ' else pref=', '
  cmd = cmd + pref + labels[i] + ': reform(data[' + strtrim(i,2) + ',*])'
endfor
cmd = cmd + ' }'

if (execute(cmd) ne 1) then $
    message, 'There was a problem executing <' + cmd + '>', /INFO



; Unwrap and quantities that may have been separately requested from object
n = (size(data))[1]
if (in_list('t',labels))    then t = object.t
if (in_list('dt',labels))   then dt = object.dt
if (in_list('dtc',labels))  then dtc = object.dtc
if (in_list('urms',labels)) then urms = object.urms
if (in_list('ekin',labels)) then ekin = object.ekin
if (in_list('eth',labels))  then eth = object.eth
if (in_list('rhom',labels)) then rhom = object.rhom
if (in_list('ssm',labels))  then ssm = object.ssm

; If requested print a summary
if keyword_set(PRINT) then begin
    print, 'For GLOBAL calculation domain:'
    print, '    NO SUMMARY INFORMATION CONFIGURED - edit pc_read_params.pro'
endif

end
;
function parse_tsheader, hline
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
;
;  Return index if label is contained in list, else 0
;
  return, (where(list eq label))[0]
end
; ---------------------------------------------------------------------- ;
function in_list, label, list
;
;  Return 1 if label is contained in list, else 0
;
  return, (list_idx(label,list)+1) ne 0
end
; ---------------------------------------------------------------------- ;
