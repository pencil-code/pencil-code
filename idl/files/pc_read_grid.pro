; $Id: pc_read_grid.pro,v 1.5 2004-03-23 12:45:19 ajohan Exp $
;
;   Read grid.dat
;
;  Author: Tony Mee (A.J.Mee@ncl.ac.uk)
;  $Date: 2004-03-23 12:45:19 $
;  $Revision: 1.5 $
;
;  27-nov-02/tony: coded 
;
;  
pro pc_read_grid,t=t,x=x,y=y,z=z,dx=dx,dy=dy,dz=dz,object=object, $
                 datadir=datadir,proc=proc,PRINT=PRINT,QUIET=QUIET,HELP=HELP
  COMMON pc_precision, zero, one
; If no meaningful parameters are given show some help!
  IF ( keyword_set(HELP) ) THEN BEGIN
    print, "Usage: "
    print, ""
    print, "pc_read_grid, t=t, x=x, y=y, z=z, dx=dx, dy=dy, dz=dz, datadir=datadir, proc=proc,          "
    print, "              /PRINT, /QUIET, /HELP                                                         "
    print, "                                                                                            "
    print, "Returns the grid position arrays and grid deltas of a Pencil-Code run. For a specific       "
    print, "processor. Returns zeros and empty in all variables on failure.                             "
    print, "                                                                                            "
    print, "  datadir: specify the root data directory. Default is './data'                    [string] "
    print, "     proc: specify a processor to get the data from. Default is 0                 [integer] "
    print, ""
    print, "        t: array of x mesh point positions in code length units                [single(mx)]"
    print, "        x: array of x mesh point positions in code length units                [single(mx)]"
    print, "        y: array of y mesh point positions in code length units                [single(my)]"
    print, "        z: array of z mesh point positions in code length units                [single(mz)]"
    print, "       dx: x mesh spacing in code length units                                     [single]"
    print, "       dy: y mesh spacing in code length units                                     [single]"
    print, "       dz: z mesh spacing in code length units                                     [single]"
    print, ""
    print, "   object: optional structure in which to return all the above as tags          [structure] "
    print, ""
    print, "   /PRINT: instruction to print all variables to standard output                            "
    print, "   /QUIET: instruction not to print any 'helpful' information                               "
    print, "    /HELP: display this usage information, and exit                                         "
    return
  ENDIF

; Default data directory

default, datadir, 'data'

; Get necessary dimensions, inheriting QUIET
pc_read_dim,mx=mx,my=my,mz=mz,datadir=datadir,proc=proc,QUIET=QUIET 
; and check pc_precision is set!
pc_set_precision,precision=precision,QUIET=QUIET

;
; Initialize / set default returns for ALL variables
;
t=zero
x=fltarr(mx)*one & y=fltarr(my)*one & z=fltarr(mz)*one
dx=zero &  dy=zero &  dz=zero & dxyz=zero

; Get a unit number
GET_LUN, file

; Build the full path and filename
default,proc,0
filename=datadir+'/proc'+str(proc)+'/grid.dat'   ; Read processor box dimensions

; Check for existance and read the data
dummy=findfile(filename, COUNT=cgrid)
if (cgrid gt 0) then begin
  IF ( not keyword_set(QUIET) ) THEN print, 'Reading ' , filename , '...'

  openr,file,filename,/F77
  readu,file, t,x,y,z
  readu,file, dx,dy,dz
  close,file 
end else begin
  message, 'ERROR: cannot find file ' + filename
end

; Build structure of all the variables
object = CREATE_STRUCT(name=filename,['t','x','y','z','dx','dy','dz'],t,x,y,z,dx,dy,dz)

; If requested print a summary
fmt = '(A,4G15.6)'
if keyword_set(PRINT) then begin
  print, FORMAT='(A,I2,A)', 'For processor ',proc,' calculation domain:'
  print, '             t = ', t
  print, 'min(x), max(x) = ',min(x),', ',max(x)
  print, 'min(y), max(y) = ',min(y),', ',max(y)
  print, 'min(z), max(z) = ',min(z),', ',max(z)
  print, '    dx, dy, dz = ' , dx , ', ' , dy , ', ' , dz
endif

end


