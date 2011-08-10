;
; $Id$
;
;   Read slice files
;
;  Author: Tony Mee (A.J.Mee@ncl.ac.uk)
;  $Date: 2007-08-03 09:53:26 $
;  $Revision: 1.7 $
;
;  28-nov-02/tony: coded 
;
;  
pro pc_read_slice,field=field,plane=plane,slice=slice,t=t,x=x,y=y,z=z,dx=dx,dy=dy,dz=dz, $
                 object=object, dim1=dim1, dim2=dim2, $
                 datadir=datadir,proc=proc,PRINT=PRINT,QUIET=QUIET,HELP=HELP
COMPILE_OPT IDL2,HIDDEN
  COMMON pc_precision, zero, one
; If no meaningful parameters are given show some help!
  IF ( keyword_set(HELP) ) THEN BEGIN
    print, "Usage: "
    print, ""
    print, "pc_read_slice, t=t, x=x, y=y, z=z, dx=dx, dy=dy, dz=dz, datadir=datadir, proc=proc,          "
    print, "              /PRINT, /QUIET, /HELP                                                         "
    print, "                                                                                            "
    print, "Returns slice data with the grid position arrays and grid deltas of a Pencil-Code run. By   "
    print, "default the first slice is loaded.  If an already open file is provided the NEXT slice is   "
    print, "read."
    print, "Returns zeros and empty in all variables on failure.                             "
    print, "                                                                                            "
    print, "  datadir: specify the root data directory. Default is './data'                    [string] "
    print, "     proc: specify a processor to get the data from. Default is 0                 [integer] "
    print, ""
    print, "    field: data field to read (lnrho,ux,uy,uz,ss,aa etc.) Default is lnrho         [string] "
    print, "    plane: selected plane (xy,xz,yz,Xy) Default is . Default is 0                  [string] "
    print, ""
    print, "    slice: returned 2-D array of slice data. Size depends on plane.        [single(m?)(m?)]"
    print, "     dim1: appropriate x,y or z array for first dimension of the slice         [single(m?)]"
    print, "     dim2: appropriate x,y or z array for second dimension of the slice        [single(m?)]"
    print, "        t: array of x mesh point positions in code length units                [single(mx)]"
    print, "        x: array of x mesh point positions in code length units                [single(mx)]"
    print, "        y: array of y mesh point positions in code length units                [single(my)]"
    print, "        z: array of z mesh point positions in code length units                [single(mz)]"
    print, "       dx: x mesh spacing in code length units                                     [single]"
    print, "       dy: y mesh spacing in code length units                                     [single]"
    print, "       dz: z mesh spacing in code length units                                     [single]"
    print, ""
    print, "     file: file logical unit number returned of of previously opened slice file   [integer]"
    print, ""
    print, "   object: optional structure in which to return all the above as tags          [structure] "
    print, ""
    print, "   /PRINT: instruction to print all variables to standard output                            "
    print, "   /QUIET: instruction not to print any 'helpful' information                               "
    print, "    /HELP: display this usage information, and exit                                         "
    return
  ENDIF

; Default data directory

if (not keyword_set(datadir)) then datadir=pc_get_datadir()

; Get necessary dimensions, inheriting QUIET
pc_read_dim,mx=mx,my=my,mz=mz,datadir=datadir,proc=proc,QUIET=QUIET 
; and check pc_precision is set!
pc_set_precision,dim=dim,QUIET=QUIET

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
if (file_test(filename)) then begin
  IF ( not keyword_set(QUIET) ) THEN print, 'Reading ' , filename , '...'

  openr,file,filename
  readu,file, t,x,y,z
  readu,file, dx,dy,dz
  close,file 
  FREE_LUN, file
end else begin
  FREE_LUN, file
  message, 'ERROR: cannot find file ' , filename
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


