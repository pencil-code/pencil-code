   function identify_varfile, filename=filename, path=file, nohdf5=nohdf5, datadir=datadir
;
;  Identifies which of allprocs/var.dat, allprocs/var.h5, proc0/var.dat is the youngest
;  If filename is provided, correspondingly which of allprocs/<filename>, allprocs/<h5filename>, proc0/<filename> 
;  is youngest where <h5filename> is <filename> with .dat replaced by .h5 or .h5 appended if <filename> is VAR[0-9]+.
;
;  Returns filename.
;
;     path: relative path to file.
;  /nohdf5: ignore .h5 files.
;
     default, filename, 'var.dat'

     if (size (datadir, /type) ne 7) then datadir=pc_get_datadir()

     dotpos=strpos(filename,'.dat')
     issnap=0
     if dotpos eq -1 then begin
       stem=strtrim(filename,2)
       issnap=stregex(stem,'VAR[0-9]+') ne -1
     endif else $
       stem=strmid(filename,0,dotpos)

     file5=stem+'.h5'
     mt=[-1L,-1L,-1L]

     if not keyword_set(nohdf5) then begin
       if file_test (datadir+'/allprocs/'+file5) then $
         mt[0]=(file_info(datadir+'/allprocs/'+file5)).mtime
     endif

     fileb=issnap ? stem : stem+'.dat'
     if file_test (datadir+'/proc0/'+fileb) then $
       mt[1]=(file_info(datadir+'/proc0/'+fileb)).mtime
     if file_test (datadir+'/allprocs/'+fileb) then $
       mt[2]=(file_info(datadir+'/allprocs/'+fileb)).mtime
     
     if max(mt,maxind) eq -1 then begin
       message, 'pc_read: ERROR: No '+file5+' or '+fileb+' found. - Please either give a filename or open an HDF5 file!'
       return, ''
     endif else begin
       file = maxind eq 0 ? datadir+'/allprocs/'+file5 : (maxind eq 1 ? datadir+'/proc0/'+fileb : datadir+'/allprocs/'+fileb)
       return, maxind eq 0 ? file5 : fileb
     endelse

   end
