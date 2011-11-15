pro pc_write_vapor,vdffile=vdffile,varfile=varfile,datadir=datadir, $
                   ivar=ivar_, $
                   variables=variables, $
                   proc=proc,varcontent=varcontent, $
                   scalar=scalar, _extra=_extra, additional=additional
;
; Program to convert a directory of pencil files into a vapor dataset.
; Arguments are:
;	topdir = absolute path to directory where pencil data is stored
;	nprocs = number of procN subdirectories (N goes from 0 to nprocs -1)
;	vdffile = absolute path, including filename, of vdf file for the result.
;		The vdf file should be created by running vdffile before this.
;		The vdf file should specify names for the variables, in the same order
;		as they appear in the pencil var.dat files
;	varfilename = a string used to identify the name of the var.dat files.  If it
;		is not specified, then the default is "var"
;	 
; 	The chunks associated with the different processors are not assumed to be
;	the same size, however, all chunks at a given x-coordinate have the same x-thickness,
;	and likewize for the chunks at a given y or z coordinate
;
    pc_read_dim,obj=dim,proc=proc
    if (n_elements(proc) eq 1L) then nprocs=1 else nprocs = dim.nprocx*dim.nprocy*dim.nprocz
    numslabs=dim.nprocz

    if (n_elements(ivar_) eq 1L) then begin
      numts=1L
      ivar=ivar_
    endif else if (n_elements(ivar_) eq 2L) then begin
   
      default,ivarmin,min(ivar_)
      numts=max(ivar_)-min(ivar_)+1
      ivar=ivarmin
    endif else begin
      print,'ivar must be either a scalar value or range [ivarmin,ivarmax]'
      return
    endelse

    default,datadir,'data'
    default,vdffile,datadir+'/var.vdf'
    ;default,vdffile,'var.vdf'

    default,varcontent,pc_varcontent(datadir=datadir,dim=dim,param=param,quiet=quiet,scalar=scalar)
 
    filevars=varcontent[where((varcontent[*].idlvar ne 'dummy'))].idlvar
    if n_elements(variables) ne 0 then begin
      if keyword_set(ADDITIONAL) then begin
        variables=[filevars,variables]
      endif
    endif else begin
      variables=filevars
    endelse

;    pc_read_var,proc=0,varfile=varfile,ivar=ivar,variables=variables,tags=tags, $
;                /quiet,varcontent=varcontent,_extra=_extra
    
 
    dimstr=strcompress(string(dim.nx),/remove_all)+'x' $
          +strcompress(string(dim.ny),/remove_all)+'x' $
          +strcompress(string(dim.nz),/remove_all)

    vdfcreate_command='vdfcreate -dimension '+dimstr $
;                              +' -bs '+dimstr $
                              +' -numts '+strcompress(string(numts)) $
                              +' -comment "Created by pc_write_vapor"' $
                              +' -gridtype regular' $
                              +' -coordsys cartesian' $
                              +' -varnames '+arraytostring(variables,list=':',/noleader) $
                              +' ' + vdffile 

    print,vdfcreate_command
    spawn,vdfcreate_command,exit_status=exit_status
    if (exit_status ne 0) then begin
      print,"Failed to create vdf file"
    endif
;Usage: vdfcreate [options] filename
;    -dimensio arg0            Volume dimensions (NXxNYxNZ)
;    -numts    arg0            Number of timesteps
;    -bs       arg0            Internal storage blocking factor (NXxNYxNZ)
;    -level    arg0            Maximum refinement level. 0 => no refinement
;    -nfilter  arg0            Number of wavelet filter coefficients
;    -nlifting arg0            Number of wavelet lifting coefficients
;    -comment  arg0            Top-level comment
;    -gridtype arg0            Data grid type (regular|streched|block_amr)
;    -coordsys arg0            Top-level comment (cartesian|spherical)
;    -extents  arg0            Domain extents in user coordinates
;    -varnames arg0            Colon delimited list of variable names
;    -mtkcompa                 Force compatibility with older mtk files
;    -help                     Print this message and exit


;
;	Read dim.dat files:
;
;	(while finding the thickest z-size)
	maxzthick = 0
	for I = 0, nprocs-1 do begin
        pc_read_dim,obj=procdim,proc=i
		maxzthick = max([maxzthick , procdim.nz])

        if (i eq 0) then begin
          procdims=[procdim]
        endif else begin
          procdims=[procdims,procdim]
        endelse
	endfor
;
;	Sort the tables on z-coordinate:
;
	fileorder = SORT(procdims[*].ipz)

;	Create the metadata from the vdffile.
;	Find the dimensions of the data.
;	These had better agree with the pencil data!
; 
print,"Call vdf_create: ", vdffile
	mfd = vdf_create(vdffile)
print,"Call vdf_getdimension"
	dim = vdf_getdimension(mfd)
print,"Call vdf_getvarnames"
	varnames = vdf_getvarnames(mfd)
	sz = size(varnames)
	numvariables = sz[1]

;
;	Determine how many chunks of pencil data are associated with a slab:
;
	chunksperslab = nprocs/numslabs
;	Allocate enough memory to hold the largest slab:

	slabdata = FLTARR(dim[0],dim[1],maxzthick)

;	Loop over each variable
	for timestep = 0L, numts -1 DO BEGIN
	for varnum = 0L, numvariables -1 DO BEGIN
        if (ivarmin ge 0) then ivar=ivarmin+timestep
		print, 'ivar: ',strcompress(string(ivar),/remove_all),' variable: ',varnames[varnum]
		dfd = vdc_bufwritecreate(mfd)
		vdc_openvarwrite, dfd, timestep, varnames[varnum], -1
;
;  	Loop over the slabs:
	
        chunknum=0
		for slab = 0L, numslabs - 1 do begin
;	For each slab, loop over the proc directories (chunks) associated with it:
			for chunk = 0L, chunksperslab - 1 do begin 		
				proc = fileorder[chunknum]
;				read the chunk into the dataarray
;				Create an array to hold a variable chunk as it is read from the proc directory 
                pc_read_var,obj=data,proc=proc,varfile=varfile,ivar=ivar, $
                            variables=[ varnames[varnum] ], $
                            /trimall,/quiet,varcontent=varcontent, $
                            additional=additional, $
                            scalar=scalar,_extra=_extra
;
;				Then copy the data chunks into the data slab:
                minx = procdims[proc].ipx * procdims[proc].nx
				maxx = FIX(minx + procdims[proc].nx - 1L)
                miny = procdims[proc].ipy * procdims[proc].ny
				maxy = FIX(miny + procdims[proc].ny - 1L)
                minz = procdims[proc].ipz * procdims[proc].nz
				maxz = FIX(minz + procdims[proc].nz - 1L)
				res=execute('slabdata[minx:maxx,miny:maxy,0:(maxz-minz)] = float(data.'+varnames[varnum]+')')
				;print,'slabdata[minx:maxx,miny:maxy,0:(maxz-minz)] = float(data.'+varnames[varnum]+')'
				chunknum = chunknum+1
			endfor
; 			now the full slab is populated, write it to the vdf, one
; 			z-slice at a time
;
			for z = 0, (maxz-minz) do begin
				vdc_bufwriteslice, dfd, slabdata[*,*,z]
			endfor
		endfor
;  	Now close the dfd
		vdc_closevar, dfd
		vdc_bufwritedestroy,dfd
	endfor
	endfor
end
