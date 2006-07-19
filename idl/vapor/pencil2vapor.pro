PRO pencil2vapor,topdir,numprocs, vdffile,timestep,varfilename
;
; Program to convert a directory of pencil files into a vapor dataset.
; Arguments are:
;	topdir = absolute path to directory where pencil data is stored
;	numprocs = number of procN subdirectories (N goes from 0 to numprocs -1)
;	vdffile = absolute path, including filename, of vdf file for the result.
;		The vdf file should be created by running vdffile before this.
;		The vdf file should specify names for the variables, in the same order
;		as they appear in the pencil var.dat files
;	timestep = an integer timestep for the data to be converted.  Must be within
;		the number of time steps specified in the vdf file 
;	varfilename = a string used to identify the name of the var.dat files.  If it
;		is not specified, then the default is "var"
;	 
; 	The chunks associated with the different processors are not assumed to be
;	the same size, however, all chunks at a given x-coordinate have the same x-thickness,
;	and likewize for the chunks at a given y or z coordinate
;
	IF (N_PARAMS() LT 5 ) THEN varfilename = 'var'
	varfilename = varfilename + '.dat' 
;
;    Create tables to hold information about data:
;
	xsize = intarr(numprocs)
	ysize = intarr(numprocs)
	zsize = intarr(numprocs)
	xposition = intarr(numprocs)
	yposition = intarr(numprocs)
	zposition = intarr(numprocs)
;
;	Read dim.dat files:
;
;	(while finding the thickest z-size)
	maxzthick = 0
	FOR I = 0, numprocs-1 DO BEGIN
		dimfile = topdir+'/proc'+STRTRIM(STRING(I),1)+'/dim.dat'
		OPENR, 1, dimfile 
;	Read the first line (contains sizes)
		READF,1,sizex,sizey,sizez,foo,bar
		xsize[I] = sizex
		ysize[I] = sizey
		zsize[I] = sizez
		IF (sizez GT maxzthick) THEN maxzthick = sizez
;	Read second line (precision, a string) 
		prec = ' '
		READF,1,prec
;	Read third line (ghost widths)
		READF,1,ghostx,ghosty,ghostz
;	Read the last line (contains positions)
		READF,1,xpos,ypos,zpos
		xposition[I] = xpos
		yposition[I] = ypos
		zposition[I] = zpos
		CLOSE, 1
	ENDFOR
	maxzthick = maxzthick - 2*ghostz
;  	Determine the x,y,and z-grid spacing

	xspacing = intarr(numprocs)
	yspacing = intarr(numprocs)
	zspacing = intarr(numprocs)
	xspacing[*] = 0
	yspacing[*] = 0
	yspacing[*] = 0

	FOR I = 0, numprocs-1 DO BEGIN
		xspacing[xposition[I]] = xsize[I] - 2*ghostx
		yspacing[yposition[I]] = ysize[I] - 2*ghosty
		zspacing[zposition[I]] = zsize[I] - 2*ghostz
	ENDFOR

;	Accumulate the spacings
	xtot = 0
	ytot = 0
	ztot = 0
	FOR I = 0, numprocs -1 DO BEGIN
		IF (xspacing[I] NE 0) THEN BEGIN
			xtot = xtot + xspacing[I]
			xspacing[I] = xtot - xspacing[I]		
		ENDIF
		IF (yspacing[I] NE 0) THEN BEGIN
			ytot = ytot + yspacing[I]
			yspacing[I] = ytot - yspacing[I]		
		ENDIF
		IF (zspacing[I] NE 0) THEN BEGIN
			numslabs = I+1
			ztot = ztot + zspacing[I]
			zspacing[I] = ztot - zspacing[I]		
		ENDIF
	ENDFOR

;	Sort the tables on z-coordinate:
;
	fileorder = SORT(zposition)

;	Create the metadata from the vdffile.
;	Find the dimensions of the data.
;	These had better agree with the pencil data!
; 
	mfd = vdf_create(vdffile)
	dim = vdf_getdimension(mfd)
	varnames = vdf_getvarnames(mfd)
	sz = size(varnames)
	numvariables = sz[1]

;
;	Determine how many chunks of pencil data are associated with a slab:
;
	chunksperslab = numprocs/numslabs
;	Allocate enough memory to hold the largest slab:

	slabdata = FLTARR(dim[0],dim[1],maxzthick)

;	Loop over each variable
	FOR varnum = 0, numvariables -1 DO BEGIN
		print, 'assembling variable ',varnames[varnum]
		dfd = vdc_bufwritecreate(mfd)
		vdc_openvarwrite, dfd, timestep, varnames[varnum], -1
;
;  	Loop over the slabs:
	
		FOR slab = 0, numslabs -1 DO BEGIN
;	For each slab, loop over the proc directories (chunks) associated with it:
			FOR chunk = 0, chunksperslab-1 DO BEGIN 		
				chunknum = slab*chunksperslab + chunk 
				dirnum = fileorder[chunknum]
;				read the chunk into the dataarray
;				Create an array to hold a variable chunk as it is read from the proc directory 
				IF (prec EQ 'D') THEN BEGIN
					dataarray = DBLARR(xsize[dirnum],ysize[dirnum],zsize[dirnum],numvariables)
				ENDIF ELSE dataarray = FLTARR(xsize[dirnum],ysize[dirnum],zsize[dirnum],numvariables)
				varfile = topdir+'/proc'+STRTRIM(STRING(dirnum),1)+'/'+varfilename
				openr,1,varfile,/f77
				readu,1,dataarray
				close,1
;
;				Then copy the data chunks into the data slab:
				minx = xspacing[xposition[dirnum]]
				maxx = FIX(minx + xsize[dirnum] - 2*ghostx-1)
				miny = yspacing[yposition[dirnum]]
				maxy = FIX(miny + ysize[dirnum] - 2*ghosty-1)
				minz = zspacing[zposition[dirnum]]
				maxz = FIX(minz + zsize[dirnum] - 2*ghostz-1)
				slabdata[minx:maxx,miny:maxy,0:(maxz-minz)] = FLOAT(dataarray[ghostx:xsize[dirnum]-ghostx-1,ghosty:ysize[dirnum]-ghosty-1,ghostz:zsize[dirnum]-ghostz-1,varnum])
			ENDFOR
; 			now the full slab is populated, write it to the vdf, one
; 			z-slice at a time
;
			FOR z = 0, (maxz-minz) DO BEGIN
				vdc_bufwriteslice, dfd, slabdata[*,*,z]
			ENDFOR
		ENDFOR
;  	Now close the dfd
		vdc_closevar, dfd
		vdc_bufwritedestroy,dfd
	ENDFOR
END
