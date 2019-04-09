; +
; NAME:
;       PC_READ
;
; PURPOSE:
;       Pencil-Code reading routine for old snapshots.
;
; CALLING:
;       pc_read_old, quantity, filename=filename, datadir=datadir, trimall=trim, processor=processor, dim=dim, start=start, count=count
;
; PARAMETERS:
;       quantity [string]: f-array component to read (mandatory).
;       filename [string]: name of the file to read. Default: last opened file
;       datadir [string]: path to the data directory. Default: 'data/'
;       trimall [boolean]: do not read ghost zones. Default: false
;       processor [integer]: number of processor subdomain to read. Default: all
;       dim [structure]: dimension structure. Default: load if needed
;       start [integer]: start reading at this grid position (includes ghost cells)
;       count [integer]: number of grid cells to read from starting position
;
; MODIFICATION HISTORY:
;       $Id$
;       09-Apr-2019/PABourdin: coded
;
function pc_read_old, quantity, filename=filename, datadir=datadir, trimall=trim, processor=processor, dim=dim, start=start, count=count

	COMPILE_OPT IDL2,HIDDEN

	if (strpos (strlowcase (quantity) ,'/') ge 0) then begin
		message, 'pc_read_old: ERROR: "'+quantity+'" can not be read from the old file format!'
	end

	; replace vector components by full vector
	load = quantity
	if (strlowcase (load) eq 'ux') then begin & load = 'uu' & component = 0 & end
	if (strlowcase (load) eq 'uy') then begin & load = 'uu' & component = 1 & end
	if (strlowcase (load) eq 'uz') then begin & load = 'uu' & component = 2 & end
	if (strlowcase (load) eq 'ax') then begin & load = 'aa' & component = 0 & end
	if (strlowcase (load) eq 'ay') then begin & load = 'aa' & component = 1 & end
	if (strlowcase (load) eq 'az') then begin & load = 'aa' & component = 2 & end
	if (strlowcase (load) eq 'bx') then begin & load = 'bb' & component = 0 & end
	if (strlowcase (load) eq 'by') then begin & load = 'bb' & component = 1 & end
	if (strlowcase (load) eq 'bz') then begin & load = 'bb' & component = 2 & end
	if (strlowcase (load) eq 'jx') then begin & load = 'jj' & component = 0 & end
	if (strlowcase (load) eq 'jy') then begin & load = 'jj' & component = 1 & end
	if (strlowcase (load) eq 'jz') then begin & load = 'jj' & component = 2 & end

	if (keyword_set (processor) or (not keyword_set (start) and not keyword_set (count))) then begin
		; read from old format
		pc_read_var_raw, obj=data, var_list=load, proc=processor, dim=dim, trimall=trim
	end else begin
		; read subvolume from old format
		if (keyword_set (start)) then begin
			xs = start[0]
			ys = start[1]
			zs = start[2]
			if (keyword_set (trim)) then begin
				if (size (dim, /type) eq 0) then pc_read_dim, obj=dim, datadir=datadir, /quiet
				if (dim.nxgrid gt 1) then xs += dim.nghostx
				if (dim.nygrid gt 1) then ys += dim.nghosty
				if (dim.nzgrid gt 1) then zs += dim.nghostz
			end
		end
		if (keyword_set (count)) then begin
			default, start, [ 0, 0, 0 ]
			xe = start[0] + count[0] - 1
			ye = start[1] + count[1] - 1
			ze = start[2] + count[2] - 1
			if (keyword_set (trim)) then begin
				if (size (dim, /type) eq 0) then pc_read_dim, obj=dim, datadir=datadir, /quiet
				if (dim.nxgrid gt 1) then xe -= dim.nghostx
				if (dim.nygrid gt 1) then ye -= dim.nghosty
				if (dim.nzgrid gt 1) then ze -= dim.nghostz
			end
		end
		pc_read_subvol_raw, obj=data, var_list=load, xs=xs, xe=xe, ys=ys, ye=ye, zs=zs, ze=ze, dim=dim
	end

	if (size (component, /type) ne 0) then begin
		; get component
		sizes = size (data, /dimensions)
		data = reform (data[*,*,*,component], sizes[0:2])
	end

	return, data
end

