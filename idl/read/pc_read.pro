;
; $Id$
;
;  Pencil-Code unified reading routine
;
;  Author: Philippe Bourdin
;  $Date: 2019-04-07 11:56:32 $
;  $Revision: 1.0 $
;
;  07-Apr-2019/PABourdin: coded
;
function pc_read, quantity, filename=filename, datadir=datadir, trim=trim, proc=proc, dim=dim, grid=grid, start_param=start_param, run_param=run_param

	COMPILE_OPT IDL2,HIDDEN

	if (keyword_set (filename)) then file = pc_get_datadir (datadir)+'/'+filename

	if (keyword_set (proc)) then begin
		if (size (dim, /type) eq 0) then pc_read_dim, obj=dim
		num_ghosts = dim.nghostx
		ipx = proc mod dim.nprocx
		ipy = proc / dim.nprocx
		ipz = proc / (dim.nprocx*dim.nprocy)
		start = [ ipx*dim.nx, ipy*dim.ny, ipz*dim.nz ]
		count = [ dim.mx, dim.my, dim.mz ]
		if (keyword_set (trim)) then begin
			first[remove] += num_ghosts
			sizes[remove] -= num_ghosts*2
		end
	end

	if (keyword_set (trim)) then begin
		if (size (dim, /type) eq 0) then pc_read_dim, obj=dim
		num_ghosts = dim.nghostx
		default, start, [ 0, 0, 0 ]
		default, count, [ dim.mxgrid, dim.mygrid, dim.mzgrid ]
		remove = where (count ge num_ghosts*2+1, num_remove)
		if (num_remove gt 0) then begin
			start[remove] += num_ghosts
			count[remove] -= num_ghosts*2
		end
	end

	return, hdf5_read (label, filename=file, start=start, count=count)
end

