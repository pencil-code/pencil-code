; Writes a data snapshot in the HDF5 format

pro pc_write_qvar, varfile, obj, tags=tags, group=group, time=time, datadir=datadir, quiet=quiet

	datadir = pc_get_datadir (datadir)
	default, group, 'points'

	; case-insensitve replacements for dataset names
	replace = { t:'', xx:'_', vv:'v_' }
	search = strlowcase (tag_names (replace))
	num_replace = n_elements (search)

	if (strmid (group, strlen (group)-1) eq '/') then group = strmid (group, 0, strlen (group)-1)
	if (strmid (varfile, strlen (varfile)-4) eq '.dat') then varfile = strmid (varfile, 0, strlen (varfile)-4)
	if (strmid (varfile, strlen (varfile)-3) ne '.h5') then varfile += '.h5'
	h5_open_file, datadir+'/allprocs/'+varfile, /write, /truncate
	h5_write, 'number', n_elements (obj.mass)
	if (keyword_set (group)) then h5_create_group, group

	; write from qvar structure (pc_read_qvar)
	labels = strlowcase (tag_names (obj))
	tags = labels

	num_content = n_elements (labels)
	for i = 0, num_content-1 do begin
		label = labels[i]
		tag_name = label
		for j = 0, num_replace-1 do begin
			if (label eq search[j]) then label = replace.(j)
		end
		if (label eq '') then continue

		pos = (where (tag_name eq strlowcase (tag_names (obj))))[0]
		if (size (obj.(pos), /n_dimensions) eq 2) then begin
			num_dims = (size (obj.(pos), /dimensions))[1]
			components = [ 'x', 'y', 'z' ]
			for comp = 0, num_dims-1 do begin
				if (num_dims eq 3) then begin
					space = strpos (label, '_')
					if (space ge 0) then begin
						comp_label = strmid (label, 0, space) + components[comp] + strmid (label, space+1)
					end else begin
						comp_label = strmid (label, 0, strlen (label)-1) + components[comp]
					end
				end else begin
					comp_label = label + str (comp+1)
				end
				h5_write, group+'/'+comp_label, reform ((obj.(pos))[*,comp])
			end
		end else begin
			h5_write, group+'/'+label, obj.(pos)
		end
	end

	if (size (time, /type) eq 0) then time = obj.t
	if (size (obj.mass, /type) eq 4) then time = float (time)
	if (size (obj.mass, /type) eq 5) then time = double (time)
	h5_write, 'time', time

	h5_close_file
end

