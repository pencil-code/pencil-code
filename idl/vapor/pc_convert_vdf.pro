; $Id$
;
; Description:
;   Converts a given set of varfiles to a VDF file for Vapor.
;
; Parameters:
;   * varfiles       varfile name (or array of names) to convert.
;
; Optional parameters:
;   * coarsening     Number of coarsening levels (Default: 0 = off).
;   * reduce         Factor for reduction of the data (Default: 1 = off).
;   * output_dir     Output directory (Default: "data/vdf").
;   * vdf_files      Output file names (Default: "var.vdf" or "VAR#.vdf").
;   * varfiles       varfile name (or array of names) to convert.
;   * quantities     Quantity name(s) to write (Default: MHD = [u,rho,Temp,B]).
;                    More quantities are listed in "pc_check_quantities.pro":
;                    IDL> help, pc_check_quantities (/all), /str
;
; Examples:
; =========
;
;   Load "var.dat" and convert to a default VDF2 set:
;   IDL> pc_convert_vdf2
;
;   Load several varfiles and convert to a default VDF2 set:
;   IDL> pc_convert_vdf2, ['var.dat', 'VAR0', 'VAR1']
;
;   Load varfile and convert to VDF2 with a given set of quantities:
;   IDL> pc_convert_vdf2, 'VAR1', ['HR_viscous', 'HR_ohm', 'B_z']

pro pc_convert_vdf, varfiles, vdf_files=vdf_files, quantities=quantities, coarsening=coarsening, reduce=reduce, output_dir=output_dir, datadir=datadir, proc=proc, allprocs=allprocs, var_list=var_list, units=units, dim=dim, grid=grid, start_param=start_param, run_param=run_param, varcontent=varcontent

	default, varfiles, 'var.dat'
	default, output_dir, 'data/vdf'
	default, quantities, ['u_x','u_y','u_z','rho','Temp','B_x','B_y','B_z']
	default, coarsening, 0

	if (not file_test (output_dir, /directory)) then file_mkdir, output_dir

	num_files = n_elements (varfiles)
	reset = 1
	for pos = 0, num_files-1 do begin
		if (varfiles[pos] eq 'var.dat') then vdf_file = 'var.vdf' else vdf_file = varfiles[pos]+'.vdf'
		if (size (vdf_files, /type) eq 0) then begin
			vdf_files = output_dir+'/'+vdf_file
		end else begin
			if (n_elements (vdf_files) le pos) then vdf_files = [ vdf_files, output_dir+'/'+vdf_file ]
		end
		pc_read_var_raw, obj=var, tags=tags, varfile=varfile, var_list=var_list, datadir=datadir, proc=proc, allprocs=allprocs, dim=dim, grid=grid, start_param=start_param, run_param=run_param, quiet=(pos gt 0)
		pc_write_vdf, vdf_files[pos], var, tags=tags, timestep=pos, max_timesteps=num_files, coarsening=coarsening, reduce=reduce, quantities=quantities, units=units, dim=dim, grid=grid, start_param=start_param, run_param=run_param, varcontent=varcontent, reset=reset
	end
end

