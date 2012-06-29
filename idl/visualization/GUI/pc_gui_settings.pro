;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_gui_settings.pro     ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  $Id$
;;;
;;;  Description:
;;;   Default settings for the Pencil Code GUI
;;;   If you like to have different settings, please make a copy of this file
;;;   and name it: "pc_gui_user_settings.pro". There, you can make your changes.

;;; Load user-defined settings, if available.
pencil_home = getenv ('PENCIL_HOME')
if (pencil_home eq "") then begin
	message, "ERROR: please 'source sourceme.sh', before using this function."
end
if (file_test (pencil_home+"/idl/visualization/GUI/pc_gui_user_settings.pro")) then begin
	@pc_gui_user_settings
end

;;; Physical quantities to be visualized
;;; Available quantities can be found and defined in 'pc_get_quantity'.
default, quantities, { $
	Temp:'temperature', $
	j_abs:'current density', $
	B_x:'magnetic field x', $
	B_y:'magnetic field y', $
	B_z:'magnetic field z', $
	u_abs:'velocity', $
	u_x:'velocity x', $
	u_y:'velocity y', $
	u_z:'velocity z', $
	P_therm:'thermal pressure', $
	grad_P_therm_abs:'thermal pressure gradient', $
	rho:'density' $
}

;;; Quantities to be overplotted
;;; Available quantities can be found and defined in 'pc_get_quantity'.
default, overplot_quantities, { $
	u:'velocities', $
	grad_Temp:'temperature gradient', $
	grad_P_therm:'thermal pressure gradient', $
;	B:'magnetic_field', $
	A_contour:'fieldlines' $
}

;;; Default filenames
default, varfile, 'var.dat'
default, crashfile, 'crash.dat'
default, pattern, 'VAR[0-9]*'

;;; Preferred units for display
default, display_units, { $
	default_length:1, default_length_str:'m', $
	default_velocity:1, default_velocity_str:'m/s', $
	default_time:1, default_time_str:'s', $
	default_temperature:1, default_temperature_str:'K', $
	default_density:1, default_density_str:'kg/m^3', $
	default_mass:1, default_mass_str:'kg', $
	default_magnetic_field:1, default_magnetic_field_str:'Tesla', $
	default_current_density:1, default_current_density_str:'A/m^2' $
}

