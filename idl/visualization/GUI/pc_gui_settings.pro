;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_gui_settings.pro     ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  $Id$
;;;
;;;  Description:
;;;   User settings for the Pencil Code GUI..

;;; Physical quantities to be visualized
;;; Available quantities can be found and defined in 'pc_get_quantity'.
default, quantities, { $
	Temp:'temperature', $
	j_abs:'current density', $
	HR_ohm:'Ohmic heating rate', $
	HR_viscous:'viscous heating rate', $
	B_z:'magnetic field z', $
	u_abs:'velocity', $
	u_z:'velocity z', $
	P_therm:'thermal pressure', $
	grad_P_therm_abs:'thermal pressure gradient', $
	Rn_visc:'viscous Rn', $
	Rn_mag:'magnetic Rn', $
	log_rho:'logarithmic density', $
	n_rho:'particle density' $
}

;;; Quantities to be overplotted
;;; Available quantities can be found and defined in 'pc_get_quantity'.
default, overplot_quantities, { $
	u:'velocities', $
	j:'current density', $
;	grad_Temp:'temperature gradient', $
;	grad_P_therm:'thermal pressure gradient', $
	B:'magnetic_field', $
	A_contour:'fieldlines' $
}

;;; Initial varfile
default, varfile, 'var.dat'
default, crashfile, 'crash.dat'
default, pattern, 'VAR[0-9]*'

;;; Preferred units for display
default, display_units, { $
	default_length:1.e6, default_length_str:'Mm', $
	default_velocity:1.e3, default_velocity_str:'km/s', $
	default_time:1, default_time_str:'s', $
	default_temperature:1, default_temperature_str:'K', $
	default_density:1, default_density_str:'kg/m^3', $
	default_mass:1, default_mass_str:'kg', $
	default_magnetic_field:1.e-4, default_magnetic_field_str:'Gauß', $
	default_current_density:1, default_current_density_str:'A/m^2' $
}

