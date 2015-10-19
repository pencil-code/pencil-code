; IDL autotest for essential IDL routines
;
; To start the script, just type at a shell prompt:
; $ idl -demo -e ".r reference"
;
; This test reads in a varfile, computes some derived physical quantity
; (Ohmic heating), and compares mean, min, and max to reference values.
;
; Tested IDL routines:
; ====================
; + lib/default
; + lib/derivatives/curlcurl
; + lib/derivatives/xder2
; + lib/derivatives/xder2_6th_ghost
; + lib/derivatives/xderyder
; + lib/derivatives/xderyder_6th_ghost
; + lib/derivatives/xderzder
; + lib/derivatives/xderzder_6th_ghost
; + lib/derivatives/yder2
; + lib/derivatives/yder2_6th_ghost
; + lib/derivatives/yderxder
; + lib/derivatives/yderxder_6th_ghost
; + lib/derivatives/yderzder
; + lib/derivatives/yderzder_6th_ghost
; + lib/derivatives/zder
; + lib/derivatives/zder_6th_ghost
; + lib/derivatives/zder2
; + lib/derivatives/zder2_6th_ghost
; + lib/derivatives/zderxder
; + lib/derivatives/zderxder_6th_ghost
; + lib/derivatives/zderyder
; + lib/derivatives/zderyder_6th_ghost
; + lib/functions/all
; + lib/my_rebin
; + lib/operators/dot
; + lib/operators/dot2
; + lib/spread
; + lib/undefine
; + pc_check_quantities
; + pc_get_quantity
; + pc_get_parameter
; + pc_set_precision
; + read/pc_get_datadir
; + read/pc_read_dim
; + read/pc_read_grid
; + read/pc_read_param
; + read/pc_read_var_raw
; + read/pc_varcontent
;
; 19-Oct-2015/PABourdin: initial version

; Settings:

default, tolerance, 1d-10
default, reference_file, 'reference.xdr'

; Code:

HR_ohm = pc_get_quantity ('HR_ohm', 'var.dat')

test_mean = mean (HR_ohm)
test_min = min (HR_ohm)
test_max = max (HR_ohm)

restore, reference_file

success = 1

if (abs (test_mean / ref_mean - 1.0) gt tolerance) then begin
	print, "ERROR: mean differs: ", test_mean, " != ", ref_mean
	success = 0
end

if (abs (test_min / ref_min - 1.0) gt tolerance) then begin
	print, "ERROR: min differs: ", test_min, " != ", ref_min
	success = 0
end

if (abs (test_max / ref_max - 1.0) gt tolerance) then begin
	print, "ERROR: max differs: ", test_max, " != ", ref_max
	success = 0
end

if (success) then print, "IDL autotest: >SUCCESS<"

END

