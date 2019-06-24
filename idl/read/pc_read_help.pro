;
; $Id$
;
;  Help on Pencil-Code read routines
;
;  Author: Tony Mee (A.J.Mee@ncl.ac.uk)
;  $Date: 2005-10-24 01:47:28 $
;  $Revision: 1.4 $
;
;  27-nov-02/tony: coded
;  07-Apr-2019/PABourdin: renamed to 'pc_read_help'
;
pro pc_read_help
COMPILE_OPT IDL2,HIDDEN
;
;  Display list of routines
;
    print, "  ------ UNIVERSAL READ ROUTINES ------                                    "
    print, "pc_read(quantity, file=file)                                               "
    print, "    Reads data from a HDF5 snapshot (like 'var.h5').                       "
    print, ""
    print, "  ------ 'BIG-DATA'-COMPATIBLE READ ROUTINES ------                        "
    print, "pc_read_var_raw                                                            "
    print, "    Reads data from a snapshot (like 'var.dat').                           "
    print, ""
    print, "pc_read_subvol_raw                                                         "
    print, "    Reads any 3D subvolume from a snapshot.                                "
    print, ""
    print, "pc_read_slice_raw                                                          "
    print, "    Reads any 2D slice from a snapshot.                                    "
    print, ""
    print, "  ------ OLD-STYLE FILE READ ROUTINES ------                               "
    print, "pc_read_var                                                                "
    print, "    Reads a snapshot and constructs an IDL data structure from it.         "
    print, ""
    print, "  ------ PENCIL-CODE AUXILIARY READ ROUTINES ------                        "
    print, "pc_read_var_time                                                           "
    print, "    Reads only the time from a snapshot.                                   "
    print, ""
    print, "pc_read_dim                                                                "
    print, "    Reads the processor parameters into a data structure.                  "
    print, ""
    print, "pc_read_grid                                                               "
    print, "    Reads the grid setup into a data structure.                            "
    print, ""
    print, "pc_read_param                                                              "
    print, "    Reads the 'start.in' or 'run.in' parameters of a Pencil-Code run.      "
    print, ""
    print, "pc_read_ts                                                                 "
    print, "    Read the time series from 'time_series.dat' into a data structure.     "
    print, ""
    print, " Type 'routine_name, /HELP' for help or check the source code file headers."
    print, ""
    print, " NOTE: All read information will be held temporarily in memory.            "
    print, "       There is NO caching - all data is always reread from file.          "

end

