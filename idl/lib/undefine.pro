;+
; NAME:
;       UNDEFINE
;
; PURPOSE:
;       The purpose of this program is to delete or undefine
;       an IDL program variable from within an IDL program or
;       at the IDL command line. It is a more powerful DELVAR.
;
; AUTHOR:
;       FANNING SOFTWARE CONSULTING
;       David Fanning, Ph.D.
;       2642 Bradbury Court
;       Fort Collins, CO 80521 USA
;       Phone: 970-221-0438
;       E-mail: davidf@dfanning.com
;       Coyote's Guide to IDL Programming: http://www.dfanning.com
;
; CATEGORY:
;       General Programming
;
; CALLING SEQUENCE:
;       UNDEFINE, variable
;
; REQUIRED INPUTS:
;       variable: The variable to be deleted.
;
; SIDE EFFECTS:
;       The variable no longer exists.
;
; EXAMPLE:
;       To delete the variable "info", type:
;
;        IDL> Undefine, info
;
; MODIFICATION HISTORY:
;       Written by David Fanning, 8 June 97, from an original program
;       given to me by Andrew Cool, DSTO, Adelaide, Australia.
;       mgs, 26 Oct 1998: - added test for undefined variable
;-
PRO UNDEFINE, varname
   On_Error, 1
   IF N_Params() EQ 0 THEN $
      Message, 'One argument required in call to UNDEFINE'
   if (n_elements(varname) eq 0) then return ; already undefined
   tempvar = SIZE(TEMPORARY(varname))
END
