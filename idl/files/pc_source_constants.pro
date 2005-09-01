;;
;;  $Id: pc_source_constants.pro,v 1.1 2005-09-01 10:48:56 ajohan Exp $
;;
;;  Read Pencil Code constants from an arbitrary directory 'datadir'.
;;  Usage: @pc_source_constants
;;
cmd = 'perl -000 -ne '+"'"+'s/[ \t]+/ /g; print join(" & ",split(/\n/,$_)),     "\n"'+"' "+datadir+'/pc_constants.pro'
spawn, cmd, result
res = flatten_strings(result)
if (execute(res) ne 1) then print, 'There was a problem with pc_constants.pro'
