;;
;; $Id$
;;
;;  Read xy-averages from file.
;;
pro pc_read_xyaver, object=object, varfile=varfile, datadir=datadir, $
    monotone=monotone, quiet=quiet, njump=njump, tmin=tmin

  pc_read_1d_aver, 'z', object=object, varfile=varfile, datadir=datadir, $
    monotone=monotone, quiet=quiet, njump=njump, tmin=tmin

end
