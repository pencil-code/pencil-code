;;
;; $Id$
;;
;;  Read xz-averages from file.
;;
pro pc_read_xzaver, object=object, varfile=varfile, datadir=datadir, $
    monotone=monotone, quiet=quiet

  pc_read_1d_aver, 'y', object=object, varfile=varfile, datadir=datadir, $
    monotone=monotone, quiet=quiet

end
