;;
;; $Id$
;;
;;  Read yz-averages from file.
;;
pro pc_read_yzaver, object=object, varfile=varfile, datadir=datadir, $
    monotone=monotone, quiet=quiet

  pc_read_1daver, 'x', object=object, varfile=varfile, datadir=datadir, $
    monotone=monotone, quiet=quiet

end
