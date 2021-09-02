;;
;; $Id$
;;
;;  Read yz-averages from file.
;;
pro pc_read_yzaver, object=object, varfile=varfile, datadir=datadir, dim=dim, grid=grid, $
    monotone=monotone, quiet=quiet, njump=njump, tmin=tmin, single=single

  pc_read_1d_aver, 'x', object=object, varfile=varfile, datadir=datadir, dim=dim, grid=grid, $
    monotone=monotone, quiet=quiet, njump=njump, tmin=tmin, single=single

end
