import pencil as pc

pc.io.sim2h5(
             newdir='/users/fagent/pencil-code/python/tutorials/mps/ism_binary_h5',
             laver2D=False,
             lvars=True,
             lvids=False,
             laver=False,
             snap_by_proc=True, 
             lremove_old_snapshots=False,
             lremove_old_slices=False,
             lremove_old_averages=False, 
             execute=False,
             lremove_deprecated_vids=True,
             lread_all_videoslices=False,
             vlarge=100000000,
             quiet=True,
            )

