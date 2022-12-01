#####################################################################################
##
##	class: IO
##
##	This scripts checks the existance of a folder an creates it if necessary
##
#####################################################################################


def mkdir(destination, rank=0, lfs=False, MB=1, count=1, comm=None):
    ##
    ## external modules
    import os
    import subprocess as sub

    ##
    ## create dir if not existing
    if not os.path.exists(destination):
        if comm:
            comm.Barrier()
        if rank == 0:
            os.makedirs(destination)
            if lfs:
                cmd = "lfs setstripe -S {}M -c {} ".format(MB, count) + destination
                process = sub.Popen(cmd.split(), stdout=sub.PIPE)
                output, error = process.communicate()
                print(cmd, output, error)
        if comm:
            comm.Barrier()
        return True
    return False
