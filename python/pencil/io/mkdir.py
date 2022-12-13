#####################################################################################
##
##	class: IO
##
##	This scripts checks the existance of a folder an creates it if necessary
##
#####################################################################################


def mkdir(destination, lfs=False, MB=1, count=1):
    ##
    ## external modules
    import os
    import subprocess as sub

    ##
    ## create dir if not existing
    if not os.path.exists(destination):
        try:
            os.makedirs(destination)
        except:
            if os.path.exists(destination):
                pass
        if lfs:
            cmd = "lfs setstripe -S {}M -c {} ".format(MB, count) + destination
            process = sub.Popen(cmd.split(), stdout=sub.PIPE)
            output, error = process.communicate()
            print(cmd, output, error)
        return True
    return False
