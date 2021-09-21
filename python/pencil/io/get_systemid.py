def get_systemid():
    """This function returnes informations over the used system.

    Please add further functionallity by yourself!

    Return Values:
        system_name:  name of computer/cluster system as identified. Default return is 'pc'
                            Known systems are: - THEO (Garchingen)
                                               - HYDRA (Garchingen)
                                               - JUQUEEN
        raw_name:        raw system name
        job_name_key:    to identify job name in submit script, e.g. '-N', 'job_name'
        submit_scriptfile: default name of submit script
        submit_line:     batch line to submit the submit.csh/.sh"""

    import socket

    system_name = False
    raw_name = socket.gethostname()  # get computer name
    job_name_key = False
    submit_scriptfile = False
    submit_line = False

    if raw_name.startswith("isaac"):  # slurm on draco
        system_name = "isaac"
        job_name_key = "#SBATCH -J"
        submit_scriptfile = "submit.csh"
        submit_line = "llsubmit submit.csh"

    elif raw_name.startswith("hy"):  # loadleveler on hydra
        system_name = "hydra"
        job_name_key = "job_name"
        submit_scriptfile = "submit.csh"
        submit_line = "llsubmit submit.csh"

    elif raw_name.startswith("juqueen"):  # loadleveler on juqueen
        system_name = "juqueen"
        job_name_key = "job_name"
        submit_scriptfile = "submit.sh"
        submit_line = "llsubmit submit.sh"

    elif raw_name.startswith("draco"):  # loadleveler on juqueen
        system_name = "draco"
        job_name_key = "#SBATCH -J"
        submit_scriptfile = "submit.csh"
        submit_line = "llsubmit submit.csh"

    else:
        system_name = "pc"
        job_name_key = False
        submit_line = False

    return system_name, raw_name, job_name_key, submit_scriptfile, submit_line
