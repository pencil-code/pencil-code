def rename_in_submit_script(
    new_name, submit_script_path=False, sim=False, old_name=False
):
    import os
    from os.path import exists, join, abspath, dirname
    from pencil.io import timestamp, get_systemid, mkdir

    if submit_script_path != False:
        path = dirname(abspath(submit_script_path))
        filename = basename(abspath(submit_script_path))

    # locate submit script file if not given
    if submit_script_path == False:
        # get path of simulation folder where submit script should be around
        if sim == False:
            path = "."
        else:
            path = sim.path

        list_of_possible_submit_scripts = [
            f
            for f in os.listdir(path)
            if f.startswith("submit") and f.split(".")[-1] in ["sh", "csh"]
        ]

        # only works if a single submit script could be identified, else prompt error and let the user do it manually
        if len(list_of_possible_submit_scripts) == 0:
            print("!! ERROR: Could not find a submit script in " + str(path))
            return False
        elif len(list_of_possible_submit_scripts) > 1:
            print(
                "!! ERROR: Could not identify submit script, please specify manually:"
                + str(list_of_possible_submit_scripts)
            )
            return False
        else:
            filename = list_of_possible_submit_scripts[0]

    # path to submit script should now be clear, but better check its a string
    if type(path) != type("STRING"):
        print(
            "!! ERROR: Could not identify submit script path, please check manually: "
            + str(path)
        )

    path_filename = join(path, filename)

    # open submit script as f and read content into s
    with open(path_filename) as f:
        s = f.read()

    # if the old name is known, we can simply replace that string with the new name
    if old_name != False and type(old_name) == type("STRING"):
        with open(path_filename) as f:
            s = f.read()
            if old_name in s:
                s = s.replace(old_name, new_name)
            else:
                print(
                    "?? ERROR: Could not find old_name "
                    + str(old_name)
                    + " in submit script "
                    + str(path_filename)
                )
                return False

    # else we need to look for specific identifiers, that differ from queue system and cluster
    else:
        # get submit name line identifier
        identify = get_systemid()[2]
        if identify == False:
            print(
                "!! ERROR: Could not identify an submit script name identifier, please update pc.io.get_systemid.py by adding your machine."
            )
        if identify in s:
            if s.count(identify) > 1:
                print(
                    "ERROR: Job name identifier has multiple appearences in submit script!"
                )
                return False
            s = s.split("\n")
            for ii, line in enumerate(s):
                if identify in line:
                    break
            s[ii] = identify + " " + new_name

            s = "\n".join(s)

        else:
            print(
                "!! ERROR: Could not find name identifier in submit script, identifier is "
                + str(identify)
            )

    # s should now be updated, we can save it in now in the file submit_script_path, but let us backup the old one
    # backup
    from shutil import copyfile

    target_dir = join(path, "pc/backups/")
    mkdir(target_dir)
    copyfile(path_filename, join(target_dir, filename + ".BAK" + str(timestamp())))

    # save new submit script
    with open(path_filename, "w") as f:
        f.write(s)

    # done
    return True
