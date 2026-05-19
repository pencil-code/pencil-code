import os
import sys
import argparse
def main():
    argparser = argparse.ArgumentParser(description="Update *.f90 files to be GPU ready",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    argparser.add_argument("-p", "--pencil-home", help="Where PC is installed", required=True)
    argparser.add_argument("-s", "--sample", help="Which sample", required=True)
    args = argparser.parse_args()
    config = vars(args)
    PC_HOME = config["pencil_home"]
    sample = config["sample"]
    os.environ["PENCIL_USER"] = "GPU_AUTOTEST"
    source_command = ""
    submodule_command = "git submodule update --init --remote"
    build_command = f"cd {PC_HOME}/samples/build-samples/{sample} && pc_setupsrc --force-astaroth && pc_build -f GNU-GCC_MPI+GNU-GCC_GPU+GNU-GCC_debug MODIFY_SOURCE_CODE=on"
    command = f"{submodule_command} && {build_command}"
    return os.system(command)

if __name__ == "__main__":
    sys.exit(main())
