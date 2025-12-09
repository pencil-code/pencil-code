import argparse
import os
from pathlib import Path
def main():
    argparser = argparse.ArgumentParser(description="Wrapper",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    argparser.add_argument("--modify_source_code",default=False,action="store_true", help="Whether to modify the fortran source code to include the pushpars etc. automatically")
    current = Path.cwd()

    args = argparser.parse_args()
    config = vars(args)
    modify_source_code = config["modify_source_code"]
    sample_dir  = current.parents[1]
    path = str(sample_dir)
    file = "equ.f90"
    PC_HOME = os.getenv("PENCIL_HOME")
    command = f"python3 {PC_HOME}/fortran-parser/parse.py -f rhs_cpu -F {path}/src/{file} -o -s --sample-dir {path} --dir {path}/src"
    if(modify_source_code): command += " --modify_source_code"
    
    sources_are_different = False
    if not os.path.exists("PC_modulesources_old.h"):
        sources_are_different = True
    else:
        sources_are_different = os.system("diff PC_modulesources.h PC_modulesources_old.h")
    skip_generation = os.path.exists("DSL/local/mhdsolver.ac") and not sources_are_different
    if skip_generation:
        return

    os.system(command)
    os.system("rm res-inlined.txt")
    rhs_already_exists = os.path.exists("DSL/local/mhdsolver.ac") and not os.system("diff mhdsolver-rhs.inc DSL/local/mhdsolver.ac")
    if (rhs_already_exists): 
        return
    os.system("mv mhdsolver-rhs.inc DSL/local/mhdsolver.ac")
    os.system("mv GW-rhs.ac DSL/local/GW-rhs.h")
    os.system("mv cparam.h DSL/local")
    os.system("mv static_var_declares.h DSL/local")
    os.system("cp DSL/solve_two.ac DSL/local")

if __name__ == "__main__":
    main()
