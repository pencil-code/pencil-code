import os
from pathlib import Path
def main():
    current = Path.cwd()
    sample_dir  = current.parents[1]
    path = str(sample_dir)
    file = "equ.f90"
    PC_HOME = os.getenv("PENCIL_HOME")
    command = f"python3 {PC_HOME}/fortran-parser/parse.py -f rhs_cpu -F {path}/src/{file} -o -s --sample-dir {path} --dir {path}/src"
    os.system(command)
    os.system("mv mhdsolver-rhs.inc DSL/local/mhdsolver.ac")
    os.system("mv cparam.h DSL/local")
    os.system("mv static_var_declares.h DSL/local")
    os.system("cp DSL/solve_two.ac DSL/local")
    os.system("rm res-inlined.txt")

if __name__ == "__main__":
    main()
