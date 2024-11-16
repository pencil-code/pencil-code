from ctypes import*
import odop
from odop.odop_obs import OdopObs
import os


def main():
    odop.start(config_file="./odop_conf.yaml",task_folder="./op-tasks")
    #odop.start(task_folder="./op-tasks")
    so_file = "./src/libPC.so"
    my_funcs = CDLL(so_file)
    my_funcs.run_start()
    odop.stop()

if __name__ == "__main__":
    main()

