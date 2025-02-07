from ctypes import*
import odop
from odop.odop_obs import OdopObs
import os
import time


def main():
    os.environ["LD_PRELOAD"] = ""
    #odop.start(config_file="./odop_conf.yaml",task_folder="./empty-op-tasks")
    odop.start(config_file="./odop_conf.yaml",task_folder="./op-tasks")

    #time.sleep(20)
    so_file = "./src/libPC.so"
    my_funcs = CDLL(so_file)
    my_funcs.run_start()
    odop.stop()

if __name__ == "__main__":
    main()

