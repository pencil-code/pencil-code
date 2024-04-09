from ctypes import*
from odop.odop_obs import OdopObs
from qoa4ml.common import ODOP_PATH

so_file = "./src/libPC.so"
my_funcs = CDLL(so_file)

def main():
    odop_obs = OdopObs(config_path=ODOP_PATH + "config/odop_obs_conf.yaml")
    odop_obs.start()

    my_funcs.run_start()
    print("End from python")

    odop_obs.stop()

if __name__ == "__main__":
    main()

