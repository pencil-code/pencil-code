# Script for plotting dynamo growth rate from specified input file(s) to a specified output file

import matplotlib.pyplot as plt
import numpy as np
import pencil as pc
import sys


def plot_dynamo_growth_rate(name_or_desc, input_dir, input_filename, label, input_dir_2, input_filename_2, label_2, output_filepath):    
    ts = pc.read.ts(file_name=input_filename, datadir=input_dir)
    plt.plot(ts.t, np.log(ts.brms), label=label, color="blue")
    
    if input_dir_2 is not None and input_filename_2 is not None and label_2 is not None:
        ts_2 = pc.read.ts(input_filename_2, datadir=input_dir_2)
        plt.plot(ts_2.t, np.log(ts_2.brms), label=label_2, color="green")

    plt.legend()

    plt.title("Dynamo Growth Rate for " + name_or_desc)
    plt.xlabel("t")
    plt.ylabel("ln(brms)")

    plt.savefig(output_filepath)


if __name__ == "__main__":
    if len(sys.argv) != 6 and len(sys.argv) != 9:
        # TODO: Make error message more descriptive to specify which parameters
        #   need to be included
        print("Please specify input parameters")
        sys.exit()

    name_or_desc = sys.argv[1]
    input_dir = sys.argv[2]
    input_filename = sys.argv[3]
    label = sys.argv[4]

    input_dir_2 = None
    input_filename_2 = None
    label_2 = None
    if len(sys.argv) == 6:
        output_filepath = sys.argv[5]
    else: # len(sys.argv) == 9:
        input_dir_2 = sys.argv[5]
        input_filename_2 = sys.argv[6]
        label_2 = sys.argv[7]
        output_filepath = sys.argv[8]

    plot_dynamo_growth_rate(name_or_desc, input_dir, input_filename, label, input_dir_2, input_filename_2, label_2, output_filepath)
