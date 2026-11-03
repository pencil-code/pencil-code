# Script for plotting training loss and brms

import matplotlib.pyplot as plt
import numpy as np
import pencil as pc
import sys


def plot_training_loss_and_brms(name_or_desc, input_dir, input_filename, output_filepath):
    ts = pc.read.ts(file_name=input_filename, datadir=input_dir)

    plt.plot(ts.t, np.log(ts.loss), label="ln(Training Loss)", color="blue")
    
    plt.plot(ts.t, np.log(ts.brms), label="ln(brms)", color="green")

    plt.legend()

    plt.title("Training Loss and brms for " + name_or_desc)
    plt.xlabel("t")
    plt.ylabel("ln(Training Loss) or ln(brms)")

    plt.savefig(output_filepath)


if __name__ == "__main__":
    if len(sys.argv) != 5:
        # TODO: Make error message more descriptive to specify which parameters
        #   need to be included
        print("Please specify input parameters")
        sys.exit()

    name_or_desc = sys.argv[1]
    input_dir = sys.argv[2]
    input_filename = sys.argv[3]
    output_filepath = sys.argv[4]

    plot_training_loss_and_brms(name_or_desc, input_dir, input_filename, output_filepath)
