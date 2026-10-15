# Script for plotting loss from a specified input file to a specified output file.


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys


def plot_loss(model_name_or_desc, input_filepath, output_filepath):
    df = pd.read_csv(input_filepath, sep=",")

    df["ln_train_loss"] = np.log(df["train_loss"])

    ax = df.plot(x="epoch", y="ln_train_loss", kind="line", title="Training Loss for "+model_name_or_desc,
        xlabel="Iteration", ylabel="Training Loss")

    ax.get_figure().savefig(output_filepath)

if __name__ == "__main__":
    if (len(sys.argv) < 4):
        print("Please specify model name or description, input file path (to file containing losses), and output file path (to file [does not need to exist yet; will be overwritten if it does] to contain plot)")
        sys.exit()

    model_name_or_desc = sys.argv[1]
    input_filepath = sys.argv[2]
    output_filepath = sys.argv[3]

    plot_loss(model_name_or_desc, input_filepath, output_filepath)
