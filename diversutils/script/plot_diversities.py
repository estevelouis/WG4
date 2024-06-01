#!/usr/bin/python3

import os
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

COLUMNS_TO_IGNORE = {
    "num_active_files",
    "num_files",
    "num_active_sentences",
    "num_sentences",
    "num_all_sentences",
    "num_documents",
    "num_discarded_types",
    "w2v",
}

#DEFAULT_COLUMN = "num_documents"
DEFAULT_COLUMN = "n"

ENABLE_LOG_X = True
ENABLE_LOG_Y = True

ENABLE_Z_NORMALISATION = True

def main(dir_input):
    subdir_list = filter(lambda x: os.path.isdir(x), os.listdir(dir_input))
    for subdir in subdir_list:
        current_abs_path_subdir = os.path.abspath(subdir)
        print(f"Working with directory {repr(current_abs_path_subdir)}")

        abs_path_output_scores = None
        abs_path_output_timing = None
        abs_path_output_memory = None

        file_list = filter(lambda x: os.path.isfile(os.sep.join([subdir, x])), os.listdir(current_abs_path_subdir))
        for file in file_list:
            current_abs_path_file = os.path.abspath(os.sep.join([subdir, file]))
            """
            if __debug__:
                print(f"Detected file: {repr(current_abs_path_file)}")
            """

            if file.endswith(".tsv"):
                if "timing" in file or "_time" in file:
                    abs_path_output_timing = current_abs_path_file
                elif "memory" in file:
                    abs_path_output_memory = current_abs_path_file
                else:
                    abs_path_output_scores = current_abs_path_file

        print(f"Inferred file for scores: {repr(abs_path_output_scores)}")
        print(f"Inferred file for timing: {repr(abs_path_output_timing)}")
        print(f"Inferred file for memory: {repr(abs_path_output_memory)}")

        if abs_path_output_scores == None or abs_path_output_timing == None or abs_path_output_memory == None:
            continue

        df_scores = pd.read_csv(abs_path_output_scores, sep='\t')
        df_timing = pd.read_csv(abs_path_output_timing, sep='\t')
        df_memory = pd.read_csv(abs_path_output_memory, sep='\t')

        print(df_scores)
        print(df_timing)
        print(df_memory)

        for df, paradigm in zip([df_scores, df_timing, df_memory], ["scores", "timing", "memory"]):
            fig, ax = plt.subplots()
            legend = []

            if paradigm in ["timing", "memory"]:
                # plotting special distributions
    
                special_distribution_args = {"linewidth": 4}
    
                """
                # plotting O(n^n)
                ax.plot(df[DEFAULT_COLUMN], np.power(df[DEFAULT_COLUMN], df[DEFAULT_COLUMN]), **special_distribution_args)
                legend.append(r"$O\left(n^n\right)$")
    
                # plotting O(2^n)
                ax.plot(df[DEFAULT_COLUMN], np.power(2.0, df[DEFAULT_COLUMN]), **special_distribution_args)
                legend.append(r"$O\left(2^n\right)$")
                """
    
                # plotting O(n^2)
                ax.plot(df[DEFAULT_COLUMN], np.power(df[DEFAULT_COLUMN], 2.0), **special_distribution_args)
                legend.append(r"$O\left(n^2\right)$")
    
                # plotting O(n ln n)
                ax.plot(df[DEFAULT_COLUMN], df[DEFAULT_COLUMN] * np.log(df[DEFAULT_COLUMN]), **special_distribution_args)
                legend.append(r"$O\left(n \ln\left(n\right)\right)$")
    
                # plotting O(n)
                ax.plot(df[DEFAULT_COLUMN], df[DEFAULT_COLUMN], **special_distribution_args)
                legend.append(r"$O\left(n\right)$")
    
                # plotting O(sqrt n)
                ax.plot(df[DEFAULT_COLUMN], np.sqrt(df[DEFAULT_COLUMN]), **special_distribution_args)
                legend.append(r"$O\left(\sqrt{n}\right)$")
    
                # plotting O(ln n)
                ax.plot(df[DEFAULT_COLUMN], np.log(df[DEFAULT_COLUMN]), **special_distribution_args)
                legend.append(r"$O\left(\ln\left(n\right)\right)$")
    
                # plotting O(1)
                ax.plot(df[DEFAULT_COLUMN], np.ones(len(df[DEFAULT_COLUMN])), **special_distribution_args)
                legend.append(r"$O\left(1\right)$")

            # plotting observed distributions

            for column in df.columns:
                if column in COLUMNS_TO_IGNORE:
                    continue
                if ENABLE_Z_NORMALISATION:
                    ax.plot(df[DEFAULT_COLUMN], (df[column] - np.mean(df[column])) / np.std(df[column]))
                else:
                    ax.plot(df[DEFAULT_COLUMN], df[column])
                legend.append(column)

            fig.legend(legend)
            plt.title(paradigm)
            plt.xlabel(DEFAULT_COLUMN)
            plt.ylabel("Other columns")
            if ENABLE_LOG_X:
                ax.set_xscale("log")
            if ENABLE_LOG_Y:
                ax.set_yscale("log")
            fig.tight_layout()
            plt.show()

    return 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir_input")

    args = parser.parse_args()

    exit(
        main(
            dir_input = args.dir_input
        )
    )
