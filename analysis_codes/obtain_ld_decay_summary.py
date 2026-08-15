import numpy as np
import pandas as pd
import pickle

from glob import glob


def ld_decay_meanbin(df, bin1=10, bin2=100, break_bp=100):
    short_mask = df["#Dist"] < break_bp
    bin_size = np.where(short_mask, bin1, bin2)

    df["bin_size"] = bin_size
    df["bin_idx"] = np.floor((df["#Dist"] - 0.1) / df["bin_size"]).astype(int)

    df["x_bp"] = (df["bin_idx"] + 1) * df["bin_size"]
    df["r2_x_count"] = df["Mean_r^2"] * df["NumberPairs"]

    binned = (
        df.groupby(["bin_size", "bin_idx", "x_bp"], as_index=False)
        .agg(
            sum_r2=("r2_x_count", "sum"),
            n_pairs=("NumberPairs", "sum"),
        )
        .sort_values("x_bp")
    )

    binned["mean_r2"] = binned["sum_r2"] / binned["n_pairs"]
    return binned[["x_bp", "mean_r2", "n_pairs"]]


def main():
    thousand_df = {}
    sim_df = {}

    dir = "/data/smew01/not-backed-up/scratch/tagami/ukb-simulation/output/slim_mu_1.7e-10_demography_relate_OutOfAfrica_4J17_ceu_slim_seed_100_chr_1p/data_process/ld_decay/"

    for pop in ["ceu", "yri", "chb", "jpt"]:
        sim_file = glob(dir + f"{pop}_sim*")[0]
        thousand_file = glob(dir + f"{pop}_thousand*")[0]
        sim_original_df = pd.read_csv(sim_file, sep="\t")
        thousand_original_df = pd.read_csv(thousand_file, sep="\t")

        thousand_df[pop] = ld_decay_meanbin(thousand_original_df)
        sim_df[pop] = ld_decay_meanbin(sim_original_df)

    with open("thousand_ld_decay_plot_chr1p.pcl", "wb") as f:
        pickle.dump(thousand_df, f)

    with open("simulation_ld_decay_plot_chr1p.pcl", "wb") as f:
        pickle.dump(sim_df, f)


if __name__ == "__main__":
    main()
