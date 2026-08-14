import pandas as pd
import pickle
import numpy as np


from glob import glob

pop_list = ["ceu", "chb", "yri", "jpt"]

sim_data = {}
thousand_data = {}

for pop in pop_list:
    sim_frame = []
    thousand_frame = []

    thousand_file_list = glob(
        f"/vols/bitbucket/tagami/ukb/raw_maf_data/thousand_pop_{pop}_chr_*_maf.txt"
    )
    sim_file_list = glob(
        f"/vols/bitbucket/tagami/ukb/raw_maf_data/sim_pop_{pop}_*_maf.csv"
    )

    for sim_file in sim_file_list:
        sim_frame.append(pd.read_csv(sim_file))

    for thousand_file in thousand_file_list:
        thousand_frame.append(
            pd.read_csv(thousand_file, header=None, names=["MAF"], sep="\t")
        )

    sim_final = pd.concat(sim_frame)
    thousand_final = pd.concat(thousand_frame)

    sim_final["plot_MAF"] = sim_final["MAF"].apply(lambda x: x if x < 0.5 else 1 - x)
    thousand_final["plot_MAF"] = thousand_final["MAF"].apply(
        lambda x: x if x < 0.5 else 1 - x
    )

    sim_data[pop] = sim_final.plot_MAF
    thousand_data[pop] = thousand_final.plot_MAF

boundary = [0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5]
groups = ["10-15%", "15-20%", "20-25%", "25-30%", "30-40%", "40-50%"]

frames = []

for pop in ["ceu", "yri", "chb", "jpt"]:
    sim_plot, _ = np.histogram(sim_data[pop], bins=boundary)
    thousand_plot, _ = np.histogram(thousand_data[pop], bins=boundary)

    afs_plot = pd.DataFrame(
        {
            "simulation": sim_plot,
            "1000 Genomes Project": thousand_plot,
            "groups": groups,
            "population": [pop] * len(sim_plot),
        },
    )
    frames.append(afs_plot)

final_data = pd.concat(frames, ignore_index=True)
final_data.to_csv("plot_data_afs.csv", index=False)
