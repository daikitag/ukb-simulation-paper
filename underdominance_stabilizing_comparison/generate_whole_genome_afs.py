import pandas as pd
import pickle

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

# Save the dictionary as pickle object
with open("simulation_maf.pcl", "w") as f:
    pickle.dump(sim_data, f)

with open("thousand_maf.pcl", "w") as f:
    pickle.dump(thousand_data, f)
