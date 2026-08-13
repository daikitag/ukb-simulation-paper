import numpy as np
import pandas as pd
import tskit
import sys


def obtain_mutation_df(ts):
    variant = tskit.Variant(ts)
    # All mutations are happening at the same time, so we can obtain the node for the 1st
    # mutation
    # Time is in forward simulation
    mutation_list = []
    for i in range(ts.num_mutations):
        mutation = ts.mutation(i)
        variant.decode(mutation.site)
        alleles = np.array(variant.alleles)
        genotype = alleles[variant.genotypes]
        allele_freq = np.count_nonzero(genotype == mutation.derived_state) / len(
            variant.genotypes
        )

        mutation_list.append(
            {
                "site_id": mutation.site,
                "selection_coeff": mutation.metadata["mutation_list"][0][
                    "selection_coeff"
                ],
                "causal_allele": mutation.derived_state,
                "trait_id": 0,
                "allele_freq": allele_freq,
            }
        )

    mutation_df = pd.DataFrame(mutation_list)
    mutation_df["selection_coeff"] = mutation_df["selection_coeff"] * 1e10 * (-2)

    return mutation_df


# v_s is the expected contribution at sites under strong selection, which is defined as the following:
def vs(N, w2=1, n=1):
    """
    This determines the expected contribution to genetic variance per site under strong selection.
    N is the population size
    w is the constant that influences the trait distribution and the default is 1.
    n is the degree of pleiotropy, and the default is 1.
    """
    return 2 * w2 / (n * N)


def single_simulate_effect(selection_coeff, w2, rng):
    effect_size = (2 * rng.binomial(1, p=0.5) - 1) * np.sqrt(w2 * selection_coeff)
    return effect_size


def single_simulate_effect_mutation_df(mutation_df, w2, seed=None):
    rng = np.random.default_rng(seed=seed)
    mutation_df["effect_size"] = mutation_df.apply(
        lambda row: single_simulate_effect(row["selection_coeff"], w2=w2, rng=rng),
        axis=1,
    )


def compute_v(mutation_df):
    mutation_df["v"] = (
        2
        * (mutation_df["effect_size"] ** 2)
        * mutation_df["allele_freq"]
        * (1 - mutation_df["allele_freq"])
    )


def single_ts_simulate_v(ts, w2, N, threshold, seed):
    mutation_df = obtain_mutation_df(ts)
    single_simulate_effect_mutation_df(mutation_df, w2=w2, seed=seed)
    compute_v(mutation_df)
    strong_selection_df = mutation_df[
        mutation_df["selection_coeff"] * 2 * N > threshold
    ]
    strong_selection_df["scaled_v"] = strong_selection_df["v"] / vs(N=N, w2=w2, n=1)
    return strong_selection_df


def compute_proportion(variance, xaxis):
    percentage = []
    variance = np.array(variance)
    total_var = np.sum(variance)

    for v in xaxis:
        variance = variance[variance > v]
        proportion = np.sum(variance) / total_var
        percentage.append(proportion)

    return percentage


def simulate_pleiotropy_effect(selection_coeff, n, w2, rng):
    effect_size = rng.normal(loc=0, scale=np.sqrt(w2 / n * selection_coeff))
    return effect_size


def sim_pleiotropy(mutation_df, n, N, w2, seed=None):
    rng = np.random.default_rng(seed=seed)
    mutation_df["effect_size"] = mutation_df.apply(
        lambda row: simulate_pleiotropy_effect(
            row["selection_coeff"], n=n, w2=w2, rng=rng
        ),
        axis=1,
    )


def pleiotropic_ts_simulate_v(ts, w2, n, N, threshold, seed):
    mutation_df = obtain_mutation_df(ts)
    sim_pleiotropy(mutation_df, n=n, N=N, w2=w2, seed=seed)
    compute_v(mutation_df)
    strong_selection_df = mutation_df[
        mutation_df["selection_coeff"] * 2 * N > threshold
    ]
    strong_selection_df["scaled_v"] = strong_selection_df["v"] / vs(N=N, w2=w2, n=n)
    return strong_selection_df


def main():
    ts = tskit.load("underdominance_popsize_1000_U_1_S_10_seed_25.tree")

    w2 = 1
    n = int(sys.argv[1])
    seed = int(sys.argv[2])
    N = 1000
    threshold = 30
    xaxis = np.linspace(0, 5, num=50)
    simulation_result = np.zeros(50)
    for _ in range(100):
        strong_selection_df = pleiotropic_ts_simulate_v(
            ts, w2=w2, n=n, N=N, threshold=threshold, seed=seed
        )
        simulation_result += np.array(
            compute_proportion(strong_selection_df["scaled_v"], xaxis)
        )
        seed += 1

    simulation_result /= 100

    result_df = pd.DataFrame({"xaxis": xaxis, "data": simulation_result})

    result_df.to_csv(
        f"var_propseed_{seed}_n_{n}_underdominance_popsize_1000_U_1_S_10_seed_25.csv",
        index=False,
    )


if __name__ == "__main__":
    main()
