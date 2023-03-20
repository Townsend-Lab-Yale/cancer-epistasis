import os
import numpy as np

from importing_data_from_tcga import filter_db_by_mutation

from count_combinations import most_commonly_mutated_genes
from count_combinations import count_pts_per_combination
from count_combinations import convert_samples_to_dict

from cancer_epistasis import estimate_lambdas
from cancer_epistasis import convert_lambdas_to_dict
from cancer_epistasis import compute_gammas
from cancer_epistasis import asymp_CI_lambdas

from mutation_rates import convert_mus_to_dict

from locations import location_results
from locations import results_to_save

from load_results import load_results

from print_results import print_table



default_genes = ['TP53', 'KRAS', 'STK11', 'LRP1B']



def produce_results(genes=default_genes, save=True):

    print("Importing data set...")
    db = filter_db_by_mutation()
    print("...done")
    print("")


    print("Computing most commonly mutated genes in data set...")
    most_commonly_mutated_genes(db)
    print("...done")
    print("")

    results = {}

    print("Counting samples on each genotype...")
    samples = count_pts_per_combination(db, mutations=genes)
    results["samples"] = convert_samples_to_dict(samples)
    print("...done")
    print("")

    print("Estimating fluxes MLE...")
    mle = estimate_lambdas(samples, draws=1)
    results["lambdas"] = convert_lambdas_to_dict(mle)
    print("...done")
    print("")

    print("Computing fluxes confidence intervals (CI)...")
    results["lambdas_cis"] = convert_lambdas_to_dict(
        asymp_CI_lambdas(mle['lambdas'],
                         samples,
                         print_progress=True))
    print("...done")

    print("Importing mutation rates...")
    results["mus"] = convert_mus_to_dict(genes)
    print("...done")
    print("")

    print("Computing selection coefficients...")
    results["gammas"] = compute_gammas(results["lambdas"],
                                       results["mus"])
    results["gammas_cis"] = compute_gammas(results["lambdas_cis"],
                                           results["mus"])
    print("...done")

    if save:
        print("")
        print("Saving results...")
        for key, value in results.items():
            np.save(os.path.join(location_results,
                                 f'{key}.npy'), value)
        print("...done")

    return results


def main(force_produce_results=False,
         genes=default_genes):

    if (not all([os.path.exists(os.path.join(location_results,
                                            f"{key}.npy"))
                 for key in results_to_save])
        or force_produce_results):
        results = produce_results(genes)
    else:
        print("Loading results previously computed...")
        print("(If you want produce the results again set force_produce_results to True)")
        results = load_results()
        print("...done")
        print("")

    print("Printing table to copy to LaTeX main file...")
    print("")
    print_table(results, genes)
    print("")
    print("done.")
    print("")

    return results

if __name__ == "__main__":
    main()
