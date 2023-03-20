import pandas as pd
import numpy as np

from importing_data_from_tcga import filter_db_by_mutation
from importing_data_from_tcga import genes

from theory import build_S_as_array

default_M = 4
"""Default number of mutations"""


def most_commonly_mutated_genes(db,
                                number_to_output=default_M,
                                print_all_results=True):
    """Give the most commonly mutated genes in the `db`.

    :type db: pandas.core.frame.DataFrame or str
    :param db: Pandas DataFrame with patient IDs and mutations. As
        returned by :func:`filter_db_by_mutation`. If a str, it is
        passed to `filter_db_by_mutation`.

    :type n: int
    :param n: Number of most commonly mutated genes to output. Default
        4.

    :rtype: list
    :return: List with the most commonly mutated genes.

    """

    mutations_per_gene = {gene:len(set(
        db[
            (db['Start_Position'] >= genes.loc[gene, 'start']) &
            (db['Start_Position'] <= genes.loc[gene, 'end']) &
            (db['Chromosome'] == str(genes.loc[gene, 'chromosome']))]
        ['Patient ID']))
                          for gene in genes.index}

    mutations_per_gene = sorted(mutations_per_gene.items(),
                                key=lambda x: x[1], reverse=True)

    if print_all_results:
        for gene, value in mutations_per_gene:
            print(f"Patients with mutations in {gene}: {value}")

    return [key for key, value in mutations_per_gene[:number_to_output]]


def count_pts_per_combination(db,
                              mutations=None,
                              tumor_col_name=None,
                              print_info=True,
                              patient_id_col_name=None):
    """Compute samples numbers for each mutation combinations in S.

    :type db: pandas.core.frame.DataFrame or str
    :param db: Pandas DataFrame with patient IDs and mutations. As
        returned by :func:`filter_db_by_mutation`. If a str, it is
        passed to `filter_db_by_mutation`.

    :type mutations: list or int or NoneType
    :param mutations: List with the mutation names to extract
        information about. This list determines the total number of
        mutations, M, and thus the space of mutation combinations,
        S. The order in the list is also important as the vectors in S
        correspond to that order. Each item on the list can be a
        string, or a list of strings. In the latter case (NOT
        IMPLEMENTED YET), those mutations will be aggregated into a
        single category. If instead of a list an int M is passed, pick
        the M most common mutations from the database, and if None
        (default), pick M to be equal to :const:`default_M`.

    :type tumor_col_name: str
    :param tumor_col_name: Name of the column in the data that
        contains the tumor allele. If None, try to infer it from the
        data.

    :type patient_id_col_name: str
    :param patient_id_col_name: Name of the column in the data that
        contains the patient identification. If None, try to infer it
        from the data.

    :rtype: tuple
    :return: A tuple with a one dimensional array of size 2^M with the
        computed sample numbers per mutation combination (ordered as
        items in S), M as in integer and S as a numpy array.

    """

    if isinstance(db, str):
        db = filter_db_by_mutation(
            db,
            tumor_col_name=tumor_col_name,
            patient_id_col_name=patient_id_col_name)

    if mutations is None:
        mutations = default_M

    if not isinstance(mutations, list):
        mutations = most_commonly_mutated_genes(db, number_to_output=mutations)

    M = len(mutations)

    S = build_S_as_array(M)

    if print_info:
        print("Gene labels")
        for k, mutation in enumerate(mutations):
            print(f"  {np.eye(1, M, k, dtype=int)[0]} : {mutation}")

    pts_per_mutation = [
        set(db[
            (db['Start_Position'] >= genes.loc[mutation, 'start']) &
            (db['Start_Position'] <= genes.loc[mutation, 'end']) &
            (db['Chromosome'] == str(genes.loc[mutation, 'chromosome']))]
            ['Patient ID'])
        if mutation in genes.index else
        set(db[db['Mutation'] == mutation]['Patient ID'])
        for mutation in mutations]

    pts_per_combination = (
        [len(db['Patient ID'].unique())]
        + [len(set.intersection(*[pts_per_mutation[i]
                                  for i in indices]))
           for indices in [[i for i in range(M)
                            if Sj[i] == 1]
                           for Sj in S[1:]]])

    pts_per_combination = np.array(pts_per_combination)

    for m in range(M, 0, -1):
        for k, y in enumerate(S):
            if np.sum(y) == m:
                in_y_indices = np.array(
                    [True
                     if np.sum(y >= S[i]) == M
                     and np.sum(y - S[i]) >= 1
                     else False
                     for i in range(len(S))])
                pts_per_combination[in_y_indices] = (
                    pts_per_combination[in_y_indices]
                    - pts_per_combination[k])

    if print_info:
        print("")
        print("Counts per genotype:")
        for x, pts in zip(S, pts_per_combination):
            print(f"  {x} : {pts}")

    return pts_per_combination


def convert_samples_to_dict(samples):
    """Convert a samples array to a dictionary.

    :type samples: list
    :param samples: Number of patients in each mutation combination
        (as returned by :func:`count_pts_per_combination`).

    :rtype: dict
    :return: Dictionary with the samples, indexed by tuples of 1's and
        0's representing whether the mutation occur fot the gene or
        not.

    """

    M = int(np.log2(len(samples))) # 2^M is the number mutation combinations

    S = build_S_as_array(M)

    results_as_dict = {tuple(x):value
                       for x, value in zip(S, samples)}

    return results_as_dict
