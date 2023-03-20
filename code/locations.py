"""Locations of directories and files for the cancer epistasis project"""

import os

from string import ascii_uppercase as alphabet


if '__file__' not in globals():
    __file__ = '.'

location_data = os.path.abspath(
    os.path.join(os.path.dirname(__file__),
                 "../"
                 "data/"))
"""Location of directory that contains data for the model."""


location_results = os.path.abspath(
    os.path.join(os.path.dirname(__file__),
                 "../"
                 "results/"))
"""Location of directory that contains the results."""


location_figures = os.path.abspath(
    os.path.join(os.path.dirname(__file__),
                 "../"
                 "figures/"))
"""Location of directory that contains the figures."""


default_db_file_name = os.path.join(location_data, "tcga.luad.maf.txt")
"""Default file name to use as data source."""


gene_coordinates_file = os.path.join(location_data, "genes.xlsx")
"""Location of the genes coordinates file."""


default_mutation_names = list(alphabet)
"""Default mutation names, in case `mutation_names` is not provided
for the plotting functions. The default is A, B, C, etc.

"""


results_to_save = ["samples", "lambdas", "lambdas_cis", "mus",
                   "gammas", "gammas_cis"]
