import os
import numpy as np

from locations import location_results
from locations import results_to_save


def load_results():
    results = {key:np.load(os.path.join(location_results,
                                        f'{key}.npy'),
                           allow_pickle=True).item()
               for key in results_to_save}

    return results
