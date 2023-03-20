# Results from cancereffectsizer for LUAD

mutation_rates = {'KRAS':1.992647e-08,
                  'TP53':1.808964e-07,
                  'LRP1B':1.559768e-07,
                  'BRAF':6.277220e-09,
                  'STK11':2.543602e-08}

def convert_mus_to_dict(genes):
    """Convert the dictionary with :const:`mutation_rates` to a dictionary
    indexed by the order of `genes`.

    :type genes: list
    :param genes: List with the names of the mutations.

    """
    M = len(genes)
    mus = {(i*(0,) + (1,) + (M-i-1)*(0,)):mutation_rates[genes[i]]
           for i in range(M)}
    return mus
