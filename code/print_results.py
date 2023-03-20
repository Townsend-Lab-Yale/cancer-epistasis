import numpy as np

from theory import human_order_single_jumps
from locations import default_mutation_names

table_formats ={
    'LaTeX':{
        'sep':" & ",
        'start_line':"",
        'end_line':"\\\\",
        'top':"\\toprule",
        'mid':"\\midrule",
        'bottom':"\\bottomrule"},
    'org-mode':{
        'sep':" | ",
        'start_line':"|",
        'end_line':"|",
        'top':"",
        'mid':"|-",
        'bottom':""}}


def print_table(results, mutation_names=None, table_format=None,
                tol=10**(-6), gammas_scale=10**6,
                include_epistatic_effect=False):
    """Print a table with the results to use in the paper.


    :type results: dict
    :param results: Dictionary with all the results (as returned by
        `produce_results` from the main module). It should have at
        least fluxes estimates (key: 'lambdas'), their confidence
        intervals (key: 'lambdas_cis'), the scale selection
        coefficients (key: 'gammas'), their confidence intervals (key:
        'gammas_cis'), and the mutation rates (key: 'mus')

    :type tol: flot
    :param tol: Assume that fluxes and selection coefficients smaller
        than this number are zero on the table.

    :rtype: dict
    :return: Dictionary with the selection coefficients, indexed by a
        pair of tuples representing the mutation combination where the
        flux is coming from and going to.

    """

    if table_format is None:
       table_format = 'org-mode'

    chars = {key:value
             for key, value in table_formats[table_format].items()}

    M = len(results["mus"])

    if mutation_names is None:
        mutation_names = default_mutation_names[:M]
        mutation_names_sep = ''
    else:
        mutation_names_sep = '+'

    rounded = {}

    ## Round and scale values and set small ones to zero
    rounded["lambdas"] = {
        key:np.round(value, 2) if value >= tol else 0
        for key, value in results["lambdas"].items()}
    rounded["lambdas_cis"] = {
        key:[np.round(value[0], 2), np.round(value[1], 2)]
        if value[0] >= tol else [0, np.round(value[1], 2)]
        for key, value in results["lambdas_cis"].items()}

    rounded["gammas"] = {
        key:np.round(value/gammas_scale, 2)
        if value/gammas_scale >= tol else 0
        for key, value in results["gammas"].items()}
    rounded["gammas_cis"] = {
        key:[np.round(value[0]/gammas_scale, 2),
             np.round(value[1]/gammas_scale, 2)]
        if value[0]/gammas_scale >= tol else
        [0, np.round(value[1]/gammas_scale, 2)]
        for key, value in results["gammas_cis"].items()}


    order_lambdas = human_order_single_jumps(M)

    # print(chars['top'])
    # print(chars['start_line']
    #       + chars['sep'].join(
    #           ["{Genotype}",
    #            "{Mutation}",
    #            "{Flux}",
    #            "{Flux CI\\tnote{$a$}}",
    #            "{Mutation rate $\\times 10^{-6}$}",
    #            "{Selection coefficient $\\times 10^{3}$}"]
    #           + (["{Epistatic effect}"] if include_epistatic_effect
    #              else []))
    #       + chars['end_line'])
    print(chars['mid'])
    for xy in order_lambdas:
        acting_in = tuple(np.array(xy[1])-np.array(xy[0]))
        mu = results["mus"][acting_in]
        if xy[0] == M*(0,):
            epistatic_effect = ""
        elif results["lambdas"][xy] > results["lambdas"][(M*(0,), acting_in)]:
            epistatic_effect = "synergistic"
        else:
            epistatic_effect = "antagonistic"
        rows = ([mutation_names_sep.join([mutation_names[i]
                                          for i in range(M) if xy[0][i] == 1]),
                 mutation_names[acting_in.index(1)],
                 ## The {\,\,\,} below fixes a siunitx bug in spacing in LaTeX
                 ("{\,\,\,}0" if rounded["lambdas"][xy] == 0
                  else "{:.2f}{{)}}".format(rounded["lambdas"][xy])),
                 ("{\,\,\,(} 0{,}" if rounded["lambdas_cis"][xy][0] == 0
                  else "{{(}} {:.2f}{{,}}".format(rounded["lambdas_cis"][xy][0])),
                 "{:.2f}{{)}}".format(rounded["lambdas_cis"][xy][1]),
                 "{:.2f}{{e-8}}".format(round(mu*10**8, 2)),
                 ("{\,\,\,}0" if rounded["gammas"][xy] == 0
                  else "{:.2f}".format(rounded["gammas"][xy])),
                 ("{\,\,\,(} 0{,}" if rounded["gammas_cis"][xy][0] == 0
                  else "{{(}} {:.2f}{{,}}".format(rounded["gammas_cis"][xy][0])),
                 "{:.2f}{{)}}".format(rounded["gammas_cis"][xy][1])]

                + ([epistatic_effect] if include_epistatic_effect
                    else []))
        rows[0] = ("normal" if rows[0] == "" else rows[0])

        print(chars['start_line'] + chars['sep'].join(rows) + chars['end_line'])
    print(chars['bottom'])
