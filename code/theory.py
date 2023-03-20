import numpy as np
from itertools import product


# Numbers of possible positive lambdas (I doubt M would ever be >20)
numbers_positive_lambdas = [2**(n-1)*n  ## because this is
                                        ## np.sum([binom(n, i)*i
                                        ##         for i in range(n+1)]).astype(int)
                                        ## (https://mathworld.wolfram.com/BinomialSums.html)
                            for n in range(0, 20+1)]
"""List with number of possible positive lambdas per mutation number
M.

"""



def build_S_with_tuples(M):
    """Build entire space of mutations combinations S.

    It will represented by a list of size 2^M with items being tuples
    of size M with 0s and 1s.

    :type M: int
    :param M: Number of mutations.

    :rtype: list
    :return: S as an ordered list of tuples.

    """
    return [x for x in product([0, 1], repeat=M)]


def build_S_as_array(M):
    """Build entire space of mutations combinations S.

    It will represented by a numpy array of shape (2^M, M).

    :type M: int
    :param M: Number of mutations.

    :rtype: numpy.ndarray
    :return: S as a numpy array of 0s and 1s.

    """
    return np.array(build_S_with_tuples(M))


def obtain_pos_lambdas_indices(S):
    """Obtain the lambda indices that could be positive.

    If Q is a matrix with entries its i-th row, j-th column being
    lambda[S[j], S[i]], then this function would produce a matrix with
    i-th row, j-th column being True if it is posible for lambda[S[j],
    S[i]] > 0, just by the definition of the lambdas, in other words
    if S[i] is equal S[j] plus one more mutation.

    :type: numpy.ndarray
    :param: S as a numpy array of 0s and 1s.

    :rtype: numpy.ndarray
    :return: A numpy array of booleans.

    """

    M = np.shape(S)[1]

    pos_lambdas_indices = np.array(
        [[True
          if np.sum(S[j] >= S[i]) == M
          and np.sum(S[j] - S[i]) == 1
          and i < len(S)-1 # for last x=(1,...,1) all are 0
          else False
          for i in range(len(S))]
         for j in range(len(S))])

    return pos_lambdas_indices


def order_pos_lambdas(S):
    """Order the tuples in S that subscript the positive lambdas.

    This function can be used to translate a flatten array of positive
    lambdas (such as the one returned by the MCMC or MAP finding
    routine) to a dictionary indexed by two elements of S.
    If pos_lambdas_ordered represents the ordered list returned by
    this function, lambda[x, y] is the flux from x to y, and Q is a
    matrix with i-th row j-th column being lambda[S[j], S[i]] then:

        Q[pos_lambdas_ordered]

    is the same as

        [lambdas[build_S_with_tuples(M).index(j),
                 build_S_with_tuples(M).index(i)]
         for i, j in positive_lambdas_ordered]


    :type: numpy.ndarray
    :param: S as a numpy array of 0s and 1s.

    :rtype: list
    :return: A list with the lambda tuples numpy array of booleans.

    """

    pos_lambdas_indices = obtain_pos_lambdas_indices(S)
    pos_lambdas_ordered = [(tuple(S[j]), tuple(S[i]))
                           for i in range(len(S))
                           for j in range(len(S))
                           if pos_lambdas_indices[i, j]]
    return pos_lambdas_ordered



def human_order_single_jumps(n, first_by_destiny=False):
    if not first_by_destiny:
        pos_lambdas = order_pos_lambdas(build_S_as_array(n))
        ordered = sorted(
            pos_lambdas,
            key=lambda xy: (np.sum(xy[0])*10**(2*n) +
                            np.dot(xy[0], [10**i for i in range(n, 2*n)]) +
                            np.dot(xy[1], [10**i for i in range(n)])))
    else:
        ordered = []
        for i in range(n):
            ordered = (
                ordered
                + [xy for xy in human_order_single_jumps(n, False)
                   if tuple(np.array(xy[1])-np.array(xy[0])).index(1) == i])
    return ordered
