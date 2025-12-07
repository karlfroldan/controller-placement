import numpy as np

def normalize_probabilities(ps, eps: float = 1e-9):
    if not isinstance(ps, np.ndarray):
        ps = np.array(ps)
    # Set small values to absolute 0
    ps[ps < eps] = 0.0

    # Renormalize
    total = np.sum(ps)
    if total > 0:
        ps = ps / total
    else:
        raise ValueError('Probabilities sum to zero')
    return ps.tolist()

def set_range(lst):
    return set(range(len(lst)))

def values_to_list(vs):
    """
      Convert an AMPL object.values() to a list containing those values.
    """
    return [x for _, x in vs.to_list()]

def to_set_ids(lst):
    """Convert the list into a dictionary of ID to item pair"""
    return { idx: x for idx, x in enumerate(lst) }

def make_ones(m, n):
    """Create a list of size M with N ones"""
    assert(m >= n)
    ones = np.ones(n, dtype='int')
    zeros = np.zeros(m - n, dtype='int')

    rng = np.random.default_rng()
    arr = np.concatenate((ones, zeros))
    rng.shuffle(arr)
    return arr

def one_indices(xs):
    # add1 = lambda x: x + 1

    if isinstance(xs, list):
        xs = np.array(xs)
    binary = (xs >= 0.5).astype(int)
    return np.where(binary == 1)[0].tolist()
