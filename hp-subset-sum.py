"""
Construct a 2d hp puzzle that encodes a subset sum problem.

"""
from __future__ import division, print_function, absolute_import

import sys


def main(sizes):

    # define the comb dimensions
    nsizes = len(sizes)
    nteeth = max(sizes)
    tooth_lengths = [sum(1 for x in sizes if x > i) for i in range(nteeth)]
    tooth_terminals = [sizes.count(i+1) for i in range(nteeth)]

    # check some invariants
    assert nsizes > 0
    assert nteeth > 0
    assert tooth_lengths[0] == nsizes
    assert tooth_terminals[-1] > 0

    # begin constructing the legs of the walk
    walk = []

    # add each tooth
    for i in range(nteeth):

        # get the tooth dimensions
        ncols = tooth_lengths[i]
        nterms = tooth_terminals[i]
        pfirst = (i == 0)
        plast = (i == nteeth - 1)

        # The top of the first tooth is fully corrugated.
        # The top of the subsequent teeth are flat.
        ncorr = ncols if pfirst else 0
        nflat = ncols - ncorr
        top = 'rurdr' * ncorr + 'rrr' * nflat

        # The bottom of some teeth may have some corrugations.
        ncorr = nterms
        nflat = ncols - nterms
        bot = 'ldlul' * ncorr + 'lll' * nflat

        # Add the leg of the walk corresponding to the tooth.
        walk.append(top + 'd' + bot)

        # Add the connection to the next tooth.
        if not plast:
            walk.append('d')

    return ''.join(walk)


if __name__ == '__main__':
    sizes = sorted(int(x) for x in sys.argv[1:])
    walk_string = main(sizes)
    print(walk_string)
