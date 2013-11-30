"""
Construct a 2d hp puzzle that encodes a subset sum problem.

"""
from __future__ import division, print_function, absolute_import

import itertools
import argparse
import sys


def rep(s, n, compress):
    if not compress:
        return s * n
    elif n == 0:
        return ''
    elif n == 1:
        return s
    elif len(s) == 1:
        return s + str(n)
    else:
        return '(' + s + ')' + str(n)


def gen_teeth(sizes, compress):

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

    for i in range(nteeth):

        tooth = []

        # get the tooth dimensions
        ncols = tooth_lengths[i]
        nterms = tooth_terminals[i]
        pfirst = (i == 0)
        plast = (i == nteeth - 1)

        # The top of the first tooth is fully corrugated.
        # The top of the subsequent teeth are flat.
        ncorr = ncols if pfirst else 0
        nflat = ncols - ncorr
        top = rep('rurdr', ncorr, compress) + rep('r', 3*nflat, compress)

        # The bottom of some teeth may have some corrugations.
        ncorr = nterms
        nflat = ncols - nterms
        bot = rep('ldlul', ncorr, compress) + rep('l', 3*nflat, compress)

        # Add the leg of the walk corresponding to the tooth.
        tooth.append(top + 'd' + bot)

        # Add the connection to the next tooth.
        if not plast:
            tooth.append('d')
        
        yield ''.join(tooth)



def main(sizes, compress):
    comb = []
    for k, g in itertools.groupby(gen_teeth(sizes, compress)):
        n = sum(1 for tooth in g)
        comb.append(rep(k, n, compress))
    return ''.join(comb)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--compress', action='store_true',
            help='use a compressed representation of the walk')
    parser.add_argument('sizes', nargs='*')
    args = parser.parse_args()
    sizes = sorted(int(x) for x in args.sizes)
    print(main(sizes, args.compress))

