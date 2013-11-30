"""
Expand an encoded HP string.

Example encoded strings:
(HP)2PH2PHP2HPH2P2HPH
HPHPHHHHHHPHPPHHPPHHPHHHHPPHHHHPPHPHHPPPPPHHH
"""

import argparse
import sys

from pyparsing import (Forward, Word, Group, OneOrMore, Literal, Optional,
        nums, alphas)

def process_integer(t):
    return int(t[0])

def process_term(tokens):
    t = tokens[0]
    if t.subgroup:
        return t.subgroup[0] * t.mult
    else:
        return t[0] * t.mult

def process_formula(tokens):
    return ''.join(tokens)

def main(args):
    s = ''.join(sys.stdin.read().split())
    lpar = Literal('(').suppress()
    rpar = Literal(')').suppress()
    integer = Word(nums)
    element = Word(alphas, exact=1)
    formula = Forward()
    term = Group((element | Group(lpar + formula + rpar)('subgroup')) +
            Optional(integer, default=1)('mult'))
    formula << OneOrMore(term)
    integer.setParseAction(process_integer)
    term.setParseAction(process_term)
    formula.setParseAction(process_formula)
    hp = formula.parseString(s)[0].lower()
    if not (set(hp) <= set('hp')):
        raise Exception('expected an encoded HP string')
    return hp

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out-format', choices=('HP', 'hp', '10'),
            help='sequence output format')
    parser.add_argument('-n', action='store_true',
            help='report sequence length n')
    parser.add_argument('-k', action='store_true',
            help='report sequence weight k')
    args = parser.parse_args()
    hp = main(args)
    if args.n:
        print 'n', len(hp)
    if args.k:
        print 'k', hp.count('h')
    if args.out_format == 'HP':
        print hp.upper()
    elif args.out_format == 'hp':
        print hp
    elif args.out_format == '10':
        hp10 = {'h':'1', 'p':'0'}
        print 'sequence', ''.join(hp10[x] for x in hp)

