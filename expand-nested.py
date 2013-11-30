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

def main(s):
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
    return formula.parseString(s)[0]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--in-format',
            choices=('generic', 'hp'),
            default='generic',
            help='sequence input format')
    parser.add_argument('--out-format',
            choices=('same', 'none', 'upper', 'lower', '10'),
            default='same',
            help='sequence output format')
    parser.add_argument('--noheader', action='store_true',
            help='suppress the output headers')
    parser.add_argument('-n', action='store_true',
            help='report sequence length n')
    parser.add_argument('-k', action='store_true',
            help='report sequence weight k')
    args = parser.parse_args()
    s = main(''.join(sys.stdin.read().split()))
    if args.in_format == 'hp':
        if not (set(s.lower()) <= set('hp')):
            raise Exception('expected an encoded hp string')
    if args.n:
        if not args.noheader:
            print 'n',
        print len(s)
    if args.k:
        if args.in_format != 'hp':
            raise Exception('k flag is only allowed for hp input format')
        if not args.noheader:
            print 'k',
        print s.lower.count('h')
    out = None
    if args.out_format == 'same':
        out = s
    elif args.out_format == 'upper':
        out = s.upper()
    elif args.out_format == 'lower':
        out = s.lower()
    elif args.out_format == '10':
        if args.in_format != 'hp':
            raise Exception('the 10 output format is only compatible '
                    'with the hp input format')
        hp10 = {'h':'1', 'p':'0'}
        out = ''.join(hp10[x] for x in s.lower())
    if out is not None:
        if not args.noheader:
            print 'sequence',
        print out
