import argparse
import itertools

def gen_bond_indicators():
    edge = 0
    while True:
        for i in range(2):
            for j in range(edge - 1):
                yield 1
            yield 0
        edge += 1

def main(args):
    if args.n is None:
        print 'n=k', 'bonds'
        ntotal = 0
        for i, bond in enumerate(gen_bond_indicators()):
            ntotal += bond
            print i, ntotal
    else:
        print sum(itertools.islice(gen_bond_indicators(), args.n+1))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', type=int)
    main(parser.parse_args())
