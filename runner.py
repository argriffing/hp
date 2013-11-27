"""
Run the 'plan c' solver for various combinations of n and k.

For each k > 7, run the search for n >= k until the score reaches k + 1.

"""

from StringIO import StringIO
import itertools
import subprocess

def main():
    filename = 'nksolutions-compact.out'
    with open(filename, 'a') as fout:
        for k in itertools.count(8):
            for n in itertools.count(k):

                # Run the solver and extract the output.
                args = ('./solve-brglez-plan-c', '-n', str(n), '-k', str(k))
                result = subprocess.check_output(args)

                # Append the output to the solutions file.
                fout.write(result + '\n')
                fout.flush()

                # If the score has attained k + 1 then proceed to the next k.
                d = dict(x.split() for x in result.splitlines())
                score = int(d['score'])
                if score > k:
                    break


if __name__ == '__main__':
    main()
