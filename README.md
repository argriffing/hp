Introduction
------------

This pile of code does things related to some variants of the 2d
[hp](http://en.wikipedia.org/wiki/Hydrophobic-polar_protein_folding_model)
protein folding model.
The variants are described [here](http://arxiv.org/abs/1309.7508).

"Plan A" variant
----------------

This is the traditional 2d hp protein folding puzzle.
Given a specific sequence of hydrophobic and polar residues,
find the 2d lattice conformation that maximizes the number of
adjacent hydrophobic residue pairs.
This problem has been proved to be NP hard.


"Plan B" variant
----------------

In this variant the conformation of the sequence of length n is fixed,
and the puzzle is to find the best positions of k hydrophobic sites
to again maximize the number of adjacent hydrophobic residue pairs.
Residues that are next to each other along the primary sequence
do not count towards this objective.

Some of the code in this pile efficiently solves this variant.

```
$ gcc -o solve-brglez-plan-b solve-brglez-plan-b.c sparsetools.c subsetsum.c \
graphgirth.c connectedcomponents.c breadthfirst.c bondsubgraph.c
$ echo dlluruluurdrdrurddl | ./solve-brglez-plan-b 10
9
10100110100101100101
```

Here's an example of constructing and solving a "plan b" 2d hp puzzle
from the example subset sum problem instance from
[here](http://www.math.sunysb.edu/~scott/blair/Subset_sum_problems_are.html).
The target sum in the example is 5842 which is multiplied by 4
because of technical details in the problem reduction.  The transformed
[problem](https://raw.github.com/argriffing/hp/master/puzzle.out)
has over 70,000 vertices, and by finding a coloring
that has one bond per sequence weight,
it solves the unary subset sum problem whose target sum
is 1/4 the sum of the requested sequence weight.

```
$ python hp-subset-sum.py 267 493 869 961 1000 1153 1246 1598 1766 1922 > puzzle.out
$ ./solve-brglez-plan-b 23368 < puzzle.out | fold -w 90 | head -n 5
23368
000000100100000010010000001001010010100100000000000000000011011011000011000011000000001100
001100001101101100000000000000110110110000110000110000000011000011000011011011000000000000
001101101100001100001100000000110000110000110110110000000000000011011011000011000011000000
001100001100001101101100000000000000110110110000110000110000000011000011000011011011000000
```

The 2d hp subset sum problem can also be constructed using compressed output,
which can be decompressed using a separate script.

```
$ python hp-subset-sum.py 267 493 869 961 1000 1153 1246 1598 1766 1922 \
--compress | fold -w 90 > puzzle.compressed
$ cat puzzle.compressed
(rurdr)10dl30d(r30dl30d)265r30dldlull27d(r27dl27d)225r27dldlull24d(r24dl24d)375r24dldlull2
1d(r21dl21d)91r21dldlull18d(r18dl18d)38r18dldlull15d(r15dl15d)152r15dldlull12d(r12dl12d)92
r12dldlull9d(r9dl9d)351r9dldlull6d(r6dl6d)167r6dldlull3d(r3dl3d)155r3dldlul
$ python expand-nested.py --noheader < puzzle.compressed > puzzle.decompressed
$ fold -w 90 puzzle.decompressed | head -n 5
rurdrrurdrrurdrrurdrrurdrrurdrrurdrrurdrrurdrrurdrdlllllllllllllllllllllllllllllldrrrrrrrr
rrrrrrrrrrrrrrrrrrrrrrdlllllllllllllllllllllllllllllldrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrdlllll
llllllllllllllllllllllllldrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrdlllllllllllllllllllllllllllllldrr
rrrrrrrrrrrrrrrrrrrrrrrrrrrrdlllllllllllllllllllllllllllllldrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
dlllllllllllllllllllllllllllllldrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrdlllllllllllllllllllllllllll
$ diff -s puzzle.out puzzle.decompressed
Files puzzle.out and puzzle.decompressed are identical
```

The decompression script can also be used for hp strings.

```
$ echo "(HP)2PH2PHP2HPH2P2HPH" | python expand-nested.py \
-nk --in-format=hp --out-format=10
n 20
k 10
sequence 10100110100101100101
$ echo "HPHPHHHHHHPHPPHHPPHHPHHHHPPHHHHPPHPHHPPPPPHHH" | \
python expand-nested.py -nk --in-format=hp --out-format=10
n 45
k 27
sequence 101011111101001100110111100111100101100000111
```

"Plan C" variant
--------------

In this variant the goal is still to maximize the number of
adjacent hydrophobic pairs, again not counting pairs that are next to each
other along the primary sequence, but we are allowed to choose
both the conformation of the sequence of length n
and also the positions of the k hydrophobic residues along the sequence.

The pile of code includes an inefficient solver for small n and k.
When n and k are both large (maybe `8 < n`)
and `k << n` (maybe `k < c*sqrt(n)` for some c)
it should be easy to efficiently construct a conformation and sequence
that meets a parity upper bound, but these constructions
are not implemented here.

```
$ gcc -o solve-brglez-plan-c solve-brglez-plan-c.c sparsetools.c subsetsum.c \
graphgirth.c connectedcomponents.c breadthfirst.c bondsubgraph.c plancgrid.c \
compactness.c
$ ./solve-brglez-plan-c -n 23 -k 9
n 23
k 9
score 10
conformation rululdldrdrrdrdldlulur
hp hpphpphpphphphpphpphpph
```

This [visualization](http://bl.ocks.org/argriffing/7619198)
shows optimal "Plan C" solutions for given sequence lengths n
and sequence weights k.
Except for the last row,
each row can be extended arbitrarily towards the right
by just increasing the lengths of polar loops
to construct optimal solutions that meet the parity upper bound.


oeis sequences
--------------

The sequence of counts of rooted self avoiding walks is
[A001411](https://oeis.org/A001411)
and related sequences can be found by searching for
[references](https://oeis.org/search?q=A001411)
to this sequence in the database.
For example the sequence
[A046171](https://oeis.org/A046171)
counts rooted self avoiding walks in a way that removes
much of the redundancy caused by symmetry.


notes on various definitions of compactness
-------------------------------------------

In
Predicting the non-compact conformation of amino acid sequence by
particle swarm optimization
by Yuzhen Guo and Yong Wang (2013)
doi 10.1109/ISB.2013.6623805
compactness is defined as follows.
Suppose that the number of amino acids in a protein sequence is n
and the number of lattice points is m.
If m=n, the conformation of the sequence in the lattice is defined as compact.

In this github repo I define compactness of a connected set of square lattice
points to mean that there are no holes.
