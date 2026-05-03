#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

import os
import sys
import pickle

 
def orderchrs(x):
    return int(''.join([c for c in x if c.isdigit()]))

 
def phasing_chr(bins, cells):
    """
    bins: list of (bin_key, {cell: (A,B)})
    """
    swapped = lambda p: (p[1], p[0])

    norm = [b[1] for b in bins]
    swap = [{e: swapped(norm[i][e]) for e in cells} for i in range(len(norm))]

    def dist(x, y):
        return sum(
            (x[e][0] != y[e][0]) + (x[e][1] != y[e][1])
            for e in cells
        )

    D = []
    B = []

    for i in range(len(norm)):
        if i == 0:
            D.append({'0': 0, '1': 0})
            B.append({'0': -1, '1': -1})
            continue

        d00 = D[i-1]['0'] + dist(norm[i-1], norm[i])
        d10 = D[i-1]['1'] + dist(swap[i-1], norm[i])
        d01 = D[i-1]['0'] + dist(norm[i-1], swap[i])
        d11 = D[i-1]['1'] + dist(swap[i-1], swap[i])

        D.append({
            '0': min(d00, d10),
            '1': min(d01, d11)
        })

        B.append({
            '0': '0' if d00 <= d10 else '1',
            '1': '0' if d01 <= d11 else '1'
        })

    # backtrack
    res = []
    prev = '0' if D[-1]['0'] < D[-1]['1'] else '1'

    for i in reversed(range(len(norm))):
        if prev == '0':
            res.insert(0, norm[i])
        else:
            res.insert(0, swap[i])
        prev = B[i][prev]

    return [ (bins[i][0], res[i]) for i in range(len(bins)) ]


def phase(cns, cells):
    # group by chromosome
    chroms = {}
    for b in cns:
        chroms.setdefault(b[0], []).append(b)

    result = {}

    for c in chroms:
        bins = sorted(chroms[c], key=lambda b: (orderchrs(b[0]), b[1], b[2]))
        data = [(b, cns[b]) for b in bins]

        phased = phasing_chr(data, cells)

        for b, val in phased:
            result[b] = val

    return result

 
def run(input_dir):
    infile = os.path.join(input_dir, "phasing_input.pkl")
    outfile = os.path.join(input_dir, "phasing_output.pkl")

    with open(infile, "rb") as f:
        data = pickle.load(f)

    cns = data["copynumbers"]
    cells = sorted(data["cells"])

    phased = phase(cns, cells)

    with open(outfile, "wb") as f:
        pickle.dump(phased, f)

    print("Done → {}".format(outfile))


if __name__ == "__main__":
    run(sys.argv[1])
