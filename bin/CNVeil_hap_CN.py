#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

import os
import sys
import csv
import pickle

def convert_csv_to_dict(csv_path):
    if not os.path.exists(csv_path):
        raise IOError("Input file not found: {}".format(csv_path))

    print("[INFO] Reading {}".format(csv_path))

    with open(csv_path, "r") as f:
        reader = csv.reader(f)
        header = next(reader)
        cells = header[3:]

        copynumbers = {}

        for row in reader:
            if len(row) < 3:
                continue

            chrom = row[0]
            start = int(row[1])
            end = int(row[2])
            key = (chrom, start, end)

            copynumbers[key] = {}

            for i, cell in enumerate(cells):
                val = row[i + 3] if i + 3 < len(row) else ""

                if "|" in val:
                    try:
                        a, b = map(int, val.split("|"))
                    except:
                        a, b = (0, 0)
                else:
                    try:
                        a = int(val)
                        b = 0
                    except:
                        a, b = (0, 0)

                copynumbers[key][cell] = (a, b)

    return copynumbers, cells


def orderchrs(x):
    return int(''.join([c for c in x if c.isdigit()]))


def phasing_chr(bins, cells):
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

        D.append({'0': min(d00, d10), '1': min(d01, d11)})
        B.append({
            '0': '0' if d00 <= d10 else '1',
            '1': '0' if d01 <= d11 else '1'
        })

    # backtrack
    res = []
    prev = '0' if D[-1]['0'] < D[-1]['1'] else '1'

    for i in reversed(range(len(norm))):
        res.insert(0, norm[i] if prev == '0' else swap[i])
        prev = B[i][prev]

    return [(bins[i][0], res[i]) for i in range(len(bins))]


def phase(cns, cells):
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


def run(csv_path, output_pkl):
    copynumbers, cells = convert_csv_to_dict(csv_path)

    print("[INFO] {} bins × {} cells".format(len(copynumbers), len(cells)))

    phased = phase(copynumbers, cells)

    with open(output_pkl, "wb") as f:
        pickle.dump(phased, f)

    print("[INFO] Saved phased result to {}".format(output_pkl))


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python2.7 CNVeil_hap_CN.py <ascn_seg.csv> <output.pkl>")
        sys.exit(1)

    csv_path = sys.argv[1]
    output_pkl = sys.argv[2]

    run(csv_path, output_pkl)
