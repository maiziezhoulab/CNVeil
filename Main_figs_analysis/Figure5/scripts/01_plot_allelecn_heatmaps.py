#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

import os, sys, glob, re, argparse, warnings
from collections import defaultdict
from itertools import cycle

import numpy as np
import scipy.cluster.hierarchy as hier
import pandas as pd

import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
 
def log(msg, level='INFO'):
    sys.stderr.write("[{}] {}\n".format(level, msg))
    sys.stderr.flush()

def save_discrete_colorbar(cmap, labels, out_path, orientation='horizontal', dpi=600):
    """
    Save a standalone discrete colorbar with readable ticks.

    labels: list of strings, length = number of discrete colors.
    """
    n = len(labels)
    if n <= 0:
        return

    boundaries = np.arange(n + 1) - 0.5
    norm = mpl.colors.BoundaryNorm(boundaries, n)
    sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    if orientation == 'horizontal':
        # Wider + taller so rotated labels are readable
        fig_w = max(6.0, 0.70 * n)
        fig_h = 2.6
        fig = plt.figure(figsize=(fig_w, fig_h))
        # left, bottom, width, height
        cax = fig.add_axes([0.06, 0.62, 0.92, 0.22])
        cb = fig.colorbar(sm, cax=cax, orientation='horizontal', ticks=np.arange(n))

        cb.ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=11)
        cb.ax.tick_params(axis='x', length=6, width=1.2, pad=8)
        cb.outline.set_linewidth(1.5)

    else:
        fig_w = 2.4
        fig_h = max(6.0, 0.45 * n)
        fig = plt.figure(figsize=(fig_w, fig_h))
        cax = fig.add_axes([0.35, 0.06, 0.25, 0.90])
        cb = fig.colorbar(sm, cax=cax, orientation='vertical', ticks=np.arange(n))

        cb.ax.set_yticklabels(labels, fontsize=11)
        cb.ax.tick_params(axis='y', length=6, width=1.2, pad=6)
        cb.outline.set_linewidth(1.5)

    fig.savefig(out_path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
 
orderchrs = (lambda x: int(''.join([l for l in x if l.isdigit()])))
order_bin = (lambda b: (orderchrs(b[0]), int(b[1]), int(b[2])))
 
# Barcode mapping
# collected_barcodes.txt: 2 cols -> SRR barcode 
_SRR_RE = re.compile(r'^SRR\d+$')
_SUFFIX_RE = re.compile(r'^(.*)-\d+$')

def read_barcode_mapping(barcode_file):
    """
    Returns:
      barcode_to_srr: dict barcode -> SRR
      srr_to_barcode: dict SRR -> barcode
    """
    if not os.path.isfile(barcode_file):
        raise ValueError("Barcode file not found: {}".format(barcode_file))

    barcode_to_srr = {}
    srr_to_barcode = {}
    with open(barcode_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            p = line.split()
            if len(p) < 2:
                continue
            srr, bc = p[0], p[1]
            barcode_to_srr[bc] = srr
            srr_to_barcode[srr] = bc
    if not barcode_to_srr:
        raise ValueError("Empty/invalid barcode mapping: {}".format(barcode_file))
    return barcode_to_srr, srr_to_barcode

def cell_to_srr(cell_id, barcode_to_srr):
    """
    Convert a CELL string to SRR if possible.
    Handles:
      - already SRR#### -> keep
      - CNRein style: BARCODE-0 / BARCODE-12 -> strip suffix then map
      - raw BARCODE -> map
    If not found -> return original string (so you can see what didn't map).
    """
    if cell_id is None:
        return cell_id
    cell_id = str(cell_id)

    if _SRR_RE.match(cell_id):
        return cell_id

    # strip trailing -digits (CNRein, sometimes others)
    m = _SUFFIX_RE.match(cell_id)
    if m:
        base = m.group(1)
        if base in barcode_to_srr:
            return barcode_to_srr[base]
        # sometimes the full id is in mapping (rare)
        if cell_id in barcode_to_srr:
            return barcode_to_srr[cell_id]
        return cell_id

    # raw barcode
    if cell_id in barcode_to_srr:
        return barcode_to_srr[cell_id]

    return cell_id

# -----------------------------
# Read your tool TSVs (not the old 12/13-field combined format)
# Expected columns (like your head screenshot):
# #CHR START END CELL NORM_COUNT COUNT RDR A_COUNT B_COUNT BAF CLUSTER HAP_CN CORRECTED_HAP_CN
# -----------------------------
def read_tool_tsv_as_bins(tsv_path):
    """
    Read tool TSV with header starting '#CHR' and convert to bins dict.
    """
    # IMPORTANT: do NOT use comment='#' because header begins with '#CHR'
    df = pd.read_csv(tsv_path, sep='\t', header=0)

    # Normalize header name
    if "#CHR" in df.columns:
        df = df.rename(columns={"#CHR": "CHR"})

    # Some tools may have stray comment rows after header; drop those if present
    if "CHR" in df.columns:
        df = df[~df["CHR"].astype(str).str.startswith("#")]

    required = ["CHR", "START", "END", "CELL", "RDR", "BAF", "CLUSTER", "HAP_CN"]
    for c in required:
        if c not in df.columns:
            raise ValueError("Missing required column '{}' in {}; columns={}".format(c, tsv_path, list(df.columns)))

    iscorr = ("CORRECTED_HAP_CN" in df.columns)

    def parse_cn(s):
        if pd.isnull(s):
            return (0, 0)
        s = str(s)
        if '|' not in s:
            return (0, 0)
        a, b = s.split('|', 1)
        try:
            return (int(float(a)), int(float(b)))
        except:
            return (0, 0)

    bins = defaultdict(lambda: dict())
    cells = set()

    for _, r in df.iterrows():
        binkey = (str(r["CHR"]), int(r["START"]), int(r["END"]))
        cell = str(r["CELL"])
        rdr = float(r["RDR"]) if not pd.isnull(r["RDR"]) else 0.0
        baf = float(r["BAF"]) if not pd.isnull(r["BAF"]) else 0.0
        clu = str(r["CLUSTER"])
        cns = parse_cn(r["HAP_CN"])

        rec = {"RDR": rdr, "BAF": baf, "Cluster": clu, "CNS": cns}
        if iscorr:
            rec["CORR-CNS"] = parse_cn(r["CORRECTED_HAP_CN"])

        bins[binkey][cell] = rec
        cells.add(cell)

    pos = sorted(bins.keys(), key=order_bin)

    for x, b in enumerate(pos):
        for e in bins[b]:
            bins[b][e]["Genome"] = x

    return bins, pos, sorted(cells), iscorr

# -----------------------------
# CNVeil order
# -----------------------------
def clustering_tot(bins, pos, cells):
    # Use only bins/cells that exist; if missing bins for a cell, fill with (0,0)
    data = []
    for e in cells:
        row = []
        for b in pos:
            if e in bins[b]:
                a,bcn = bins[b][e]['CNS']
            else:
                a,bcn = (0,0)
            row.extend([a,bcn])
        data.append(row)

    linkage = hier.linkage(data, method='average', metric='hamming', optimal_ordering=True)
    clus = hier.fcluster(linkage, t=len(cells), criterion='maxclust')
    mapc = [(clus[i], e) for i, e in enumerate(cells)]
    index = [pp[1] for pp in sorted(mapc, key=(lambda x: x[0]))]
    clones = {e: int(clus[i]) for i, e in enumerate(cells)}
    return index, clones

def make_index_from_cnveil(cnveil_order, tool_cells):
    tool_set = set(tool_cells)
    kept = [c for c in cnveil_order if c in tool_set]   # CNVeil-only missing in tool -> skip
    extra = [c for c in tool_cells if c not in set(cnveil_order)]  # tool-only -> append
    return kept + extra
 
def set_style(gridsize):
    plt.style.use('ggplot')
    sns.set_style("whitegrid")
    plt.rcParams["axes.grid"] = True
    plt.rcParams["axes.edgecolor"] = "k"
    plt.rcParams["axes.linewidth"] = 1.5

def addchr(g, pos):
    corners = []
    prev = 0
    for x, b in enumerate(pos):
        if x != 0 and pos[x-1][0] != pos[x][0]:
            corners.append((prev, x))
            prev = x
    corners.append((prev, x))

    ax = g.ax_heatmap
    ticks = []
    for o in corners:
        ax.set_xticks(np.append(ax.get_xticks(), int(float(o[1] + o[0] + 1) / 2.0)))
        ticks.append(pos[o[0]][0])
    ax.set_xticklabels(ticks, rotation=45, ha='center')
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)

def draw(table, bins, pos, palette, center, title, out_path, gridsize):
    chr_palette = cycle(['#525252', '#969696', '#cccccc'])
    chr_colors = {c: next(chr_palette) for c in sorted(set(b[0] for b in bins), key=orderchrs)}

    para = {}
    para['data'] = table
    para['cmap'] = palette
    if center is not None:
        para['center'] = center
    para['yticklabels'] = False
    para['row_cluster'] = False
    para['xticklabels'] = False
    para['col_cluster'] = False
    para['figsize'] = gridsize
    para['rasterized'] = True
 
    para['col_colors'] = pd.DataFrame(
        [{'index': s, 'chromosomes': chr_colors[pos[x][0]]} for x, s in enumerate(table.columns)]
    ).set_index('index')

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        g = sns.clustermap(**para)

        # --- remove the heatmap colorbar (old seaborn compatible) ---
        # seaborn clustermap usually creates g.cax for the colorbar axis
        try:
            if hasattr(g, 'cax') and g.cax is not None:
                g.cax.set_visible(False)
                # also remove it from the figure so it doesn't take space
                try:
                    g.fig.delaxes(g.cax)
                except:
                    pass
        except:
            pass

        addchr(g, pos)
        g.fig.suptitle(title)

    plt.savefig(out_path, bbox_inches='tight', dpi=600)
    plt.close()

def plot_allelecn_states(bins, pos, index, out_path, gridsize, val='CNS', title='Copy-number states'):
    avail = [(t - i, i) for t in xrange(7) for i in reversed(xrange(t + 1)) if i <= t - i]
    _ord = (lambda p: (max(p), min(p)))
    convert = (lambda p: _ord(p) if sum(p) <= 6 else min(avail, key=(lambda x: abs(p[0] - x[0]) + abs(p[1] - x[1]))))

    df = []
    for x, e in enumerate(index):
        for b in pos:
            if e not in bins[b]:
                continue
            df.append({'Cell': x, 'Genome': bins[b][e]['Genome'], 'Value': convert(bins[b][e][val])})

    if not df:
        # nothing to plot
        log("No data to plot for {} (val={})".format(out_path, val), level="WARN")
        return

    df = pd.DataFrame(df)

    found = [v for v in avail if v in set(df['Value'])]
    if not found:
        log("No CN states found for {} (val={})".format(out_path, val), level="WARN")
        return

    smap = {v: i for i, v in enumerate(found)}
    df['CN states'] = df['Value'].map(smap)

    table = pd.pivot_table(df, values='CN states', columns=['Genome'], index=['Cell'], aggfunc='first')

    palette = {}
    palette.update({(0, 0): 'darkblue'})
    palette.update({(1, 0): 'lightblue'})
    palette.update({(1, 1): 'lightgray', (2, 0): 'dimgray'})
    palette.update({(2, 1): 'lightgoldenrodyellow', (3, 0): 'gold'})
    palette.update({(2, 2): 'navajowhite', (3, 1): 'orange', (4, 0): 'darkorange'})
    palette.update({(3, 2): 'salmon', (4, 1): 'red', (5, 0): 'darkred'})
    palette.update({(3, 3): 'plum', (4, 2): 'orchid', (5, 1): 'purple', (6, 0): 'indigo'})

    colors = [palette[c] for c in found]
    cmap = LinearSegmentedColormap.from_list('multi-level', colors, len(colors))

    # standalone colorbar
    cb_labels = ["{}|{}".format(v[0], v[1]) for v in found]
    root, ext = os.path.splitext(out_path)
    cb_path = root + ".colorbar" + ext
    save_discrete_colorbar(cmap, cb_labels, cb_path, orientation='horizontal')

    # heatmap (no colorbar, no clone label) — relies on updated draw()
    draw(table, bins, pos, palette=cmap, center=None, title=title, out_path=out_path, gridsize=gridsize)
 
def rename_tsv_cells_to_srr(in_tsv, out_tsv, barcode_to_srr):
    df = pd.read_csv(in_tsv, sep='\t', comment=None, header=0)
    if "CELL" not in df.columns:
        raise ValueError("No CELL column in {}".format(in_tsv))
    df["CELL"] = df["CELL"].apply(lambda x: cell_to_srr(x, barcode_to_srr))
    df.to_csv(out_tsv, sep='\t', index=False)
 
def tool_name_from_path(path):
    base = os.path.basename(path)
    # "CNRein.tsv" -> "CNRein"
    if base.endswith(".tsv"):
        base = base[:-4]
    return base

def parse_args():
    ap = argparse.ArgumentParser(description="Rename CELL to SRR using collected_barcodes.txt, then plot allelecn using CNVeil order.")
    ap.add_argument("--folder", required=True, help="Folder containing tool .tsv files (TN1 etc).")
    ap.add_argument("--barcode_file", required=True, help="collected_barcodes.txt (SRR barcode).")
    ap.add_argument("--cnveil_tsv", required=True, help="CNVeil.tsv inside that folder (or full path).")
    ap.add_argument("--pattern", default="*.tsv", help="Which TSVs to process (default: *.tsv).")
    ap.add_argument("--format", default="png", choices=["png","pdf"], help="Output format.")
    ap.add_argument("--gridsize", default="12,6", help="Heatmap figsize, e.g. 12,6")
    args = ap.parse_args()

    def get_size(s):
        p = s.split(',')
        if len(p) != 2:
            raise ValueError("Bad gridsize: {}".format(s))
        return (float(p[0]), float(p[1]))

    return {
        "folder": os.path.abspath(args.folder),
        "barcode_file": os.path.abspath(args.barcode_file),
        "cnveil_tsv": os.path.abspath(args.cnveil_tsv),
        "pattern": args.pattern,
        "format": args.format,
        "gridsize": get_size(args.gridsize),
    }

def main():
    args = parse_args()
    set_style(args["gridsize"])

    folder = args["folder"]
    barcode_to_srr, _ = read_barcode_mapping(args["barcode_file"])
    log("Loaded mapping: {} barcodes".format(len(barcode_to_srr)))

    # 1) Rename all TSVs to SRR (write *.srr.tsv)
    files = sorted(glob.glob(os.path.join(folder, args["pattern"])))
    if not files:
        raise ValueError("No files matched: {}".format(os.path.join(folder, args["pattern"])))

    renamed = []
    for f in files:
        out = f[:-4] + ".srr.tsv" if f.endswith(".tsv") else f + ".srr.tsv"
        rename_tsv_cells_to_srr(f, out, barcode_to_srr)
        renamed.append(out)
    log("Renamed CELL -> SRR for {} TSVs".format(len(renamed)))

    # 2) Build CNVeil order from CNVeil.srr.tsv
    cnveil_srr = args["cnveil_tsv"]
    if cnveil_srr.endswith(".tsv"):
        cnveil_srr = cnveil_srr[:-4] + ".srr.tsv"
    if not os.path.isfile(cnveil_srr):
        raise ValueError("Expected CNVeil.srr.tsv not found: {}".format(cnveil_srr))

    log("Reading CNVeil SRR TSV: {}".format(cnveil_srr))
    cn_bins, cn_pos, cn_cells, cn_iscorr = read_tool_tsv_as_bins(cnveil_srr)
    cn_order, _ = clustering_tot(cn_bins, cn_pos, cn_cells)

    order_path = os.path.join(folder, "CNVeil.cell_order.srr.txt")
    with open(order_path, "w") as o:
        for c in cn_order:
            o.write(c + "\n")
    log("Saved CNVeil SRR order: {}".format(order_path))

    # Plot CNVeil itself
    out_cn = os.path.join(folder, "allelecn.CNVeil.{}".format(args["format"]))
    plot_allelecn_states(cn_bins, cn_pos, cn_order, out_cn, args["gridsize"], val="CNS", title="CNVeil Copy-number states")
    log("Wrote {}".format(out_cn))

    if cn_iscorr:
        out_cnc = os.path.join(folder, "allelecn-corrected.CNVeil.{}".format(args["format"]))
        plot_allelecn_states(cn_bins, cn_pos, cn_order, out_cnc, args["gridsize"], val="CORR-CNS", title="CNVeil Corrected copy-number states")
        log("Wrote {}".format(out_cnc))

    # 3) For each tool, reorder using CNVeil order and plot
    for f in renamed:
        if os.path.abspath(f) == os.path.abspath(cnveil_srr):
            continue
        tname = tool_name_from_path(f.replace(".srr.tsv", ".tsv"))
        log("Processing {}".format(tname))

        bins, pos, cells, iscorr = read_tool_tsv_as_bins(f)
        idx = make_index_from_cnveil(cn_order, cells)

        out1 = os.path.join(folder, "allelecn.{}.{}".format(tname, args["format"]))
        plot_allelecn_states(bins, pos, idx, out1, args["gridsize"], val="CNS", title="{} Copy-number states".format(tname))
        log("Wrote {}".format(out1))

        if iscorr:
            out2 = os.path.join(folder, "allelecn-corrected.{}.{}".format(tname, args["format"]))
            plot_allelecn_states(bins, pos, idx, out2, args["gridsize"], val="CORR-CNS", title="{} Corrected copy-number states".format(tname))
            log("Wrote {}".format(out2))

    log("Done.")

if __name__ == "__main__":
    main()
