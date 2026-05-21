# %%
import os
import sys
import pickle
import matplotlib.pyplot as plt
from collections import Counter
from tqdm import tqdm
from joblib import Parallel, delayed
from scipy.stats import shapiro, ttest_1samp
from sklearn.mixture import GaussianMixture
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import leaves_list
from sklearn.cluster import AgglomerativeClustering
import logging

# Set up logger
logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
    filename="run.log",
    filemode="w"
)

logger = logging.getLogger(__name__)


def _to_linkage(model):
    children, dist = model.children_, model.distances_
    n = model.n_leaves_
    counts = np.zeros(len(children), dtype=float)

    for i, (l, r) in enumerate(children):
        cl = 1 if l < n else counts[l - n]
        cr = 1 if r < n else counts[r - n]
        counts[i] = cl + cr

    return np.column_stack([children, dist, counts]).astype(float)


def hclust_categorical_sorted(X, linkage="average"):
    A = pd.DataFrame(X).astype(str)
    Z = np.column_stack([pd.factorize(A[c])[0] for c in A.columns])

    D = squareform(pdist(Z, metric="hamming"))

    ac = AgglomerativeClustering(
        n_clusters=None,
        distance_threshold=0.0,
        linkage=linkage,
        compute_distances=True,
    ).fit(D)

    L = _to_linkage(ac)
    order = leaves_list(L)

    return np.asarray(X, dtype=object)[order], order


def is_single_normal_mu_in_range(x, lo=0.4, hi=0.6):
    x = np.asarray(x, float)
    x = x[~np.isnan(x)]

    if x.size < 8:
        return False

    X = x.reshape(-1, 1)

    g1 = GaussianMixture(1, random_state=0).fit(X)
    g2 = GaussianMixture(2, random_state=0).fit(X)

    # print(float(g1.means_[0, 0]), float(np.sqrt(g1.covariances_[0, 0, 0])))

    return (g1.bic(X) <= g2.bic(X)) and (lo <= float(g1.means_[0, 0]) <= hi)


def is_normal_centered_at_half(x, alpha=0.05):
    x = np.asarray(x, float)
    x = x[~np.isnan(x)]

    return (
        (len(x) >= 3)
        and (shapiro(x).pvalue < alpha)
        and (ttest_1samp(x, 0.5).pvalue < alpha)
    )


def M_step(df, df1):
    df2 = df.merge(df1[["cell", "log_theta_A", "log_theta_B"]])

    df2["probA"] = df2["ref"] * df2["log_theta_A"] + df2["alt"] * df2["log_theta_B"]
    df2["probB"] = df2["ref"] * df2["log_theta_B"] + df2["alt"] * df2["log_theta_A"]

    df3 = df2.groupby(["chrom", "pos"])[["probA", "probB"]].sum().reset_index()

    df3["ref_is_A"] = np.exp(df3["probA"]) / (
        np.exp(df3["probA"]) + np.exp(df3["probB"])
    )

    df3["ref_is_A"] = df3["ref_is_A"].fillna(0)
    df3["F"] = df3["probA"] * df3["ref_is_A"] + df3["probB"] * (1 - df3["ref_is_A"])

    return df3


def E_step(df, snp_flip):
    df2 = df.merge(snp_flip)

    df2["alt_is_A"] = 1 - df2["ref_is_A"]
    df2["A_read"] = df2["ref"] * df2["ref_is_A"] + df2["alt"] * df2["alt_is_A"]
    df2["cov"] = df2["ref"] + df2["alt"]

    df3 = df2.groupby("cell")[["A_read", "cov"]].sum().reset_index()

    df3["theta_A"] = df3["A_read"] / df3["cov"]
    df3["theta_A"] = df3["theta_A"].clip(0.001, 0.999)
    df3["theta_B"] = 1 - df3["theta_A"]
    df3["log_theta_A"] = np.log(df3["theta_A"])
    df3["log_theta_B"] = np.log(df3["theta_B"])

    return df3


def EM(df, max_iter=50):
    initial_theta = df.groupby("cell")[["ref", "alt"]].sum().reset_index()
    initial_theta["theta_A"] = initial_theta["alt"] / (
        initial_theta["ref"] + initial_theta["alt"]
    )
    initial_theta["theta_A"] = initial_theta["theta_A"].clip(0.001, 0.999)
    initial_theta["theta_B"] = 1 - initial_theta["theta_A"]
    initial_theta["log_theta_A"] = np.log(initial_theta["theta_A"])
    initial_theta["log_theta_B"] = np.log(initial_theta["theta_B"])

    last_snp_flip = M_step(df, initial_theta)

    for i in range(max_iter):
        updated_theta = E_step(df, last_snp_flip)
        cur_snp_flip = M_step(df, updated_theta)

        cur_F = cur_snp_flip["F"].sum()

        if i:
            # print(last_F, cur_F, abs(cur_F - last_F))
            if abs(cur_F - last_F) < 0.001:
                break

        last_F = cur_F
        last_snp_flip = cur_snp_flip

    updated_theta = E_step(df, last_snp_flip)

    return updated_theta, cur_snp_flip


def expand2d(a, rr, cr):
    a = np.asarray(a)
    rr = np.full(a.shape[0], rr, int) if np.isscalar(rr) else np.asarray(rr, int)
    cr = np.full(a.shape[1], cr, int) if np.isscalar(cr) else np.asarray(cr, int)

    return np.repeat(np.repeat(a, rr, axis=0), cr, axis=1)


def expand_rows(a, rr):
    a = np.asarray(a)
    rr = np.full(a.shape[0], rr, int) if np.isscalar(rr) else np.asarray(rr, int)

    logger.info(str(a.shape) + str(rr.shape))

    return np.repeat(a, rr, axis=0)


def get_allele_cn(T, baf):
    a_b = np.char.add(
        (T - np.rint(T * baf)).astype(int).astype(str),
        np.char.add("|", np.rint(T * baf).astype(int).astype(str))
    )

    return a_b


def flexible_load(tsv):
    has_header = False

    with open(tsv, "r") as f:
        for line in f:
            if line.startswith("chrom"):
                has_header = True
                break

    if has_header:
        df = pd.read_csv(tsv, sep="\t")
    else:
        df = pd.read_csv(tsv, sep="\t", header=None)
        df.columns = ["chrom", "pos", "cell", "ref", "alt"]

    return df


# ============================================================
# Threshold-model helper functions
# ============================================================

def load_threshold_model(threshold_file, feature="median_mBAF"):
    """
    threshold_file example:
        best_mBAF_threshold_by_CN.tsv

    Expected columns:
        feature
        rounded_total_cn
        threshold

    Rule:
        predict LOH if feature <= threshold
    """
    if threshold_file is None:
        return None

    th = pd.read_csv(threshold_file, sep="\t")

    if "feature" in th.columns:
        th = th[th["feature"] == feature].copy()

    if len(th) == 0:
        raise ValueError(f"No threshold rows found for feature={feature}")

    model = {}

    for _, row in th.iterrows():
        cn = int(row["rounded_total_cn"])
        model[cn] = float(row["threshold"])

    return model


def load_pos_baf_file(baf_dir, sample, chrom):
    """
    Load precomputed position-level BAF.

    Expected files:
        {baf_dir}/{sample}_pos_baf_chr1.tsv.gz
        or
        {baf_dir}/{sample}_pos_baf_chr1.tsv

    Expected columns:
        chrom
        pos
        depth
        baf
        mbaf

    If mbaf is missing, it is computed as min(baf, 1 - baf).
    """
    f1 = os.path.join(baf_dir, f"{sample}_pos_baf_{chrom}.tsv.gz")
    f2 = os.path.join(baf_dir, f"{sample}_pos_baf_{chrom}.tsv")

    if os.path.exists(f1):
        pos_baf_file = f1
    elif os.path.exists(f2):
        pos_baf_file = f2
    else:
        logger.warning(f"No precomputed position-level BAF file found for {sample} {chrom}")
        return None

    df_pos = pd.read_csv(pos_baf_file, sep="\t")

    if "mbaf" not in df_pos.columns:
        if "baf" not in df_pos.columns:
            raise ValueError(f"{pos_baf_file} needs either mbaf or baf column")
        df_pos["mbaf"] = np.minimum(df_pos["baf"], 1 - df_pos["baf"])

    if "depth" not in df_pos.columns:
        if ("ref_sum" in df_pos.columns) and ("alt_sum" in df_pos.columns):
            df_pos["depth"] = df_pos["ref_sum"] + df_pos["alt_sum"]
        else:
            raise ValueError(f"{pos_baf_file} needs depth column or ref_sum/alt_sum columns")

    df_pos = df_pos[df_pos["chrom"] == chrom].copy()
    df_pos["pos"] = df_pos["pos"].astype(int)

    return df_pos


def predict_loh_by_threshold_model(
    df_pos_baf,
    seg_start,
    seg_end,
    avg_total_cn,
    threshold_model,
    min_depth=10,
    min_snps=10,
    feature="median_mBAF"
):
    """
    Return:
        is_loh_pred
        feature_value
        n_snps
        rounded_total_cn
        threshold
    """
    rounded_cn = int(np.rint(avg_total_cn))

    if threshold_model is None:
        return False, np.nan, 0, rounded_cn, np.nan

    if rounded_cn not in threshold_model:
        logger.warning(f"No threshold for rounded CN={rounded_cn}; predict nonLOH")
        return False, np.nan, 0, rounded_cn, np.nan

    threshold = threshold_model[rounded_cn]

    sub = df_pos_baf[
        (df_pos_baf["pos"] > seg_start)
        & (df_pos_baf["pos"] < seg_end)
        & (df_pos_baf["depth"] >= min_depth)
    ].copy()

    n_snps = len(sub)

    if n_snps < min_snps:
        return False, np.nan, n_snps, rounded_cn, threshold

    if feature == "mean_mBAF":
        val = sub["mbaf"].mean()
    elif feature == "median_mBAF":
        val = sub["mbaf"].median()
    else:
        raise ValueError(f"Unsupported feature: {feature}")

    is_loh = val <= threshold

    return is_loh, val, n_snps, rounded_cn, threshold


def identify_loh_segments_one_chrom(
    df,
    chrom,
    sample,
    baf_dir,
    threshold_model,
    threshold_feature="median_mBAF",
    min_depth=10,
    min_snps=10
):
    """
    Fast LOH-only mode.

    Uses:
        total CN segments
        precomputed position-level BAF
        CN-specific mBAF threshold model

    Does NOT run EM.
    Does NOT load raw cell-level BAF.
    """
    df_chr = df[df["chrom"] == chrom].reset_index(drop=True)

    if len(df_chr) == 0:
        return pd.DataFrame()

    df_pos_baf = load_pos_baf_file(
        baf_dir=baf_dir,
        sample=sample,
        chrom=chrom
    )

    if df_pos_baf is None:
        logger.warning(f"{sample} {chrom}: no position-level BAF file; skip")
        return pd.DataFrame()

    cnv = df_chr.iloc[:, 3:].values

    bnds = list(np.where(((cnv[1:, :] - cnv[:-1, :]) != 0).sum(1) > 0)[0] + 1)
    bnds = [0] + bnds

    if bnds[-1] != df_chr.shape[0]:
        bnds.append(df_chr.shape[0])

    rows = []

    for seg_id in range(len(bnds) - 1):
        seg_start = int(df_chr["start"].iloc[bnds[seg_id]])
        seg_end = int(df_chr["end"].iloc[bnds[seg_id + 1] - 1])

        cn_vals = df_chr.iloc[bnds[seg_id]:bnds[seg_id + 1], 3:].values.ravel()
        cn_vals = pd.to_numeric(pd.Series(cn_vals), errors="coerce").dropna()

        avg_total_cn = float(np.nanmean(cn_vals)) if len(cn_vals) > 0 else np.nan

        pred_loh, mbaf_val, n_pos_snps, rounded_cn, cn_threshold = predict_loh_by_threshold_model(
            df_pos_baf=df_pos_baf,
            seg_start=seg_start,
            seg_end=seg_end,
            avg_total_cn=avg_total_cn,
            threshold_model=threshold_model,
            min_depth=min_depth,
            min_snps=min_snps,
            feature=threshold_feature
        )

        rows.append({
            "sample": sample,
            "chrom": chrom,
            "loh_seg_id": seg_id,
            "start": seg_start,
            "end": seg_end,
            "avg_total_cn": avg_total_cn,
            "rounded_total_cn": rounded_cn,
            "threshold_feature": threshold_feature,
            "feature_value": mbaf_val,
            "threshold": cn_threshold,
            "n_snps": n_pos_snps,
            "min_depth": min_depth,
            "min_snps": min_snps,
            "is_loh": "LOH" if pred_loh else "nonLOH"
        })

    return pd.DataFrame(rows)


# ============================================================
# Full EM mode
# ============================================================

def phase_one_chrom(
    df,
    count_file,
    chrom,
    cell_cluster,
    snp_source,
    normal_cells,
    convert=None,
    snp_file=None,
    logger=None,
    sample=None,
    baf_dir=None,
    threshold_model=None,
    threshold_feature="median_mBAF",
    min_depth=10,
    min_snps=10
):
    if convert:
        try:
            df_conv = pd.read_csv(convert, header=None, sep=" ")
            name_to_bc = dict(zip(df_conv.iloc[:, 0], df_conv.iloc[:, 1]))
            bc_to_name = dict(zip(df_conv.iloc[:, 1], df_conv.iloc[:, 0]))
        except Exception:
            df_conv = pd.read_csv(convert, header=None, sep="\t")
            name_to_bc = dict(zip(df_conv.iloc[:, 0], df_conv.iloc[:, 1]))
            bc_to_name = dict(zip(df_conv.iloc[:, 1], df_conv.iloc[:, 0]))

    secter2cell = pickle.load(open(cell_cluster, "rb"))

    cell2sector = {}
    ordered_cells = []
    n = 0

    for x, cells in secter2cell.items():
        cells = [c.split("-")[0] for c in cells]

        for cell in cells:
            cell2sector[cell] = x
            ordered_cells.append(cell)

        secter2cell[x] = cells
        n += len(cells)

    ordered_subclones = list(secter2cell.keys())

    df_baf = flexible_load(count_file)
    df_baf["cell"] = [x.split("-")[0] for x in df_baf["cell"]]
    df_baf["cov"] = df_baf["ref"] + df_baf["alt"]

    df_pos_baf = None

    if snp_source == "tumor":
        df_pos_baf = load_pos_baf_file(
            baf_dir=baf_dir,
            sample=sample,
            chrom=chrom
        )

    if snp_file is not None:
        df_snp = pd.read_csv(snp_file, sep="\t", header=None)
        df_snp.columns = ["chrom", "pos", "gt"]
        df_snp = df_snp[df_snp["chrom"] == chrom]

    df = df[df["chrom"] == chrom].reset_index(drop=True)

    if convert:
        if len(set(df_baf["cell"]) & set(df.columns[3:])) == 0:
            df_baf["cell"] = df_baf["cell"].map(bc_to_name)

    df_baf["subclone"] = df_baf["cell"].map(cell2sector)

    olp_cells = set(ordered_cells) & set(df.columns[3:]) & set(df_baf["cell"])
    ordered_cells = [x for x in ordered_cells if x in olp_cells]

    df = df[["chrom", "start", "end"] + ordered_cells]
    df.iloc[:, 3:] = df.iloc[:, 3:].clip(upper=6).astype(int)

    df_baf = df_baf[df_baf["cell"].isin(olp_cells)].reset_index(drop=True)

    for x in secter2cell:
        secter2cell[x] = list(set(secter2cell[x]) & olp_cells)

    for x in list(secter2cell.keys()):
        if len(secter2cell[x]) == 0:
            del secter2cell[x]

    ordered_subclones = [x for x in ordered_subclones if x in secter2cell.keys()]
    cell2sector = {x: cell2sector[x] for x in cell2sector if x in olp_cells}

    cnv = df.iloc[:, 3:].values

    bnds = list(np.where(((cnv[1:, :] - cnv[:-1, :]) != 0).sum(1) > 0)[0] + 1)
    bnds = [0] + bnds

    if bnds[-1] != df.shape[0]:
        bnds.append(df.shape[0])

    seg_boundary = []

    for i in range(len(bnds) - 1):
        seg_boundary.append((df["start"][bnds[i]], df["end"][bnds[i + 1] - 1]))

    print(chrom, ": ", len(seg_boundary), " segments")

    phase_segs = []
    loh_seg_ids = []
    balance_seg_ids = []

    for i in range(len(seg_boundary)):
        seg = seg_boundary[i]

        print("seg ", i + 1, seg)

        df_baf_seg = df_baf[(df_baf["pos"] > seg[0]) & (df_baf["pos"] < seg[1])]

        if snp_file is not None:
            df_snp_seg = df_snp[(df_snp["pos"] > seg[0]) & (df_snp["pos"] < seg[1])]
            homo_ratio = (df_snp_seg["gt"] == "1|1").mean()

            if homo_ratio > 0.8:
                loh_seg_ids.append(i)
            elif homo_ratio < 0.6:
                df1 = df_baf_seg.groupby("pos")[["ref", "alt"]].sum().reset_index()
                df1["baf"] = df1["alt"] / (df1["ref"] + df1["alt"])
                df1["cov"] = df1["ref"] + df1["alt"]

                use_pos = set(
                    df1[
                        (df1["baf"] > 0.1)
                        & (df1["baf"] < 0.9)
                        & (df1["cov"] > 5)
                    ]["pos"]
                )

                df_baf_seg = df_baf_seg[df_baf_seg["pos"].isin(use_pos)].reset_index()

        if True:
            logger.info("snp_source:" + snp_source)
            logger.info("Perform filter...")
            logger.info(f"Before filter: {df_baf_seg.shape[0]} snps")

            df1 = df_baf_seg.groupby("pos")[["ref", "alt"]].sum().reset_index()
            df1["baf"] = df1["alt"] / (df1["ref"] + df1["alt"])
            df1["cov"] = df1["ref"] + df1["alt"]

            homo_ratio = ((df1["baf"] < 0.2) | (df1["baf"] > 0.8)).mean()
            avg_total_cn = np.nanmean(
                    pd.to_numeric(
                        pd.Series(df.iloc[bnds[i]:bnds[i + 1], 3:].values.ravel()),
                        errors="coerce"
                    )
                )
            if snp_source == "tumor":
                

                if df_pos_baf is not None:
                    pred_loh, mbaf_val, n_pos_snps, rounded_cn, cn_threshold = predict_loh_by_threshold_model(
                        df_pos_baf=df_pos_baf,
                        seg_start=seg[0],
                        seg_end=seg[1],
                        avg_total_cn=avg_total_cn,
                        threshold_model=threshold_model,
                        min_depth=min_depth,
                        min_snps=min_snps,
                        feature=threshold_feature
                    )

                    logger.info(
                        f"{chrom} seg={i} start={seg[0]} end={seg[1]} "
                        f"avgCN={avg_total_cn:.4f} roundedCN={rounded_cn} "
                        f"{threshold_feature}={mbaf_val} n_pos_snps={n_pos_snps} "
                        f"threshold={cn_threshold} pred_LOH={pred_loh}"
                    )

                    print(f"{chrom} seg={i} start={seg[0]} end={seg[1]} "
                        f"avgCN={avg_total_cn:.4f} roundedCN={rounded_cn} "
                        f"{threshold_feature}={mbaf_val} n_pos_snps={n_pos_snps} "
                        f"threshold={cn_threshold} pred_LOH={pred_loh}")

                    if pred_loh:
                        loh_seg_ids.append(i)
                else:
                    logger.warning(
                        f"{chrom} seg={i}: no position-level BAF available; fallback to homo_ratio"
                    )
                    if homo_ratio > 0.4:
                        loh_seg_ids.append(i)
            if snp_source == 'tumor':
                if (i not in loh_seg_ids) and (round(avg_total_cn)<4):
                    # for non-loh segs
                    #   snp source is tumor: filter out SNPs with extreme baf for CN<5 non-loh seg
                    #   snp source is tumor-normal: filter out SNPs with extreme baf for all non-loh seg
                    use_pos = set(
                        df1[
                            (df1["baf"] > 0.3)
                            & (df1["baf"] < 0.7)
                            & (df1["cov"] > 2)
                        ]["pos"]
                    )

                    df_baf_seg = df_baf_seg[df_baf_seg["pos"].isin(use_pos)].reset_index()

                    logger.info(f"After filter: {df_baf_seg.shape[0]} snps")
            else:
                if (i not in loh_seg_ids) :
                    # for non-loh segs
                    #   snp source is tumor: filter out SNPs with extreme baf for CN<5 non-loh seg
                    #   snp source is tumor-normal: filter out SNPs with extreme baf for all non-loh seg
                    use_pos = set(
                        df1[
                            (df1["baf"] > 0.3)
                            & (df1["baf"] < 0.7)
                            & (df1["cov"] > 2)
                        ]["pos"]
                    )

                    df_baf_seg = df_baf_seg[df_baf_seg["pos"].isin(use_pos)].reset_index()

                    logger.info(f"After filter: {df_baf_seg.shape[0]} snps")
                

        updated_theta, df_phase_seg = EM(df_baf_seg)

        try:
            df_phase_seg["ref_is_A"] = df_phase_seg["ref_is_A"].round().astype(int)
        except Exception:
            print(set(df_phase_seg["ref_is_A"]))
            exit()

        s = pd.to_numeric(df.iloc[bnds[i], 3:], errors="coerce")
        even = (s.dropna().astype(int) % 2 == 0)
        frac_even = even.mean()

        most_cn_evn = frac_even > 0.95

        df1 = df_baf_seg.groupby("pos")[["ref", "alt"]].sum().reset_index()
        df1["baf"] = df1["alt"] / (df1["ref"] + df1["alt"])

        og_vaf_balance = is_normal_centered_at_half(df1["baf"])
        em_baf_balance = is_single_normal_mu_in_range(updated_theta["theta_A"])

        if most_cn_evn & og_vaf_balance & em_baf_balance:
            balance_seg_ids.append(i)

        phase_segs.append(df_phase_seg)

    df_phase = pd.concat(phase_segs)
    df_phase["z"] = df_phase["ref_is_A"]
    df_phase = df_phase[["chrom", "pos", "z"]]

    df_baf = df_phase.merge(df_baf)

    df_baf["bin_id"] = df_baf["pos"] // 5_000_000
    bin_id_cn = (df["start"] // 5_000_000).tolist()

    seg_map = {i: bin_id_cn[bnds[i]:bnds[i + 1]] for i in range(len(bnds) - 1)}

    bin_to_seg = {}

    for seg_id, bin_ids in seg_map.items():
        for bin_id in bin_ids:
            bin_to_seg[bin_id] = seg_id

    df_baf["seg_id"] = df_baf["bin_id"].map(bin_to_seg)

    df_baf["B"] = df_baf["ref"].where(df_baf["z"].eq(0), df_baf["alt"])

    if snp_file is not None:
        mask = df_baf["seg_id"].isin(set(loh_seg_ids))
        df_baf.loc[mask, "B"] = df_baf.loc[mask, "cov"].to_numpy()

    if snp_source == "tumor":
        # baf based LOH decision overides EM and balance seg
        # print(set(df_baf["seg_id"]), set(loh_seg_ids))
        mask = df_baf["seg_id"].isin(set(loh_seg_ids))
        df_baf.loc[mask, "B"] = df_baf.loc[mask, "cov"].to_numpy()
        # print('=======')
        # print(df_baf.head())

    balance_seg_ids = set(balance_seg_ids) - set(loh_seg_ids)
    # print(loh_seg_ids)
    # print(balance_seg_ids)
    df_seg = df_baf.groupby(["seg_id", "cell"])[["B", "cov"]].sum().reset_index()

    mask_balance = (
        df_seg["seg_id"].isin(set(balance_seg_ids))
        | df_seg["cell"].isin(set(normal_cells))
    )

    df_seg.loc[mask_balance, "B"] = (
        df_seg.loc[mask_balance, "cov"] / 2
    ).to_numpy().round()


    
    df_seg["baf_seg"] = df_seg["B"] / df_seg["cov"]
    print(df_seg.head())

    df_sub = df_baf.groupby(["seg_id", "subclone"])[["B", "cov"]].sum().reset_index()

    mask_balance = df_sub["seg_id"].isin(set(balance_seg_ids))

    df_sub.loc[mask_balance, "B"] = (
        df_sub.loc[mask_balance, "cov"] / 2
    ).to_numpy().round()

    df_sub["baf_sub"] = df_sub["B"] / df_sub["cov"]

    baf_seg = df_seg.pivot(index="seg_id", columns="cell", values="baf_seg")
    baf_seg = baf_seg.reindex(columns=ordered_cells, fill_value=0.5)

    baf_seg.index = baf_seg.index.astype(int)

    print(baf_seg)

    full_index = np.arange(baf_seg.index.min(), len(bnds) - 1)
    baf_seg = baf_seg.reindex(full_index, fill_value=0.5)

    logger.info(str(baf_seg.index) + str(baf_seg.shape))

    baf_sub = df_sub.pivot(index="seg_id", columns="subclone", values="baf_sub")[ordered_subclones].fillna(0.5)
    baf_sub.index = baf_sub.index.astype(int)

    full_index = np.arange(baf_sub.index.min(), len(bnds) - 1)
    baf_sub = baf_sub.reindex(full_index, fill_value=0.5)

    row_repeats = [bnds[i + 1] - bnds[i] for i in baf_seg.index.tolist()]

    logger.info(f"{len(bnds) - 1}, {len(row_repeats)} " + str(row_repeats))

    col_repeats = [len(cells) for x, cells in secter2cell.items()]

    baf_seg_expand = expand_rows(baf_seg, row_repeats)
    baf_sub_expand = expand2d(baf_sub, row_repeats, col_repeats)

    baf_seg_expand = np.where(np.isnan(baf_seg_expand), baf_sub_expand, baf_seg_expand)

    logger.info(str(bnds))

    ascn_seg = get_allele_cn(cnv, baf_seg_expand)
    ascn_sub = get_allele_cn(cnv, baf_sub_expand)

    ascn_seg_df = pd.concat(
        [df.iloc[:, :3], pd.DataFrame(ascn_seg, columns=ordered_cells)],
        axis=1
    )

    ascn_sub_df = pd.concat(
        [df.iloc[:, :3], pd.DataFrame(ascn_sub, columns=ordered_cells)],
        axis=1
    )

    baf_seg_df = pd.concat(
        [df.iloc[:, :3], pd.DataFrame(baf_seg_expand, columns=ordered_cells)],
        axis=1
    )

    baf_sub_df = pd.concat(
        [df.iloc[:, :3], pd.DataFrame(baf_sub_expand, columns=ordered_cells)],
        axis=1
    )

    return ascn_seg_df, ascn_sub_df, baf_seg_df, baf_sub_df


# ============================================================
# Plotting
# ============================================================

import argparse
import re
from matplotlib.colors import ListedColormap, BoundaryNorm


def plot_phased_heatmap_unordered(
    arr_str,
    default_color="white",
    row_labels=None,
    col_labels=None,
    chroms=None,
    out=None,
    figsize=(12, 6),
    draw_chrom_boundaries=True,
    title=None,
    sort=True
):
    if sort:
        arr_str, labels = hclust_categorical_sorted(arr_str)

    arr_str = np.array([
        ["|".join(sorted(s.split("|"))) for s in row]
        for row in arr_str
    ])

    palette = {}

    palette.update({(0, 0): "darkblue"})
    palette.update({(1, 0): "lightblue"})
    palette.update({(1, 1): "lightgray", (2, 0): "dimgray"})
    palette.update({(2, 1): "lightgoldenrodyellow", (3, 0): "gold"})
    palette.update({(2, 2): "navajowhite", (3, 1): "orange", (4, 0): "darkorange"})
    palette.update({(3, 2): "salmon", (4, 1): "red", (5, 0): "darkred"})
    palette.update({(3, 3): "plum", (4, 2): "orchid", (5, 1): "purple", (6, 0): "indigo"})

    palette_str = {}

    for (a, b), color in palette.items():
        a2, b2 = (a, b) if a <= b else (b, a)
        palette_str[f"{a2}|{b2}"] = color

    if not isinstance(arr_str, np.ndarray) or arr_str.ndim != 2:
        raise ValueError("arr_str must be a 2D NumPy array of strings.")

    def norm_key(s):
        if s is None or s != s:
            return None

        if "|" not in s:
            return None

        a, b = s.split("|", 1)

        try:
            a, b = int(a), int(b)
        except Exception:
            return None

        a2, b2 = (a, b) if a <= b else (b, a)

        return f"{a2}|{b2}"

    norm = np.vectorize(norm_key, otypes=[object])(arr_str)

    flat = norm.ravel()
    mask = np.array([x is not None for x in flat], dtype=bool)
    present = np.unique(flat[mask])

    def parse_pair(s):
        a, b = s.split("|", 1)
        return int(a), int(b)

    present_sorted = sorted(
        present.tolist(),
        key=lambda s: (lambda a, b: (a + b, a))(*parse_pair(s))
    )

    colors = [palette_str.get(k, default_color) for k in present_sorted]
    cat2int = {k: i for i, k in enumerate(present_sorted)}

    code = np.vectorize(lambda x: cat2int.get(x, -1), otypes=[int])(norm)

    if len(colors) == 0:
        cmap = ListedColormap([default_color])
        boundaries = np.array([-0.5, 0.5])
    else:
        cmap = ListedColormap(colors)
        boundaries = np.arange(-0.5, len(colors) + 0.5, 1)

    norm_bins = BoundaryNorm(boundaries, cmap.N, clip=False)
    cmap.set_under(default_color)

    plt.figure(figsize=figsize)

    im = plt.imshow(
        code,
        aspect="auto",
        cmap=cmap,
        norm=norm_bins,
        interpolation="nearest"
    )

    if len(colors) > 0:
        cbar = plt.colorbar(im, ticks=np.arange(len(colors)))
        cbar.ax.set_yticklabels(present_sorted)

    r, c = arr_str.shape

    if row_labels is None:
        row_labels = [f"row{i}" for i in range(r)]

    if chroms is not None:
        if len(chroms) != c:
            raise ValueError("Length of chroms must equal number of columns in arr_str.")

        xticks, xticklabels, boundaries_x = [], [], []
        start = 0

        for j in range(1, c + 1):
            if j == c or chroms[j] != chroms[j - 1]:
                end = j - 1
                mid = (start + end) / 2.0

                xticks.append(mid)
                xticklabels.append(str(chroms[start]))

                if j < c:
                    boundaries_x.append(j - 0.5)

                start = j

        plt.xticks(xticks, xticklabels, rotation=90)

        if draw_chrom_boundaries:
            for x in boundaries_x:
                plt.axvline(x=x, linestyle=":", linewidth=0.8, alpha=0.5)
    else:
        if col_labels is None:
            col_labels = [f"col{j}" for j in range(c)]

        plt.xticks(range(c), col_labels, rotation=90)

    if title is None:
        plt.title("ASCN heatmap (unordered a|b, a≤b)")
    else:
        plt.title(title)

    plt.xlabel("Segments")
    plt.ylabel("Cells")
    plt.tight_layout()

    if out:
        outdir = os.path.dirname(out)

        if outdir:
            os.makedirs(outdir, exist_ok=True)

        plt.savefig(out + ".pdf")
        plt.savefig(out + ".png")

    plt.show()

    return code, cat2int, colors


def run_one_chrom(
    i,
    baf_dir,
    sample,
    snp_source,
    df,
    convert,
    snp_file,
    cell_cluster_file,
    normal_cells,
    logger,
    threshold_model=None,
    threshold_feature="median_mBAF",
    min_depth=10,
    min_snps=10
):
    chrom = f"chr{i}"

    count_file = f"{baf_dir}/{sample}_baf_{chrom}.tsv"

    out = phase_one_chrom(
        df=df,
        count_file=count_file,
        chrom=chrom,
        cell_cluster=cell_cluster_file,
        snp_source=snp_source,
        normal_cells=normal_cells,
        convert=convert,
        snp_file=snp_file,
        logger=logger,
        sample=sample,
        baf_dir=baf_dir,
        threshold_model=threshold_model,
        threshold_feature=threshold_feature,
        min_depth=min_depth,
        min_snps=min_snps
    )

    return out


def allele_CN(
    total_cn_dir,
    outdir,
    baf_dir,
    sample,
    convert,
    snp_source,
    snp_file,
    n_thread=10,
    threshold_file=None,
    threshold_feature="median_mBAF",
    min_depth=10,
    min_snps=10,
    loh_only=False,
    only_chrom=None
):
    cell_cluster_file = f"{total_cn_dir}/fine_cluster_result.pkl"
    cn_file = f"{total_cn_dir}/{sample}_smooth_seg.csv"

    if not os.path.exists(outdir):
        os.system("mkdir -p " + outdir)

    df = pd.read_csv(cn_file)
    df.columns = ["chrom", "start", "end"] + [x.split("-")[0] for x in df.columns[3:]]

    cell_cols = df.columns[3:]
    means = df[cell_cols].mean()

    tol = 0.1
    normal_cells = means.index[abs(means - 2) <= tol]

    print("Total cells:", len(cell_cols))
    print("Normal cells:", len(normal_cells))
    print("Percentage normal: {:.2f}%".format(100 * len(normal_cells) / len(cell_cols)))

    chroms = set(df["chrom"])
    
    useful_is = []

    for i in tqdm(range(1, 23)):
        chrom = f"chr{i}"

        if chrom in chroms:
            useful_is.append(i)

    if only_chrom is not None:
        only_chrom = str(only_chrom)

        if only_chrom.startswith("chr"):
            only_i = int(only_chrom.replace("chr", ""))
        else:
            only_i = int(only_chrom)

        if only_i not in useful_is:
            raise ValueError(f"{only_chrom} is not present in CN file.")

        useful_is = [only_i]

        # Usually no need for many jobs if running one chromosome
        n_thread = 1

        print("Chromosomes to run:", [f"chr{i}" for i in useful_is])

    print(snp_source)

    threshold_model = None

    if snp_source == "tumor":
        threshold_model = load_threshold_model(
            threshold_file=threshold_file,
            feature=threshold_feature
        )

        print("Loaded tumor LOH threshold model:")
        print(threshold_model)

    # ========================================================
    # Fast LOH-only mode
    # ========================================================
    if loh_only:
        # if snp_source != "tumor":
        #     raise ValueError("--loh_only currently requires --snp_source tumor")

        if threshold_model is None:
            raise ValueError("--loh_only requires --threshold_file")

        print("Running LOH-only mode: skip EM, only identify LOH segments.")

        loh_parts = Parallel(n_jobs=n_thread)(
            delayed(identify_loh_segments_one_chrom)(
                df=df,
                chrom=f"chr{i}",
                sample=sample,
                baf_dir=baf_dir,
                threshold_model=threshold_model,
                threshold_feature=threshold_feature,
                min_depth=min_depth,
                min_snps=min_snps
            )
            for i in tqdm(useful_is)
        )

        loh_parts = [
            x for x in loh_parts
            if x is not None and len(x) > 0
        ]

        if len(loh_parts) == 0:
            raise ValueError("No LOH segments generated. Check BAF files and threshold model.")

        loh_df = pd.concat(loh_parts, ignore_index=True)

        out_tsv = f"{outdir}/quick_loh_segments.tsv"
        out_bed = f"{outdir}/quick_loh_segments.bed"
        out_loh_bed = f"{outdir}/quick_loh_only_segments.bed"

        loh_df.to_csv(out_tsv, sep="\t", index=False)

        loh_df[
            [
                "chrom", "start", "end",
                "avg_total_cn", "rounded_total_cn",
                "threshold_feature", "feature_value",
                "threshold", "n_snps", "is_loh"
            ]
        ].to_csv(
            out_bed,
            sep="\t",
            index=False,
            header=False
        )

        loh_df[loh_df["is_loh"] == "LOH"][
            [
                "chrom", "start", "end",
                "avg_total_cn", "rounded_total_cn",
                "threshold_feature", "feature_value",
                "threshold", "n_snps", "is_loh"
            ]
        ].to_csv(
            out_loh_bed,
            sep="\t",
            index=False,
            header=False
        )

        print(f"[OUT] {out_tsv}")
        print(f"[OUT] {out_bed}")
        print(f"[OUT] {out_loh_bed}")
        print("Done LOH-only mode.")

        return loh_df

    # ========================================================
    # Full EM mode
    # ========================================================
    if n_thread == 1:
        outs = []

        for i in tqdm(useful_is):
            logger.info(f"============={i}")

            out = run_one_chrom(
                i,
                baf_dir,
                sample,
                snp_source,
                df,
                convert,
                snp_file,
                cell_cluster_file,
                normal_cells,
                logger,
                threshold_model,
                threshold_feature,
                min_depth,
                min_snps
            )

            outs.append(out)
    else:
        outs = Parallel(n_jobs=n_thread)(
            delayed(run_one_chrom)(
                i,
                baf_dir,
                sample,
                snp_source,
                df,
                convert,
                snp_file,
                cell_cluster_file,
                normal_cells,
                logger,
                threshold_model,
                threshold_feature,
                min_depth,
                min_snps
            )
            for i in tqdm(useful_is)
        )

    df = pd.concat([out[0] for out in outs], ignore_index=True)
    df.to_csv(f"{outdir}/ascn_seg.csv", index=False)

    df = pd.concat([out[1] for out in outs], ignore_index=True)
    df.to_csv(f"{outdir}/ascn_sub.csv", index=False)

    df = pd.concat([out[2] for out in outs], ignore_index=True)
    dfs = [out[2] for out in outs]
    pickle.dump(dfs, open(f"{outdir}/baf_seg.pkl", "wb"))
    df.to_csv(f"{outdir}/baf_seg.csv", index=False)

    df = pd.concat([out[3] for out in outs], ignore_index=True)
    dfs = [out[3] for out in outs]
    pickle.dump(dfs, open(f"{outdir}/baf_sub.pkl", "wb"))
    df.to_csv(f"{outdir}/baf_sub.csv", index=False)

    df = pd.read_csv(f"{outdir}/ascn_seg.csv")
    df = df.fillna("0|6")
    arr = df.iloc[:, 3:].values.T

    plot_phased_heatmap_unordered(
        arr,
        default_color="white",
        chroms=df["chrom"],
        figsize=(12, 6),
        draw_chrom_boundaries=True,
        title=f"{sample} ASCN smooth by seg",
        out=f"{outdir}/{sample}_ASCN_smooth_seg"
    )

    df = pd.read_csv(f"{outdir}/ascn_sub.csv")
    df = df.fillna("0|6")
    arr = df.iloc[:, 3:].values.T

    plot_phased_heatmap_unordered(
        arr,
        default_color="white",
        chroms=df["chrom"],
        figsize=(12, 6),
        draw_chrom_boundaries=True,
        title=f"{sample} ASCN smooth by subclone",
        out=f"{outdir}/{sample}_ASCN_smooth_sub"
    )


if __name__ == "__main__":

    def parse_args():
        parser = argparse.ArgumentParser(
            description="Run allele-specific copy number inference."
        )

        parser.add_argument(
            "--total_cn_dir",
            required=True,
            help="Directory containing total copy-number input files."
        )

        parser.add_argument(
            "--outdir",
            required=True,
            help="Output directory."
        )

        parser.add_argument(
            "--baf_dir",
            required=True,
            help="Directory containing per-chromosome BAF files."
        )

        parser.add_argument(
            "--sample",
            required=True,
            help="Sample name."
        )

        parser.add_argument(
            "--convert",
            default=None,
            help="Optional cell-name / barcode converter file."
        )

        parser.add_argument(
            "--snp_source",
            default="tumor-normal",
            choices=["tumor", "tumor-normal", "matched-normal"],
            help="SNP source type."
        )

        parser.add_argument(
            "--snp_file",
            default=None,
            help="Optional SNP genotype file."
        )

        parser.add_argument(
            "--n_thread",
            "-t",
            type=int,
            default=10,
            help="Number of threads."
        )

        parser.add_argument(
            "--threshold_file",
            default=None,
            help="Threshold model TSV, e.g. best_mBAF_threshold_by_CN.tsv"
        )

        parser.add_argument(
            "--threshold_feature",
            default="median_mBAF",
            choices=["median_mBAF", "mean_mBAF"],
            help="Feature used by threshold model"
        )

        parser.add_argument(
            "--min_depth",
            type=int,
            default=10,
            help="Minimum depth for position-level BAF SNPs"
        )

        parser.add_argument(
            "--min_snps",
            type=int,
            default=10,
            help="Minimum number of SNP positions required in a segment"
        )

        parser.add_argument(
            "--loh_only",
            action="store_true",
            help="Only identify LOH/nonLOH segments using precomputed position-level BAF and threshold model; skip EM/phasing."
        )

        parser.add_argument(
            "--chrom",
            default=None,
            help="Optional chromosome to run only one chromosome, e.g. chr20 or 20"
        )
        
        return parser.parse_args()

    args = parse_args()

    allele_CN(
        total_cn_dir=args.total_cn_dir,
        outdir=args.outdir,
        baf_dir=args.baf_dir,
        sample=args.sample,
        convert=args.convert,
        snp_source=args.snp_source,
        snp_file=args.snp_file,
        n_thread=args.n_thread,
        threshold_file=args.threshold_file,
        threshold_feature=args.threshold_feature,
        min_depth=args.min_depth,
        min_snps=args.min_snps,
        loh_only=args.loh_only,
        only_chrom=args.chrom
    )