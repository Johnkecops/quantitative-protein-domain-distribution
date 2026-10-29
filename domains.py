#!/usr/bin/env python3
"""
Module: domains
Purpose: Domain annotation parsing, counting and co-occurrence statistics after
         Parikesit, Stadler & Prohaska (2010).
Author: Dr. Arli Aditya Parikesit
Date: 2026
Parameters:
    - evalue: per-domain (independent) E-value cut-off, paper default 1e-4
    - n_proteins: size of the peptide set searched (denominator of expectations)
References:
    - Parikesit AA, Stadler PF, Prohaska SJ (2010). Quantitative comparison of
      genomic-wide protein domain distributions. GCB 2010, LNI P-173: 93-102.
    - Eddy SR (2011). Accelerated profile HMM searches. PLoS Comput Biol 7:e1002195.
    - Larralde M, Zeller G (2023). PyHMMER. Bioinformatics 39:btad214.
"""
import io
import math
import os
import zlib

import numpy as np
import pandas as pd
from scipy import stats

COLS = ["protein", "domain", "start", "end", "score", "evalue"]
PAPER_EVALUE = 1e-4
# Resource limits for uploads (shared deployments). 150k seqs covers a human genscan set (118,894).
MAX_BYTES = 500 * 2**20
MAX_SEQS = 150_000
MAX_HMMS = 200
SEARCH_CPUS = min(4, os.cpu_count() or 1)


def _gunzip(data, limit):
    """Bounded gunzip, multi-member aware; refuses output above `limit` bytes (gzip bombs)."""
    out = bytearray()
    while data:
        d = zlib.decompressobj(wbits=31)
        buf = data
        while buf and not d.eof:
            out += d.decompress(buf, limit + 1 - len(out))
            if len(out) > limit:
                raise ValueError(f"Decompressed input exceeds {limit // 2**20} MB limit")
            buf = d.unconsumed_tail
        data = d.unused_data
    return bytes(out)


def _text(src, limit=MAX_BYTES):
    """Bytes / str / file-like -> str, transparently gunzipping, size-capped."""
    data = src.read() if hasattr(src, "read") else src
    if isinstance(data, str):
        return data
    if data[:2] == b"\x1f\x8b":
        data = _gunzip(data, limit)
    if len(data) > limit:
        raise ValueError(f"Input exceeds {limit // 2**20} MB limit")
    return data.decode()


def defang(s):
    """Neutralise spreadsheet formula prefixes in IDs from uploads (CSV export from tables)."""
    return s.map(lambda v: "'" + v if isinstance(v, str) and v[:1] in ("=", "+", "-", "@", "\t", "\r") else v)


def resolve_overlaps(df):
    """Per protein keep the best-scoring hit of every overlapping group (paper rule)."""
    kept = []
    for _, g in df.sort_values("score", ascending=False).groupby("protein", sort=False):
        taken = []
        for row in g.itertuples(index=False):
            # ponytail: O(k^2) per protein, fine for k < few hundred hits; use an interval tree beyond
            if all(row.end < s or row.start > e for s, e in taken):
                taken.append((row.start, row.end))
                kept.append(row)
    out = pd.DataFrame(kept, columns=df.columns)
    return out.sort_values(["protein", "start"]).reset_index(drop=True)


def _finish(df, evalue):
    df = df[df.evalue <= evalue].astype({"start": int, "end": int})
    df = df.assign(protein=defang(df.protein), domain=defang(df.domain))
    return resolve_overlaps(df) if len(df) else pd.DataFrame(columns=COLS)


def read_domtblout(src, evalue=PAPER_EVALUE, program="hmmsearch"):
    """HMMER3 --domtblout. Sequence coords are ali/env columns in both programs;
    only target/query roles swap between hmmsearch and hmmscan."""
    rows = []
    for line in _text(src).splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        f = line.split(maxsplit=22)
        prot, dom = (f[0], f[3]) if program == "hmmsearch" else (f[3], f[0])
        rows.append((prot, dom, f[19], f[20], float(f[13]), float(f[12])))
    return _finish(pd.DataFrame(rows, columns=COLS), evalue)


def read_interproscan(src, evalue=PAPER_EVALUE, analysis=None):
    """InterProScan TSV. Column 9 is the E-value; score = -log10(E) for ranking."""
    df = pd.read_csv(io.StringIO(_text(src)), sep="\t", header=None, usecols=range(9),
                     names=["protein", "md5", "len", "analysis", "domain", "desc",
                            "start", "end", "evalue"], na_values="-")
    if analysis:
        df = df[df.analysis.str.lower() == analysis.lower()]
    df["evalue"] = df.evalue.fillna(0.0)
    df["score"] = -np.log10(df.evalue.clip(lower=1e-300))
    return _finish(df[COLS], evalue)


def read_pairs(src):
    """Two-column TSV protein<TAB>domain; no coordinates, so no overlap rule."""
    df = pd.read_csv(io.StringIO(_text(src)), sep="\t", header=None, comment="#",
                     usecols=[0, 1], names=["protein", "domain"], dtype=str)
    return df.assign(protein=defang(df.protein), domain=defang(df.domain),
                     start=0, end=0, score=np.nan, evalue=np.nan)[COLS]


def _s(x):
    return x.decode() if isinstance(x, bytes) else x


def run_hmmsearch(fasta, hmm, evalue=PAPER_EVALUE, cpus=SEARCH_CPUS):
    """In-process hmmsearch with pyhmmer. Returns (hits, n_sequences).
    Refuses jobs above MAX_SEQS sequences or MAX_HMMS models."""
    import pyhmmer

    abc = pyhmmer.easel.Alphabet.amino()
    try:
        with pyhmmer.plan7.HMMFile(io.BytesIO(_text(hmm).encode())) as hf:
            hmms = []
            for h in hf:
                hmms.append(h)
                if len(hmms) > MAX_HMMS:
                    raise ValueError(f"More than {MAX_HMMS} HMMs; run large libraries offline and upload domtblout")
        with pyhmmer.easel.SequenceFile(io.BytesIO(_text(fasta).encode()), format="fasta",
                                        digital=True, alphabet=abc) as sf:
            seqs = sf.read_block(sequences=MAX_SEQS + 1)
    except EOFError as e:  # pyhmmer signals empty/truncated input this way
        raise ValueError(str(e)) from e
    if len(seqs) > MAX_SEQS:
        raise ValueError(f"More than {MAX_SEQS} sequences; run large proteomes offline and upload domtblout")
    rows = []
    for hits in pyhmmer.hmmer.hmmsearch(hmms, seqs, cpus=cpus, domE=evalue):
        dom_name = _s(hits.query.name)
        for hit in hits:
            for d in hit.domains:
                rows.append((_s(hit.name), dom_name, d.env_from, d.env_to, d.score, d.i_evalue))
    return _finish(pd.DataFrame(rows, columns=COLS), evalue), len(seqs)


def apply_mapping(df, mapping):
    """Collapse domain IDs to families (e.g. several Pfam -> one SCOP superfamily)."""
    return df.assign(domain=defang(df.domain.map(mapping).fillna(df.domain))) if mapping else df


def counts(df):
    """Per domain: number of domain copies and number of proteins carrying it."""
    return (df.groupby("domain")
              .agg(domains=("protein", "size"), genes=("protein", "nunique"))
              .sort_values("domains", ascending=False))


def compare_counts(a, b):
    """Fig. 1 statistics on two count series sharing an index.
    Correlations and OLS use log10 of entries > 0 in both."""
    m = pd.concat({"x": a, "y": b}, axis=1).fillna(0)
    pos = m[(m.x > 0) & (m.y > 0)]
    lx, ly = np.log10(pos.x), np.log10(pos.y)
    res = {"n_shared": len(pos), "n_x_only": int(((m.x > 0) & (m.y == 0)).sum()),
           "n_y_only": int(((m.y > 0) & (m.x == 0)).sum())}
    if len(pos) >= 3:
        lr = stats.linregress(lx, ly)
        res.update(pearson_r=lr.rvalue, pearson_p=lr.pvalue,
                   spearman_rho=stats.spearmanr(pos.x, pos.y).statistic,
                   slope=lr.slope, intercept=lr.intercept,
                   median_log2_ratio=float(np.median(np.log2(pos.y / pos.x))))
    m["log2_ratio"] = np.log2((m.y + 0.5) / (m.x + 0.5))
    return m, res


def cooccurrence_test(k, n_a, n_b, n):
    """Observed k proteins with both A and B among n proteins.
    Hypergeometric null; returns expected, both tails and log2(obs/exp)."""
    exp = n_a * n_b / n
    hg = stats.hypergeom(n, n_a, n_b)
    return {"observed": k, "expected": exp,
            "log2_obs_exp": math.log2((k + 0.5) / (exp + 0.5)),
            "p_avoid": float(hg.cdf(k)), "p_attract": float(hg.sf(k - 1))}


def paper_expectation(d_a, d_b, n):
    """Paper's estimate from domain copy counts (e.g. 5090*274/118894 = 11.7)
    with a Poisson probability of seeing zero co-occurrences."""
    exp = d_a * d_b / n
    return {"expected": exp, "p_zero": math.exp(-exp)}


def pair_detail(df, a, b, n):
    """Gene-level test plus the paper's domain-level expectation for one pair."""
    ga, gb = set(df.protein[df.domain == a]), set(df.protein[df.domain == b])
    both = sorted(ga & gb)
    res = cooccurrence_test(len(both), len(ga), len(gb), n)
    res.update(genes_a=len(ga), genes_b=len(gb),
               domains_a=int((df.domain == a).sum()), domains_b=int((df.domain == b).sum()))
    res["paper"] = paper_expectation(res["domains_a"], res["domains_b"], n)
    return res, df[df.protein.isin(both) & df.domain.isin([a, b])]


def cooccurrence_matrix(df, domains, n):
    """All pairs among `domains`: gene-level test, two-sided p, BH-FDR."""
    sets = {d: set(df.protein[df.domain == d]) for d in domains}
    rows = []
    for i, a in enumerate(domains):
        for b in domains[i + 1:]:
            r = cooccurrence_test(len(sets[a] & sets[b]), len(sets[a]), len(sets[b]), n)
            r.update(a=a, b=b, p=min(1.0, 2 * min(r["p_avoid"], r["p_attract"])))
            rows.append(r)
    out = pd.DataFrame(rows)
    if len(out):
        out["q_bh"] = stats.false_discovery_control(out.p)
    return out
