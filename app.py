#!/usr/bin/env python3
"""
Module: app
Purpose: Streamlit front end for quantitative protein domain distribution and
         co-occurrence analysis, after Parikesit, Stadler & Prohaska (2010).
Author: Dr. Arli Aditya Parikesit
Date: 2026
Parameters:
    - E-value cut-off (sidebar), default 1e-4 as in the paper
References:
    - Parikesit AA, Stadler PF, Prohaska SJ (2010). GCB 2010, LNI P-173: 93-102.
Run:
    streamlit run app.py
"""
from pathlib import Path

import altair as alt
import numpy as np
import pandas as pd
import streamlit as st

import domains as D

DATA = Path(__file__).parent / "data"
st.set_page_config(page_title="Protein Domain Distribution", layout="wide")


@st.cache_data(show_spinner="Running hmmsearch (pyhmmer)...", max_entries=8, ttl=3600)
def _hmmsearch_cached(fasta: bytes, hmm: bytes, evalue: float):
    return D.run_hmmsearch(fasta, hmm, evalue)


def hmmsearch_cached(fasta, hmm, evalue):
    try:
        return _hmmsearch_cached(fasta, hmm, evalue)
    except ImportError:
        st.sidebar.error("pyhmmer is not installed in the Python running Streamlit. "
                         "Run `python -m pip install pyhmmer` with that interpreter, or start the "
                         "app with `.venv/bin/streamlit run app.py`. Tabs 1-2 and annotation "
                         "uploads work without it.")
    except ValueError as e:
        st.sidebar.error(f"hmmsearch refused or failed: {e}")
    return None, None


def parse(file, fmt, evalue, program):
    data = file.getvalue()
    try:
        if fmt == "HMMER domtblout":
            return D.read_domtblout(data, evalue, program)
        if fmt == "InterProScan TSV":
            return D.read_interproscan(data, evalue)
        return D.read_pairs(data)
    except ValueError as e:
        st.error(f"Could not read {file.name}: {e}")
        st.stop()


def scatter(m, xlab, ylab, res):
    pos = m[(m.x > 0) & (m.y > 0)].reset_index(names="label")
    lo, hi = 0.8, max(pos.x.max(), pos.y.max()) * 1.25
    scale = alt.Scale(type="log", domain=[lo, hi])
    pts = alt.Chart(pos).mark_circle(size=60).encode(
        x=alt.X("x", scale=scale, title=xlab), y=alt.Y("y", scale=scale, title=ylab),
        tooltip=["label", "x", "y", alt.Tooltip("log2_ratio", format=".2f")])
    xs = np.geomspace(lo, hi, 50)
    lines = [pd.DataFrame({"x": xs, "y": xs, "line": "identity"})]
    if "slope" in res:
        lines.append(pd.DataFrame({"x": xs, "y": 10 ** (res["intercept"] + res["slope"] * np.log10(xs)),
                                   "line": "OLS on log10"}))
    ln = alt.Chart(pd.concat(lines)).mark_line(strokeDash=[4, 3]).encode(
        x="x", y="y", color=alt.Color("line", title=None))
    return (ln + pts).properties(height=480).interactive()


# ---------------- sidebar: annotation source for tabs 3 and 4 ----------------
sb = st.sidebar
sb.header("Annotation source")
evalue = sb.select_slider("Per-domain E-value cut-off", [1e-10, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2],
                          value=D.PAPER_EVALUE, format_func=lambda v: f"{v:g}")
sb.caption("Keeps domain hits with independent E-value <= cut-off (paper: 1e-4), then re-applies the "
           "overlap rule. Affects tabs 3-4 and uploads in tab 2, not the published paper tables. "
           "Uploaded hit tables can only be filtered more strictly than the threshold they were made "
           "with; 2-column files carry no E-values and are not filtered.")
src = sb.radio("Source", ["Demo: yeast proteome x Pfam DNA-binding HMMs",
                          "Run hmmsearch (upload FASTA + HMM)", "Upload annotation file"])
fmt = sb.selectbox("Annotation format", ["HMMER domtblout", "InterProScan TSV", "protein<TAB>domain"])
program = sb.radio("domtblout produced by", ["hmmsearch", "hmmscan"], horizontal=True)

hits, n_prot, mapping = None, None, None
if src.startswith("Demo"):
    hits, n_prot = hmmsearch_cached((DATA / "raw/yeast_UP000002311.fasta.gz").read_bytes(),
                                    (DATA / "raw/pfam_demo.hmm.gz").read_bytes(), evalue)
    fam = pd.read_csv(DATA / "demo/pfam_families.csv")
    if sb.checkbox("Collapse Pfam to paper domain families", value=True):
        mapping = dict(zip(fam.domain, fam.family))
elif src.startswith("Run"):
    fa = sb.file_uploader("Protein FASTA (.fa, .fa.gz)")
    hm = sb.file_uploader("HMM file(s), HMMER3 format (.hmm, .hmm.gz)")
    if fa and hm:
        hits, n_prot = hmmsearch_cached(fa.getvalue(), hm.getvalue(), evalue)
else:
    up = sb.file_uploader("Annotation file")
    n_prot = sb.number_input("Number of proteins searched (N)", min_value=1, value=20000,
                             help="Size of the peptide set, e.g. genscan predictions per genome. "
                                  "All expectations scale with 1/N; it cannot be read from hit tables.")
    if up:
        hits = parse(up, fmt, evalue, program)
map_up = sb.file_uploader("Optional mapping CSV (domain,family)")
if map_up:
    mp = pd.read_csv(map_up)
    mapping = dict(zip(mp.iloc[:, 0], mp.iloc[:, 1]))
if hits is not None:
    hits = D.apply_mapping(hits, mapping)
    sb.success(f"{len(hits)} domain hits in {hits.protein.nunique()} proteins; N = {n_prot}")

st.title("Quantitative comparison of genome-wide protein domain distributions")
st.caption("After Parikesit, Stadler & Prohaska (2010), German Conference on Bioinformatics, "
           "LNI P-173: 93-102. Domain counts and co-occurrences are computational estimates; "
           "every co-occurrence hit needs inspection for gene-prediction artifacts (fusion "
           "proteins, polycistrons, retro-pseudogenes) and experimental support.")

t1, t2, t3, t4 = st.tabs(["1. Paper data", "2. Domain count comparison",
                          "3. Co-occurrence test", "4. Co-occurrence matrix"])

# ---------------- tab 1 ----------------
with t1:
    t1d = pd.read_csv(DATA / "paper_table1.csv")
    t2d = pd.read_csv(DATA / "paper_table2.csv")
    st.subheader("Table 1: source data (100 random SUPERFAMILY 1.73 models, HMMER 3.0rc1, E <= 1e-4)")
    st.dataframe(t1d, hide_index=True)
    st.markdown("Annotated transcripts give nearly 3x more peptides in human than in gorilla, "
                "yet genscan predictions give similar domain counts across the three apes.")
    st.subheader("Table 2: zinc finger [SF 57667] and winged helix [SF 46785]")
    st.dataframe(t2d, hide_index=True)
    st.caption("sf_multi_strain = 1: several strains/variants in SUPERFAMILY, SF tends to over-count. "
               "Dom = domain copies, Genes = genscan predictions with the domain, co_* = co-occurrences.")
    long = t2d.melt(id_vars=["group", "species"], value_vars=["znf_genes", "wh_genes"],
                    var_name="domain", value_name="genes")
    st.altair_chart(alt.Chart(long).mark_bar().encode(
        y=alt.Y("species", sort=list(t2d.species), title=None),
        x=alt.X("genes", scale=alt.Scale(type="symlog"), title="genscan genes carrying domain"),
        color="domain", yOffset="domain", tooltip=["group", "species", "domain", "genes"]),
        width="stretch")

    st.subheader("Expected co-occurrence, recomputed")
    h = t2d.set_index("species").loc["Homo sapiens"]
    n_gs = int(t1d.set_index("species").loc["Human", "peptides_genscan"])
    pd_dom = D.paper_expectation(h.znf_dom, h.wh_dom, n_gs)
    gene = D.cooccurrence_test(0, h.znf_genes, h.wh_genes, n_gs)
    sf = D.cooccurrence_test(0, 1324, 414, 30712)
    st.dataframe(pd.DataFrame([
        ["Human genscan, domain copies (paper)", f"{h.znf_dom} x {h.wh_dom} / {n_gs}",
         pd_dom["expected"], pd_dom["p_zero"]],
        ["Human genscan, genes (hypergeometric)", f"{h.znf_genes} x {h.wh_genes} / {n_gs}",
         gene["expected"], gene["p_avoid"]],
        ["SUPERFAMILY human transcripts, genes", "1324 x 414 / 30712", sf["expected"], sf["p_avoid"]],
    ], columns=["Estimate", "Formula", "Expected", "P(observed = 0)"]), hide_index=True)
    st.info("The paper's p < 1e-5 uses domain copy counts. ZNF occurs in tandem arrays "
            "(5090 copies in 1048 genes), so copy counts inflate the expectation. At gene level "
            "the single-genome genscan test alone is weak; the SUPERFAMILY gene-level test and the "
            "absence of co-occurrence across all 18 species carry the evidence for avoidance.")

# ---------------- tab 2 ----------------
with t2:
    st.markdown("Figure 1 style comparison: per-domain counts from two peptide sources "
                "(transcripts vs gene predictions) or two species. Points on the identity "
                "line mean consistent estimates.")
    mode = st.radio("Data", ["Paper Table 2: genscan genes vs SUPERFAMILY genes",
                             "Upload two annotation files"], horizontal=True)
    if mode.startswith("Paper"):
        dom = st.radio("Domain", ["znf", "wh"], horizontal=True,
                       format_func={"znf": "zinc finger", "wh": "winged helix"}.get)
        t2d = pd.read_csv(DATA / "paper_table2.csv").set_index("species")
        a, b, xl, yl = t2d[f"{dom}_genes"], t2d[f"{dom}_sf"], "genscan genes", "SUPERFAMILY genes"
    else:
        c1, c2 = st.columns(2)
        f1 = c1.file_uploader("Source X (format from sidebar)")
        f2 = c2.file_uploader("Source Y (format from sidebar)")
        xl, yl = c1.text_input("X label", "source X"), c2.text_input("Y label", "source Y")
        level = st.radio("Count", ["domains", "genes"], horizontal=True)
        a = b = None
        if f1 and f2:
            a = D.counts(D.apply_mapping(parse(f1, fmt, evalue, program), mapping))[level]
            b = D.counts(D.apply_mapping(parse(f2, fmt, evalue, program), mapping))[level]
    if a is not None:
        m, res = D.compare_counts(a, b)
        cols = st.columns(4)
        cols[0].metric("Shared entries", res["n_shared"])
        cols[1].metric("Only in X / only in Y", f"{res['n_x_only']} / {res['n_y_only']}")
        if "slope" in res:
            cols[2].metric("Pearson r (log10)", f"{res['pearson_r']:.3f}")
            cols[3].metric("Spearman rho", f"{res['spearman_rho']:.3f}")
            st.caption(f"OLS log10(Y) = {res['intercept']:.3f} + {res['slope']:.3f} log10(X); "
                       f"median log2(Y/X) = {res['median_log2_ratio']:.2f}; Pearson p = {res['pearson_p']:.2g}. "
                       "Zero counts are excluded from the log-log fit.")
        st.altair_chart(scatter(m, xl, yl, res), width="stretch")
        st.markdown("Largest disagreements, log2((Y+0.5)/(X+0.5))")
        st.dataframe(m.reindex(m.log2_ratio.abs().sort_values(ascending=False).index).head(20))

# ---------------- tab 3 ----------------
with t3:
    if hits is None or hits.empty:
        st.info("Choose an annotation source in the sidebar.")
    else:
        cnt = D.counts(hits)
        st.dataframe(cnt.T)
        opts = sorted(cnt.index)  # stable order, so a cut-off change does not swap the pair
        c1, c2 = st.columns(2)
        a = c1.selectbox("Domain A", opts, index=0, key="dom_a")
        b = c2.selectbox("Domain B", opts, index=min(1, len(opts) - 1), key="dom_b")
        if a == b:
            st.warning("Pick two different domains.")
        else:
            r, both = D.pair_detail(hits, a, b, n_prot)
            cols = st.columns(4)
            cols[0].metric("Proteins with A / B", f"{r['genes_a']} / {r['genes_b']}")
            cols[1].metric("Observed together", r["observed"])
            cols[2].metric("Expected (genes)", f"{r['expected']:.2f}")
            cols[3].metric("log2(obs/exp)", f"{r['log2_obs_exp']:.2f}")
            st.dataframe(pd.DataFrame({
                "test": ["Hypergeometric, gene level: P(X <= obs) avoidance",
                         "Hypergeometric, gene level: P(X >= obs) attraction",
                         "Paper estimate, domain copies: P(0) Poisson"],
                "expected": [r["expected"], r["expected"], r["paper"]["expected"]],
                "p": [r["p_avoid"], r["p_attract"], r["paper"]["p_zero"]]}), hide_index=True)
            st.caption(f"Domain copies: A = {r['domains_a']}, B = {r['domains_b']}, N = {n_prot}. "
                       "The Poisson P(0) row is only meaningful when observed = 0.")
            if len(both):
                st.markdown("Proteins carrying both domains; check each for fusion predictions, "
                            "polycistronic artifacts or pseudogene-derived exons.")
                st.dataframe(both, hide_index=True)

# ---------------- tab 4 ----------------
with t4:
    if hits is None or hits.empty:
        st.info("Choose an annotation source in the sidebar.")
    else:
        opts = sorted(D.counts(hits).index)
        sel = st.multiselect("Domains", opts, default=opts[:8], key="mat_sel")
        if len(sel) >= 2:
            mat = D.cooccurrence_matrix(hits, sel, n_prot)
            sym = pd.concat([mat, mat.rename(columns={"a": "b", "b": "a"})])
            st.altair_chart((alt.Chart(sym).mark_rect().encode(
                x=alt.X("a", sort=sel, title=None), y=alt.Y("b", sort=sel, title=None),
                color=alt.Color("log2_obs_exp", scale=alt.Scale(scheme="redblue", domainMid=0, reverse=True),
                                title="log2(obs/exp)"),
                tooltip=["a", "b", "observed", alt.Tooltip("expected", format=".2f"),
                         alt.Tooltip("p", format=".2g"), alt.Tooltip("q_bh", format=".2g")])
                + alt.Chart(sym).mark_text(fontSize=11).encode(
                    x=alt.X("a", sort=sel), y=alt.Y("b", sort=sel), text="observed"))
                .properties(height=60 * len(sel) + 40), width="stretch")
            st.caption("Cell text = observed co-occurring proteins. Red = attraction, blue = avoidance. "
                       "Two-sided hypergeometric p, Benjamini-Hochberg q over the pairs shown. "
                       "Avoidance needs large expected counts to reach significance; in a single small "
                       "genome, zero observations are usually not informative on their own.")
            st.dataframe(mat[["a", "b", "observed", "expected", "log2_obs_exp", "p_avoid",
                              "p_attract", "p", "q_bh"]].sort_values("p"), hide_index=True)
