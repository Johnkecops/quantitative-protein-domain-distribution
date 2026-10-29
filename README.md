# Protein Domain Distribution (Streamlit)

Interactive re-implementation of the analysis in Parikesit, Stadler & Prohaska (2010):
compare genome-wide protein domain counts between peptide sources or species, and test
whether two domain types co-occur in the same protein more or less often than chance.

## Run

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/streamlit run app.py
.venv/bin/python -m pytest tests     # 5 self-checks, includes the paper's 11.7 and 17.8
```

## Tabs

1. **Paper data**: Tables 1 and 2 transcribed from the paper, with the expected
   ZNF / winged-helix co-occurrence recomputed three ways (domain copies as in the paper,
   genes with a hypergeometric test, SUPERFAMILY genes).
2. **Domain count comparison**: Fig. 1 style log-log scatter, identity line, OLS on log10
   counts, Pearson and Spearman. Demo uses Table 2 (genscan genes vs SUPERFAMILY genes);
   or upload two annotation files.
3. **Co-occurrence test**: one domain pair, observed vs expected, both hypergeometric
   tails, the paper's Poisson estimate, and the list of proteins carrying both (inspect
   these for fusion predictions and other gene-prediction artifacts).
4. **Co-occurrence matrix**: all pairs among selected domains, log2(obs/exp) heatmap,
   two-sided hypergeometric p, Benjamini-Hochberg q.

## Inputs (sidebar)

- **Demo**: *S. cerevisiae* S288C reference proteome (UniProt UP000002311, 6066 proteins,
  downloaded 2026-10-29) searched in-process with pyhmmer against 12 Pfam HMMs
  (InterPro API, 2026-10-29), collapsed to the paper's domain families via
  `data/demo/pfam_families.csv`.
- **Run hmmsearch**: upload a protein FASTA and HMMER3 HMM file(s); gz accepted.
- **Upload annotation**: HMMER `--domtblout` (hmmsearch or hmmscan), InterProScan TSV,
  or a 2-column `protein<TAB>domain` file. Enter N, the number of proteins searched.

Paper rules applied to all hit tables: per-domain (i-)E-value <= 1e-4 (adjustable) and,
for overlapping hits on one protein, only the best-scoring hit is kept.

Exact tested versions: `pip install -r requirements.lock` (Python 3.11, streamlit 1.66.1,
pandas 3.1.0, scipy 1.17.1, pyhmmer 0.12.3).

Demo data SHA-256 (verify with `shasum -a 256 data/raw/*.gz`):

```
8ed21b6c6ce28bc0b74eac295aab65ec2136286616769757b39758427642f768  data/raw/yeast_UP000002311.fasta.gz
d2b3a4c956162447d14e163e1e366528219100eea3691f8b3ae9e98412e116ad  data/raw/pfam_demo.hmm.gz
```

## Upload limits (for shared deployments)

Security: a focused `/security-audit` review (2026-10-29, source-only, not the full 6-phase audit) was run and its
upload-handling risks were patched: gzip-bomb cap, hmmsearch job limits, bounded cache, CSV formula escaping.

Set in `domains.py`: decompressed input <= 500 MB (gzip-bomb guard), in-app hmmsearch
<= 150,000 sequences and <= 200 HMMs, pyhmmer limited to 4 CPU threads. Search results
are cached for at most 8 inputs for 1 hour. Protein and domain IDs that start with
`= + - @` are prefixed with `'` so CSV exports cannot run as spreadsheet formulas.
Larger jobs (e.g. full Pfam-A against a vertebrate proteome) belong offline; upload the
resulting `--domtblout`. `.streamlit/config.toml` keeps XSRF protection on and the
200 MB upload cap explicit.

## Deviations from the paper and caveats

- The paper used SUPERFAMILY 1.73 (SCOP) HMMs. SUPERFAMILY models need a licence, so the
  demo uses Pfam families as proxies: zf-C2H2 (SF 57667), Forkhead / HSF_DNA-bind /
  RFX_DNA_binding / WHD_E2F_TDP (winged helix, SF 46785), bZIP_1 / bZIP_2 (SF 57979),
  HMG_box / HMG_box_2 (SF 47095). Zn_clus + Fungal_trans(_2) are a fungal positive control
  for co-occurrence. Pfam covers far less of the winged-helix superfamily than SCOP does.
- genscan is not run in the app (academic licence, chromosome chunking). Bring your own
  gene predictions as FASTA.
- The paper's p < 1e-5 for human uses domain copy counts (5090 x 274 / 118894 = 11.7).
  ZNF occurs in tandem arrays, so at gene level the expectation is 1048 x 256 / 118894 = 2.3
  and P(0) about 0.10. The app shows both. Avoidance in one small genome (e.g. yeast,
  expected < 1) is not detectable; the evidence comes from large counts and from consistent
  absence across many lineages.
- All results are computational hypotheses. Co-occurrence calls need manual curation and,
  where it matters, transcript or protein-level validation (RT-PCR across the predicted
  domain junction, mass spectrometry, or orthologous gene models in related species).

## Files

```
app.py                     Streamlit UI
domains.py                 parsers, overlap rule, counts, statistics, pyhmmer wrapper
data/paper_table1.csv      Table 1 of the paper
data/paper_table2.csv      Table 2 of the paper
data/demo/pfam_families.csv
data/raw/yeast_UP000002311.fasta.gz
data/raw/pfam_demo.hmm.gz  12 concatenated Pfam HMMs
tests/test_domains.py
```

## References

- Parikesit AA, Stadler PF, Prohaska SJ (2010). Quantitative comparison of genomic-wide
  protein domain distributions. German Conference on Bioinformatics 2010, Lecture Notes in
  Informatics P-173: 93-102.
- Eddy SR (2011). Accelerated profile HMM searches. PLoS Comput Biol 7:e1002195.
  doi:10.1371/journal.pcbi.1002195
- Larralde M, Zeller G (2023). PyHMMER: a Python library binding to HMMER for efficient
  sequence analysis. Bioinformatics 39:btad214. doi:10.1093/bioinformatics/btad214

AI Assistance Disclaimer: This codebase was developed with the assistance of Claude Code. While the AI provided code generation, debugging, and structural support, the human developer maintains full responsibility for reviewing, testing, and maintaining all content and functionality
