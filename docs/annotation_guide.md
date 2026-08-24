# Cluster Annotation Guide

Clustering gives you numbered groups. Turning those numbers into cell types is
the step where most of the error in a single-cell paper gets introduced, and
it is the step that is hardest to check afterwards.

This guide describes annotation by **triangulation**: run several independent
methods, put their answers side by side in one table, and have a human make the
call from a dotplot. No method is treated as the answer. Where methods
disagree, the disagreement is recorded rather than resolved silently — because
the disagreement is usually informative about the cluster.

---

## 1. Why not just use an automated annotator

Automated methods are fast, reproducible, and biased in ways that do not
announce themselves.

- **Reference mapping** (CellTypist, Azimuth, scANVI) returns a label for every
  cell with a confidence score, including for cell types absent from the
  reference. A cluster that is genuinely a doublet, a technical artifact, or a
  cell type the reference never saw still gets a confident-looking name.
- **Reference granularity leaks into your results.** Azimuth's `predicted.celltype.l2`
  systematically over-calls CD4 central memory: in one PBMC cohort it put CD4
  TCM at ~36% of cells against 10–15% in matched CyTOF and matched scRNA-seq
  data from other cohorts. Adopt those labels and the distortion becomes a
  finding.
- **Marker panels** are interpretable but only see what you put in them. A
  cluster of mast cells in a PBMC dataset is invisible to a panel without
  tryptase genes, and will be assigned to whatever it is nearest.

Each of these fails differently. That is the argument for running several.

---

## 2. The evidence sources

| Source | What it is good for | How it fails |
|---|---|---|
| CellTypist `Immune_All_Low` | Fine-grained immune subsets | Confident on absent types; misses non-immune cells |
| CellTypist `Immune_All_High` | Lineage-level sanity check | Too coarse to resolve subsets |
| Same two models on **full gene set** | Usually the more confident run | Slower; can disagree with the HVG run |
| Marker panel z-scores | Interpretable, hypothesis-driven | Blind to anything not in the panel |
| Unbiased Wilcoxon DE | Panel-free; catches technical splits | Unlabelled; ribosomal/mito genes dominate |
| Original authors' labels | Free comparison, catches metadata errors | Carries the original paper's biases |
| **The dotplot** | **Decisive** | Requires you to look at it |

Running CellTypist **twice** — once on the HVG matrix, once on the full gene
set — is not redundancy. When the two runs flip a cluster between
central-memory and effector-memory, that flip is usually telling you the
cluster contains both, and is a better description of the cluster than either
label alone.

---

## 3. Running it

Four steps, after clustering. See `pipeline/README.md` for the full sequence.

```bash
# 1e — CellTypist, twice: HVG and full-gene, two models each = 4 label columns
python pipeline/step1e_celltypist.py --file work/step1b_pca.h5ad \
    --out-prefix work/data/celltypist --mode hvg
python pipeline/step1e_celltypist.py --file work/step1b_pca.h5ad \
    --out-prefix work/data/celltypist --mode fullgene

# 1g — marker z-scores, raw means, unbiased DE, and the dotplot
python pipeline/step1g_marker_evidence.py --file work/step1b_pca.h5ad \
    --panel pbmc --out-dir figures

# 1h — assemble everything into one CSV with the label columns left blank
python pipeline/step1h_build_annotation_template.py --file work/step1b_pca.h5ad \
    --celltypist-dir work/data --evidence-dir figures \
    --out work/data/cluster_annotation_template.csv

# --- you now open the dotplot and fill in manual_label_final ---

# 1i — apply the filled-in labels back onto the cells
python pipeline/step1i_apply_labels.py --file work/step1b_pca.h5ad \
    --template work/data/cluster_annotation_template.csv --out work/annotated.h5ad
```

Step 1h writes a 24-column table: cluster, size, four CellTypist calls with
confidences and top-3 breakdowns, top panel markers by z-score, top-10 unbiased
DE genes, and three **blank** columns — `user_label`, `manual_label_final`,
`notes`. The blanks are the point. The table gathers evidence; it does not
decide.

### Using the authors' labels as a column

If the GEO submission ships labels (a Seurat `.rds`, a metadata CSV), carry
them into the template with `--author-labels`:

```bash
python pipeline/step1h_build_annotation_template.py ... \
    --author-labels work/data/authors_label_by_cluster.csv   # cluster,author_label
```

Carry them, compare against them, and do not adopt them. They come from a
different annotation philosophy at a different granularity, and importing them
imports that distortion into your dataset. What they are genuinely good for is
catching metadata errors: if the authors' `diagnosis` field disagrees with the
GEO-derived `disease_state` for some samples, one of the two is wrong, and it is
worth finding out which before running anything downstream.

---

## 4. Interpretation traps

### The z-score trap

`cluster_means_zscore.csv` is z-scored **across clusters**. A high value means
"higher here than in the other clusters" — not "expressed here". For sparse
genes those come apart badly:

| Gene | Cluster mean (raw) | z-score | Reading |
|---|---|---|---|
| `PDCD1` | 0.10 | +2.8 | **Not** exhaustion. It is near-zero everywhere. |
| `IL10` | 0.01 | +3.0 | **Not** a Tr1 cell. Noise. |
| `IL17F` | 0.00 | +5.0 | Pure noise. Not Tc17. |

This is why step 1g writes `cluster_means_raw.csv` next to the z-scores.
**Check the raw mean before acting on any high z-score.** A marker at raw 0.01
is not a marker.

### Low-complexity clusters

Clusters with markedly lower library complexity than the rest of the dataset
(e.g. median 944 genes / 1,280 counts against a typical ~1,300 / ~2,100) and no
lineage-defining positive marker are a recurring problem. Their unbiased DE is
characteristically long, intron-rich, nuclear-enriched transcripts: `MALAT1`,
`ARHGAP15`, `FOXP1`, `CAMK4`, `INPP4B`, `MAML2`, `MBNL1`, `PDE3B`.

Two readings are usually both live: genuine resting cells with low RNA content,
or degraded/nuclear-heavy cells that Leiden separated on a technical axis.
There is often no way to settle it from the data alone. Decide explicitly,
write the reasoning into `notes`, and check how much the decision moves your
results — a low-complexity cluster at 9% of the dataset shifts every downstream
frequency.

### Doublet clusters

Scrublet (step 1a) works per cell and misses **clusters** of doublets. Those
show up here as clusters co-expressing two lineages:

- platelet + monocyte: `PPBP`, `PF4` alongside `CD14`, `VCAN`
- B + monocyte: `MS4A1`, `CD79A` alongside `CD14`, `CD68`
- T + monocyte: `CD3D` alongside `LYZ`, `S100A8`

Mark these `doublet` in `notes`. Step 1i will flag them in
`obs['cell_type_flag']`, and `--drop-flagged` removes them.

The reverse also happens: Scrublet's automatic threshold is unreliable on small
or low-complexity samples and can flag an implausible fraction of cells. This is
why the compiler flags rather than drops by default. Check the distribution of
`obs['doublet_score']` per sample before trusting `predicted_doublet`.

### Mitochondria-heavy, mixed-lineage clusters

A cluster whose top DE genes are mostly `MT-` and whose markers span unrelated
lineages is usually dying cells rather than a cell type. Mark it `artifact`.

### Cell types absent from your schema

If you are matching a published paper's cell-type list, you will find real
populations that the list does not contain — MAIT cells (`SLC4A10`, `TRAV1-2`,
`KLRB1`), iNKT cells (`TRAJ18`), mast cells/basophils (`TPSAB1`, `TPSB2`,
`CPA3`, `MS4A2`, `HDC`). These signatures are categorical; do not force them
into the nearest schema category. Add the label and note the deviation.

---

## 5. Marker panel — PBMC

The built-in panel (`--panel pbmc`, or `anndata_compiler.PBMC_MARKER_PANEL`).
Human symbols.

### T cells (all: `CD3D`, `CD3E`, `CD3G`, `TRAC`)

| Cell type | Positive | Negative | Note |
|---|---|---|---|
| Naive CD4 T | `CD4`, `CCR7`, `SELL`, `LEF1`, `TCF7` | `IL7R` low, granzymes off | DE is ribosomal-dominated — the classic naive signature |
| Memory CD4 T | `CD4`, `IL7R`, `S100A4`, `IL32` | `CCR7`/`SELL` off | `IL7R` (CD127) is the key |
| Th1 | `CXCR3`, `TBX21`, `CCR5`, `IFNG` | `CCR6`, `RORC` off | `TBX21` is sparse in scRNA; `CXCR3` more reliable |
| Th17 | `RORC`, `CCR6`, `IL17A`, `IL17F`, `IL23R` | `CXCR3` off | `IL17A`+`IL17F` together is far more specific than `IL17A` alone |
| Treg | `FOXP3`, `IL2RA`, `CTLA4`, `TIGIT` | `IL7R` low | `FOXP3` is definitive |
| Naive CD8 T | `CD8A`, `CD8B`, `CCR7`, `SELL`, `LEF1`, `TCF7` | granzymes off | `CD8B` is naive-enriched vs memory |
| Effector memory CD8 | `CD8A`, `GZMK`, `CCL5` | `CCR7` off | `GZMK`-high is the EM marker |
| TEMRA / terminal effector CD8 | `GZMH`, `FGFBP2`, `KLRG1`, `EOMES`, `CX3CR1` | `CCR7` off | `FGFBP2`+`GZMH`+`KLRG1` reads terminal, not plain EM |
| Exhausted CD8 | `PDCD1`, `LAG3`, `HAVCR2`, `TOX` | `GZMB` lower than effector | Needs **multiple** checkpoints co-expressed at real raw values |
| MAIT | `SLC4A10`, `TRAV1-2`, `KLRB1`, `IL18R1` | — | Categorical signature |
| γδ T | `TRDC`, `TRGC1`, `XCL1`, `XCL2` | `CD4`, `CD8A` off/low | `TRDC` is the decider. CellTypist routinely calls these NK |

### NK / ILC (`CD3D`, `CD3E` negative)

| Cell type | Positive | Note |
|---|---|---|
| CD56bright NK | `NCAM1`, `XCL1`, `XCL2`, `SELL`, `CD27` | `FCGR3A` (CD16) **low** — this is what separates it |
| Mature CD16+ NK | `FCGR3A`, `FGFBP2`, `PRF1`, `GNLY`, `KLRF1` | `NCAM1` dim |
| Activated NK | above plus `GZMB`, `NKG7`, `IFNG` | Whether this is a separate call or a state of mature NK is a schema decision — make it once |

### B cells (`CD19`, `MS4A1`, `CD79A`, `CD79B`; except plasma)

| Cell type | Positive | Negative | Note |
|---|---|---|---|
| Naive B | `IGHD`, `IGHM`, `TCL1A`, `FCER2`, `CR2` | `CD27` off | `IGHD`+`IGHM` = unswitched |
| Memory B | `CD27`, `AIM2`, `BANK1`, `TNFRSF13B` | `IGHD` off | Class-switched |
| Breg | `CD24` hi, `CD38` hi, `IL10` | — | Very hard in scRNA — requires real `IL10`. Check the raw mean, not the z-score |
| Plasma | `JCHAIN`, `XBP1`, `MZB1`, `PRDM1`, `SDC1`, `TNFRSF17` | `MS4A1` off | `MS4A1` is downregulated in plasma cells |

### Myeloid

| Cell type | Positive | Negative | Note |
|---|---|---|---|
| Classical monocytes | `CD14`, `VCAN`, `S100A8`, `S100A9`, `S100A12`, `LYZ`, `FCN1` | `FCGR3A` low | Textbook S100A signature |
| Non-classical monocytes | `FCGR3A`, `CDKN1C`, `MS4A7`, `LST1`, `HES4`, `CX3CR1` | `CD14` low | `CDKN1C` is very specific |
| cDC1 | `CLEC9A`, `XCR1`, `BATF3`, `IRF8` | `CD14` off | `BATF3` alone is **not** enough — it is also high in non-classical monocytes. Require `CLEC9A`+`XCR1` co-expression |
| cDC2 | `CD1C`, `CLEC10A`, `FCER1A`, HLA class II | `CD14` off | |
| pDC | `TCF4`, `IRF7`, `LILRA4`, `CLEC4C`, `PLD4`, `BCL11A` | `CD14` off | `TCF4`+`IRF7` is definitive |
| Macrophage | `CD68`, `CD163`, `MARCO`, `MSR1`, `MRC1`, `MERTK` | — | Uncommon in PBMC — suspect M-MDSC or tissue contamination |

### Other

| Cell type | Positive |
|---|---|
| Platelets | `PPBP`, `PF4`, `GP1BB`, `TUBB1`, `NRGN`, `ITGA2B` |
| Mast / basophil | `TPSAB1`, `TPSB2`, `CPA3`, `MS4A2`, `HDC` |
| HSPC | `CD34`, `SPINK2`, `PRSS57` |
| Proliferating | `MKI67`, `TOP2A`, `STMN1`, `PCLAF` |

---

## 6. Marker panel — microglia and CNS myeloid

The PBMC panel does not cover brain myeloid cells, and microglia should not be
annotated with a hand-written panel. Use the maintained atlas:

### → **https://ldhawka.github.io/microglia-subtype-markers/**

**9 umbrella families over 23 sub-states**, parallel human (HGNC) and mouse
(MGI) symbols, up- **and** down-regulated markers, curated protein/surface
markers, gene sets **denoised** of dissociation-response artifacts, and every
sub-state linked to its primary reference DOI. The page is searchable and
sortable; [`LITERATURE_VALIDATION.md`](https://github.com/ldhawka/microglia-subtype-markers/blob/main/LITERATURE_VALIDATION.md)
carries the validation.

Load it directly:

```bash
# sub-state level (23 panels)
python pipeline/step1g_marker_evidence.py --file work/step1b_pca.h5ad --panel microglia

# umbrella-family level (9 panels) — usually the level you should call at
python pipeline/step1g_marker_evidence.py --file work/step1b_pca.h5ad \
    --panel microglia --level umbrella

# mouse
python pipeline/step1g_marker_evidence.py --file work/step1b_pca.h5ad \
    --panel microglia --species mouse
```

```python
from anndata_compiler import load_microglia_panel
panel = load_microglia_panel(species='human', level='umbrella', denoised=True)
```

**The key point from that atlas: disease-associated states are a continuum, not
discrete clusters.** Validation against annotated microglia in three ALS
datasets found the disease-associated signatures overlap heavily (max Jaccard
0.42). Scoring per sub-state and rolling up to the umbrella family reaches ~73%
family-level accuracy; forcing a one-of-23 call does not. Prefer
`--level umbrella` for the label and keep the sub-state scores as a continuous
readout.

Use `denoised=True` (the default) for anything that went through enzymatic
dissociation — otherwise the dissociation-response genes will make every
activated-looking cluster look more activated.

---

## 7. Granularity, and why it matters beyond this dataset

Annotation granularity is not a presentation choice. It changes results.

If you are comparing cohorts — across datasets, across platforms, or against a
published result — the vocabularies must match. The arithmetic is unforgiving:
24 cell types give 276 possible cell-type pairs; 13 give 78. A cohort annotated
at 13 lumped categories cannot reproduce a pairwise result found at 24, because
the two top-ranked pairs may have been merged into a single bucket.

Practical consequences:

1. **Fix the vocabulary before annotating**, not after seeing results. Decide
   once whether TEMRA is its own label or a state of effector memory, whether
   activated NK is separate from mature NK, and apply it to every cohort.
2. **Re-annotating a cohort after seeing it fail is not defensible.** If a
   granularity mismatch is the suspected explanation for a negative result, the
   honest routes are to demonstrate the claim positively in a *different* cohort
   at full granularity, or to pre-register the re-annotation as a sensitivity
   analysis and report **both** label sets.
3. **Reference-mapped labels are not neutral.** See §1.

Relevant literature: Hao et al., *Cell* 2021 (PMID 34062119) — Azimuth's three
levels of 8/30/57 types, where CD4 memory only appears at level 2; Abdelaal
et al., *Genome Biol* 2019 (PMID 31500660); omnideconv, *Genome Biol* 2026
(PMID 41582216).

---

## 8. Downstream caveat: cell-type frequencies are compositional

Once labels are applied, per-sample cell-type frequencies are the usual next
output. They are **compositional** — forced to sum to 100% per sample — and that
closure manufactures correlations that are not biology. A variable dominant
population squeezes the minor populations together, inducing spurious *positive*
correlations among them.

This matters if you correlate cell-type frequencies with each other:

- Stability selection and leave-one-out **do not catch it**. Closure is
  deterministic and reproducible, which is exactly why resampling-based checks
  pass. Only a log-ratio transform catches it.
- Test with a **centered log-ratio (CLR)** transform: per sample, `log(count + 0.5)`
  per cell type minus the per-sample mean across cell types. The `+0.5` Haldane
  pseudocount handles zeros. Report which findings survive.
- Empirically, the strongest and FDR-significant pairs — usually between
  abundant populations — are the most closure-robust. The fragile layer is
  network *topology*: hub claims and the long tail of low-abundance positive
  pairs. In one PBMC analysis ~68% of frequency-selected edges survived CLR;
  restricted to a granulocyte-inclusive convention it fell to ~37%.
- CLR is not free. Its per-sample geometric-mean reference is platform-specific,
  so it can break comparability when transferring between platforms (e.g. CyTOF
  to scRNA-seq).

The short version: if a frequency-correlation result is load-bearing, run the
CLR version and report it as a sensitivity analysis.

---

## 9. Checklist

- [ ] Clustered on the **Harmony** embedding, not raw PCA
- [ ] One resolution chosen for annotation (not annotating several)
- [ ] CellTypist run **four** ways (2 models × HVG/full-gene)
- [ ] CellTypist given log1p(CP10K) — **not** scaled data, not raw counts
- [ ] Marker z-scores computed on the **full gene set**, not just HVGs
- [ ] Raw means checked for every marker driving a call
- [ ] Unbiased DE inspected for ribosomal / mito / long-nuclear-transcript signatures
- [ ] Doublet and artifact clusters flagged in `notes`
- [ ] Label strings exactly consistent so sister clusters merge
- [ ] Deviations from any target schema written down
- [ ] Vocabulary matched across cohorts **before** looking at results
- [ ] Template CSV kept with the analysis — it is the record of why each call was made
