# Literature Review — Cell–Cell Communication & Spatial Signaling Methods

Structured summaries of the papers in `alarmist/literature/`. All are cell–cell communication (CCC) / ligand–receptor (LR) inference methods — the direct method family around **ALARMIST** (Bayesian tensor factorization of spatial LR interactions into latent multicellular **motifs/programs**, characterized by cell-type senders/receivers and downstream gene impact).

PDFs are named `Method_Author_Year_Journal.pdf`.

---

## Quick reference

| # | Method | Venue (year) | Spatial? | Unit of inference | Core model | Factorization of communication? | Role vs ALARMIST |
|---|--------|--------------|:--------:|-------------------|-----------|:-------------------------------:|------------------|
| 1 | **COMPOTES** | bioRxiv (2025) | ✅ spot | per-spot LRI → cohort programs | diffusion+competition LR score → **non-negative CP/PARAFAC (≈NMF)** | ✅ | **Closest sibling / primary competitor** |
| 2 | **Tensor-cell2cell** | Nat Commun (2022) | ❌ (context axis) | context×LR×sender×receiver | **4D non-negative CP/PARAFAC** | ✅ | Direct conceptual predecessor / baseline |
| 3 | **LIANA+** | Nat Cell Biol (2024) | ✅ + ❌ | all-in-one | consensus LR + **NMF / Tensor-cell2cell / MOFA+** + causal nets | ✅ | Closest *published* relative / main baseline |
| 4 | **CytoSignal** | Nat Genet (2026) | ✅ cell | per-cell, per-LR | LRscore `L×R`, contact vs diffusion, permutation | ❌ | Per-LR spatial scorer; PLA-validated; input/baseline |
| 5 | **SpatialDM** | Nat Commun (2023) | ✅ spot/cell | per-LR + per-spot | **bivariate Moran's R** + analytical null | ❌ (patterns via SpatialDE) | Input scorer / per-pair baseline |
| 6 | **COMMOT** | Nat Methods (2023) | ✅ spot/cell | per-LR, spot→spot | **collective optimal transport** | ❌ | Competition-aware input scorer |
| 7 | **stLearn** | Nat Commun (2023) | ✅ spot/cell | per-LR + trajectory + imputation | **SCTP** 2-level permutation; PSTS; stSME | ❌ | Per-LR spatial scorer / baseline |
| 8 | **NICHES** | Bioinformatics (2023) | ✅ (NN) | per cell-pair / per-spot niche | bilinear `L×R`, then cluster | ❌ | Lightweight input scorer / predecessor |
| 9 | **CellChat v2** | Nat Protoc (2024) | ❌ | cell-group → cell-group | **mass-action + Hill** + permutation; NMF patterns | ❌ (group-level NMF) | Non-spatial aggregator; DB & baseline |

**The three tiers below:** (A) methods that **factorize communication into latent programs** — ALARMIST's own class and its true peers; (B) **spatial per-LR-pair scorers** — the front-end signal ALARMIST factorizes, and per-interaction baselines; (C) a **non-spatial cell-type aggregator** — an upstream scorer / database.

---

# Tier A — Communication-factorization / program-discovery methods (ALARMIST's class)

## 1. COMPOTES — Herpin et al., bioRxiv 2025 (Owkin / MOSAIC)
**File:** `COMPOTES_Herpin_2025_bioRxiv.pdf` · *Deciphering Cellular Ecosystems Driving Tumor Progression and Immune Escape from Spatial Transcriptomics and Single-Cell with COMPOTES*

**One-liner:** Computes spatially-resolved, diffusion- and competition-aware LR interaction scores per Visium spot, then applies non-negative CP/PARAFAC (≈NMF) across an entire multi-sample cohort to extract recurrent multicellular communication **programs**.

- **Background:** TME CCC is driven by LR interactions governing survival, proliferation, immune evasion. Prior scorers use curated DBs (CellPhoneDB, CellChatDB) via mean/product of L,R expression; multi-sample methods increasingly use matrix/tensor factorization to find recurrent programs.
- **Motivation:** No prior method combined **local, spatially-aware** LRI modeling **with** cohort-scale program extraction. Single-sample spatial tools (COMMOT, BATCOM) model diffusion/competition but find no cross-sample programs; cohort program methods (Tensor-cell2cell, MOFAcell, DIALOGUE) use non-spatial inputs → false positives from non-proximal pairs.
- **Data:** **MOSAIC** — 10x Visium V2 CytAssist FFPE + matched snRNA-seq, bulk RNA, WES, H&E, clinical. 3.5M spots, 1,084 samples, 7 indications (NSCLC 248, ovarian 240, DLBCL 154, MIBC 146, breast 118, GBM 112, meso 66). Deep dive: **146 MIBC** samples (~500k spots). LR prior: CellChatDB. Plus pseudo-Visium simulations.
- **Method:** (1) *Per-spot LR scoring* — diffuse ligand over neighbor rings by CellChatDB range annotation (short/medium/long), weight `w_dist = d/(2πr)`; **competition** step normalizes each ligand by the sum of all ligands binding the same receptor; `LRI = (diffused, competed ligand) × (receptor in the spot)`; complexes use **min over subunits**. (2) *Cohort factorization* — stack all spots × all LRIs into `X`, apply **non-negative CP** via TensorLy HALS (2D case = NMF); each of K programs = LR-loading vector × spot-loading vector. Downstream: cell-type sender/receiver maps from matched snRNA pseudobulk; associations via decoupleR ULM (TF/signatures), PROGENy, InferCNVpy CNV, linear mixed-effects, Mann-Whitney+BH.
- **Model:** `X ≈ Σ_r λ_r · w_r ⊗ h_r`, `w,h ≥ 0`, minimize Frobenius residual. Program similarity = cosine of `h_r` (match threshold 0.8). No generative count model / no permutation null.
- **Output:** K programs, each = (weighted LR-pair list `h_r`) + (spatial activity map `w_r`); per-program cell-type sender→receiver matrix, sample-level activity, LR-mode counts (contact/secreted/ECM), complexity class, and associations to subtype/clinical/mutation/CNV/TF/pathway.
- **Results:** On MIBC, **K=45** programs in 3 classes (immune / ECM / growth-factor); no center batch effect; 19 programs shared across 5 cancers (GBM, DLBCL most singular). **Program 15** (6/7 cancers): "mounted-yet-inefficient" anti-tumor immunity (IL16-CD4, CCL5-CCR5, NECTIN2-TIGIT), MoMac receivers, tumor–stroma interface, dual IFN-γ + exhaustion signature, Basal subtype. **Program 7** (MIBC-specific): high-proliferation "engine" (Semaphorin-Plexin, Ephrin-Eph, Netrin-UNC5, WNT5A), autocrine malignant loops, T2/Luminal/tumor-core, +NACT response (p=0.03), **KMT2D LoF** (p=0.0091, q=0.047), CNV-driven.
- **Benchmark:** vs **Tensor-cell2cell** + two ablations (no-competition, no-diffusion) + permutation control. On simulations with deliberately **non-proximal** interacting pairs, COMPOTES recovered only the true spatially-proximal programs and **rejected the false pairs**, while Tensor-cell2cell included them; on a 50-slide sim it recovered true interactions with far fewer false positives. No experimental gold standard (synthetic + biological plausibility).
- **Use cases:** Recurrent spatially-grounded CCC programs across large multi-indication cohorts; pan-cancer vs indication-specific programs; linking programs to clinical/mutation/CNV/TF/pathway; immunosuppression axes & vulnerabilities. Visium spot-res; GPU-memory-bound (~16GB at 146 samples).
- **Relevance to ALARMIST:** **The single most conceptually related paper — a near-sibling / concurrent competitor.** Identical recipe: score spatial LRIs per spot across a cohort → (spots × LR-pairs) matrix → non-negative factorization → latent programs (ALARMIST's "motifs") → cell-type sender/receiver + covariate characterization. Contrasts to exploit: COMPOTES uses **deterministic** CP/HALS vs ALARMIST's **Bayesian** tensor factorization (uncertainty, ON/OFF motif states); COMPOTES's LR score has no statistical count model, no permutation null, and **no built-in downstream target-gene impact** (it cites this as future work) — both are ALARMIST differentiators. Its diffusion/competition scorer could be an alternative input; its Tensor-cell2cell false-positive comparison and cosine cross-indication matching are ready templates.

## 2. Tensor-cell2cell — Armingol et al., Nature Communications 2022
**File:** `Tensor-cell2cell_Armingol_2022_NatCommun.pdf` · *Context-aware deconvolution of cell–cell communication with Tensor-cell2cell* (13:3665)

**One-liner:** Unsupervised non-negative CANDECOMP/PARAFAC decomposition that factorizes a 4D (context × LR-pair × sender × receiver) communication tensor into a few latent context-driven communication **modules**.

- **Background:** CCC is shaped by "cellular context" (severity, life stage, tissue, subject) treated as a non-binary variable. Tensor Component Analysis preserves cross-context correlation structure better than matrix methods (PCA/NMF/UMAP) that collapse it.
- **Motivation:** Prior tools ignore context (CellPhoneDB, NicheNet, scTensor) or handle only pairwise/two-context comparisons (CellChat, NATMI, iTALK); per-context-then-stitch loses correlation structure, power, and scalability (n-choose-2). Gap: one scalable unsupervised analysis jointly modeling many contexts.
- **Data:** Non-spatial sc/snRNA-seq. Simulation (3 cell types, 300 LRs, 12 time points, 4 embedded temporal patterns); **COVID-19 BALF** (Liao, 12 samples, 65,813 cells → tensor 12×189×6×6); **COVID PBMC** (Ren, 1.46M cells → 60×1639×6×6, efficiency test); **Autism PFC** snRNA (Velmeshev, 104,559 nuclei, 23 subjects → 23×749×16×16). LR DB: CellChat 2,005 pairs. Built in `cell2cell` on TensorLy.
- **Method:** (1) per context, compute communication score (default **mean L×mean R**, min-subunit for complexes) for each LR × sender × receiver; (2) stack into 3D per-context tensor; (3) stack contexts → 4D tensor; apply non-negative TCA. Rank R by reconstruction-error elbow. Downstream: rank LRs by loadings, pre-ranked GSEA over 22 KEGG "LR pathway sets"; Spearman of context loadings vs severity; group t-tests; Ward clustering; sender⊗receiver outer-product networks with **Gini** imbalance. Can ingest external LIANA edge-lists.
- **Model:** `χ ≈ Σ_r c_r ⊗ p_r ⊗ s_r ⊗ t_r`, all loadings ≥ 0, minimize Frobenius residual by ALS. Rank by normalized reconstruction error (best: 4 sim, 10 BALF, 6 ASD). Decomposition consistency via **CorrIndex** (permutation/scale-invariant). Gini on edge weights.
- **Output:** Per factor, four loading vectors (context dynamics, LR pairs, senders, receivers — directional); reconstruction-elbow plots; factor cell-cell networks + Gini; GSEA per factor; Spearman severity; loading heatmaps; RF disease classifier on context loadings.
- **Results:** Simulation: 4-factor decomposition **exactly recovered** the 4 embedded temporal patterns + assigned LRs and cells (Jaccard/Pearson vs truth); noise-robust. Scoring-robust (mean 1−CorrIndex 0.82 across CellChat/CellPhoneDB/NATMI/SingleCellSignalR); preprocessing-robust (0.86). BALF: factor 8 most severity-correlated (Spearman 0.92, macrophage CCL→CCR); factor 9 homeostatic (Gini 0.09, MIF-CD74). ASD: factors 3/4 separate ASD vs control (interneurons; axon-guidance/ECM/ERBB down).
- **Benchmark:** vs **CellChat** (only prior multi-context tool). Faster + less memory on PBMC (GPU V100). Downstream RF (ROC-AUC) predicting COVID status from loadings: **~10–20% higher** than CellChat's manifold embeddings. Consistency vs 5 scorers via CorrIndex. Table 1 positions it as the only tool supporting **unlimited simultaneous contexts**. Ground-truth-limited (hence synthetic tensor).
- **Use cases:** Deciphering CCC across ordered/unordered context sets (severity, time, subjects, tissues) in one unsupervised run; ranking LRs and senders/receivers per module; linking modules to phenotype (correlation, tests, GSEA, classifiers). Works on sc/sn/bulk; can wrap any external scorer; spatial only by pre-defining niches as contexts.
- **Relevance to ALARMIST:** **Closest conceptual predecessor.** Same idea — factorize LR communication into latent modules via non-negative tensor decomposition, interpret each by top LRs + sender/receiver cell types. Differences: **deterministic CP/ALS** vs ALARMIST **Bayesian**; **fundamentally non-spatial** (context axis = sample/time), spatial only by pretending niches are contexts. It consumes LR scores (not a scorer itself), is not a spatial tool. For ALARMIST it's the direct baseline for the factorization step and a **template for the downstream loading-interpretation toolkit** (rank LRs, GSEA, sender/receiver outer-product networks, Gini, loading clustering/classification). ALARMIST's edge: native spatial motif discovery + probabilistic model.

## 3. LIANA+ — Dimitrov et al., Nature Cell Biology 2024
**File:** `LIANA+_Dimitrov_2024_NatCellBiol.pdf` · *LIANA+ provides an all-in-one framework for cell–cell communication inference* (26:1613–1622)

**One-liner:** A scalable, modular scverse/AnnData framework unifying LR scoring, spatial co-localization metrics, cross-condition factorization, and causal intracellular-network inference for CCC from single-cell & spatial multi-omics.

- **Background:** 100+ CCC tools each capture a partial view; single-cell methods assume co-expression = signaling, spatial methods quantify co-localization globally or locally. Most handle only protein-mediated interactions, one task, one data type, heterogeneous priors.
- **Motivation:** No single framework spans tasks/data types/mediators. LIANA+ (1) homogenizes many LR methods under one consensus, (2) adds global (MISTy-style) + local (8 spatial metrics) modeling, (3) extends to multi-omics mediators (metabolites), (4) adds cross-condition **program** discovery + differential CCC, (5) uniquely links extracellular CCC → intracellular signaling via sign-coherent **causal networks**.
- **Data:** Murine **6-OHDA Parkinson's** (paired MALDI-MSI metabolomics + Visium, 3 mice); human **myocardial infarction** (Kuppe: 29 snRNA + 28 Visium); 5 slide-tags datasets (co-localization benchmark); 4 breast-cancer Visium (spot classification); cross-condition atlases. KBs: OmniPath/LIANA consensus, MetalinksDB, PROGENy, CollecTRI, BioCypher.
- **Method:** (a) re-implements 8–9 LR methods + rank-aggregate consensus; (b) **global** MISTy multi-view learning (intraview + spatially-weighted extra views → CV linear meta-model → per-predictor importances, R²); (c) **local** 8 spatially-weighted bivariate metrics (cosine default, Pearson, Spearman, Jaccard, Moran's R from SpatialDM) with permutation p-values + pos/neg/neither categories; (d) **cross-condition programs**: concatenate local LR scores across slides → **sklearn NMF** (k via kneedle) = intercellular programs (per-spot scores + per-LR loadings); for dissociated data → **Tensor-cell2cell** or **MOFA+**; (e) differential CCC via pseudobulk → PyDESeq2; (f) causal net via **CORNETO ILP** (CARNIVAL) from LR nodes to deregulated TFs over OmniPath.
- **Model:** RBF spatial kernels `w_ij = exp(−d²/2l²)` (bandwidth 150 local); global bivariate Moran's R for significance; NMF `X ≈ WH` (k by elbow); Tensor-cell2cell CP + MOFA+ variational Bayes into 10–20 factors; DESeq2 NB-GLM; ILP with sign-coherence penalties.
- **Output:** Aggregated LR statistics (magnitude + specificity); per-spot local bivariate scores + p-values + categories; multi-view importances/R²; **NMF factors = spatial intercellular programs**; MOFA+/Tensor-cell2cell factors (per-sample scores, per-cell-type-pair LR loadings, variance per view); pathway enrichment; differential LR; causal LR→TF networks. Standardized AnnData/MuData.
- **Results:** Parkinson's: multi-view explained dopamine in intact (R²=0.535) vs lesioned (~0) striatum; localized dopamine–D2R to intact striatum. MI: NMF → 5 programs with condition-specific enrichment (ischaemic/fibrotic/myogenic; hypoxia/EGFR/TGF-β); recovered integrin–FN1/SPP1 pro-fibrotic axis. MOFA+ factor 1 separated conditions (ANOVA P=9.9e-12), pinpointed FB→myeloid axis. Causal search linked FN1/SPP1–ITGA5_ITGB1 → MAPK/SMAD.
- **Benchmark:** (1) co-localization "assumed truth" on 5 slide-tags (global Moran's R = indirect truth): spatially-agnostic LR methods scored **near chance (AUROC ~0.45–0.60)**. (2) spot-calling (malignant breast spots; heart cell-type proportions) → cosine chosen default. (3) condition separation: Tensor-cell2cell vs MOFA+ (RF on factor scores). No gold standard — orthogonal proxies (flagged).
- **Use cases:** Consensus LR scoring; spatial CCC (global + local) incl. metabolite mediators; cross-condition programs (NMF / Tensor-cell2cell / MOFA+); differential CCC; LR→TF causal networks. A hypothesis-generation platform.
- **Relevance to ALARMIST:** **Closest *published* relative and main multi-axis baseline.** Its **NMF-on-local-spatial-LR-scores** module is almost exactly ALARMIST's problem: concatenate spatially-weighted per-spot LR scores across slides → factor into k latent spatial "intercellular programs" (per-spot scores + per-LR loadings). The dissociated arm (Tensor-cell2cell / MOFA+) is ALARMIST's tensor/matrix-factorization class. ALARMIST's advances: purpose-built **Bayesian tensor model** vs off-the-shelf NMF/MOFA+; native **ON/OFF motif states** and **motif→gene-impact** (LIANA+ only post-hoc PROGENy/TF-causal). Practical roles: (1) LIANA+ local metrics (cosine, Moran's R) as **input scorers**; (2) its NMF / Tensor-cell2cell / MOFA+ factorizations as **primary baselines**; (3) its co-localization + spot-calling eval designs and heart/Parkinson's data as reusable scaffolds.

---

# Tier B — Spatial per-LR-pair scorers (input signal & per-interaction baselines)

## 4. CytoSignal — Liu et al., Nature Genetics 2026
**File:** `CytoSignal_Liu_2026_NatGenet.pdf` · *CytoSignal detects locations and dynamics of ligand–receptor signaling at cellular resolution from spatial transcriptomic data* (58:1396–1408) · Welch Lab · code: github.com/welch-lab/CytoSignal

**One-liner:** Infers LR signaling at **single-cell resolution** per LR pair per position, distinguishing **contact- vs diffusion-dependent** mechanisms, with a permutation-tested LRscore and a dynamic (RNA-velocity) extension.

- **Background:** LR signaling is contact-dependent (touching cells) or diffusion-dependent (secreted ligands to nearby cells). scRNA CCC methods ignore spatial proximity and work at cell-type level; spatial methods mostly work at group level or don't score each LR pair, and those that do treat contact/diffusion identically and infer only a static snapshot.
- **Motivation:** Need cellular-resolution, per-LR-pair, per-position, mechanism-aware, **temporally dynamic**, and **in-situ-validatable** inference. Prior methods rely on computational validation only.
- **Data:** Slide-seqV2 & Slide-tags (embryonic mouse brain/cortex); **Stereo-seq** whole mouse embryo E9.5–E16.5 (MOSTA); **new Visium HD** (2 sections, GEO **GSE291088**) paired with **PLA** ground truth (5 LR pairs); Parkinson's mouse model (4 young + 5 aged); scDesign3 simulations.
- **Method:** `LRscore = L × R` in a cell's neighborhood; complexes `(L1+L2)×(R1+R2)`. **Diffusion**: ε-ball (200 µm), ligand Gaussian-weighted + density-normalized within sender's ε-ball, receptor raw (membrane-bound). **Contact**: Delaunay neighbors, uniform weight. Gaussian bandwidth derived analytically `σ ≤ t/√(−2 log ε)`. Mean over Delaunay neighbors + Gaussian smoothing. Significance via **spatial permutation** (100k) + spatial FDR; **SPARK-X** for spatial variability. LR DB: CellPhoneDB v2. Signaling-associated genes via elastic-net (LRscore ~ genes + cluster labels, cell-type confounder controlled) → GO. **VeloCytoSignal**: `dS/dt` via product rule using RNA velocity (VeloVAE). Multisample differential via **NEBULA** negative-binomial regression (library-size offset).
- **Model:** `S = L×R`; analytic Gaussian null; NB regression (FDR<0.05, |log2FC|>0.1); VeloVAE ODE velocity `vu=α−βu, vs=βu−γs`.
- **Output:** Per-interaction LRscore maps + significant-location maps; 3D sender→receiver edge plots; cluster Circos (contact vs diffusion); spatial gradients; signaling-gene heatmaps + GO; differential volcanoes; VeloCytoSignal 3D velocity arrows.
- **Results:** Recovered canonical developmental interactions (Sema3a–PlexinA4/Nrp1, Wnt5a–Epha7, Efnb1–Epha4) with correct GO; diffusible signaling near fibroblasts, contact near radial glia; choroid-plexus contact signaling peaks E14.5. PD model: 42 diffusion + 3 contact interactions changed young→old (immune, e.g. Spp1–Cd44, up in aged, matching microglia increase). VeloCytoSignal correctly predicted (blinded) temporal trends validated on consecutive time points.
- **Benchmark:** Literature review of >100 methods → 8 categories; **Category 1** (spatial + single-spot + per-LR) = CytoSignal, stLearn, SpatialDM, COMMOT, LIANA+, NICHES. **PLA-anchored AUC/AUPRC**: CytoSignal higher than the five baselines across most regions (gain from spatial smoothing). More physically plausible sender–receiver distances vs CellChat/CellPhoneDB. Markedly faster/lighter (72-run efficiency). scDesign3: low FPR/FNR (real-replicate FPR median 0).
- **Use cases:** TME signaling; tissue-repair dynamics (VeloCytoSignal on time series); disease/age/stage/genotype differential signaling. Needs cell×gene + cell×position; VeloCytoSignal needs spliced/unspliced.
- **Relevance to ALARMIST:** A **per-LR-pair, per-cell spatial scorer** (not a factorizer) — the kind of front-end signal ALARMIST factorizes into motifs, and a strong per-interaction baseline. Its diffusion-vs-contact neighborhood scoring, Gaussian ligand imputation, and especially its **PLA in-situ validation** strategy are directly comparable to and could inform ALARMIST's LRI scoring and evaluation. Also from the **Welch Lab** (same VeloVAE lineage).

## 5. SpatialDM — Li et al., Nature Communications 2023
**File:** `SpatialDM_Li_2023_NatCommun.pdf` · *SpatialDM for rapid identification of spatially co-expressed ligand–receptor and revealing cell–cell communication patterns* (14:3995)

**One-liner:** Detects spatially co-expressed LR pairs via a **bivariate Moran's R** with a closed-form analytical null, plus single-spot interaction hits and shared communication patterns.

- **Background:** scRNA CCC methods lack coordinates → false positives from non-proximal cell types. Early ST methods (Giotto, SVCA, SpaOTsc, SpaTalk) test *which cell types* interact for all LRIs rather than first selecting *which LR pairs* are genuinely spatially co-expressed.
- **Motivation:** Reframe the primary task as **feature selection of spatially co-expressed LR pairs** + single-spot localization; remove the slow permutation bottleneck via an analytical null.
- **Data:** Thrane melanoma (ST, 293 spots); Fawkner-Corbett human intestine (Visium, 8 slices, GSE158328); mouse SVZ (seqFISH+). LR DB: CellChatDB v1.1.3. Deconvolution via RCTD. SVCA-based simulations.
- **Method:** Global bivariate Moran's R per LR (RBF-kernel weights; complexes = mean over subunits) → z-score/permutation test (FDR<0.1) → local bivariate Moran's R per spot (sender + receiver components) for single-spot hits → cluster binary local-hit matrix into "patterns" via **SpatialDE.aeh** → per-pattern pathway enrichment (Fisher) → cell-type chord diagrams → differential interactions via LM on z-scores (LR test). Secreted: RBF kernel; ECM/contact: k-NN weights.
- **Model:** Global `R = Σᵢⱼ w_ij (x_i−x̄)(y_j−ȳ) / [√Σ(x−x̄)²·√Σ(y−ȳ)²]` (x=ligand, y=receptor); **analytical** `Var(R)` closed-form → z, survival-fn p (>100× faster than permutation). Local Moran's R per spot = sender + receiver terms. Differential: χ² likelihood-ratio on LM of z-scores.
- **Output:** Ranked spatially co-expressed LR pairs (global R, z, FDR); per-spot interaction hit maps; pattern assignments; per-pattern pathway dot plots; sender→receiver chords; differential volcanoes/heatmaps. Python/Scanpy.
- **Results:** Analytical ≈ permutation (Spearman R>0.9); scales to **1M spots in ~12 min, 1 CPU** (>100×). Sim AUROC 0.912. Melanoma: 103 pairs (FDR<0.1), 3 patterns (lymphoid/melanoma/CAF); recovered VEGF/SPP1/CSF1; CD23 validated. Intestine: recovered 326/414 prior + 220 novel; high cross-replicate reproducibility; EGF in adult crypt-top; 146 adult- vs 97 fetus-specific pairs.
- **Benchmark:** vs CellChat, Giotto, SpaTalk, SpatialCorr. Led AUROC (0.912 vs 0.570–0.723 at 25%), only method scalable to 1M spots, best negative-control FPR (shuffled LR list). No experimental gold standard.
- **Use cases:** Prioritize spatially co-expressed LR pairs before interpreting CCC; single-spot hotspots; regional pattern annotation; differential interactions across conditions/stages; million-spot scale. ST/Visium/seqFISH+ and spatial proteomics.
- **Relevance to ALARMIST:** A **per-LR-pair spatial co-expression scorer** at single-spot resolution — **not** a factorizer (its "patterns" come from SpatialDE clustering of the binary hit matrix, not a joint decomposition). The paper explicitly flags pair-independent analysis as a limitation — the multi-pair coupling ALARMIST solves. Best framed as (1) a **strong input scorer** (its local sender/receiver Moran interaction scores are ideal tensor entries), (2) a **per-pair baseline**, and (3) a shared-goal predecessor (cell-type attribution via chords; the open LR→downstream-gene challenge ALARMIST tackles).

## 6. COMMOT — Cang et al., Nature Methods 2023
**File:** `COMMOT_Cang_2023_NatMethods.pdf` · *Screening cell–cell communication in spatial transcriptomics via collective optimal transport* (20:218–228)

**One-liner:** Spatially-constrained, **competition-aware** LR CCC inference that solves a "collective optimal transport" problem to infer per-LR signaling flux between every spot/cell pair.

- **Background:** scRNA CCC tools lack spatial range → false positives; early spatial tools (Giotto, CellPhoneDB v3, stLearn, SVCA, NCEM) treat cell pairs/neighborhoods independently, ignoring **collective competition** among species and cells.
- **Motivation:** Two unmodeled realities: most ligands/receptors bind **multiple partners** (72%/60% in FANTOM5 → competition), and signaling has **finite range**. Classical OT forces equal mass (breaks unit comparability); unbalanced/partial OT have flaws. **Collective OT** preserves species comparability, caps total signal at available amounts, enforces range, handles many competing species.
- **Data:** 8 spatial datasets across 5 technologies + PDE simulations: MERFISH hypothalamus, STARmap placenta/cortex, seqFISH+ cortex, Slide-seqV2 hippocampus, Visium breast cancer & mouse brain; SpaOTsc-integrated epidermis & Drosophila. LR DB: CellChatDB (secreted). Target-gene truth from scSeqComm (TRRUSTv2/RegNetwork).
- **Method:** Solve one global collective-OT problem → 4D coupling tensor `P[i,j,k,l]` (ligand i, receptor j, sender spot k, receiver spot l). Inequality marginals cap total ligand sent ≤ expression and receptor received ≤ expression (competition/capacity); costs infinite beyond per-pair range `T(i,j)`; complexes = min over subunits. Downstream: interpolate to sender/receiver **vector fields**; cluster×cluster CCC with permutation p-values; per-cell CCC profile for clustering; signaling-DE genes vs received signal via **tradeSeq**; RF impact scoring (Gini of received-signal covariate).
- **Model:** `min_P Σ⟨P[i,j],C(i,j)⟩_F + Σ F(μ_i) + Σ F(ν_j)` with non-binding zeros + two inequality marginals; entropy-regularized, solved by **stabilized log-domain Sinkhorn**. Linear scaling in #locations (range constraint bounds nonzeros). Eval: cosine (vector fields), Jaccard (binarized networks), Spearman (received signal vs target-gene activity).
- **Output:** Per-LR spot×spot directed matrix; sender/receiver vector fields; cluster CCC networks + permutation p; per-spot CCC profiles + clusters; ranked signaling-DE genes + RF impact. Python: github.com/zcang/COMMOT.
- **Results:** On PDE sims with competition, outperforms pairwise/unbalanced/partial OT. Epidermis: GAS6/PROS1-TYRO3 (validated by IF/RNAscope). MERFISH: OXT most active. STARmap placenta: midkine vs IGF opposing directions. Visium breast: midkine most active; COL1A1/S100G signaling-DE.
- **Benchmark:** (1) vs pairwise/unbalanced/partial OT on PDE truth (more accurate). (2) vs CellChat, Giotto, CellPhoneDB v3 on real data via Spearman(inferred CCC, scSeqComm target activity): medians 0.237/0.180/0.230; stronger than baselines on most, resolves localized hotspots cluster methods miss.
- **Use cases:** Competition-aware per-pair spatial LR scoring on any ST platform; signaling-direction vector fields; cluster sender/receiver summaries; nominating signal-regulated genes.
- **Relevance to ALARMIST:** **Complementary upstream input scorer**, not a factorizer — a per-LR, per-spot-pair OT scorer + cluster aggregator producing exactly the spot×LR CCC signal ALARMIST decomposes. No latent motifs / rank / cross-cohort programs. Its **competition-aware spatial-range OT** could improve the fidelity of per-LRI values ALARMIST factorizes; its downstream philosophy (senders/receivers, direction, CCC→gene impact) mirrors ALARMIST's characterization.

## 7. stLearn — Pham et al., Nature Communications 2023
**File:** `stLearn_Pham_2023_NatCommun.pdf` · *Robust mapping of spatiotemporal trajectories and cell–cell interactions in healthy and diseased tissues* (14:7739)

**One-liner:** A Python ST toolkit fusing expression + distance + H&E morphology for (1) spatial trajectories (**PSTS**), (2) low-FDR spatial LR CCC via a two-level permutation test (**SCTP**), and (3) dropout imputation (**stSME**).

- **Background:** ST yields expression + coordinates + (Visium) H&E, but most tools use expression alone. stLearn fuses all three via `f(G,I,D)`.
- **Motivation:** scRNA-derived CCI methods predict interactions between spatially **distant** cell types (true interactions within ~0–200 µm) → high FDR; trajectory tools ignore space; PCR-based ST has dropout. stLearn encodes distance + morphology to fix all three.
- **Data:** In-house Visium mouse TBI brain + human BCC (RNAscope IL34–CSF1R validation); public Visium breast cancer, sci-Space embryo, seqFISH+ cortex, Slide-seqV2 hippocampus, Visium human brain. Simulated ST via per-gene negative-binomial fits to 11 scRNA cell types.
- **Method:** **PSTS** — spatial layer on Diffusion Pseudotime; combine pseudo-temporal distance (cosine) + spatial distance `dPTS = ω·dPT + (1−ω)·dS`; spatial-PAGA graph → **Chu-Liu/Edmonds** min directed spanning arborescence; ω via Laplacian spectral-distance balancing (~0.46–0.51). **SCTP CCI** — neighborhoods via cKDTree; `LRscore = ½[mean(Expr_L·1[R>0]) + mean(Expr_R·1[L>0])]` × optional heterogeneity; **level-1** permutation (expression-matched background gene pairs, Canberra-quantile, k=1000, BH) → per-spot LR p-values; **level-2** permutation (cell-type label shuffles) → directional sender→receiver p. Numba-parallel + binning for millions of cells. **stSME** — combine distance, PCA-expression correlation, and ResNet50 H&E morphology features → weighted imputation.
- **Model:** Modified DPT; `dPTS` convex weighting; graph Laplacian spectral distance for ω; two nested permutation tests + optional NB tail test; variogram/Matheron semi-variance for benchmarking spatial continuity.
- **Output:** PSTS pseudo-time-space + rooted directed trajectory + transition genes; SCTP ranked LR pairs + per-spot significance maps + sender→receiver networks; stSME imputed matrix + gap-filling + improved clusters. Python (AnnData) + no-code i-stLearn.
- **Results:** PSTS: dorsoventral microglia-activation trajectory in TBI (validated across 6 time points); novel inside-out cortical-migration branch; DCIS→IDC progression. SCTP: Gas6–Axl (SVZ), Apoe–Lrp1 (hippocampus), GPC3–IGF1R (breast DCIS); IL34–CSF1R confirmed by RNAscope. stSME: rescued dropouts, separated CA1/CA3, improved ARI.
- **Benchmark:** PSTS vs Slingshot/Monocle3/SpaceFlow → lowest variogram semi-variance + unique correct branching. **SCTP vs 8 methods** (Squidpy, CellPhoneDB, NATMI, SingleCellSignalR, CellChat, NCEM, SpaTalk, spaOTsc) on NB-ground-truth sim + Visium breast: **only method with zero false positives** and no distal-cluster interactions. stSME via leave-out ARI + Moran's I reproducibility.
- **Use cases:** Spatial/spatio-temporal trajectories (activation, development, cancer progression, cross-section integration); low-FDR spatial LR CCI (scalable via binning); dropout imputation + clustering when H&E available.
- **Relevance to ALARMIST:** A **per-LR-pair, per-location spatial scorer + cell-type-pair aggregator** — not a factorizer. SCTP produces exactly the per-spot per-LR scores that could be ALARMIST tensor entries, and its directional sender→receiver matrices mirror ALARMIST's motif characterization — but it ranks LRs/cell-type pairs one at a time and never decomposes into shared programs. The paper itself contrasts LR-testing methods with program-discovery (DIALOGUE); **ALARMIST sits on the program-discovery side but built on LR interactions.** Best positioned as a per-LR spatial scorer/baseline feeding ALARMIST.

## 8. NICHES — Raredon et al., Bioinformatics 2023
**File:** `NICHES_Raredon_2023_Bioinformatics.pdf` · *Comprehensive visualization of cell–cell interactions in single-cell and spatial transcriptomics with NICHES* (39(1):btac775)

**One-liner:** An R tool that scores LR signaling at single-cell/single-spot resolution via `ligand(sender) × receptor(receiver)`, yielding an LR-by-interaction matrix that is embedded and clustered like a standard single-cell object.

- **Background:** Applications Note (formalism in supplement). A crowded CCC field (CellPhoneDB, NicheNet, CellChat, Connectome…) already exists.
- **Motivation:** All prior tools use **cluster/cell-type mean** expression, hiding single-cell signaling heterogeneity. NICHES scores individual **cell–cell pairs** and **spatial spots** so intra/inter-cluster heterogeneity and archetype shifts become observable.
- **Data:** Demonstrations only: rat pulmonary alveolus scRNA (3 co-localized types); brain ST (oligodendrocyte spots). Built-in LR from OmniPath/FANTOM5; custom lists accepted. Outputs feed Seurat/Scanpy/Monocle3.
- **Method:** Build a matrix (rows = LR mechanisms, cols = interactions). **CellToCell**: per sampled cell pair, entry = `L(sender) × R(receiver)`. **CellToNiche**: entry = `mean L over neighbors × R(cell)`; spatial mode restricts to nearest neighbors → per-spot niche vector. Then treat as ordinary features → UMAP/clustering/subclustering + Seurat marker DE.
- **Model:** Deliberately simple bilinear product `L×R` per mechanism; **no** generative model, OT, factorization, or permutation significance in main text (significance via downstream Seurat DE).
- **Output:** LR-mechanism × interaction matrix (cell-pair "signaling" + spot "niche" objects), embeddable/clusterable; per-cross marker LR mechanisms; spot niche signatures plotted in situ; differential signaling across conditions.
- **Results:** Qualitative only. Cross-level signaling signatures + intra-relationship heterogeneity invisible to mean-based tools; Fgf1–Fgfr2 recovers oligodendrocyte spots (sanity check). Claims: archetype shifts with conserved means; system-level effects of adding/removing populations.
- **Benchmark:** **None quantitative** — Applications Note; comparison to prior methods is conceptual (single-cell/spot resolution vs cluster means).
- **Use cases:** Per-cell-pair or per-spot LR scoring; embeddable signaling/niche atlases; nearest-neighbor niche mapping; differential signaling; generating a per-observation LR feature matrix for Seurat/Scanpy.
- **Relevance to ALARMIST:** A **per-LR-pair / per-interaction scorer**, the kind of upstream primitive ALARMIST factorizes — **not** a competitor to the factorization. Its spatial niche mode (NN-restricted, mean-ligand) is essentially a per-spot LR-interaction tensor slice. It has no motif/program layer, no tensor rank, no cross-sample latent structure — it discovers niches by empirical clustering, not by factorizing communication. Best framed as an **input-scoring/preprocessing step** and a lightweight conceptual predecessor/baseline.

---

# Tier C — Non-spatial cell-type-level aggregator (upstream scorer / database)

## 9. CellChat v2 — Jin et al., Nature Protocols 2024
**File:** `CellChat_Jin_2024_NatProtoc.pdf` · *CellChat for systematic analysis of cell–cell communication from single-cell transcriptomics* (20:180–219)

**One-liner:** An LR database + mass-action interaction-scoring & network-analysis toolkit inferring cell-type→cell-type communication from scRNA-seq (this paper is the **v2 step-by-step protocol**; method from Jin et al., Nat Commun 2021).

- **Background:** scRNA-seq probes signaling between annotated groups using LR databases. Prior tools (CellPhoneDB) score LR pairs but a versatile toolkit combining network analysis + visualization + cross-condition comparison was lacking. CellChat ranks among top performers in independent benchmarks of >15 methods.
- **Motivation:** v2 adds (1) expanded **CellChatDB v2** (~3,300 interactions incl. non-protein/metabolic/synaptic from NeuronChatDB + CellPhoneDB), (2) systematic **multi-condition comparison**, (3) an **interactive Explorer** (Shiny/Plotly). Explicitly non-spatial, hypothesis-generating.
- **Data:** Human skin scRNA (atopic dermatitis NL vs LS; demo 5,011 cells / 12 groups; GSE147424); differing-composition example (GSM3453535–38); ~300k-cell skin atlas (E-MTAB-8142) for runtime. Cited apps: hair-follicle development, wound-healing aging, COVID brain/lung, kidney injury, CRC + PD-1.
- **Method:** Inputs: normalized gene×cell + group labels. Per group, average L/R expression via robust statistic (**triMean** default — fewer/stronger interactions; or truncatedMean to recover weak signaling). `computeCommunProb` scores each (sender, receiver, LR) via **law-of-mass-action + Hill saturation + cofactor** (agonist/antagonist) modulation, handling complexes. Permutation test (shuffle group labels) → p. Sum LR probs into pathways. Systems analysis: network centrality (sender/receiver/mediator/influencer); **NMF pattern recognition** (K via Cophenetic/Silhouette) → outgoing/incoming patterns (river plots); manifold learning (functional/structural similarity → UMAP). Cross-condition: merge, compare interaction number/strength, joint manifold + rankSimilarity, rankNet + paired Wilcoxon, DE-driven up/down LR pairs.
- **Model:** 3D array `P(sender, receiver, LR)` (`@net$prob`) + p-value array; mass-action/Hill `~ K1·[L][R]` modulated by agonist/antagonist cofactors; pathway prob sums LRs; **NMF** for patterns (Cophenetic/Silhouette rank); paired Wilcoxon + logFC thresholds for cross-condition.
- **Output:** 3D sender×receiver×(LR or pathway) probability + permutation p; ranked senders/receivers/mediators/influencers; outgoing/incoming patterns; pathway similarity groupings + UMAP; cross-condition differential counts/strengths, rewired pathways, up/down LR pairs. Circle/chord/hierarchy/heatmap/bubble plots + interactive Explorer. ~5 min (~15 min at 300k cells).
- **Results:** Reproduces biology: on atopic-dermatitis skin, **CCL19–CCR7** top event in lesional skin (inflammatory fibroblasts → dendritic cells). Interaction counts depend strongly on averaging statistic (triMean 19 vs 10%-trim 82 vs 5%-trim 140). Near-linear runtime scaling.
- **Benchmark:** **No new head-to-head accuracy benchmark** — cites independent evaluations (Dimitrov 2022, Liu 2022, ESICCC 2023) ranking CellChat among top performers of >15 methods. Internal quantitative study = trim-parameter sensitivity + runtime scaling.
- **Use cases:** Infer/rank cell-type→cell-type signaling from annotated scRNA; aggregate LR into pathway networks; identify dominant senders/receivers + coordinated patterns; **compare two conditions** (disease vs normal, young vs aged, treated vs control). Non-spatial; points to other tools for proximity.
- **Relevance to ALARMIST:** A **non-spatial, cell-type-level LR aggregator and canonical input scorer/database** — an upstream/baseline component, not a motif-discovery competitor. Its output is literally a sender×receiver×LR tensor (the shape ALARMIST decomposes), and **CellChatDB v2** is a commonly borrowed curated LR prior. Its NMF pattern-recognition is the closest CellChat-internal analog to ALARMIST motifs but operates on **non-spatial group-level** networks with no spatial location or downstream-gene-impact modeling. R-only (vs ALARMIST Python/Torch).

---

## Synthesis — how these map onto ALARMIST

**ALARMIST's true peers (communication factorization):** **COMPOTES** (spatial LR + non-negative CP across a cohort) is the closest sibling and the primary head-to-head competitor; **Tensor-cell2cell** is the direct non-spatial predecessor of the factorization idea; **LIANA+** is the closest *published* relative because it bundles NMF-on-spatial-LR-scores *and* Tensor-cell2cell/MOFA+ dissociated factorization. ALARMIST differentiates via a **purpose-built Bayesian tensor model** (uncertainty, **ON/OFF motif states**) and a **motif→downstream-gene-impact** characterization that none of the three build in.

**The front-end scorers ALARMIST factorizes / benchmarks against (per-LR spatial):** **CytoSignal** (PLA-validated, contact-vs-diffusion, cellular resolution), **SpatialDM** (bivariate Moran's R, million-spot scale), **COMMOT** (competition-aware optimal transport), **stLearn/SCTP** (low-FDR two-level permutation), **NICHES** (bilinear per-pair). Any of these can supply the per-spot × per-LR interaction signal that ALARMIST decomposes into motifs, and each is a per-interaction baseline.

**Upstream aggregator / prior:** **CellChat v2** — non-spatial cell-type-level scoring and the widely reused **CellChatDB** LR database.

**Recurring gap ALARMIST fills:** every per-LR scorer here (CytoSignal, SpatialDM, COMMOT, stLearn, NICHES) explicitly analyzes LR pairs **independently** and stops at ranking; every factorizer here is either **non-spatial** (Tensor-cell2cell, CellChat-NMF, MOFA+) or uses **off-the-shelf deterministic NMF/CP** (COMPOTES, LIANA+). ALARMIST's contribution — a **Bayesian, spatial, motif-state-aware** factorization of LR interactions with downstream gene-impact — is precisely the intersection none of these occupy alone.

---

# Part 2 — Fine-grained dimension comparison (for single-cell spatial / Xenium)

How each method resolves three orthogonal axes — **(A) how sender/receiver cell types enter**, **(B) how ligand & receptor are modeled**, **(C) the spatial grain** (strict single-cell vs neighborhood/region vs cell-type-pooled) — plus **the single-cell-resolution (Xenium) path specifically** and **the multisample/differential mode**. Read this when deciding what to run on Xenium (also CosMx/MERFISH).

## Axis A × B × C — the master matrix

| Method | (A) Cell-type sender/receiver handling | (B) Ligand/Receptor modeling | (C) Spatial grain |
|--------|----------------------------------------|------------------------------|-------------------|
| **CytoSignal** | **Optional, post-hoc.** Score never uses cell type; each receiving cell's ligand sum *names its sender cells* → directionality is intrinsic. Cell type only as GLM covariate, differential covariate, Circos/3D coloring. Not pooled at scoring. | `S = L×R`. **Receptor = cell's own raw expr (no smoothing).** Ligand = neighbor sum: Gaussian-weighted ε-ball (diffusion) or Delaunay uniform (contact), sender-density-normalized. Complexes **summed** `(L1+L2)(R1+R2)`. Analytic σ. + Delaunay-mean + 200 µm smoothing. No competition. | **Strict per-cell** (per receiving cell × LR). Neighborhood collapsed to one number on index cell; index included. |
| **COMPOTES** | **Not in scoring; post-hoc** from matched snRNA pseudobulk (cell-type-pooled): L1-normalized sending/receiving capacity → per-program sender→receiver map. Visium cell types via **deconvolution** (downstream only). | `L×R`. Receptor = spot's own; ligand = diffusion sum `w=d/(2πr)`, 3 hardwired Visium rings (membrane/medium/long). Complexes **min** over subunits. **Competition** (normalize ligand by all ligands binding that receptor). | **Per-Visium-spot × LRI**; matrix = (all spots, all samples) × LRI. |
| **Tensor-cell2cell** | **Mandatory — cell type IS the unit** (tensor dims 3&4 = sender/receiver types). **Cell-type-pooled** aggregation (fraction-nonzero default, or mean). Must be shared across all contexts. | Cell-type-level score: **Expression Mean** (default), product, or geometric. Complexes **min** over subunits. **No spatial term, no competition.** | **Aspatial / cell-type-pooled.** Score = (context, LR, sender-type, receiver-type). Coordinates never used. |
| **LIANA+** | **Component-dependent.** Dissociated LR methods: mandatory, cell-type-pooled, per-type-pair. Local metrics: none (post-hoc). Multi-view: cell-type proportions (deconvolution) as a predictor. | (A) 9 dissociated scorers; (B) 8 **local bivariate** metrics (**cosine default**, bivariate Moran's R…). Complexes **min**. RBF kernels (Gaussian default, bandwidth 150, cutoff 0.1); local self-weight=1, multi-view=0. | **All three:** per-type-pair (dissociated), per-slide-global (multi-view), **per-spot/cell (local)**. |
| **CytoSignal / SpatialDM / COMMOT / stLearn / NICHES** | *(see individual rows)* | *(see individual rows)* | *(spatial scorers — per cell/spot or per cell-pair)* |
| **SpatialDM** | **Not in core test** (feature-selects LR pairs). Post-hoc: LM vs composition; chord weight `n_AB` = per-spot local sender/receiver score × composition. Spot: **deconvolution** (RCTD) downstream. Directionality at molecule level (ligand=sender). | **Bivariate Moran's R** (spatial cross-correlation, *not* a product); z-scored expr for local. Complexes **arithmetic mean** (geometric optional). **No competition** (pair-independent — stated limitation). RBF kernel bandwidth `l` (coordinate units!); secreted=RBF, contact/ECM=**kNN (default 6)**. | **Two:** global (per-LR per-sample) + **local (per-spot per-LR)**. Soft neighborhood; index included **unless single-cell → set w_ii=0**. |
| **COMMOT** | **Not in core** (per-cell-pair OT). Post-hoc: cluster→cluster directed `S^cl` + permutation p; per-type enrichment. No deconvolution for spots. | **Collective optimal transport** (*not* a product); raw expr as mass. **Competition** via inequality marginals (coupled signal ≤ available L/R). Complexes **min**. **Hard distance cutoff `T`** (not a smooth kernel) + φ scaling + entropy ε; Sinkhorn. | **Per-cell-pair** `P_{i,j,k,l}` → per-cell summaries (received signal, CCC profile). Neighborhood implicit via `T`. |
| **stLearn (SCTP)** | **Optional 2nd stage.** Base LR score type-agnostic. Cell type → optional diversity weight `HET`, and directional **per-type-pair CCI count matrix** (permutation). Discrete or **mixture** (deconvolution/label-transfer for spots, C=0.2). | Symmetrised half-thresholded neighborhood co-expression (one side = continuous neighborhood mean, other = **binary detected** in index). **Hard uniform binary neighborhood** (no distance decay), cKDTree. Complexes/competition **not modeled** (connectomeDB). | **Per-unit (spot/bin/cell) × LR** neighborhood score; downstream **per-type-pair** CCI. |
| **NICHES** | **Used, but score stays per-cell-pair** (not cluster mean — its headline claim). Cell type defines which crosses to sample (dissociated) / labels edges (spatial). Directionality explicit. | **Product** `L(sender)×R(receiver)`; niche modes = mean ligand over neighbors × receptor on index. Complexes supported (subunit rule in supplement). **Hard binary neighbor mask** (no kernel). | **Per-cell-pair** (cell-cell modes) or **per-cell/spot niche** (niche modes). |
| **CellChat v2** | **Mandatory, atomic unit. Cell-type-pooled** (triMean default → zeros genes in <25% of a group; or truncatedMean/mean). Per-type-pair. Directionality explicit; autocrine included; `min.cells=10`. No deconvolution (scRNA only). | **Mass-action + Hill** on group-mean expr; complexes **geometric mean** (from 2021 paper); agonist/antagonist cofactors; permutation (label shuffle); pathway = sum. **No spatial** (optional PPI-network smoothing ≠ spatial). | **Aspatial / cell-type-pooled.** Per-type-pair × LR (or pathway). |

**One-line grain summary:** strict single-cell = **CytoSignal**; per-cell-pair = **COMMOT, NICHES**; per-cell/spot neighborhood = **SpatialDM (local), LIANA+ (local), stLearn**; per-spot (multi-cell) = **COMPOTES**; cell-type-pooled/aspatial = **Tensor-cell2cell, CellChat, LIANA+ (dissociated)**.

## The single-cell (Xenium) path

| Method | Native at single-cell? | Neighborhood at single-cell res | Must-set / key params | Deconv? | Precomputed cell types? |
|--------|:----------------------:|--------------------------------|-----------------------|:-------:|:-----------------------:|
| **CytoSignal** | ✅ **designed for it** | Delaunay (contact) / 200 µm ε-ball (diffusion); **index included, self-weight ≠ 0** | radius in µm (200 default; shrink for dense tissue); contact-vs-diffusion per pair; **imputeLR essential** | ❌ | Only for GLM/Circos/differential (not the score) |
| **COMMOT** | ✅ **native "high-res" path** (=MERFISH/seqFISH+ class) | implicit: all cell pairs with dist ≤ `T`; autocrine included | distance cutoff `T` (µm), φ, entropy ε, direction top-k | ❌ | Only for cluster-level/enrichment output |
| **SpatialDM** | ✅ (SVZ seqFISH+ done as single-cell) | RBF over µm dist; **set w_ii = 0**; kNN=6 for contact/ECM | **bandwidth `l` in µm** (paper's `l` not µm — re-derive, e.g. ~50–75 µm); use **z-score null** for scale | ❌ | Only for chord/interpretation |
| **stLearn (SCTP)** | ✅ (seqFISH+ single-cell; bin for scale) | cKDTree, **between-spot** mode, hard uniform radius | **µm distance (don't use "2× spot" default)**; grid-bin for 10⁵–10⁶ cells | ❌ (discrete labels) | ✅ required (discrete mode) |
| **NICHES** | ✅ **in its wheelhouse** | spatial nearest-neighbor mask (CellToCellSpatial / NeighborhoodToCell) | neighbor rule (k/radius) + self-inclusion **live in supplement/vignette — set them** | ❌ | ✅ required (ident column) |
| **LIANA+ (local)** | ✅ resolution-agnostic | RBF kernel over µm coords; local self-weight=1 (multi-view=0) | **recalibrate bandwidth `l` to µm** (default 150 = Visium-tuned); cosine default | ❌ | Optional for local; ✅ for dissociated/factorization |
| **Tensor-cell2cell** | ⚠️ **not spatial** | none — you must supply **niches/regions as "contexts"** (external clustering) | context definition; shared cell types across contexts | ❌ | ✅ required (the unit) |
| **COMPOTES** | ❌ **Visium-only; needs reimplementation** | (A) bin→pseudospots ~50–100 µm, or (B) per-cell + rebuilt kernel | rebuild neighbor rings + `w`; scale likely **infeasible** at cell res | (Visium) | Post-hoc only |
| **CellChat v2** | ❌ aspatial in this protocol | none (coordinates ignored) | triMean→truncatedMean + PPI smoothing for sparsity | ❌ | ✅ required (the unit) |

> **For a Xenium single-cell workflow, the "just works" tier is: CytoSignal, COMMOT, SpatialDM, stLearn, NICHES, LIANA+ (local metrics).** COMPOTES and Tensor-cell2cell are not single-cell-spatial tools; CellChat (this protocol) is non-spatial. In **all** spatial tools, remember Xenium coords are already µm — but every neighborhood parameter (`l`, `T`, radius) whose published default was tuned on Visium/ST **must be reset to a micron-appropriate value**, and SpatialDM specifically needs **w_ii = 0** at single-cell resolution.

## Multisample / cross-condition (differential-LRI) mode

| Method | Native? | Statistical model | Unit tested | Significance / notes |
|--------|:-------:|-------------------|-------------|----------------------|
| **CytoSignal** | ✅ | **Negative-binomial mixed model** (nebula) on floored, *unnormalized* LRscores; covariates age/phenotype/**cell type**; library-size offset | **per-LR-pair** (optionally attributed to a cell type) | Wald χ²; **FDR<0.05 & \|log2FC\|>0.1**; low FPR (real-replicate median 0); run contact/diffusion separately |
| **LIANA+** | ✅ (3 ways) | (1) pseudobulk→**PyDESeq2** per cell type → join to LR; (2) **Tensor-cell2cell** 4D; (3) **MOFA+** multi-view; (+spatial NMF across slides) | (1) **per-type-pair LR**; (2/3) **per-factor** | sample = replicate throughout; DESeq2 handles design covariates; ANOVA/Fisher on factor scores |
| **SpatialDM** | ✅ | Per-LR per-sample **global Moran z** → LM full vs reduced → **likelihood-ratio χ²** | **per-LR-pair** (also pathway) | FDR<0.1; supports **continuous covariate/time**; fit each sample independently first; small-n is hard |
| **Tensor-cell2cell** | ✅ (central) | Context axis = samples/conditions/timepoints/niches; **non-negative CP**; post-hoc on context loadings | **per-factor** (module) | t-test / Spearman / clustering / RF on loadings; **no covariate regression** |
| **CellChat v2** | ⚠️ pairwise only | mergeCellChat; `compareInteractions`; `rankNet` **paired Wilcoxon**; DE-based up/down LR; `liftCellChat` for differing compositions | per-LR / pathway / per-type-pair | **Not replicate-aware, pairwise only**; >2 conditions = repeated pairwise; for replicates → MultiNicheNet/Tensor-cell2cell |
| **COMPOTES** | ⚠️ pooled, post-hoc | Concatenate all samples → shared **NMF/CP** (no condition axis); differential is post-hoc on per-sample **program-activity descriptor** | **per-program**, unit = sample/patient | Mann-Whitney / linear mixed-effects vs clinical/mutation/CNV; cosine≥0.8 to match programs across runs |
| **COMMOT** | ❌ none | all stats within one sample (permutation, tradeSeq) | — | workaround: fit per sample, compare summaries externally |
| **stLearn (SCTP)** | ❌ none | single-sample permutation vs within-section background | — | workaround: per-sample → external test (sample = replicate) |
| **NICHES** | ❌ no engine | hand off matrix to Seurat/Scanpy DE, tag columns by condition | per cell-pair/niche column | **pseudoreplicated at sample level** (unaddressed); good at "archetype shifts with conserved means" |

**Multisample takeaways:** for a *statistically principled per-LR-pair* differential across a Xenium cohort, **CytoSignal (NB mixed model)** and **SpatialDM (LRT on per-sample Moran-z)** are the cleanest, and **LIANA+'s pseudobulk→DESeq2** path is the most standard (sample = replicate). **Tensor-cell2cell / LIANA+-factorization / COMPOTES** give *program-level* condition differences (closest to how you'd compare ALARMIST motifs across conditions). **CellChat** is pairwise and not replicate-aware; **COMMOT, stLearn, NICHES** have **no** native cross-condition test — you fit per sample and compare downstream, treating the sample (not the cell) as the replicate.

## Practical recommendation for single-cell (Xenium) LRI

- **Per-cell, per-LR score with a significance test, mechanism-aware, PLA-validated:** **CytoSignal** — closest to ALARMIST's per-cell interaction unit; use `imputeLR` (do not skip — Xenium is sparse), keep µm radius, supply a panel-matched DB, and it has a native cohort NB differential mode.
- **Competition/diffusion-physics-aware per-cell-pair fluxes + signaling direction fields:** **COMMOT** — native single-cell path; tune `T` in µm; no built-in differential.
- **"Which LR pairs are spatially co-localized, and where," at massive scale + clean per-pair differential:** **SpatialDM** — set **w_ii=0**, z-score null (scales to 1M cells), LRT multisample.
- **Signaling as an embeddable single-cell feature space for clustering/heterogeneity:** **NICHES** (per-cell-pair, spatial modes) or **stLearn SCTP** (low-FDR, bin for scale).
- **One framework spanning local scoring + program factorization + differential:** **LIANA+** — recalibrate the kernel bandwidth to µm; its local cosine feeds a factorization and its pseudobulk→DESeq2 gives replicate-aware differential.
- **Not for native Xenium single-cell:** **COMPOTES** (Visium-spot; would need bin-to-pseudospot and is memory-bound), **Tensor-cell2cell** (aspatial; needs externally-defined niches as contexts), **CellChat v2 protocol** (aspatial cell-type-pooled — gives the same answer regardless of where cells sit).

---

# Part 3 — What questions can each method answer? (ALARMIST included)

## Can a method return a value for *any* (sender cell type × receiver cell type × ligand × receptor) tuple?

**Not all of them, and not in the same sense.** Universal constraints first: (i) every method is limited to **known L–R pairs in its database** — you cannot freely pair an arbitrary ligand with an arbitrary receptor; (ii) on Xenium a gene must be **on the panel**; (iii) "can return a value" ≠ "meaningful" — sparse/low counts give zeros or noise; (iv) directionality is explicit sender→receiver in most, but **SpatialDM** only at the molecule level (ligand end = sender).

Legend: ✅ native (returns it directly) · ◐ only by aggregation / within a program / post-hoc · ✗ no

| Method | Native output unit | Value for any (senderCT×receiverCT×LR)? | Spatially localized? | Directional? | Significance for that unit? | Key limitation |
|--------|--------------------|:--------------------------------------:|:--------------------:|:------------:|:---------------------------:|----------------|
| **CellChat v2** | CT-pair × LR (3D prob array) | ✅ | ✗ | ✅ | ✅ permutation | Aspatial — one value per whole slide |
| **Tensor-cell2cell** | (context × LR × senderCT × receiverCT) tensor | ✅ | ✗ | ✅ | ✗ (it factorizes) | Aspatial; cell types must be shared across contexts |
| **LIANA+ (dissociated)** | CT-pair × LR | ✅ | ✗ | ✅ | ✅ | Dissociated arm has no space; local metrics have no CT |
| **stLearn / SCTP** | per spot/cell × LR, + CT-pair CCI | ✅ | ✅ | ✅ | ✅ 2-level permutation | No cross-condition test; complexes/competition not modeled |
| **COMMOT** | cell-pair × LR | ✅ (aggregated to cluster level) | ✅ | ✅ | ◐ cluster-level | No cross-condition test; must set distance `T` |
| **CytoSignal** | per receiving cell × LR | ◐ aggregate edges → Circos | ✅ | ✅ | ◐ per receiving cell | No cross-sample program layer |
| **NICHES** | cell-pair × LR | ◐ aggregate per CT cross | ✅ | ✅ | ✗ none built-in | No significance / no downstream genes; pseudoreplication risk |
| **SpatialDM** | per-LR (global) + per-spot (local) | ◐ chord `n_AB`, composition-weighted | ✅ | ◐ molecule-level | ✗ not per-tuple | Pair-independent; cell types only post-hoc |
| **COMPOTES** | per spot × LRI | ◐ only within a program (from snRNA) | ✅ | ◐ post-hoc | ✗ not per-tuple | Visium spot; program-level |
| **ALARMIST** | patch × LRI → per-cell motif loadings | ◐ only within a motif (motif's sender/receiver CTs × defining LRIs) | ✅ | ✅ motif direction | ◐ motif-level ON/OFF + condition tests | Returns motifs, not standalone per-tuple scores; but adds downstream-gene impact + ON/OFF + Bayesian model |

**Read:** methods that return a value for an arbitrary tuple directly = **CellChat, Tensor-cell2cell, LIANA+ (dissociated), stLearn, COMMOT**; methods that only reach it by aggregation / within a program = **CytoSignal, NICHES, SpatialDM, COMPOTES, ALARMIST**. Key divide: the **pooled** tools (CellChat, Tensor-cell2cell, LIANA+ dissociated) give **one number per slide, no location**; the **spatial** tools give a tuple value that **varies by position**; the **program/motif** tools (COMPOTES, ALARMIST) don't expose a standalone tuple score at all — the tuple lives inside each program/motif as "sender/receiver cell types × defining LRIs".

## Capability matrix — what question does it answer?

Legend: ✅ yes (direct) · ◐ partial / indirect / needs downstream / not at this resolution · ✗ no

Questions: **Q1** explicit cell→cell edge for an LR · **Q2** is that edge/location significant · **Q3** for a given cell-type pair, which LRIs · **Q4** spatial hotspot: which LR is strong and where · **Q5** per-single-cell send/receive score for an LR · **Q6** mechanism: contact vs diffusion / direction / range · **Q7** differential across conditions · **Q8** downstream genes / pathways · **Q9** recurrent multi-LR programs/motifs across samples · **Q10** temporal dynamics (up/down) · **Q11** global sender→receiver cell-type pairs · **Q12** spatial gradient

| Method | Q1 | Q2 | Q3 | Q4 | Q5 | Q6 | Q7 | Q8 | Q9 | Q10 | Q11 | Q12 |
|--------|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| **CytoSignal** | ✅ | ◐ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✗ | ◐* | ✅ | ✅ |
| **COMMOT** | ✅ | ◐ | ✅ | ✅ | ✅ | ◐ | ✗ | ✅ | ✗ | ✗ | ✅ | ◐ |
| **SpatialDM** | ✗ | ◐ | ◐ | ✅ | ✅ | ◐ | ✅ | ✗ | ◐ | ◐ | ◐ | ◐ |
| **stLearn / SCTP** | ✗ | ◐ | ✅ | ✅ | ✅ | ◐ | ✗ | ✗ | ✗ | ◐† | ✅ | ◐† |
| **NICHES** | ✅ | ✗ | ✅ | ✅ | ✅ | ◐ | ◐ | ✗ | ◐ | ✗ | ✅ | ◐ |
| **LIANA+** | ✗ | ◐ | ✅ | ✅ | ✅ | ◐ | ✅ | ✅ | ✅ | ✗ | ✅ | ◐ |
| **Tensor-cell2cell** | ✗ | ✗ | ✅ | ✗ | ✗ | ✗ | ✅ | ◐ | ✅ | ◐‡ | ✅ | ✗ |
| **COMPOTES** | ✗ | ✗ | ◐ | ✅ | ✗ | ◐ | ✅ | ◐ | ✅ | ✗ | ✅ | ◐ |
| **CellChat v2** | ✗ | ◐ | ✅ | ✗ | ✗ | ◐ | ◐ | ◐ | ◐ | ✗ | ✅ | ✗ |
| **ALARMIST** | ✗ | ◐ | ◐ | ✅ | ◐ | ◐ | ✅ | ✅ | ✅ | ◐ | ✅ | ◐ |

\* CytoSignal dynamics = VeloCytoSignal, which needs spliced/unspliced counts — **not available on Xenium**.
† stLearn's temporal/gradient answers come from the PSTS trajectory module, a separate analysis from LR scoring.
‡ Tensor-cell2cell's "time" is only time points supplied as the context axis, not predicted dynamics.

**ALARMIST note:** its distinctive strength is holding **Q7 (cross-condition) + Q8 (downstream genes) + Q9 (motifs) simultaneously** — no other method here does. It does *not* output single-cell edges (Q1) or a per-cell single-LR value (Q5); it factorizes into motifs first, projects motif loadings back to each cell, then calls motif ON/OFF states, tests condition enrichment, and quantifies motif→downstream-gene impact.

---

# Part 4 — Signaling-type treatment: paracrine (diffusion) vs juxtacrine (contact) vs autocrine

How each method (a) *decides* whether a pair is diffusible vs contact, (b) models each, (c) computes them **together or separately**, and (d) handles **autocrine**. ALARMIST rows are read from its source (`src/alarmist/core/lri.py`); comparator rows from the papers' Methods.

| Method | Splits contact vs diffusion? (source) | Diffusion / paracrine model | Contact / juxtacrine model | Together or separate | Autocrine handling |
|--------|---------------------------------------|-----------------------------|----------------------------|----------------------|--------------------|
| **CytoSignal** | ✅ from DB label (CellPhoneDB v2 category) | ε-ball radius **200 µm**, Gaussian distance kernel, sender-density-normalized; receptor kept local | **Delaunay** adjacency graph, **uniform** weights; receptor local | **Separate** at every stage — different neighborhood, and separate Circos / ranking / differential per type | Silently **included**; index cell in its own neighborhood, self-weight ≠ 0 (contact w_ii=1; diffusion self-dist≈1e-9). Not a named category |
| **ALARMIST** | ✅ from DB label (`signaling_type`) | within-**patch** cross-cell co-occurrence **count** `count_lig × count_rec` (patch = one scale, e.g. 50 µm) | same within-patch count **minus** same-cell autocrine (`juxta = count_lig×count_rec − auto`) | Each mode is its **own column**, but all share **one patch scale** and are **factorized together** (a motif can mix modes) | **Explicit own category** — a separate `autocrine` column for non-contact same-cell-type pairs; **juxtacrine explicitly subtracts** it. Autocrine = same-cell-type *and* same physical cell co-expressing L+R. **Unique among all methods here** |
| **COMPOTES** | ✅ from DB label (CellChatDB 4 cats → 3 tiers short/med/long) | ligand summed over neighbor **rings** with kernel `w=d/(2πr)`; 1 layer (~155–200 µm) or 2 layers (~255–400 µm); receptor spot-local | "short" = **home spot only** (55 µm, no diffusion) | Different neighborhood radius per pair, but **pooled** into one spot×LRI matrix → one factorization | **Included**; self-spot always in ligand neighborhood (spot-level auto/juxta conflated); **named + retained** at cell-type level (P7 autocrine loops; autocrine as sim ground truth) |
| **SpatialDM** | ✅ from DB label (CellChatDB: 1199 secreted / 421 ECM / 319 contact) | **RBF Gaussian** kernel, bandwidth `l` (coordinate units) + cutoff | **6-nearest-neighbor** hard adjacency (ECM + contact **lumped together**) | Mechanically separate (different weight matrix), but downstream **pooled**; in the actual runs one RBF was often applied to all pairs | Same-**spot** diagonal `w_ii`; **included by default**, optionally **excluded** (`w_ii=0`) only at single-cell resolution. Treated as nuisance auto-correlation, not a class |
| **COMMOT** | ✗ no internal split (type = upstream pair-inclusion filter; they used secreted only) | one collective-OT with per-pair distance limit `T`; cost `φ(D)` continuous within `T`, ∞ beyond | not modeled (could be a small `T`, but not done) | **Together** (one OT couples all pairs; competition across species) | Silently **included** — `D_kk=0 ≤ T` is stored as an explicit zero and is the **lowest-cost** edge. But the solved problem is *entropy-regularized unbalanced* OT (`cot_eps_p=1e-1`) on a cost rescaled by the max within-cutoff distance, so the plan is **not** near-diagonal: at `dis_thr=365 µm` only **~3–10 %** of mass lands on the diagonal, ≈ the heterophilic rate. Self mass is real but minor — see Part 5 |
| **stLearn (SCTP)** | ✗ no DB label; type = user radius (within- vs between-spot) | between-spot: **flat/uniform** average within a radius (default 2× spot diameter); no decay | within-spot (dist 0) — **lumps juxtacrine + autocrine + paracrine** | **Together** (one radius → one score → one test) | Silently **included** (index spot's own L/R contributes; cell-type CCI permits sender==receiver) |
| **NICHES** | ✗ (main text; details in supplement) | plain `L×R` **product**; optional binary spatial-neighbor gate, **no kernel** | same product; only a "spatial nearest-neighbor" restriction | **Together** (one product, distance-agnostic) | Not stated in main text; self cell-type cross sampled + "full connectivity" ⇒ **implicitly included**; explicit self-exclusion deferred to supplement |
| **LIANA+** | ✗ (words paracrine/juxtacrine never appear) | one **RBF kernel**, bandwidth (default Gaussian, l=150, cutoff 0.1); "contact" = just a small bandwidth | no dedicated contact graph | **Together** (one kernel; range is one continuous knob) | **Corrected 2026-08-07:** the *package* default is `set_diag=False` → diagonal **0, self excluded, in every branch** (one builder, `spatial_neighbors.py:53`, `if not set_diag: dist.setdiag(0)` at `:156`). `li.mt.bivariate` has **no** `set_diag` argument at all — it consumes a user-precomputed `obsp` matrix. The "diagonal = 1 for local scores" statement is from the **paper text and tutorial**, not the code default; our runs set it to 1 by following the tutorial |
| **CellChat v2** | ✅ label only (CellChatDB v2 4 cats), used **for filtering, not math** | **no distance** (aspatial); mass-action + Hill on cell-group means | same mass-action model; contact is only a subset label | **Together** (identical score for all; type only subsets the DB) | Sender==receiver **group diagonal retained**, shown as "autocrine/paracrine" in plots; not a computational category |
| **Tensor-cell2cell** | ✗ no signaling type at all (aspatial) | **no distance**; mean(L in sender CT) × (R in receiver CT) | not modeled | **Together** (one score function; all pairs one tensor axis) | Same-cell-type diagonal **silently included**; "autocrine/paracrine" only appear as post-hoc factor interpretation |

### "Together vs separate," decomposed

"Separate" is not one property — a method can separate the two signaling types on some axes and pool them on others. The axis that matters most operationally is **whether paracrine and juxtacrine get a *different distance limit / neighborhood radius*** — and, as the example shows, a method can label and even column-split the two modes while still applying the *same* spatial extent to both.

Legend: ✅ yes · ◐ by design but collapsed in practice / partial · ✗ no

| Method | Assigns a per-pair type label? | **Different distance limit / radius per type?** | Different kernel / weighting per type? | Separate score or column per type? | Separate downstream (rank / significance / differential)? | Aggregated jointly? |
|--------|:------------------------------:|:-----------------------------------------------:|:--------------------------------------:|:----------------------------------:|:---------------------------------------------------------:|:-------------------:|
| **CytoSignal** | ✅ DB | ✅ 200 µm ε-ball **vs** Delaunay adjacency | ✅ Gaussian + density-norm **vs** uniform | ✅ (each pair scored by its type) | ✅ separate Circos / ranking / significance / differential | n/a (no factorization) |
| **ALARMIST** | ✅ DB | ✗ **same patch scale for both** | ✗ **identical count formula**; only autocrine bookkeeping differs | ✅ distinct `juxtacrine` / `paracrine` / `autocrine` columns | ✗ joint factorization (group by mode only post-hoc) | joint |
| **COMPOTES** | ✅ DB → 3 tiers | ✅ 0 / 1 / 2 rings (≈55 / ≈200 / ≈400 µm) | ◐ same `d/2πr` kernel, different extent | ◐ one column per pair (type just sets its radius) | ✗ pooled into one matrix; categories only a post-hoc descriptor | joint |
| **SpatialDM** | ✅ DB | ◐ RBF **vs** 6-NN by design, but one RBF often used for all pairs in practice | ◐ RBF **vs** kNN by design | one score per pair (by its weight matrix) | ✗ pooled (one ranked list / pattern clustering / enrichment) | grouping not by type |
| **COMMOT** | ✗ (type only as upstream pair filter) | ✗ per-pair `T` exists but is **not** set by contact/diffusion class (uniform `T` recommended) | ✗ one OT cost for all | ✗ | ✗ | joint (single OT) |
| **stLearn** | ✗ | ✗ one radius per run; within-spot mode lumps all three | ✗ flat/uniform | ✗ | ✗ | joint |
| **NICHES** | ✗ (main text) | ✗ binary neighbor gate only | ✗ | ✗ | ✗ | joint |
| **LIANA+** | ✗ | ✗ one kernel bandwidth for all ("contact" = just a smaller `l` the user picks) | ✗ | ✗ | ✗ | joint |
| **CellChat v2** | ✅ label, **for filtering only** | ✗ aspatial — no distance at all | ✗ identical mass-action for all types | ◐ only if user subsets the DB by type | ◐ optional post-hoc `signaling.type` filter | joint |
| **Tensor-cell2cell** | ✗ (no signaling type) | ✗ aspatial | ✗ | ✗ | ✗ | joint |

**Reading the distance-limit column (the key one):**
- **Only CytoSignal, COMPOTES, and SpatialDM give paracrine vs juxtacrine a genuinely different spatial extent.** CytoSignal (200 µm Gaussian ε-ball vs Delaunay adjacency) and COMPOTES (0 vs 1 vs 2 neighbor rings) actually use it; SpatialDM defines RBF-vs-6-NN by design but the reported runs often applied a single RBF to all pairs.
- **ALARMIST distinguishes the modes but NOT by distance** — juxtacrine and paracrine share one patch scale and an identical count formula (`count_lig × count_rec − auto`); the only real computational difference ALARMIST introduces is the **autocrine** column (and subtracting autocrine from juxtacrine). Its contact/diffusion split is a DB-inherited label, not a spatial-range distinction.
- **COMMOT, stLearn, NICHES, LIANA+, CellChat, Tensor-cell2cell** apply one distance treatment (or none) to every pair regardless of type.

**Reading the downstream column:** only **CytoSignal** keeps contact vs diffusion separate all the way through (separate Circos, ranking, significance, and differential testing). ALARMIST keeps separate columns but factorizes them jointly (a motif can mix modes). Everyone else pools at or before scoring.

**Autocrine, precisely:**
- **Explicit, first-class, and subtracted from contact:** **ALARMIST only** (own `autocrine` column; `juxtacrine = count_lig×count_rec − auto`).
- **Named/visualized at cell-type level but not a distinct score:** COMPOTES, CellChat.
- **Silently included via self-weight / same-spot / same-cell-type diagonal, not named:** CytoSignal, Tensor-cell2cell, SpatialDM (spot), COMMOT, NICHES.
- **Optionally excluded:** SpatialDM (`single_cell=True` → `w_ii=0`; **our GBM runs use this**), LIANA+ (`set_diag`, default **False = excluded** in the package; our runs overrode it to 1 per the tutorial), stLearn (self is **off** by default — `distance=0` is what turns it *on*).
- **Same physical cell vs same cell type:** spatial per-cell/spot methods (CytoSignal, COMMOT, SpatialDM, stLearn, LIANA-local) mean same *location*; aspatial methods (CellChat, Tensor-cell2cell) mean same *cell type*; ALARMIST requires **both** (a cell of that type co-expressing L and R).

**Counting vs scoring layer:** all comparators produce a continuous score (expression product/mean, Moran's R, OT mass, mass-action probability). **ALARMIST is the only one that counts co-occurrences** — how many ligand+ cells × how many receptor+ cells fall in a patch (binarized), minus same-cell autocrine — rather than a distance-weighted expression score.

---

# Part 5 — The receptor–receptor / homophilic axis (added 2026-08-07)

Part 4 asked *where* a method looks (paracrine vs juxtacrine vs autocrine). This part asks a
different question: **what is in the pair to begin with**, and what happens when the two sides
of the "interaction" are the same molecule. Verified by reading each package's source at the
version we actually ran, plus the papers, plus our own saved outputs.

## The distinction that matters

L–R databases are not all ligand→receptor. Three tiers:

1. **True L–R** — secreted ligand binds a receptor. Directional, asymmetric, mechanistically a
   signalling event. (`GRN→SORT1`, `WNT3→FZD1_LRP6`.)
2. **Heterophilic R–R** — two membrane proteins binding *in trans*. Still a real interaction and
   still directionally meaningful-ish, but it is **adhesion or checkpoint**, not diffusible
   signalling. (`CD99→PILRA`, `NECTIN2→TIGIT`, `DLL3→NOTCH1`, `CADM3→CADM4`.)
3. **Homophilic** — `ligand == receptor`, the *same gene* on both sides. (`CD99–CD99`,
   `NCAM1–NCAM1`, `CDH1–CDH1`, `GJA1–GJA1`, `OCLN–OCLN`, `PCDHA1–PCDHA1`.)

In human CellChatDB v2 tier 3 is **88 / 3,233 rows (2.7 %)** — but **16.4 % of all
`Cell-Cell Contact`**, and all 88 are contact-typed. 53 of the 88 are clustered protocadherins.

**Tier 3 is where the methods break, and it breaks them the same way everywhere:** the score
function is some `f(L) ⊗ g(R)`; setting `L = R = x` collapses it to a function of `x` alone.
The output is then a **marginal expression statistic of one ubiquitous adhesion gene**, carrying
**zero interaction information**, competing on the same ranked list as real signalling.

## Axis 1 × Axis 2 — the master table

"Same cell" = are `L` and `R` scored when both sit in **one** cell/spot.

| Method | Homophilic rows in **its own default** DB | Same-cell scored? | Default? | Any `ligand == receptor` guard? | What homophilic degenerates to |
|---|---|---|---|---|---|
| **ALARMIST** | **88** (CellChatDB v2, mandatory) | **NO — subtracted** | hardwired | ✗ none in `src/` | `n(n−1)`, n = #cells expressing the one gene |
| **CellChat v2** | **88** (it *is* the DB) | **YES**, max weight | no opt-out | ✗ none in 13,956 lines | `Hill(x_s · x_t)`; on diagonal `Hill(x_s²)` |
| **CytoSignal** | **0** (CellPhoneDB v2) → **88** with our CellChatDB | **YES**, by design | no opt-out in contact slot | ✗ none in 183 objects | `x̃ᵢ² + x̃ᵢ·mean_nb(x̃) ≈ 2x̃ᵢ²` |
| **SpatialDM** | **23** (CellChatDB v1) | optional | **`single_cell=False`** → ON | ✗ none | **exactly univariate Moran's I** of that gene |
| **NICHES** | **0** (FANTOM5) → **88** with our CellChatDB | **YES**, 1 self-edge per cell | no opt-out | ✗ none in 1,998 lines | `Sᵢᵢ = xᵢ²` |
| **COMMOT** | **0** (contact rows dropped by the default `Secreted Signaling` filter) | **YES**, cheapest edge | no opt-out | ✗ none | marginals `a == b`; 100 % of mass matches |
| **stLearn** | **21 / 2,293** (connectomeDB2020) | **NO** by default | `distance=None` | ✗ (the one check *worsens* the null) | `xᵢ + x̄_N`; null forbids `l==r` → mis-specified |
| **LIANA+** | **0** in *all 17* bundled resources | package **NO** (`set_diag=False`) | our runs overrode to 1 | ✅ **`remove_self_interactions=True`** — bivariate branch only | local cosine ≡ **1.0 for every cell** |
| **Tensor-cell2cell** | 88 via CellChatDB (23 in v1) | **YES**, CT diagonal | no opt-out | ✗ none | `expression_mean` → **`a_x[A]` exactly** |
| **COMPOTES** | 88 via CellChatDB | **YES**, home spot always in ring | no opt-out | ✗ none | home-spot expression squared |

**Only LIANA+ has a guard**, `remove_self_interactions=True` in
`_spatial_bivariate.py:59` — undocumented in the paper, and it protects only the bivariate
branch. `li.mt.inflow` is unprotected and shows the worst inflation of any method we ran.

## The load-bearing finding: the two axes are *separable*, and Axis 1 dominates

The natural intuition — "the problem is that a cell is interacting with itself" — is **wrong for
most methods**. The homophilic degeneracy survives removing the diagonal:

- **SpatialDM** runs with `single_cell=True` (diagonal zeroed) in our GBM runs, and `CD99_CD99`
  is **still rank 1 of 1,661** (z = 42.15, above `COL4A1_ITGA1_ITGB1` at 37.56). 10 of the top
  25 are homophilic; pooled over 13 cores homophilic pairs are 1.8 % of tests but 14.3 % of the
  top 50 (**7.8×**, hypergeometric).
- **ALARMIST** subtracts same-cell exactly and correctly, and **still** puts 26 % of its total
  count mass on homophilic columns.

Where Axis 2 *does* dominate is the methods that keep an explicit self-edge:

- **NICHES**: self-edges are 24–25 % of edges but carry **75–95 %** of score mass for homophilic
  mechanisms vs 36–45 % for heterophilic. `CD99—CD99` draws **87.8 %** of its entire mass from
  cells paired with themselves; `MPZ—MPZ` 91.6 %.
- **CytoSignal** contact slot: `corr(Lᵢ, Rᵢ) = 0.874` for CD99 — the score is essentially
  smoothed-CD99².

So the correct decomposition is: **Axis 1 (L == R) is the primary defect; Axis 2 amplifies it.**

## Measured inflation, on our own GBM data

| Method | Ranked output | Homophilic share of the top |
|---|---|---|
| **CytoSignal** `contact-Raw_smooth` | 169 pairs | **9 of top 10**, 16 of top 20 (homophilic = 16 % of pairs). `CD99-CD99` rank 12, significant in **42 % of all cells** |
| **NICHES** CellToCellSpatial | 1,088 mech. | **9 of top 15** (homophilic = 2.5 %). `NCAM1—NCAM1` #1, `CD99—CD99` #6. Median rank **55 vs 558** |
| **SpatialDM** | 1,661 pairs | `CD99_CD99` **#1**; 10 of top 25 |
| **LIANA+ `inflow`** | 4,608 features | **26 of top 50 (52 %)** from 5.1 % of features. Ranks 1,2,3,4,6,8 all homophilic |
| **CellChat v2** high-grade | 1,638 edges | 18.7 % of edges, **29.0 % of probability mass** |
| **stLearn** | 526 pairs | 8 of top 20 after permutation; significance rate **2.3× enriched** (p = 8.4e-5) |
| **ALARMIST** | 25,271 columns | **13 of top 20**; 26.05 % of all count mass from 5.3 % of columns |
| **LIANA+ bivariate** | 131 pairs | **none** — 27 self-interactions removed. `CD99^PILRA` rank 59/131 |
| **COMMOT** | 1,468 pair-rows | **none** — 0 CD99 rows, 0 `L==R` rows survive the secreted-only default |

Two methods come out clean, and **both are clean by an upstream filter, not by understanding the
problem**: COMMOT drops all `Cell-Cell Contact` rows by default, LIANA+ bivariate drops `L==R`.
Injecting 7 homophilic rows into COMMOT's real GBM core-14 run puts **5 of the top 10**
homophilic, with `NCAM1-NCAM1` outranking the best genuine secreted pair by **4.9×**.

## Why CD99 specifically

Four compounding reasons, none of them biological:

1. **It is homophilic** — `CD99–CD99` is one of the 88, so it hits the degeneracy above.
2. **It is near-ubiquitous** — 33.5 % of cells positive in our GBM TMA (vs `FPR1` 0.7 %,
   `ANXA1` 3.9 %). Every `f(x)·g(x)` score scales with prevalence.
3. **It gets 4 DB rows** (`CD99–CD99`, `–PILRA`, `–PILRB`, `–CD99L2`), so it enters the ranked
   list four times.
4. **Symmetric pairs get emitted twice** by any method that loops over ordered cell-type pairs.

The most striking external confirmation is **Tensor-cell2cell's own flagship figure**: in
Armingol et al. Fig. 4b (COVID-19 BALF), **`CD99-CD99` is in the top-5 of 5 of 10 factors and is
the single highest-loading pair in 4 of them** (F5 0.213, F6 0.333, F7 0.307). Under
`expression_mean` — their default — the homophilic same-cell-type score is *identically*
`a_x[A]`, the fraction of type-A cells expressing CD99. A published factor is being led by a
detection rate.

## Nobody acknowledges it

Exhaustive full-text scans: **zero** occurrences of `homophilic` / `homotypic` /
`self-interaction` in the CellChat v2, COMMOT, NICHES, stLearn, or Tensor-cell2cell papers.
SpatialDM acknowledges the *diagonal* in one sentence (p8) but never the `L == R` degeneracy —
even though it writes both halves of the proof itself (its eq. 1 with `y = x` is its own eq. 10,
character for character). LIANA+ acknowledges the diagonal in Methods but never documents
`remove_self_interactions`. The ALARMIST manuscript states the `C_juxta` formula
(`methods.tex:26-35`) without noting that `L == R` makes it `n(n−1)`.

## What this means for ALARMIST

**Axis 2: clean, and better than every comparator.** `auto` is subtracted from every contact and
same-cell-type paracrine column (`lri.py:875-905`, `:914-919`); it is *same physical cell*, not
same cell type; no `autocrine` column is ever created for a contact row (verified: 100 % of
homophilic columns are `mode=juxtacrine`). Your instinct that same-cell R–R is meaningless is
right, and the code already acts on it.

**Axis 1: exposed, and quantified.** No `ligand == receptor` check exists anywhere in `src/`.
Measured on saved matrices (GBM / COAD-Xenium / AIS-LUAD):

- Homophilic = **3.6–5.3 % of columns** but **16–26 % of all count mass** (4.4–4.9× over-rep);
  **34–49 % of the juxtacrine block**.
- `n(n−1)` confirmed on **48,532 / 48,532 (100 %)** nonzero GBM entries — the observed value set
  is literally 2, 6, 12, 20, 30, 42, 56, 72, 90, 110.
- **Duplication:** every cross-cell-type homophilic column has a bit-identical mirror —
  **599/599 (GBM), 1,333/1,333 (COAD), 3,108/3,108 (LUAD)**, adding 8.6 / 1.8 / 2.8 % redundant
  mass. Negative control on heterophilic pairs: **0/2000, 1/2000, 0/2000**.
- Highest-mass feature of each run is a gene times itself: `NCAM1–NCAM1` 10.77 % of the GBM
  tensor, `CDH1–CDH1` 5.80 % of COAD, `F11R–F11R` 2.78 % of LUAD.
- In **~45 % of all patches** the single largest "interaction" is homophilic.
- **Downstream:** 24.8 % of BPTF LRI-factor mass is homophilic; **13 of 20 GBM motifs** have a
  homophilic feature as their top loading; motif 0's top feature is `CD99–CD99` (10.4 % of that
  motif's mass).
- **Motif 1 (mGAM) is the least affected of all 20** (12.8 % homophilic, the minimum). Its top 7
  loadings are genuine directed L→R (`APP→SORL1` ×4, `C3→C3AR1`, `C4A→C3AR1`, `LGALS9→CD44`).
  `CD99–CD99` first appears at ranks **8 and 9 — as the duplicate pair**. The published mGAM
  story is not driven by this artifact.

**Minimum defensible fixes, cheapest first:**

1. **Collapse `(A,B)`/`(B,A)` for `ligand == receptor` rows to one column.** Pure win — removes
   1.8–8.6 % redundant mass and stops BPTF seeing a perfectly correlated feature twice. No
   biology is lost because the two columns are bit-identical.
2. **Flag the 88 homophilic rows** and report motif compositions with and without them, or drop
   them behind a flag. This is a *sensitivity analysis*, not a claim that adhesion is
   uninteresting.
3. **State the `n(n−1)` degeneracy in Methods**, next to the `C_juxta` equation.

Note the scope limit: all of the above is tier 3 (`ligand == receptor`). **Tier 2 heterophilic
R–R is untouched by any of these fixes** and is a larger share of the DB — `DLL3→NOTCH1`,
`CADM3→CADM4`, `APP→SORL1`, `NCAM1→PTPRZ1` dominate several methods' top lists and are real
interactions, just adhesion/contact rather than diffusible signalling.

---

# Part 6 — Output data structures (what each method returns), with sender vs receiver

> (This is the "Part 5" requested for the output-object comparison; the number 5 was already taken by the receptor–receptor section above, so it lives here as Part 6.)

**Reference — ALARMIST's outputs** (read from `src/alarmist/core/`):
- **① `lri_factors` = motif × LRI** (`extract_factors`); each LRI column key is `senderCT|receiverCT|ligand|receptor|mode`, so sender/receiver cell type is baked into the column identity.
- **② `cell_loadings` = cell × motif** (`project_cell_loadings`); a **role-agnostic** per-cell motif-usage value — *not* split into a sending vs a receiving score.
- **③ GLM `de_results` = gene × motif × celltype** downstream impact.

Object names/shapes below are read from each package's actual source/API (documented) unless marked inferred.

## Master table

| Method | Primary score object (name · shape) | ① motif/factor decomposition? | ② per-cell/spot side | ③ gene×motif×celltype downstream? |
|--------|-------------------------------------|:-----------------------------:|----------------------|:---------------------------------:|
| **ALARMIST** | `lri_factors` (motif×LRI) + `cell_loadings` (cell×motif) | ✅ BPTF | **role-agnostic** (direction is in the LRI column key) | ✅ `de_results` |
| **CytoSignal** | `@lrscore[[t]]@score` (cell × interaction) | ✗ | **receiver-only** | ◐ gene×interaction (celltype collapsed) |
| **COMMOT** | `obsp['commot-…-L-R']` (cell×cell) + `obsm['…sum-sender']`/`['…sum-receiver']` (cell×LR) | ✗ | **two matrices: sender + receiver** | ◐ `df_impact` (2·LR × gene, **no celltype**) |
| **SpatialDM** | `local_stat['local_I']`(spot×pair, sender) / `['local_I_R']`(spot×pair, receiver) | ✗ | **two matrices: sender + receiver** (significance merged) | ✗ |
| **stLearn** | `obsm['lr_scores']` (spot × LR) | ✗ | **symmetric / undirected** | ✗ |
| **NICHES** | one `assay` per mode (LR-mechanism × observation) | ✗ (user clusters) | **both sides, split across different mode objects** | ✗ |
| **Tensor-cell2cell** | `.factors` = {`Contexts`, `Ligand-Receptor Pairs`, `Sender Cells`, `Receiver Cells`} each ×R | ✅ CP | **cell-type-level, two matrices; no per-cell** | ✗ |
| **LIANA+** | `uns['liana_res']` (directed rows) · `bivariate` lrdata (spot×LR) · factors | ✅ (optional) | dissociated = directed type-pair; local = **symmetric**; factor = separate S/R | ✗ |
| **COMPOTES** | `Contexts`=W (spot×program) + `Ligand-Receptor Pairs`=H (LRpair×program) | ✅ NMF/CP | **spot, receiver-anchored**; post-hoc CT×CT | ✗ |
| **CellChat v2** | `@net$prob` (**sender × receiver × LR** array) | ✅ (outgoing/incoming NMF patterns) | **cell-type-level, two axes; no per-cell** | ✗ |

## Per-cell sender-vs-receiver taxonomy (the crux)

1. **Role-agnostic single value** — **ALARMIST** `cell_loadings` (cell×motif). Direction lives *only* in the LRI column identity (`ligand_ct`=sender, `receptor_ct`=receiver); the cell's loading itself carries no send/receive role.
2. **Two separate per-cell/spot matrices (sender + receiver)** — **COMMOT** `obsm['commot-…-sum-sender']` (row-sum = signal sent) vs `['…sum-receiver']` (col-sum = signal received), each cell×LR; **SpatialDM** `local_stat['local_I']` (own ligand × neighbours' receptor = sender) vs `['local_I_R']` (own receptor × neighbours' ligand = receiver), each spot×pair; **NICHES** `CellToNeighborhood` (cell as sender) vs `NeighborhoodToCell` (cell as receiver), each LR-mech×cell.
3. **Receiver-side only** — **CytoSignal** `@score` (value on the receptor-expressing cell; `res.list` vectors literally named `receiver`); **COMPOTES** spot `Contexts`/W (receptor at the spot, ligand from neighbours).
4. **Undirected / symmetric per-spot** — **stLearn** `obsm['lr_scores']` (ligand↔receptor orientations summed, `base.py:362-364`); **LIANA+** local `bivariate` (default cosine is symmetric).
5. **Directionality only at cell-type level, no per-cell object** — **CellChat** `@net$prob` (sender×receiver×LR array; `outdeg`=sender / `indeg`=receiver centrality); **Tensor-cell2cell** `factors['Sender Cells']` vs `['Receiver Cells']`; **LIANA+** dissociated `liana_res` (`source`=sender / `target`=receiver columns).

Two traps: **stLearn**'s cell-type CCI (`uns['lr_cci_*]`, `per_lr_cci_*`) is drawn as directed arrows but counted with an OR over both orientations (`het.py:260-263`) → `M[i,j]==M[j,i]`, effectively **symmetric**. **SpatialDM** transposes axes: `local_I`/`local_I_R` are (spot×pair) but `local_z`/`selected_spots` are (pair×spot).

## Per-method exact output objects

- **CytoSignal** — `object@lrscore[[type]]@score` (n_cells × n_interactions, **receiver-indexed**); `@score.null` (permuted null); per-cell `pval.mtx`/`padj.mtx`; `@res.list` (per-interaction lists of significant **receiver** barcodes); `inferIntrDEG()` → per-interaction DE-gene sets (receiver cells, celltype-stratified then collapsed); on-the-fly **sender-cluster × receiver-cluster** edge-count matrix (only place a sender cell type is materialized). No factorization, no motif axis.
- **COMMOT** — `obsp['commot-<db>-<L>-<R>']` (n_cells × n_cells directed, one per LR pair; ~800 matrices in the GBM run); `obsm['commot-<db>-sum-sender']` / `['…sum-receiver']` (n_cells × n_LR); `uns['commot_cluster-…']` = {`communication_matrix`, `communication_pvalue`} (senderCT × receiverCT per pair); `obsm['commot_sender_vf-…']`/`['commot_receiver_vf-…']` (vector fields); `df_impact`/`df_deg` (downstream, no celltype axis).
- **SpatialDM** — `uns['global_res']` (n_pairs × stats, whole-slide, **symmetric**); `local_stat['local_I']` (sender) + `['local_I_R']` (receiver), each spot×pair; `local_z`/`selected_spots` (pair×spot, **merged** significance); multi-sample `p_df`/`tf_df`/`zscore_df` (pair×sample); `differential_test` → per-pair `p_val`/`diff_fdr`/`diff`.
- **stLearn** — `obsm['lr_scores']`/`['lr_sig_scores']` (spot×LR, **symmetric**); `obsm['p_vals']`/`['p_adjs']`; `uns['lr_summary']` (per-LR ranking); `uns['lr_cci_<label>']` (CT×CT, **symmetric** despite arrows); `uns['per_lr_cci_<label>']` (dict LR→CT×CT); `uns['lr_go']` (GO of co-expressed genes).
- **NICHES** — `RunNICHES()` returns a named list, one assay per mode: `CellToCell`/`CellToCellSpatial` (LR-mech × directed `SendingCell—ReceivingCell` pairs); `NeighborhoodToCell`/`SystemToCell` (LR-mech × cell, **receiver**); `CellToNeighborhood`/`CellToSystem` (LR-mech × cell, **sender**). "Motifs" = user-run PCA/cluster + `FindAllMarkers` (VectorType archetype × LR-mechanism).
- **Tensor-cell2cell** — input `interaction_tensor.tensor` (context × LR × senderCT × receiverCT); `interaction_tensor.factors` = OrderedDict of 4 DataFrames (`Contexts` n_ctx×R, `Ligand-Receptor Pairs` n_LR×R, `Sender Cells` n_senderCT×R, `Receiver Cells` n_receiverCT×R); `get_factor_specific_ccc_networks()` → per-factor senderCT×receiverCT outer-product. No per-cell object; downstream = GSEA on LR loadings only.
- **LIANA+** — `uns['liana_res']` (directed rows: `source`,`target`,`ligand_complex`,`receptor_complex`,scores); spatial `bivariate` → `lrdata.X` (spot×LR) + `layers['pvals']`/`['cats']`; `to_tensor_c2c` → 4 factor matrices (incl. separate `Sender Cells`/`Receiver Cells`); `lrs_to_views`→MOFA `varm['LFs']` (LR-in-`source&target`-view × factor) + `obsm['X_mofa']` (sample×factor); spatial NMF `obsm['NMF_W']` (spot×factor) + `varm['NMF_H']` (LR×factor).
- **COMPOTES** — `lr_products` (spot×LRpair) + decomposed `ligand_scores` (sender side) / `receptor_scores` (receiver side); `Contexts`/W (spot×program) + `Ligand-Receptor Pairs`/H (LRpair×program, **no celltype/direction in the row label**); `calculate_factor_cell_communication()` → per-program senderCT×receiverCT matrix (post-hoc, needs external scRNA means).
- **CellChat v2** — `@net$prob` + `@net$pval` (senderCT × receiverCT × LR); `@net$count`/`$weight` (CT×CT aggregate); `@netP$prob` (CT×CT×pathway); `subsetCommunication()` (tidy edge list, `source`/`target`); `@netP$centr` (per-pathway `outdeg`=sender / `indeg`=receiver); outgoing vs incoming role matrices (CT×pathway); `netP$pattern$outgoing`/`$incoming` (each a W: CT×K + H: K×pathway, **two separate NMF decompositions**).

## Which of ALARMIST's three outputs each method produces

- **① motif × LRI (direction in the column):** only ALARMIST fuses `senderCT|receiverCT|ligand|receptor|mode` into one column. Factorizers that come closest **split** it: Tensor-cell2cell / LIANA-c2c into 4 separate factor matrices; CellChat into two role-specific NMF pattern pairs; COMPOTES's H is direction-*agnostic* (celltype recovered post-hoc). CytoSignal/SpatialDM/COMMOT/stLearn/NICHES have **no factor/motif axis at all**.
- **② cell × motif (role-agnostic):** ALARMIST's role-agnostic per-cell loading is unusual. The per-cell/spot analogs that exist are either receiver-only (CytoSignal, COMPOTES-spot), symmetric (stLearn, LIANA-local `NMF_W`), or explicitly role-split (COMMOT, SpatialDM, NICHES). The cell-type factorizers (Tensor-c2c, CellChat, LIANA-dissociated) have **no per-cell loading** (their factor "cell" mode is cell-*type*, and split by role).
- **③ gene × motif × celltype (downstream impact):** **no comparator produces this tensor.** Closest partials: **CytoSignal** (gene×interaction, celltype collapsed) and **COMMOT** (`df_impact`, 2·LR×gene, no celltype). Everyone else has no downstream-gene-impact object — see Part 7.

---

# Part 7 — Downstream target-gene impact: signaling → which *other* genes go up/down (NOT pathway annotation)

**The strict question** (ALARMIST's `gene × motif × celltype` GLM): given a signaling program/interaction, which **other, non-ligand/receptor genes** are up- or down-regulated *as a function of the signaling activity*? This must be separated from three easily-confused things:
- **(A)** annotating the LR pair with its curated pathway, or GSEA/GO on the significant **LR pairs / LR-pair loadings / ligand-receptor partner genes** — the common but *wrong* thing.
- **(B)** modelling **downstream target-gene expression as a function of the signaling score** — ALARMIST's question.
- **(C)** linking signaling to downstream **TF / pathway ACTIVITY** via prior-knowledge networks (NicheNet/CARNIVAL/decoupleR/PROGENy) — an indirect middle ground.

**ALARMIST reference** (`src/alarmist/core/glm.py::run_univariate_de_sklearn_by_celltype`): a **per-cell-type Poisson GLM** with gene expression as the **response** and the **motif loading** as the **predictor** → signed **log2FC per gene × motif × celltype**. (Signal is the predictor, gene is the response; direction is reported; resolved per cell type.)

## Verdicts

| Method | Verdict | What it actually does (exact) | Regression direction vs ALARMIST | What's missing vs ALARMIST |
|--------|:-------:|-------------------------------|----------------------------------|----------------------------|
| **ALARMIST** | **B** (ref) | Poisson GLM: `gene ~ motif loading`, per cell type | signal→gene, up/down | — |
| **COMMOT** | **B** (+A) | `communication_deg_detection`: tradeSeq NB-GAM `gene ~ received-signal r`, per LR pair/pathway; `communication_impact`: RF Gini of `r` | **same** (signal→gene), reports up/down | per-**LR**, not per-motif; **no cell-type axis** (pooled over all spots) |
| **CytoSignal** | **B** (+A) | `inferIntrDEG`/`refine_score`: elastic net `LRscore ~ genes + clusters`, per LR pair | **reverse** (gene→signal); signed weight, **no log2FC** | per-**LR**; cell type only a **covariate**, not per-celltype fit; one-sided pre-filter (up-biased) |
| **LIANA+** | **C** | `find_causalnet` (CORNETO/CARNIVAL): LR → deregulated **TFs**; `dc.mt.ulm` on CollecTRI = TF activity | TF-level up/down via **prior network**; condition-driven (DESeq2) | not raw target genes; endpoint is **TF activity**, not gene expression |
| **COMPOTES** | **C** (+A) | decoupleR ULM **TF/signature/PROGENy activity**, then **cosine co-localization** with the program (liana+ bivariate) | correlational co-localization, per **program** | not per-gene; **no up/down**; **no cell-type axis** |
| **Tensor-cell2cell** | **A** | GSEA (`run_gsea`) + PROGENy (`run_mlm`) on the **LR-pair loadings** of each factor | — | **no target-gene model at all** |
| **SpatialDM** | **A** | Fisher's-exact pathway enrichment of the **significant LR pairs** | — | Discussion: *"an open challenge is to identify the downstream targets of LR interactions"* |
| **stLearn** | **A** | `run_lr_go` = GO of the **ligand/receptor partner genes**; hotspot-vs-rest Wilcoxon → GO | — | no signaling→target-gene model |
| **CellChat v2** | **A** | LR → curated pathway aggregation; `net.up`/`net.down` = up/down of the **LR pair's own** strength | — | Limitations: *"how to better validate … their downstream gene outputs remains to be answered"* |
| **NICHES** | **none** | `FindAllMarkers` on the signaling matrix returns differentially-active **LR mechanisms**, not genes | — | no downstream analysis; not even pathway annotation |

## The two that genuinely answer your question (B), and how they differ from ALARMIST

- **COMMOT — closest in spirit, same regression direction.** `commot.tl.communication_deg_detection` puts the per-spot **received signal `r` = Σⱼ S[j,i]** as the predictor into a **tradeSeq NB-GAM** and finds genes whose expression rises/falls with `r` ("positive/negative DE gene"); `commot.tl.communication_impact` adds a random-forest Gini importance of `r`. Documented up/down: *"as the received midkine signaling increases, COL1A1 expression increases while S100G decreases"*; *"increasing [WNT] signal → higher KRT15/KRT5, lower LOR/FLG"* (Cang 2023, Results; Methods "Downstream gene analysis" lines 1703-1718). **Differences from ALARMIST:** per-**LR-pair/pathway** (no motif/factorization), and the tradeSeq fit is **global across all spots — no cell-type stratification**, so it yields `gene × LR`, not `gene × motif × celltype`. It also validates against scSeqComm's curated target-gene DB (benchmark only).
- **CytoSignal — same question, but reversed regression.** `inferIntrDEG` fits, per LR pair, an elastic net where **the LRscore is the RESPONSE and gene expression + cluster labels are the PREDICTORS** — the paper is explicit: *"a sparse regression analysis using the LRscore from CytoSignal as the response variable … to predict the LRscore from gene expression and cluster labels"* (Liu 2026, Methods "Identifying signaling-associated genes"; source `refine_score`: `y <- scores[, test_intr]`). So it finds **"genes predictive of the signaling,"** not "expression change caused by the signaling." It gives signed coefficients (positive = up-associated) but **no log2FC**, is **per-LR** (no motif), enters **cell type as a single covariate** (not a per-celltype fit), and its Wilcoxon candidate pre-filter is one-sided toward up genes. The GO/GOrilla step layered on top is (A).

## The two "indirect" ones (C) — downstream **activity**, not raw target genes

- **LIANA+**: `find_causalnet` (CORNETO, modified CARNIVAL ILP) traces a sign-coherent OmniPath subnetwork from a chosen LR pair down to **deregulated TFs** (TF activity via `decoupler ULM` on CollecTRI regulons, seeded by PyDESeq2 Wald stats). It reports up/down **TFs** (e.g. SMAD1/3 up, FOXO3 down in myeloid cells), and log2FC only for the **ligand/receptor genes** — never a de-novo regression of downstream target-gene expression on a signaling score. Its official tutorial states the pipeline *"terminates at TF activity prediction; it does not connect TF activity back to target gene expression."* Its `Fig 5d` "pathway enrichment of NMF ligand–receptor loadings" is (A).
- **COMPOTES**: per program, infers **TF/signature/PROGENy activity** (decoupleR ULM on OmniPath/CollecTRI; hallmark & IFN-γ/exhaustion/TAM signatures) and then **cosine-co-localizes** that activity map with the program's spatial activity (liana+ bivariate). This is correlational co-localization at the **activity** level, **per program, no cell-type axis, no per-gene up/down** — it answers "which TF/pathway activities co-localize with this program," not "which genes does this program up/down-regulate."

## The four that only annotate the LR pair (A) — exactly what you said you don't want

**Tensor-cell2cell** (GSEA + PROGENy on LR-pair *loadings*), **SpatialDM** (Fisher enrichment of significant LR pairs; authors call downstream targets an *open problem*), **stLearn** (`run_lr_go` = GO of the L/R *partner* genes), **CellChat v2** (LR → curated pathway; `net.up/down` is the LR pair's *own* strength; authors call downstream gene outputs *unsolved*). **NICHES** does even less — its markers are LR mechanisms, not genes.

## Bottom line

- **Only COMMOT and CytoSignal actually ask "which other genes go up/down with the signaling" (B).** COMMOT matches ALARMIST's regression direction (signal→gene, up/down) but is **per-LR-pair and pooled over cell types**; CytoSignal **reverses** it (genes→signal) and keeps cell type only as a covariate.
- **LIANA+ and COMPOTES reach "downstream" only as TF/pathway ACTIVITY via prior networks (C)** — not raw target genes.
- **Tensor-cell2cell, SpatialDM, stLearn, CellChat give only pathway annotation (A); NICHES gives none.**
- **ALARMIST's exact deliverable — target-gene up/down as a function of a signaling *program/motif*, resolved *per cell type* (`gene × motif × celltype`) — is matched by no method here.** COMMOT is the nearest, but at per-LR-pair granularity and without the cell-type axis.
