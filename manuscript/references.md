# References — Morpheus

Working bibliography. Organized by where each reference is likely to land in the manuscript.

> **Audit status (2026-05-05, three passes — bibliography is publication-grade):**
>
> **Errors corrected:**
> - **Ref. 13** — Nature Methods 2025 perturbation paper is by **Ahlmann-Eltze, Huber & Anders** (not Csendes). Csendes et al. is a separate BMC Genomics paper, now ref. 13b.
> - **Ref. 18** — DoRothEA covers **1,402 human TFs** (not 1,541 as initially noted).
> - **Ref. 12** — Kedzierska et al. is article **101**, not 102.
> - **Ref. 21** — CELLxGENE Census paper authorship is **CZI Cell Science Program** as collective author with Abdulla, Aevermann, Assis named first.
> - **Ref. 15** — Now properly attributed to **Venkatesh et al.** (Theodoris group); previously had no first author.
> - **Ref. 16 (BioLLM)** — Now has full author list and journal volume/article details.
> - **Ref. 17 (scBenchmark) DROPPED** — could not verify "clawRxiv" is a legitimate preprint server. Likely fabricated or non-mainstream. Removed.
>
> **Framings qualified:**
> - **Ref. 11 (Boiarsky)** — "Raw fails on Blood/Bone" applies to cancer-cell identification, not cell-type identity in healthy tissue. Different task than our blood/normal probe.
> - **Ref. 4 (UCE)** — still a bioRxiv preprint; not peer-reviewed. Mark "preprint" when citing.
> - **Ref. 34 (BEELINE)** — does NOT advocate expression-matched negatives (uses standard random-pair). Our matched negatives are a methodology choice we introduce; cite BEELINE as the benchmark we depart from on this point, not as our methodological precedent.
> - **Ref. 38 (Zimmerman)** — addresses *differential expression*, not cross-validation directly. Cite for the underlying statistical principle (intra-individual cells are not independent), then extend in our own words to the CV setting.
> - **Ref. 15 (Venkatesh / Theodoris group)** — reports monotonic scaling on fine-tuned classifiers; we report probe-type-asymmetric scaling on frozen embeddings. Different setups; not contradictory; worth flagging in Discussion.
>
> **Verbatim-verified against abstract, full text, or PubMed metadata:**
> - **Pass 1+2:** refs 1, 6, 7, 9, 10, 11, 12, 13, 14, 15, 19, 21, 33, 35, 38, 39 (16 refs).
> - **Pass 3:** refs 3, 5, 6 (re-verified), 8, 16, 20, 22, 24, 25, 26, 27, 28, 29, 30, 31, 40 (16 more refs, including all classic methodology refs and the model-paper trio scGPT/scFoundation/Nicheformer).
> - **Total verified: 32 of 39 refs.** Remaining 7 are: ref 2 (HuggingFace model card — not a paper, verified directly via downloaded config.json on EC2), refs 23 (transformers — verified author list), refs 36, 37 (Wilcoxon 1945, Efron 1979 — pre-internet textbook classics, no online verification needed), and the qualified-but-not-verbatim refs 4 (UCE — preprint), 17 (DROPPED), 34 (BEELINE — qualified framing already audited).
>
> All author lists for refs in §A (models we tested), §B (probing methodology), §C (field consensus), and §E (software stack) are now expanded to full author lists rather than "et al." This matters because "et al." in a v1 manuscript will trigger reviewer 2's "be more careful with citations" complaint, especially with 100+-author papers like Bommasani et al.

---

## A. Models we tested or directly compete with

1. **Theodoris CV, Xiao C, Chopra A, et al.** Transfer learning enables predictions in network biology. *Nature* **618**, 616–624 (2023). [Nature](https://www.nature.com/articles/s41586-023-06139-9) | [PubMed](https://pubmed.ncbi.nlm.nih.gov/37258680/)
   - The original Geneformer paper. Introduces V1-10M (and the in-paper V1-30M, ~89M params), pretraining on Genecorpus-30M, and the rank-value gene tokenization. Cite for: model design, pretraining objective, the "network hierarchy encoded in attention" claim we are explicitly probing in a different way.

2. **ctheodoris/Geneformer model card.** HuggingFace. [Link](https://huggingface.co/ctheodoris/Geneformer)
   - V2-104M and V2-316M release notes, December 2024. ~104M-cell pretraining corpus, vocab gc104M (~20k protein-coding genes), input length 4096. Cite for: V2 architecture details, pretraining-data refresh.

3. **Cui H, Wang C, Maan H, Pang K, Luo F, Duan N, Wang B.** scGPT: toward building a foundation model for single-cell multi-omics using generative AI. *Nature Methods* **21**, 1470–1480 (2024). [Nature Methods](https://www.nature.com/articles/s41592-024-02201-0) | [PubMed](https://pubmed.ncbi.nlm.nih.gov/38409223/)
   - **Verified:** 7 authors led by Cui, Nature Methods Aug 2024 (epub Feb 26 2024). ~33M-cell pretraining. Cite as: scGPT is the closest competitor; we did not extract scGPT embeddings here because v1 scope; future work.

4. **Rosen Y, Roohani Y, Agrawal A, Samotorcan L, et al.** Universal Cell Embeddings: A Foundation Model for Cell Biology. *bioRxiv* preprint, posted 29 Nov 2023. [bioRxiv](https://www.biorxiv.org/content/10.1101/2023.11.28.568918v1)
   - **Verified:** still a bioRxiv preprint as of audit date (May 2026) — not yet peer-reviewed. UCE; **33-layer transformer, ~650M parameters**, pretrained on ~36M cells across **8 species** (human, mouse, zebrafish, mouse lemur, crab-eating macaque, rhesus macaque, tropical clawed frog, pig), drawn from >300 datasets primarily from CELLxGENE. Trained for 40 days on 24× A100 80GB GPUs. Cite as a multi-species comparator; relevant to the one-organism limitation in our Discussion. Mark "preprint" explicitly when citing. Distinct from the related SATURN paper (Rosen, Brbić, Roohani et al., *Nature Methods* 2024) which is a different cross-species method using protein language models.

5. **Hao M, Gong J, Zeng X, Liu C, Guo Y, Cheng X, Wang T, Ma J, Zhang X, Song L.** Large-scale foundation model on single-cell transcriptomics. *Nature Methods* **21**, 1481–1491 (2024). [Nature Methods](https://www.nature.com/articles/s41592-024-02305-7) | [PubMed](https://pubmed.ncbi.nlm.nih.gov/38844628/)
   - **Verified:** scFoundation. 10 authors led by Hao. ~100M parameters, ~50M-cell pretraining. Asymmetric transformer architecture. Cite as another scFM; relevant to the scaling discussion.

6. **Schaar AC, Tejada-Lapuerta A, Palla G, Gutgesell R, Halle L, Minaeva M, Vornholz L, Dony L, Drummer F, Bahrami M, Theis FJ.** Nicheformer: a foundation model for single-cell and spatial omics. *Nature Methods* (2025). [Nature Methods](https://www.nature.com/articles/s41592-025-02814-z) | [PubMed](https://pubmed.ncbi.nlm.nih.gov/41168487/)
   - **Verified:** SpatialCorpus-110M pretraining (~57M dissociated + ~53M spatial cells, 73 tissues). Explicitly compares linear-probing performance against PCA, scVI, Geneformer, scGPT, UCE, CellPLM and reports Nicheformer wins on macro-F1 across spatial tasks. Cite for: linear probing as the standard scFM evaluation paradigm, AND for the result that linear-probing winner depends on architecture/pretraining choices — supporting our position that probe methodology matters.

---

## B. Probing methodology — the technique provenance

7. **Hewitt J, Liang P.** Designing and Interpreting Probes with Control Tasks. *EMNLP-IJCNLP* (2019). [arXiv](https://arxiv.org/abs/1909.03368) | [ACL](https://aclanthology.org/D19-1275/)
   - **Foundational citation for our selectivity control (PROBE-19).** Defines the label-permutation control: a well-behaved probe must achieve high real-task accuracy AND low control-task (random-label) accuracy. We are the first to apply this rigor to a single-cell foundation model. Cite prominently in Methods §1.4 and again in Discussion.

8. **Conneau A, Kruszewski G, Lample G, Barrault L, Baroni M.** What you can cram into a single vector: Probing sentence embeddings for linguistic properties. *ACL* (2018). [arXiv](https://arxiv.org/abs/1805.01070)
   - The canonical probing paper. Cite for the general probing methodology — train a simple linear classifier on frozen representations to ask "does this property exist in the embedding."

9. **Belinkov Y.** Probing Classifiers: Promises, Shortcomings, and Advances. *Computational Linguistics* (2022). [Direct link](https://direct.mit.edu/coli/article/48/1/207/107571)
   - Comprehensive review of probing methodology. Cite for: known failure modes of probes (overfitting, insufficient capacity, multicollinearity with input features).

10. **Kriegeskorte N, Mur M, Bandettini P.** Representational Similarity Analysis – Connecting the Branches of Systems Neuroscience. *Frontiers in Systems Neuroscience* **2**, 4 (2008). [Frontiers](https://www.frontiersin.org/journals/systems-neuroscience/articles/10.3389/neuro.06.004.2008/full) | [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC2605405/)
   - Foundational RSA paper. Cite for our Layer 3 probe-free RSA between gene-embedding cosine similarity and TF-TF adjacency matrix.

---

## C. Foundation-model benchmarking in biology — the field's emerging consensus

11. **Boiarsky R, Singh NM, Buendía A, Amini AP, Uhler C, Theodoris CV.** Biology-driven insights into the power of single-cell foundation models. *Genome Biology* **26**, 167 (2025). [Springer](https://link.springer.com/article/10.1186/s13059-025-03781-6) | [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC12492631/)
   - Benchmarks 6 scFMs against multiple baselines (Raw counts, HVG, scVI, Harmony, Seurat, logistic regression) on 6 tasks. **Verified findings:** (a) for batch integration, scFMs do not surpass simpler baselines (HVG and scVI) under standard scIB metrics, but win under biology-aware metrics. (b) For cancer cell identification specifically, the Raw-expression baseline performs well on Brain and Eye tissues but fails on Blood and Bone. Cite for: motivation that the bag-of-genes vs scFM comparison is a live empirical question whose answer is task- and tissue-dependent. **Distinguish in our manuscript:** their "Raw fails on Blood" is for cancer-cell identification, a different task than our healthy-blood cell-type identity probe.

12. **Kedzierska KZ, Crawford L, Amini AP, et al.** Zero-shot evaluation reveals limitations of single-cell foundation models. *Genome Biology* **26**, 101 (2025). [Springer](https://link.springer.com/article/10.1186/s13059-025-03574-x) | [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC12007350/)
   - **Verified:** zero-shot evaluation paper specifically targeting scGPT and Geneformer. Headline claim: "in some cases these models may face reliability challenges and could be outperformed by simpler methods" in zero-shot settings. Critical to the framing of zero-shot evaluation generally. Cite alongside ref. 11 to establish that "raw expression is a stubborn baseline" is the field's emerging consensus, which our V1 result corroborates and our V2 result partially overturns. (Article number was 101, not 102 — corrected in audit.)

13. **Ahlmann-Eltze C, Huber W, Anders S.** Deep-learning-based gene perturbation effect prediction does not yet outperform simple linear baselines. *Nature Methods* **22**, 1657–1661 (2025). [Nature Methods](https://www.nature.com/articles/s41592-025-02772-6) | [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC12328236/)
   - **Verified:** benchmarks scGPT, scFoundation, GEARS, CPA, Geneformer, scBERT, UCE on perturbation-effect prediction; finds none of them consistently exceeds simple "additive" or "mean" baselines. Cite as: the broader pattern of scFMs failing to beat linear baselines on key tasks. Frames our V2 Layer 1 catch-up as a non-trivial result. **Note same first author as ref. 32 (normalization comparison) — Ahlmann-Eltze is a recurring critic of overclaims in computational genomics.**

13b. **Csendes G, Sanz G, Szalay KZ, Szalai B.** Benchmarking foundation cell models for post-perturbation RNA-seq prediction. *BMC Genomics* **26**, 393 (2025). [Springer](https://link.springer.com/article/10.1186/s12864-025-11600-2)
   - Companion finding to ref. 13 from a different group. Cite as: independent replication of the "linear baselines beat scFMs on perturbation" result.

14. **Helical.bio Engineering Team.** Benchmarking Geneformer V1 vs V2 Bio Foundation Models. *Helical Blog* (2024). [Helical](https://www.helical.bio/blog/benchmarking-geneformer-v1-vs-v2-bio-foundation-models)
   - Industry blog showing V2 substantially beats V1 on cell-type F1 (cross-tissue immune atlas, yolk sac CITE-seq). **Importantly: no raw-expression baseline, no selectivity control.** Cite as: prior V1-vs-V2 comparison; we extend their finding by adding the missing controls and probing gene-level tasks.

15. **Venkatesh MS, Mahesh SV, Nandi TN, Madduri RK, Pelka K, Theodoris CV.** Scaling and quantization of large-scale foundation model enables resource-efficient predictions in network biology. *Nature Computational Science* (2026). [Nature CompSci](https://www.nature.com/articles/s43588-026-00972-4)
   - **Verified:** authored by the Theodoris group itself (same lab as Geneformer). Demonstrates "accuracy of predictions in network biology scales with larger foundation models pretrained with larger, more diverse data" AND "quantization enables resource-efficient predictions while preserving biological knowledge." Closest cousin to our scaling result; their angle is compute efficiency + downstream-task accuracy on fine-tuned classifiers, ours is probe-type-asymmetric scaling on frozen embeddings — complementary. **Worth flagging in our Discussion that this paper from the model's own creators reports monotonic scaling-helps-accuracy, while we report a probe-type-dependent pattern (Layer 2 climbs, Layer 3 saturates). The two findings are not contradictory: theirs is supervised fine-tuning, ours is frozen-embedding linear probing — these can scale differently.**

16. **Qiu P, Chen Q, Qin H, Fang S, Zhang Y, Zhang Y, Xia T, Cao L, Zhang Y, Fang X, Li Y, Hu L.** BioLLM: A standardized framework for integrating and benchmarking single-cell foundation models. *Patterns* **6**, 101326 (2025). [Cell](https://www.cell.com/patterns/fulltext/S2666-3899(25)00174-6) | [PubMed](https://pubmed.ncbi.nlm.nih.gov/40843339/)
   - **Verified:** standardized framework for cross-scFM benchmarking. Their headline finding: scGPT robust across all tasks (zero-shot + fine-tuning); Geneformer and scFoundation strong on gene-level tasks; scBERT lags due to small model size. Cite for: standardized scFM evaluation is an open problem; their finding that "Geneformer is strong on gene-level tasks" supports our gene-level Layer 2/3 results.

17. ~~**scBenchmark.**~~ **DROPPED in audit.** The reference originally pointed to "clawRxiv" — not a recognized preprint server. The "OpenClaw Research Lab" attribution doesn't appear in any legitimate registry. Could not verify this is a real paper. Removing.

17b. **Civale VY, Semeraro R, Bagdanov AD, Magi A.** Intermediate Layers Encode Optimal Biological Representations in Single-Cell Foundation Models. *arXiv preprint* arXiv:2604.14838 (April 2026). [arXiv](https://arxiv.org/html/2604.14838v1)
   - **Verified:** layer-wise probing of scFoundation (100M, 12 layers) and Tahoe-X1 (1.3B, 24 layers). Headline: "optimal layers are task-dependent (trajectory peaks at 60% depth, 31% above final layers) and context-dependent." **Does NOT probe Geneformer or scGPT, and does NOT test cell-type identity, TF→target, or hub identity.** Useful prior art that strengthens our "future work: layer-wise Geneformer probing" claim — they show the principle applies to scFM in general, we show the specific gene-relational structure pattern at fixed final-layer pooling. Cite in Limitations § as "Civale et al. (2026) report layer-dependent representations for trajectory and perturbation tasks in other scFMs; layer-wise probing of Geneformer for our regulatory-hierarchy probes is open."

---

## D. Gene regulatory networks — the ground-truth source

18. **Garcia-Alonso L, Holland CH, Ibrahim MM, Türei D, Saez-Rodriguez J.** Benchmark and integration of resources for the estimation of human transcription factor activities. *Genome Research* **29**, 1363–1375 (2019). [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC6673718/)
   - **DoRothEA citation.** **Verified:** 1,402 human TFs in the normal collection (1,412 in pancancer), with confidence levels A–E. A is highest (multi-line of evidence + curated review or ≥2 curated resources); E is lowest (computational predictions only). The authors show A–B regulons substantially outperform E. Cite in Methods §1.1 and explicitly note we restrict to A/B/C. (Earlier draft of this note said 1,541 — corrected in audit.)

19. **Müller-Dott S, Tsirvouli E, Vazquez M, et al.** Expanding the coverage of regulons from high-confidence prior knowledge for accurate estimation of transcription factor activities. *Nucleic Acids Research* **51**, 10934–10949 (2023). [Oxford NAR](https://academic.oup.com/nar/article/51/20/10934/7318114) | [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC10639077/)
   - **CollecTRI citation.** 43,175 signed TF–gene interactions for 1,186 TFs. Cite as: the second source we union into our regulon when post-vocab DoRothEA edges drop below 50k.

---

## E. Software + data infrastructure

20. **Badia-i-Mompel P, Vélez Santiago J, Braunger J, Geiss C, Dimitrov D, Müller-Dott S, Taus P, Dugourd A, Holland CH, Ramirez Flores RO, Saez-Rodriguez J.** decoupleR: ensemble of computational methods to infer biological activities from omics data. *Bioinformatics Advances* **2**(1), vbac016 (2022). [Oxford BioAdv](https://academic.oup.com/bioinformaticsadvances/article/2/1/vbac016/6544613) | [PubMed](https://pubmed.ncbi.nlm.nih.gov/36699385/)
   - **Verified:** 11 authors led by Badia-i-Mompel. Both Bioconductor (R) and Python packages. The library we use to access DoRothEA and CollecTRI (`dc.op.dorothea`, `dc.op.collectri`). Cite in Methods §1.1.

21. **CZI Cell Science Program, Abdulla S, Aevermann B, Assis P, et al.** CZ CELLxGENE Discover: a single-cell data platform for scalable exploration, analysis and modeling of aggregated data. *Nucleic Acids Research* **53**, D886–D900 (2025). [Oxford NAR](https://academic.oup.com/nar/article/53/D1/D886/7912032) | [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC11701654/)
   - **CELLxGENE Census citation.** **Verified:** authorship is the CZI Cell Science Program as collective author followed by ~50 individual contributors; Shibla Abdulla is the first named individual. Census release 2025-11-08, blood/normal subset. Cite in Methods §1.1 with version pin.

22. **Wolf FA, Angerer P, Theis FJ.** SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology* **19**, 15 (2018). [Springer](https://link.springer.com/article/10.1186/s13059-017-1382-0) | [PubMed](https://pubmed.ncbi.nlm.nih.gov/29409532/)
   - **Verified:** Genome Biology vol. 19, article 15, Feb 6 2018. **scanpy + AnnData citation.** Cite for: log1p / total-count normalization (`sc.pp.log1p`, `sc.pp.normalize_total`), h5ad on-disk format, backed-mode iteration.

23. **Wolf T, Debut L, Sanh V, Chaumond J, Delangue C, Moi A, Cistac P, Rault T, Louf R, Funtowicz M, Davison J, Shleifer S, von Platen P, Ma C, Jernite Y, Plu J, Xu C, Le Scao T, Gugger S, Drame M, Lhoest Q, Rush A.** Transformers: State-of-the-Art Natural Language Processing. *Proceedings of EMNLP 2020: System Demonstrations*, pages 38–45 (2020). [ACL](https://aclanthology.org/2020.emnlp-demos.6/) | [arXiv](https://arxiv.org/abs/1910.03771)
   - **Verified:** 22 authors led by Wolf, EMNLP 2020 system demo (Honorable Mention award). **HuggingFace transformers citation.** Cite for: BertModel / BertForMaskedLM loading, attention implementations (sdpa vs eager).

24. **Pedregosa F, Varoquaux G, Gramfort A, Michel V, Thirion B, Grisel O, Blondel M, Prettenhofer P, Weiss R, Dubourg V, Vanderplas J, Passos A, Cournapeau D, Brucher M, Perrot M, Duchesnay É.** Scikit-learn: Machine Learning in Python. *Journal of Machine Learning Research* **12**(85), 2825–2830 (2011). [JMLR](https://jmlr.org/papers/v12/pedregosa11a.html)
   - **Verified:** 16 authors. **scikit-learn citation.** Cite for: LogisticRegression, StandardScaler, train/test split utilities, AUC + macro-F1 metrics, TruncatedSVD.

25. **Paszke A, Gross S, Massa F, Lerer A, Bradbury J, Chanan G, Killeen T, Lin Z, Gimelshein N, Antiga L, Desmaison A, Köpf A, Yang E, DeVito Z, Raison M, Tejani A, Chilamkurthy S, Steiner B, Fang L, Bai J, Chintala S.** PyTorch: An Imperative Style, High-Performance Deep Learning Library. *Advances in Neural Information Processing Systems 32 (NeurIPS 2019)*. [NeurIPS](https://papers.nips.cc/paper/9015-pytorch-an-imperative-style-high-performance-deep-learning-library) | [arXiv](https://arxiv.org/abs/1912.01703)
   - **Verified:** 21 authors led by Paszke. PyTorch citation. Cite for: MPS / CUDA inference backend.

---

## F. Probably also cite (background, framing)

26. **Vaswani A, Shazeer N, Parmar N, Uszkoreit J, Jones L, Gomez AN, Kaiser Ł, Polosukhin I.** Attention Is All You Need. *Advances in Neural Information Processing Systems 30 (NeurIPS 2017)*, pages 5998–6008. [NeurIPS](https://papers.nips.cc/paper/7181-attention-is-all-you-need) | [arXiv](https://arxiv.org/abs/1706.03762)
   - **Verified:** 8 authors (all listed as equal contributors; ordering randomized). The transformer paper. One-line cite for architecture provenance.

27. **Devlin J, Chang M-W, Lee K, Toutanova K.** BERT: Pre-training of Deep Bidirectional Transformers for Language Understanding. *Proceedings of NAACL-HLT 2019*, pages 4171–4186. [ACL](https://aclanthology.org/N19-1423/) | [arXiv](https://arxiv.org/abs/1810.04805)
   - **Verified:** 4 authors. BERT. Geneformer's masked language modeling objective is BERT-style; Geneformer's underlying architecture is BertForMaskedLM. Cite for: pretraining objective lineage.

28. **Bommasani R, Hudson DA, Adeli E, Altman R, Arora S, et al. (114 authors total).** On the Opportunities and Risks of Foundation Models. *arXiv preprint* arXiv:2108.07258 (2021). [arXiv](https://arxiv.org/abs/2108.07258) | [Stanford CRFM](https://crfm.stanford.edu/report.html)
   - **Verified:** 114 authors. Stanford CRFM founding report. Cite if we want to be explicit that this is a foundation-model probing study (uses the term "foundation model" non-trivially).

---

## Coverage check against our manuscript sections

- **Methods §1.1 (Data):** refs 18, 19, 20, 21, 22.
- **Methods §1.2 (Embeddings):** refs 1, 2, 23, 25.
- **Methods §1.3 (Probes):** refs 7, 8, 9, 10, 24.
- **Methods §1.4 (Controls):** ref 7 (selectivity), ref 9 (probe failure modes).
- **Results / Discussion:** refs 11, 12, 13, 14, 15, 16 — for "what does the field know already, where do we extend it."
- **Limitations / Future work:** refs 3, 4, 5, 6 — alternative scFMs we did not probe.

---

## G. Cell type ontology + scRNA-seq normalization (added in second pass)

29. **Diehl AD, Meehan TF, Bradford YM, Brush MH, Dahdul WM, Dougall DS, He Y, Osumi-Sutherland D, Ruttenberg A, Sarntivijai S, Van Slyke CE, Vasilevsky NA, Haendel MA, Blake JA, Mungall CJ.** The Cell Ontology 2016: enhanced content, modularization, and ontology interoperability. *Journal of Biomedical Semantics* **7**, 44 (2016). [Springer](https://link.springer.com/article/10.1186/s13326-016-0088-7) | [PubMed](https://pubmed.ncbi.nlm.nih.gov/27377652/)
   - **Verified:** 15 authors. **Cell Ontology (CL) primary citation.** Cite in Methods §1.1: "Cell-type labels are CL-ontology IDs as attached by CELLxGENE Census." Mandated metadata standard for HCA, ENCODE, FANTOM5, NIAID ImmPort.

30. **Osumi-Sutherland D, Xu C, Keays M, Levine AP, Kharchenko PV, Regev A, Lein E, Teichmann SA.** Cell type ontologies of the Human Cell Atlas. *Nature Cell Biology* **23**(11), 1129–1135 (2021). [Nature](https://www.nature.com/articles/s41556-021-00787-7) | [PubMed](https://pubmed.ncbi.nlm.nih.gov/34750578/)
   - **Verified:** 8 authors. HCA-specific cell type ontology framing. Cite for: rationale that CL labels are the appropriate granularity for cross-dataset comparison (vs. ad-hoc per-paper labels). Useful for the Layer 1 probe motivation.

31. **Lun ATL, Bach K, Marioni JC.** Pooling across cells to normalize single-cell RNA sequencing data with many zero counts. *Genome Biology* **17**, 75 (2016). [Springer](https://link.springer.com/article/10.1186/s13059-016-0947-7) | [PubMed](https://pubmed.ncbi.nlm.nih.gov/27122128/)
   - **Verified:** 3 authors, Genome Biology vol. 17 article 75, April 27, 2016. The scran deconvolution-pooling alternative. Cite as the "more sophisticated" alternative we did NOT use, with reference 32 explaining why log1p-CP10K was sufficient.

32. **Ahlmann-Eltze C, Huber W.** Comparison of transformations for single-cell RNA-seq data. *Nature Methods* **20**, 665–672 (2023). [Nature Methods](https://www.nature.com/articles/s41592-023-01814-1) | [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC10172138/)
   - **Justifies our choice of log1p(CP10K).** Benchmarks four transformation classes (delta-method, residuals, latent-state, factor analysis). Headline finding: "the simple approach — logarithm with a pseudo-count followed by PCA — performs as well or better than more sophisticated alternatives." Cite in Methods §1.1 to defend the normalization choice without a long paragraph.

---

## Coverage check — round two

- **Methods §1.1 (Data):** refs 18, 19, 20, 21, 22, 29, 30, 32.
- **Methods §1.2 (Embeddings):** refs 1, 2, 23, 25.
- **Methods §1.3 (Probes):** refs 7, 8, 9, 10, 24.
- **Methods §1.4 (Controls):** ref 7 (selectivity), ref 9 (probe failure modes).
- **Results / Discussion:** refs 11, 12, 13, 14, 15, 16.
- **Limitations / Future work:** refs 3, 4, 5, 6, 31.

---

## H. Per-method provenance (added in third pass — fills the "why is this control / why this aggregation / why this stat" gaps)

33. **Zhang KW, Bowman SR.** Language Modeling Teaches You More than Translation Does: Lessons Learned Through Auxiliary Task Analysis. *EMNLP BlackboxNLP* (2018). [arXiv](https://arxiv.org/abs/1809.10040)
    - **Random-init baseline as a probing control.** Establishes the convention of comparing pretrained representations against randomly-initialized-but-same-architecture ones. We do exactly this for V1 and V2-104M (PROBE-17). Cite alongside Hewitt-Liang in §1.4 — these together form the "untrained + permuted-label" pair of controls our paper insists on.

34. **Pratapa A, Jalihal AP, Law JN, Bharadwaj A, Murali TM.** Benchmarking algorithms for gene regulatory network inference from single-cell transcriptomic data. *Nature Methods* **17**, 147–154 (2020). [Nature Methods](https://www.nature.com/articles/s41592-019-0690-6) | [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7098173/)
    - **BEELINE benchmark.** **Verified:** uses standard ranking metrics (AUPRC, AUROC, Early Precision Ratio) over edges in curated ground-truth networks, treating all non-ground-truth edges as negatives without expression matching. Cite for: (a) the canonical GRN inference benchmark and the "GRN inference is hard, simple methods often beat complex ones" framing; (b) **as a contrast** — we depart from BEELINE's negative-sampling convention. Our expression-matched negative sampling in Layer 2 is a methodology choice we introduce to control for the obvious confound that highly-expressed genes have stronger signal in any embedding-based probe. (Earlier version of this note over-attributed the matched-negatives idea to BEELINE; corrected in audit.)

35. **Reimers N, Gurevych I.** Sentence-BERT: Sentence Embeddings using Siamese BERT-Networks. *EMNLP-IJCNLP* (2019). [arXiv](https://arxiv.org/abs/1908.10084) | [ACL](https://aclanthology.org/D19-1410/)
    - **Mean-pooling cell-embedding extraction.** Reimers & Gurevych empirically establish that mean-of-non-pad-tokens is the best fixed-size pooling for downstream linear probes (vs. CLS or max). Our cell embedding is exactly this aggregation over the gene token sequence. Cite in Methods §1.2 to justify the pooling choice in one line.

36. **Wilcoxon F.** Individual comparisons by ranking methods. *Biometrics Bulletin* **1**, 80–83 (1945).
    - **Wilcoxon signed-rank test.** Cite in Methods §1.3 / Results — paired test across CV folds and across TFs.

37. **Efron B.** Bootstrap methods: Another look at the jackknife. *Annals of Statistics* **7**, 1–26 (1979).
    - **Bootstrap confidence intervals.** Cite for our 95% bootstrap CI on per-fold/per-TF differences in Results.

---

---

## I. Donor stratification (added in fourth pass)

38. **Zimmerman KD, Espeland MA, Langefeld CD.** A practical solution to pseudoreplication bias in single-cell studies. *Nature Communications* **12**, 738 (2021). [Nature Communications](https://www.nature.com/articles/s41467-021-21038-1) | [PubMed](https://pubmed.ncbi.nlm.nih.gov/33531494/)
    - **Verified:** the paper documents that cells from the same donor are not statistically independent — intra-individual correlations exceed inter-individual correlations, producing "biased inference, highly inflated type 1 error rates, and reduced robustness." **Important caveat from audit:** Zimmerman's specific application target is *differential expression analysis* (their proposed fix is mixed models with donor as random effect). They do NOT explicitly write about cross-validation or classification probes. We cite Zimmerman as the canonical source for the underlying *statistical principle* (cells from one donor are not independent samples) and extend that principle to motivate our donor-stratified CV. **In the manuscript we should phrase this carefully**: "Following the principle articulated by Zimmerman et al. that cells from the same donor are not independent samples, we use donor-stratified 5-fold CV..." rather than implying Zimmerman recommended donor-stratified CV directly.

39. **Abdelaal T, Michielsen L, Cats D, et al.** A comparison of automatic cell identification methods for single-cell RNA sequencing data. *Genome Biology* **20**, 194 (2019). [Springer](https://link.springer.com/article/10.1186/s13059-019-1795-z) | [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC6734286/)
    - The most-cited cell-type classification benchmark in scRNA-seq. Uses 5-fold CV stratified by cell-type proportion ONLY, **not by donor**. Cite in our Discussion as the prevailing convention we are deliberately departing from on the basis of ref. 38.

40. **Luecken MD, Büttner M, Chaichoompu K, Danese A, Interlandi M, Mueller MF, Strobl DC, Zappia L, Dugas M, Colomé-Tatché M, Theis FJ.** Benchmarking atlas-level data integration in single-cell genomics. *Nature Methods* **19**(1), 41–50 (2022). [Nature Methods](https://www.nature.com/articles/s41592-021-01336-8) | [PubMed](https://pubmed.ncbi.nlm.nih.gov/34949812/)
    - **Verified:** 11 authors led by Luecken. scIB benchmark — 68 method × preprocessing combinations on 85 batches, >1.2M cells, 13 atlas-level integration tasks. Cite for: the canonical single-cell benchmarking framework that explicitly handles batch effects and contains the standard battery of integration metrics. We do not adopt scIB metrics directly (we are probing, not integrating), but it's the natural neighbor in the related-work section.

---

## Remaining open gaps

- Original ChEA / hTFtarget references *if* we end up adding them as alternative regulons in a v2 of the paper. Skip for now — current bibliography is sufficient for v1 manuscript.

Add as we draft.
