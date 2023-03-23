# Single-Cell Analysis of BPDCN Evolution

This repository contains single-cell analysis on bone marrow samples from healthy donors and blastic plasmacytoid dendritic cell neoplasm (BPDCN) patients. We performed single-cell RNA-sequencing (scRNA-seq) using Seq-Well S^3 ([Hughes 2020](http://dx.doi.org/10.1016/j.immuni.2020.09.015)) for healthy donor 6 and Patient 9 and 10x Genomics 3' gene expression analysis for the remaining 15 samples. In addition, we performed single-cell genotyping by enriching mutation sites using the eXpressed Variant sequencing (XV-seq) protocol (based on [Van Galen 2019](http://dx.doi.org/10.1016/j.cell.2019.01.031) and [Nam 2019](http://dx.doi.org/10.1038/s41586-019-1367-0)). Single-cell transcriptomes were aligned to hg38 using custom scripts (Seq-Well data) or Cell Ranger (10x data). Genotyping results were analyzed using [IronThrone-GoT](https://github.com/dan-landau/IronThrone-GoT) and filtered according to parameters described in the Methods of the accompanying paper: Griffin et al, 2023. The analysis is broken up into 13 semi-chronological parts as follows.

* Load gene expression data from healthy donors, integrate, reduce dimensionality, cluster *(01_Seurat_Harmony)*

* Annotate healthy donor cell clusters using canonical marker genes *(02_Annotate)*

* Annotate patient cells using healthy donors as a reference, including calculation of projected UMAP coordinates *(03_RandomForest)*

* Add single-cell genotyping metadata to the Seurat objects containing gene expression and cell type annotations *(04_XV-seq)*

* Evaluate XV-seq efficiency, compare to bulk sequencing VAFs, relate to gene expression levels *(05_Stats)*

* Project patient cells onto healthy donor UMAP colored by cell type annotation or mutation status *(06_UMAP_projections)*

* Generate heatmaps of cell type marker genes with annotation bars indicating mutation status *(07_DEG_Heatmaps)*

* Generate heatmaps of mutated cell proportions within different lineages *(08_Mutation_ratios)*

* Evaluate co-occurrence of single-nucleotide variants (SNVs) and haplotype B to infer somatic evolution trajectory *(09_Pt9_haplotype)*

* Classify host/donor origin of Pt10 relapse cells (post-transplant) and intersect with mutation status to confirm malignant cell identity and somatic evolution  *(10_Pt10_donor-host)*

* Evaluate stepwise gene expression changes from healthy pDCs to pre-malignant pDCs to malignant BPDCN cells  *(11_pDC_expr)*

* Visualize relationships between single-cell Random Forest pDC prediction scores, malignant BPDCN signature scores, and mutation status *(12_Mut_vs_scores)*

* Miscellaneous scripts for future follow-up studies *(13_Misc)*
