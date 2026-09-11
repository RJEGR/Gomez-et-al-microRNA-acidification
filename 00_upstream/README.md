# 00_upstream — read QC, ShortStack annotation, and the reference transcriptome

This stage runs on an HPC environment and produces the ShortStack outputs and the reference
transcriptome that feed the downstream modules. It is **documented** here rather than re-run locally.

## A. Small-RNA QC and locus annotation
1. **Read QC / preprocessing** — adapter and quality trimming; retain 18–45 nt, Phred > 30; contaminant
   screening with **miRTrace v1.0.1** (rRNA/tRNA and clade-specific databases).
2. **Locus annotation & quantification** — **ShortStack v4.0** on the concatenated nuclear
   (GCA_023055435.1) + mitochondrial (JALGQA010000616.1) reference, with `--dn_mirna`, `--pad 1`,
   `--mincov 0.8`, `--mmap u`, `--dicermax 30`, `--strand_cutoff 0.8`, `--known_miRNAs` (MirGeneDB 2.1).

**Outputs** (archived on Zenodo, placed in `data/`): `Results.txt` (count matrix), `Results.gff3` /
`knownRNAs.gff3` (coordinates), `mir.fasta`, `MajorRNA.fasta`, `alignment_details.tsv`.

The original command-line documentation lives in the working repository under `A_UPSTREAM/`
(`PREPROCESSING.md`, `DATABASES.md`, `SHORTSTACKS4.md`, `FUNCTIONAL_PREDICTION.md`).

## B. Reference transcriptome (target-expression check)
`REFBASED_TRANSCRIPTOME_ASSEMBLY.md` — full HISAT2 + StringTie reference-guided assembly protocol
(Pertea et al. 2016) applied to the public *H. rufescens* RNA-seq libraries
(BioProject **PRJNA488641**; Masonbrink et al. 2019). It produces the StringTie count matrices
`gene_count_matrix.csv` / `transcript_count_matrix.csv` used by `02_targets/` and `04_annotation_go/`
to confirm that predicted targets are expressed during larval development.

> These are provenance protocols (run on an HPC). Their products are provided as deposited inputs in
> `data/`, so the downstream R modules reproduce without re-running the upstream stage.
