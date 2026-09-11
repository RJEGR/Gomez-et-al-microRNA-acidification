# data/

This folder holds the **input data objects** required to run the workflow. Large
files are **not** stored in git; they are distributed with the Zenodo archive
(DOI in the top-level README) and downloaded here.

Expected contents (populated during the input-reorganisation step):

| File | Description | Feeds |
|---|---|---|
| `Results.txt` / `Counts.txt` | ShortStack count matrix | Fig 1/2, DESeq2 |
| `mir.fasta` | miRNA mature/star sequences | Fig 1 |
| `Results.gff3`, `knownRNAs.gff3` | miRNA genomic coordinates | Fig 1D, loci |
| `IDENTICAL_SEQUENCES_MERGED_COUNT.rds` | processed miRNA matrix (117 × 12) | Fig 2B PCA, DESeq2 |
| `RNA_LOCATION_DB.tsv`, `RNA_LOCATION_MIR_DB.rds` | locus location database | Fig 1D, Fig 2C |
| RNAhybrid / TargetScan raw outputs | pre-intersection target tables | targets |
| `SRNA_REGULATORY_FUNCTION_DB.tsv` | post-intersection targets | targets, Table 3 |
| `gene_count_matrix.csv` | StringTie mRNA counts (PRJNA488641) | target expression |
| `SEQUENCES_MERGED_DESEQ_RES.tsv` | DESeq2 results, 4 contrasts | Table 3, Fig 3 |
| `LONGER_RELATIONAL_DB.tsv`, `Gene_ontologies_DB.tsv` | annotation / GO map | Fig 2D, Fig 3, Table 3 |
| `prof_by_read_length_summary.rds` | read-length profiling | Fig 1A, Fig 2, Table 2 |
| `Abulon_*.xlsx`, `MTD_RESPIROMETRY.csv` | raw respirometry | Table 1 |

Do not commit these files to the repository; keep them referenced through `config.R`.

## Already included in this bundle

- `prof_by_read_length_summary.tsv` — read-length × first-nucleotide profiling (converted from the original `.rds` for interoperability; 4,621 rows × 5 cols: `rnatype, Length, first_nuc, n, sample_id`). Feeds Fig 1A / Fig 2 / Table 2.
- `water_chemistry/` — `quimicas*.tsv/csv`, `CO2_sys_inputs.tsv` (Tables S1–S2; see `06_water_chemistry/`).
- `target_prediction/mature_star_mir_vs_mir_vs_utr_rmdup_RNAhybrid.out.psig_targetscan.out` — TargetScan overlap on the RNAhybrid significant hits (post step). The raw RNAhybrid output (`mir_vs_utr_rmdup_RNAhybrid.out.tsv`, ≈2.6 GB) is **not** bundled — archive it gzip-compressed on Zenodo.

Still to add here from the connected folder (large; Zenodo-fed): `Results.txt`/`Counts.txt`, `Results.gff3`, `RNA_LOCATION_DB.tsv`, `gene_count_matrix.csv`, `SEQUENCES_MERGED_DESEQ_RES.tsv`, `LONGER_RELATIONAL_DB.tsv`, `Gene_ontologies_DB.tsv`, `IDENTICAL_SEQUENCES_MERGED_COUNT.rds`, `SRNA_REGULATORY_FUNCTION_DB.tsv`, `mir.fasta`, `MIRGENEDB_2.1.tsv`, metadata.
