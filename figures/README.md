# figures/ — manuscript figures

| Script | Produces | Key inputs |
|---|---|---|
| `FIGURE_1.R` | **Fig 1** panels A–C (miRNA identification) | `prof_by_read_length_summary.rds`, `MIRGENEDB_2.1.tsv`, `RNA_LOCATION_MIR_DB.rds`, `METADATA.tsv` |
| `LOCI_PLOT.R` | **Fig 1D** chromosomal map + inter-loci distances | `RNA_LOCATION_DB.tsv`, ENSEMBL genome (`ANNOT_DIR`) |
| `ALIGNMENT_DETAILS_EXPLORATORY.R` | Fig 1B precision / alignment details | ShortStack `alignment_details.tsv` |
| `FIGURE_2_READ_PROFILING.R` | **Fig 2** read-length profiling; clade panel | `prof_by_read_length_summary.rds`; miRTrace output (`ANNOT_DIR`, FAIR-external) |
| `PCA.R`, `PCA_variance_stabilized.R` | **Fig 2B** PCA of VST counts | `IDENTICAL_SEQUENCES_MERGED_COUNT.rds` |
| `INTRAGENIC_ANALYSIS.R` | **Fig 2C** intergenic/intragenic proportion | `RNA_LOCATION_DB.tsv`, `LONGER_RELATIONAL_DB.tsv` |
| `Figure_3.R` | **Fig 3** degree & density | `LONGER_RELATIONAL_DB.tsv`, `SEQUENCES_MERGED_DESEQ_RES.tsv`, `Biological_themes.tsv` |
| `OA_MIR_degree_analysis.R` | **Fig 3** degree analysis (+ STRING network) | `LONGER_RELATIONAL_DB.tsv`, `protein_links_full_v12.rds` (`ANNOT_DIR`, FAIR-external) |
| `FIG_S1_redesign.R` | **Fig S1** (composition / MA / concordance) | `ALL_117_miRNAs_4contrasts.tsv`, `IDENTICAL_SEQUENCES_MERGED_COUNT.rds`, `METADATA_MICRORNAS.tsv` |

## Paths
All scripts resolve paths through `config.R` (`DATA_DIR` / `ANNOT_DIR` / `RESULTS_DIR`). **No absolute
`~/Documents/…` or `/Users/cigom/…` paths remain.** `FIG_S1_redesign.R` additionally honours the env
vars `ALL117` / `COUNTS` / `META` / `OUTDIR`, which now default to the `config.R` paths.

## Inputs
Bundled in `data/`: `prof_by_read_length_summary.{tsv,rds}`, `MIRGENEDB_2.1.tsv`,
`RNA_LOCATION_MIR_DB.rds`, `IDENTICAL_SEQUENCES_MERGED_COUNT.rds`, `ALL_117_miRNAs_4contrasts.tsv`,
`METADATA.tsv`, `METADATA_MICRORNAS.tsv`, `LONGER_RELATIONAL_DB.tsv`, `Gene_ontologies_DB.tsv`.
Large — copy from your folder or Zenodo: `RNA_LOCATION_DB.tsv` (~15 MB).
FAIR-external (`data/annotation/`, see README §6–§7): ENSEMBL genome (`LOCI_PLOT.R`), miRTrace clade
output (`FIGURE_2_READ_PROFILING.R` clade panel), STRING `protein_links_full_v12.rds`
(`OA_MIR_degree_analysis.R`).

> Figures are written to `RESULTS_DIR` (i.e. `outputs/`); every `ggsave` in this folder now points there.
