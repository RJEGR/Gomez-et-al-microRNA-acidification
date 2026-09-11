# 04_annotation_go — relational annotation DB & GO enrichment

Feeds **Figure 2D** (GO enrichment) and builds the annotation table (`LONGER_RELATIONAL_DB.tsv`)
used by **Table 3** and **Figure 3**.

## Scripts
- `JOIN_RELATIONAL_DB.R` — **provenance/build**: joins eggNOG-mapper annotation, STRING best-hit homology (human 9606, *C. elegans* 6239) and the miRNA→target DB into `LONGER_RELATIONAL_DB.tsv` and `Gene_ontologies_DB.tsv`.
- `3_STRING-db.R` — **provenance/build**: parses the DIAMOND-blastx-vs-STRING hits and assembles the STRING interaction graph (`protein_links_full_v12.rds`, `gene2stringid_diamond_blastx.tsv`).
- `TopGO_by_CONTRAST.R` — original topGO enrichment run (needs `org.Hs.eg.db` GOSemSim data).
- `TopGO_by_CONTRAST_viz_RECONSTRUCTED.R` — **current Fig 2D visualisation** (WGCNA-free). Reads the reconstructed enrichment table and writes the figure. **Run this to regenerate the figure.**
- `TopGO_by_CONTRAST_viz.R` — earlier visualisation, superseded by the reconstructed one.

## Run (current figure)
```
Rscript 04_annotation_go/TopGO_by_CONTRAST_viz_RECONSTRUCTED.R
```

## Paths
All scripts now resolve paths through `config.R`: `DATA_DIR` (deposited inputs),
`ANNOT_DIR = data/annotation` (annotation sources), `RESULTS_DIR` (outputs). **No absolute
`~/Documents/…` or `/Users/cigom/…` paths remain.**

## Inputs
Deposited in `data/` (bundled): `LONGER_RELATIONAL_DB.tsv`, `Gene_ontologies_DB.tsv`,
`Boostrap_topGO_intra_inter_reconstructed.tsv`. Also needed (large — copy from your folder or Zenodo):
`RNA_LOCATION_DB.tsv`, `RNA_LOCATION_MIR_DB.rds`.

`data/annotation/` — **FAIR-external / provenance inputs, NOT redeposited** (see README §6–§7):
eggNOG-mapper outputs, DIAMOND-blastx-vs-STRING `.outfmt6`, STRING `protein.links.full/info.v12.0`,
`STRING_DB_species.v12.0.txt`, ENSEMBL `genome_features.rds` / `cogs.rds`,
`SRNA_FUNCTION_PREDICTED_LONG_EXPRESSED.rds`, and `org.Hs.eg.db` GOSemSim data.
`JOIN_RELATIONAL_DB.R`, `3_STRING-db.R` and `TopGO_by_CONTRAST.R` are **build/provenance scripts** and
require these external inputs to re-run; the manuscript figures are reproducible **downstream** from the
deposited tables above.

## Note (WGCNA)
Per reviewer comment 2, WGCNA was dropped from the manuscript. The current GO figure comes from
`TopGO_by_CONTRAST_viz_RECONSTRUCTED.R`; the older scripts still reference a
`SEQUENCES_MERGED_DESEQ_RES_WGCNA.rds` object and are kept only for provenance.
