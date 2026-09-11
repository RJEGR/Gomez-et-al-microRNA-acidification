#!/usr/bin/env python3
"""
Tabla 3 — consolidacion de contrastes DESeq2 + origen genomico + blancos destacados.

Entradas
  SEQUENCES_MERGED_DESEQ_RES.tsv          (llave: MajorRNA + CONTRAST)
  LONGER_RELATIONAL_DB.tsv                (llave: MajorRNA ; aporta biotype_best_rank y gene_id)
  Supplementary_tables_Biological_themes  (llave: gene_id ; aporta STRINGID y tema biologico)

Convencion de signo (heredada de DESeq2, sampleA -> sampleB):
  CONTRAST_A  Control -> Low        @ 24 hpf    negativo = mayor a pH 7.6
  CONTRAST_B  Control -> Low        @ 110 hpf   negativo = mayor a pH 7.6
  CONTRAST_C  Early -> Competent    @ pH 8.0    negativo = mayor a 110 hpf
  CONTRAST_D  Early -> Competent    @ pH 7.6    negativo = mayor a 110 hpf
"""
import pandas as pd, numpy as np, re

DESEQ = "29723010-1787333112900_SEQUENCES_MERGED_DESEQ_RES.tsv"
RELDB = "5032f9a2-1787334047108_LONGER_RELATIONAL_DB.tsv"
THEME = "b15d18a2-1787334001651_Supplementary_tables__Biological_themes.tsv"

# ---------------------------------------------------------------- 1. DESeq
de = pd.read_csv(DESEQ, sep="\t", dtype=str)
for c in ["log2FoldChange", "padj", "baseMean", "baseMeanA", "baseMeanB"]:
    de[c] = pd.to_numeric(de[c], errors="coerce")

# Etiqueta canonica: el ID de MirGeneDB manda; 'Name' es poco fiable
# (p. ej. LET-7-5p viene etiquetado como Cluster_88 y BANTAM-3p como Cluster_86).
# En los loci novel, MirGeneDB_ID repite la secuencia MajorRNA en lugar de un ID;
# ese es el discriminante known/novel, no la columna 'Name'.
de["known"] = (de.MirGeneDB_ID.notna() & (de.MirGeneDB_ID != de.MajorRNA)
               & ~de.MirGeneDB_ID.isin(["NA", ""]))
de["mirna"] = np.where(de.known, de.MirGeneDB_ID, de.Name)

wide = de.pivot_table(index=["mirna", "MajorRNA", "known"], columns="CONTRAST",
                      values=["log2FoldChange", "padj"], aggfunc="first").reset_index()
wide.columns = [c[0] if not c[1] else f"{'lfc' if c[0]=='log2FoldChange' else 'padj'}_{c[1][-1]}"
                for c in wide.columns]
bm = de.groupby(["mirna", "MajorRNA"], as_index=False).baseMean.mean()
wide = wide.merge(bm, on=["mirna", "MajorRNA"])
# Fragilidad: conteo normalizado del grupo MAS PEQUENO en cada contraste.
# Un log2FC grande sostenido por 1-2 lecturas no es un efecto, es un cero.
frag = (de.assign(minbm=de[["baseMeanA", "baseMeanB"]].min(axis=1))
          .pivot_table(index="MajorRNA", columns="CONTRAST", values="minbm", aggfunc="first"))
frag.columns = [f"minCount_{c[-1]}" for c in frag.columns]
wide = wide.merge(frag, left_on="MajorRNA", right_index=True, how="left")

# ---------------------------------------------------------------- 2. Origen genomico
rel = pd.read_csv(RELDB, sep="\t", dtype=str)
# Un mismo MajorRNA puede tener >1 fila de origen: (a) mismo biotipo con dos loci
# hospederos (miR-190-5p), o (b) biotipo ambiguo entre loci (miR-278-3p, no reportado).
# Se colapsa a un valor por locus y se deja constancia de los casos ambiguos.
origin = (rel.groupby("MajorRNA")
             .agg(biotype_best_rank=("biotype_best_rank", lambda s: "/".join(sorted(set(s.dropna())))),
                  n_host_loci=("Loci_geneid", lambda s: s.dropna().nunique()))
             .reset_index())
amb = origin[origin.biotype_best_rank.str.contains("/")]
if len(amb):
    print("[aviso] biotipo ambiguo en:", ", ".join(amb.MajorRNA), flush=True)
wide = wide.merge(origin, on="MajorRNA", how="left")

# ---------------------------------------------------------------- 3. Blancos + temas
th = pd.read_csv(THEME, sep="\t", dtype=str).rename(columns={"Bilogical_pathway": "theme"})
th["theme"] = th.theme.fillna("").str.strip()
th["STRINGID"] = th.STRINGID.replace({"NA": np.nan})

tg = (rel[["MajorRNA", "gene_id"]].drop_duplicates()
        .merge(th[["gene_id", "STRINGID", "theme", "ENSEMBL description"]], on="gene_id", how="left"))
# etiqueta legible del blanco: proteina STRING si existe, si no el locus
tg["label"] = tg.STRINGID.fillna(tg.gene_id)
# una fila por (locus, blanco, tema)  -> explota los temas multiples
tg["theme"] = tg.theme.fillna("")
tg = tg.assign(theme=tg.theme.str.split(",")).explode("theme")
tg["theme"] = tg.theme.str.strip()

n_targets = tg.groupby("MajorRNA").gene_id.nunique().rename("n_targets")
n_named   = (tg[tg.STRINGID.notna()].groupby("MajorRNA").gene_id.nunique()
               .rename("n_named"))

# Blancos que el manuscrito ya nombra: tienen prioridad en la celda
MS_TARGETS = set("""RNF166 CENPM CDK10 TFB1M C2CD2L MAP11 ABCG2 HM13 PARP14 ZNFX1 ANKIB1
SLC46A1 SDHC CTBP1 TTLL7 IQCM DLD PCSK6 GCNT1 RBCK1 H2AZ2 SBDS PPP2R5D ANKRD55 RAB32 KIF12
TIGD4 SLC10A5 SPOPL RPS12""".split())
MS_LOCI = {"LOC125383552": "Ank (LOC125383552)",
           "LOC124144562": "glutamic acid-rich protein (LOC124144562)",
           "LOC124113664": "zinc finger CCHC (LOC124113664)"}

THEME_ORDER = ["Development", "Growth", "Biomineralization", "Respiratory metabolism",
               "Immune system", "Genome maintenance", "Epigenetics"]

def targets_cell(major, max_themes=4, max_genes=6):
    """Celda 'Blancos destacados': blancos con proteina STRING, agrupados por tema.
    Cada gen aparece una sola vez (se asigna al tema de mayor prioridad) para que la
    celda no repita SBDS/HM13 en dos temas."""
    sub = tg[(tg.MajorRNA == major) & (tg.theme != "") & tg.theme.notna()
             & tg.STRINGID.notna()].copy()
    if sub.empty:
        return "—"
    sub["label"] = sub.label.str.replace(";", "/", regex=False)
    counts = sub.groupby("theme").gene_id.nunique().to_dict()
    themes = sorted(counts, key=lambda t: (-counts[t],
                    THEME_ORDER.index(t) if t in THEME_ORDER else 99))
    used, parts, shown_n = set(), [], 0
    for t in themes:
        labs = [l for l in sub[sub.theme == t].label.unique() if l not in used]
        labs = sorted(labs, key=lambda x: (x not in MS_TARGETS, x))
        if not labs or len(parts) >= max_themes:
            continue
        take = labs[:max(1, max_genes - shown_n)]
        used.update(take); shown_n += len(take)
        parts.append(f"{t}: {', '.join(take)}")
        if shown_n >= max_genes:
            break
    total_named = sub.gene_id.nunique()
    rest = total_named - len(used)
    tail = f" (+{rest} more)" if rest > 0 else ""
    return "; ".join(parts) + tail

wide["targets"] = wide.MajorRNA.map(targets_cell)
wide = wide.merge(n_targets, left_on="MajorRNA", right_index=True, how="left") \
           .merge(n_named,   left_on="MajorRNA", right_index=True, how="left")

# ---------------------------------------------------------------- 4. Filtro: miRNAs del manuscrito
REPORTED = ["MIR-133-3p", "LET-7-5p", "MIR-10-5p", "MIR-2-3p", "MIR-153-3p",
            "MIR-2722", "MIR-1988-5p", "BANTAM-3p", "MIR-315-5p", "MIR-190-5p",
            "Cluster_71", "Cluster_74"]
# Locus cuyos numeros cita el manuscrito, cuando el miRNA tiene mas de uno
MS_LOCUS = {"MIR-10-5p": "UACCCUGUAGAUCCGAAUUUGU",
            "MIR-2-3p":  "UAUCACAGCCAGCUUUGAUGAGCU"}

t3 = wide[wide.mirna.isin(REPORTED)].copy()
t3["in_manuscript_text"] = [
    (m not in MS_LOCUS) or (seq == MS_LOCUS[m])
    for m, seq in zip(t3.mirna, t3.MajorRNA)]
t3["_ord"] = t3.mirna.map({m: i for i, m in enumerate(REPORTED)})
t3 = t3.sort_values(["_ord", "in_manuscript_text", "baseMean"],
                    ascending=[True, False, False]).drop(columns="_ord")

# etiqueta de locus para los multi-locus (a, b, c por expresion media decreciente)
t3["locus_tag"] = ""
for m, grp in t3.groupby("mirna"):
    if len(grp) > 1:
        for rank, i in enumerate(grp.sort_values("baseMean", ascending=False).index):
            t3.loc[i, "locus_tag"] = "abcdef"[rank]

# magnitud en veces, para redactar el texto
for c in "ABCD":
    t3[f"fold_{c}"] = 2 ** t3[f"lfc_{c}"].abs()

cols = ["mirna", "locus_tag", "MajorRNA", "known", "biotype_best_rank", "n_host_loci",
        "baseMean", "lfc_C", "padj_C", "lfc_D", "padj_D", "lfc_A", "padj_A",
        "lfc_B", "padj_B", "fold_C", "fold_D", "fold_A", "fold_B",
        "minCount_A", "minCount_B", "minCount_C", "minCount_D",
        "n_targets", "n_named", "targets", "in_manuscript_text"]
t3[cols].to_csv("TABLE3_backing_data.tsv", sep="\t", index=False, float_format="%.4g")
wide.to_csv("ALL_117_miRNAs_4contrasts.tsv", sep="\t", index=False, float_format="%.4g")

# ---------------------------------------------------------------- 5. Auditoria de blancos citados
CLAIMS = {
 "MIR-133-3p":  [],
 "LET-7-5p":    ["RNF166", "CENPM", "CDK10"],
 "MIR-10-5p":   ["TFB1M", "C2CD2L", "RNF166", "MAP11", "ABCG2",
                 "LOC125383552", "LOC124144562", "LOC124113664"],
 "MIR-2-3p":    ["HM13", "PARP14", "ZNFX1", "RNF166"],
 "MIR-153-3p":  ["ANKIB1", "SLC46A1"],
 "MIR-2722":    ["RNF166", "LOC125383552"],
 "MIR-1988-5p": ["PPP2R5D", "ANKRD55", "RNF166", "SDHC", "CTBP1"],
 "BANTAM-3p":   ["TTLL7", "IQCM", "DLD"],
 "MIR-315-5p":  ["PCSK6", "GCNT1", "RBCK1", "CTBP1", "H2AZ2"],
 "MIR-190-5p":  ["CTBP1", "RNF166", "SBDS"],
}
seq2mir = dict(zip(t3.MajorRNA, t3.mirna))
tg2 = tg.assign(mirna=tg.MajorRNA.map(seq2mir)).dropna(subset=["mirna"])
audit = []
for mir, claimed in CLAIMS.items():
    have = set(tg2[tg2.mirna == mir].STRINGID.dropna()) | set(tg2[tg2.mirna == mir].gene_id)
    for c in claimed:
        audit.append({"miRNA": mir, "blanco_citado": c,
                      "en_base": "si" if c in have else "NO ENCONTRADO"})
audit = pd.DataFrame(audit)
# Promiscuidad del blanco: a cuantos de los 117 loci esta ligado en la base.
promis = rel.groupby("gene_id").MajorRNA.nunique()
# OJO: zip(th.STRINGID.dropna(), th.gene_id) desalinea (una serie filtrada contra
# una completa se aparea por posicion). Se construye desde las filas ya filtradas.
_m = th.dropna(subset=["STRINGID"])
name2gene = dict(zip(_m.STRINGID, _m.gene_id))
audit["n_miRNAs_que_lo_targetean"] = [
    promis.get(name2gene.get(b, b), np.nan) for b in audit.blanco_citado]
audit["de_117_loci"] = 117
audit.to_csv("AUDIT_targets_claimed_vs_db.tsv", sep="\t", index=False)
print("\n=== promiscuidad de los blancos citados en el texto ===")
print(audit.sort_values("n_miRNAs_que_lo_targetean", ascending=False)
          .drop_duplicates("blanco_citado")[["blanco_citado","n_miRNAs_que_lo_targetean"]]
          .head(12).to_string(index=False))

print(t3[["mirna","locus_tag","biotype_best_rank","baseMean","lfc_C","lfc_D","lfc_A","lfc_B",
          "n_targets","n_named","in_manuscript_text"]].to_string(index=False))
print()
print("=== blancos citados en el texto que NO estan en la base ===")
print(audit[audit.en_base != "si"].to_string(index=False))


# ---------------------------------------------------------------- 6. Render Tabla 3
SUP = str.maketrans("0123456789-", "\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079\u207b")

def fmt_p(p):
    if pd.isna(p):
        return "n.d."
    if p < 1e-6:
        return "< 1 \u00d7 10\u207b\u2076"
    m, e = f"{p:.1e}".split("e")
    return f"{m} \u00d7 10{str(int(e)).translate(SUP)}"

def cell(lfc, p):
    if pd.isna(lfc):
        return "n.d."
    v = f"{lfc:+.2f}".replace("+", "+").replace("-", "\u2212")
    return f"**{v}** ({fmt_p(p)})" if (pd.notna(p) and p < 0.05) else f"{v} (n.s., {fmt_p(p)})"

pretty = {"MIR-133-3p":"miR-133-3p","LET-7-5p":"let-7-5p","MIR-10-5p":"miR-10-5p",
          "MIR-2-3p":"miR-2-3p","MIR-153-3p":"miR-153-3p","MIR-2722":"miR-2722",
          "MIR-1988-5p":"miR-1988-5p","BANTAM-3p":"bantam-3p","MIR-315-5p":"miR-315-5p",
          "MIR-190-5p":"miR-190-5p","Cluster_71":"Cluster_71 (novel)",
          "Cluster_74":"Cluster_74 (novel)"}

rows = []
for _, r in t3.iterrows():
    name = pretty.get(r.mirna, r.mirna)
    if r.locus_tag:
        name += f" (locus {r.locus_tag})"
    rows.append([name, r.biotype_best_rank,
                 cell(r.lfc_C, r.padj_C), cell(r.lfc_D, r.padj_D),
                 cell(r.lfc_A, r.padj_A), cell(r.lfc_B, r.padj_B),
                 f"{r.targets} [n = {int(r.n_targets)}]"])

hdr = ["miRNA", "Origin",
       "log\u2082FC development @ pH 8.0 (p adj)", "log\u2082FC development @ pH 7.6 (p adj)",
       "log\u2082FC pH @ 24 hpf (p adj)", "log\u2082FC pH @ 110 hpf (p adj)",
       "Highlighted targets [n predicted]"]
md = ["| " + " | ".join(hdr) + " |", "|" + "---|" * len(hdr)]
md += ["| " + " | ".join(x) + " |" for x in rows]
open("TABLE3.md", "w").write("\n".join(md) + "\n")
print("\n".join(md[:5]))
