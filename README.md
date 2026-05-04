# GO Enrichment Galaxy Tool

Per-module Gene Ontology enrichment analysis for the Galaxy bioinformatics
platform. Companion tool to the
[WGCNA Galaxy Tool](https://github.com/Luiza-Romao/wgalaxy_complete).

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

---

## Table of contents

- [Overview](#overview)
- [File structure](#file-structure)
- [Requirements](#requirements)
- [Installation](#installation)
- [Input formats](#input-formats)
- [Parameters](#parameters)
- [Outputs](#outputs)
- [Transcript ID mapping guide](#transcript-id-mapping-guide)
- [Workflow with the WGCNA tool](#workflow-with-the-wgcna-tool)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [Authors and funding](#authors-and-funding)
- [License](#license)

---

## Overview

This tool performs per-module Gene Ontology (GO) enrichment analysis using
clusterProfiler. It takes as primary input the module gene assignment table
produced by the WGCNA Galaxy Tool and runs GO enrichment independently for
each selected module across the BP, MF, and CC ontologies.

A central feature is built-in **transcript ID mapping** for non-model
organisms. Species such as sugarcane (*Saccharum officinarum* hybrid) produce
transcript identifiers like `transcript_002741` that are not present in any GO
annotation database. This tool accepts a user-provided two-column mapping file
that converts those identifiers to ortholog IDs from a model organism (for
example, Arabidopsis TAIR IDs), enabling GO enrichment for genomes that lack
dedicated annotation databases.

---

## File structure

| File | Description |
|---|---|
| `go_enrichment.xml` | Galaxy tool wrapper — inputs, outputs, UI, help text |
| `go_enrichment_macros.xml` | Reusable XML macros (requirements, shared parameters) |
| `go_enrichment.R` | Main R analysis script |
| `conda_environment.yml` | Reproducible conda environment with pinned dependencies |
| `BUGS.md` | Development bug log with root causes and fixes |
| `CHANGELOG.md` | Version history |

---

## Requirements

All Bioconductor packages are pinned to release 3.20 (R 4.4). Mixing packages
from different Bioconductor releases causes dependency resolution failures at
environment build time.

| Package | Version | Source |
|---|---|---|
| r-base | 4.4 | conda-forge |
| bioconductor-clusterprofiler | 4.14.0 | Bioconductor 3.20 |
| bioconductor-enrichplot | 1.26.0 | Bioconductor 3.20 |
| bioconductor-dose | 3.32.0 | Bioconductor 3.20 |
| bioconductor-annotationdbi | 1.68.0 | Bioconductor 3.20 |
| bioconductor-org.at.tair.db | 3.20.0 | Bioconductor 3.20 |
| bioconductor-org.hs.eg.db | 3.20.0 | Bioconductor 3.20 |
| bioconductor-org.mm.eg.db | 3.20.0 | Bioconductor 3.20 |
| bioconductor-org.os.eg.db | 3.20.0 | Bioconductor 3.20 |
| bioconductor-org.sc.sgd.db | 3.20.0 | Bioconductor 3.20 |
| r-ggplot2 | 3.5.1 | conda-forge |
| r-ggrepel | 0.9.6 | conda-forge |
| r-optparse | 1.7.5 | conda-forge |

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/Luiza-Romao/go-enrichment-galaxy.git
cd go-enrichment-galaxy
```

### 2. Create the conda environment

```bash
conda env create -f conda_environment.yml
conda activate go_enrichment_env
```

### 3. Test locally with Planemo

```bash
conda activate planemo_env
planemo serve --port 8081
```

Open `http://localhost:8081` in a browser and submit a job using the test data.

### 4. Validate syntax and run automated tests

```bash
conda activate go_enrichment_env
Rscript -e "parse('go_enrichment.R'); cat('Syntax OK\n')"
planemo test go_enrichment.xml
```

### 5. Deploy to a Galaxy instance

```bash
cp -r go-enrichment-galaxy/ /path/to/galaxy/tools/go_enrichment/
```

Register in `tool_conf.xml`:

```xml
<section id="transcriptomics" name="Transcriptomics">
    <tool file="go_enrichment/go_enrichment.xml" />
</section>
```

Restart Galaxy after registering the tool.

---

## Input formats

### Module gene assignments table (required)

The `Module gene assignments` output from the WGCNA Galaxy Tool, or any
tab-separated file with a gene ID column and a module color column.

```
Gene              Module
transcript_002741 turquoise
transcript_018843 turquoise
transcript_031021 blue
```

Accepted column names for the gene column: `Gene`, `GeneID`, `gene_id`,
`transcript_id`. Accepted column names for the module column: `Module`,
`Color`, `color`. The `grey` module (unassigned genes) is automatically
excluded from enrichment.

### ID mapping file (optional)

A tab-separated file with exactly two named columns: one containing your
organism's transcript IDs and one containing the model organism's gene IDs.
A header row is required.

```
transcript_id     target_id
transcript_002741 AT1G01010
transcript_018843 AT3G54440
transcript_031021 AT5G20790
```

Column names are configurable in the tool parameters. One-to-many mappings
are accepted: if a single transcript maps to multiple target IDs (for example
from OrthoFinder orthogroups), all target IDs are included in the enrichment.
Unmapped transcripts are excluded and recorded in the mapping report output.

See [Transcript ID mapping guide](#transcript-id-mapping-guide) for
instructions on generating this file from BLAST output.

### Selected modules table (optional)

The `Selected modules summary` output from the WGCNA Galaxy Tool. When
provided, only the modules listed in this table are analyzed. This restricts
enrichment to modules that passed the correlation threshold in the WGCNA run
and is the recommended mode when running both tools in sequence.

---

## Parameters

### Module selection

| Parameter | Default | Description |
|---|---|---|
| Filter mode | all | `all` analyzes every module except grey; `sel_table` restricts to the WGCNA selected-modules table; `manual` accepts a comma-separated list of module names |
| Module names | turquoise,blue | Active only when mode is `manual`; names must match module colors exactly |
| Minimum genes per module | 10 | Modules with fewer mapped genes than this threshold are skipped |

### ID mapping

| Parameter | Default | Description |
|---|---|---|
| Use mapping | no | Whether to apply transcript ID conversion before enrichment |
| Mapping file | — | Two-column TSV with source and target IDs |
| Source column name | transcript_id | Column in the mapping file containing your transcript IDs |
| Target column name | target_id | Column in the mapping file containing model organism IDs |

### Organism database

| Parameter | Default | Description |
|---|---|---|
| OrgDb | org.At.tair.db | Annotation database for the target organism |
| Key type | TAIR | ID format expected by the selected OrgDb; must match the IDs in the mapped gene list |
| GO ontology | ALL | `BP`, `MF`, `CC`, or `ALL` (runs all three) |

### Enrichment thresholds

| Parameter | Default | Description |
|---|---|---|
| Adjusted p-value cutoff | 0.05 | BH-adjusted p-value threshold |
| q-value cutoff | 0.2 | Additional FDR control via q-value |
| Minimum gene-set size | 10 | GO terms annotating fewer genes are excluded |
| Maximum gene-set size | 500 | GO terms annotating more genes are excluded |
| Top terms per module | 15 | Number of most significant terms shown in plots |

---

## Outputs

| Output | Format | Description |
|---|---|---|
| GO enrichment results | TSV | Combined table across all modules and ontologies; columns: Module, Ontology, ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count |
| ID mapping report | TSV | Per-transcript mapping result: source ID, target ID, module; records unmapped transcripts for diagnostic purposes |
| Dotplot — all modules | PNG | Overview dotplot faceted by module and ontology; point size encodes gene count, color encodes adjusted p-value |
| Barplot — all modules | PNG | Overview barplot faceted by module and ontology |
| Per-module plots | PNG collection | Individual dotplot and barplot for each module with significant enrichment |
| Analysis log | TXT | Timestamped execution log with mapping rate, gene counts, and significant term counts per module |

### Combined table column reference

| Column | Description |
|---|---|
| Module | WGCNA module color |
| Ontology | GO sub-ontology: BP, MF, or CC |
| ID | GO term accession (e.g. GO:0009909) |
| Description | GO term name |
| GeneRatio | Fraction of module genes annotated with this term |
| BgRatio | Fraction of all background genes annotated with this term |
| pvalue | Hypergeometric test p-value |
| p.adjust | BH-corrected p-value |
| qvalue | q-value (Storey method) |
| geneID | Slash-separated list of mapped gene IDs driving the enrichment |
| Count | Number of genes driving the enrichment |

---

## Transcript ID mapping guide

Sugarcane does not have a dedicated GO annotation database. The standard
approach is to use *Arabidopsis thaliana* as a proxy by identifying sequence
homologs and transferring annotations. This requires three steps: downloading
the Arabidopsis proteome, running a BLAST search, and formatting the output
as a two-column mapping file.

### Step 1 — Download the Arabidopsis proteome

```bash
wget https://www.arabidopsis.org/download_files/Proteins/\
TAIR10_genome_release/TAIR10_pep_20110103_representative_gene_model \
-O tair10_proteins.fa
```

### Step 2 — Build a BLAST database and search

```bash
makeblastdb \
  -in    tair10_proteins.fa \
  -dbtype prot \
  -out   tair10_db

blastx \
  -query           sugarcane_transcripts.fa \
  -db              tair10_db \
  -out             blast_results.txt \
  -outfmt          6 \
  -evalue          1e-5 \
  -max_target_seqs 1 \
  -num_threads     8
```

### Step 3 — Extract the two-column mapping file

```bash
awk '{print $1 "\t" $2}' blast_results.txt \
  | sort -k1,1 -k11,11g \
  | awk '!seen[$1]++' \
  > mapping_noheader.tsv

echo -e "transcript_id\ttarget_id" \
  | cat - mapping_noheader.tsv \
  > id_mapping.tsv
```

The `target_id` column will contain Arabidopsis identifiers such as
`AT1G01010.1`. The tool accepts the isoform suffix; `org.At.tair.db` resolves
both `AT1G01010` and `AT1G01010.1` to the same gene entry.

### Mapping from OrthoFinder output

If you have OrthoFinder orthogroups, one-to-many mappings can be extracted
with the following script, adjusting the column names to match the species
headers in your output:

```python
import csv

with open("Orthogroups/Orthogroups.tsv") as f:
    reader = csv.DictReader(f, delimiter="\t")
    print("transcript_id\ttarget_id")
    for row in reader:
        src = [g.strip() for g in row["Saccharum"].split(",") if g.strip()]
        tgt = [g.strip() for g in row["Arabidopsis"].split(",") if g.strip()]
        for s in src:
            for t in tgt:
                print(f"{s}\t{t}")
```

---

## Workflow with the WGCNA tool

This tool connects directly to two outputs of the WGCNA Galaxy Tool.

```
VST matrix + sample metadata
        |
        v
 +--------------+
 | WGCNA tool   |--- out_table_modules --------+
 +--------------+--- out_table_sel_modules --+ |
                                             | |
                    ID mapping file          | |
                    (BLAST / OrthoFinder)    | |
                              |              | |
                              v              v v
                    +------------------------------+
                    |  GO Enrichment tool (this)   |
                    +------------------------------+
                              |
               +--------------+--------------+
               v              v              v
         GO table         Dotplot       Per-module
        (combined)        Barplot         plots
```

Using `out_table_sel_modules` as the module filter restricts enrichment to
modules with significant trait correlation, reducing the number of tests and
focusing results on biologically relevant modules.

---

## Supported organisms

The following OrgDb packages are included in the default installation.
Additional packages can be installed with
`BiocManager::install("org.XX.eg.db")` and used by entering the package name
in the OrgDb parameter.

| OrgDb | Organism | Recommended key type |
|---|---|---|
| org.At.tair.db | Arabidopsis thaliana | TAIR |
| org.Hs.eg.db | Homo sapiens | ENTREZID or SYMBOL |
| org.Mm.eg.db | Mus musculus | ENTREZID or SYMBOL |
| org.Os.eg.db | Oryza sativa | ENTREZID |
| org.Sc.sgd.db | Saccharomyces cerevisiae | SGD |

For sugarcane analyses, use `org.At.tair.db` with key type `TAIR` in
combination with the ID mapping feature.

---

## Troubleshooting

**No significant GO terms found in any module**

Check the mapping rate in the ID mapping report. If it is below 30%, the BLAST
search may have been too stringent — try relaxing the e-value threshold to
`1e-3`. Alternatively, relax `pvalue_cutoff` to `0.1` and `qvalue_cutoff` to
`0.5` as a diagnostic test to confirm whether any signal exists.

**Error: column not found in mapping file**

The `Source column name` and `Target column name` parameters must match the
column headers in your mapping file exactly, including capitalization. Verify
with:

```bash
head -2 id_mapping.tsv
```

**Key type mismatch error from clusterProfiler**

The key type must match the ID format in your (mapped) gene list and be
recognized by the selected OrgDb. For `org.At.tair.db`, use `TAIR` with IDs
in the format `AT1G01010`. For `org.Hs.eg.db`, use `ENTREZID` with numeric
IDs. Mixing key types causes silent failures where all genes are dropped before
enrichment.

**OrgDb package not installed**

If the tool fails with a package not found error, install the missing OrgDb:

```r
BiocManager::install("org.XX.eg.db")
```

Then add it to the `<requirements>` block in `go_enrichment_macros.xml` so it
is available in the conda environment during Galaxy jobs.

**Plot output is very large or unreadable**

When many modules are analyzed simultaneously, the combined faceted plots can
become dense. Use the per-module plots collection for individual inspection, or
restrict the analysis to a smaller set of modules using the `manual` filter
mode or the selected-modules table from WGCNA.

---

## Citation

If you use this tool in your research, please cite:

**This tool:**
> Romão, L. O. & Vicentini, R. (2026). GO Enrichment Galaxy Tool.
> GitHub: https://github.com/Luiza-Romao/go-enrichment-galaxy
> DOI: [Zenodo DOI — to be added]

**clusterProfiler:**
> Wu T, Hu E, Xu S, Chen M, Guo P, Dai Z, Feng T, Zhou L, Tang W, Zhan L,
> Fu X, Liu S, Bo X, Yu G (2021). clusterProfiler 4.0: A universal enrichment
> tool for interpreting omics data. *The Innovation*, 2(3), 100141.
> https://doi.org/10.1016/j.xinn.2021.100141

> Yu G, Wang L, Han Y, He Q (2012). clusterProfiler: an R Package for
> Comparing Biological Themes Among Gene Clusters. *OMICS: A Journal of
> Integrative Biology*, 16(5), 284–287.
> https://doi.org/10.1089/omi.2011.0118

**Galaxy platform:**
> The Galaxy Community. The Galaxy platform for accessible, reproducible, and
> collaborative data analyses: 2024 update. *Nucleic Acids Research* 2024,
> **52**(W1):W83–W94. https://doi.org/10.1093/nar/gkac995

---

## Authors and funding

Developed by **Luiza Oliveira Romão** under FAPESP doctoral fellowship
2025/08740-0, supervised by **Prof. Dr. Renato Vicentini**.

Laboratório de Bioinformática e Biologia de Sistemas
Instituto de Biologia — Universidade Estadual de Campinas (UNICAMP)

Developed within the **Centro de Melhoramento Molecular de Plantas (CeM²P)**,
funded by FAPESP grant 2022/04006-2.

---

## License

MIT License. See [LICENSE](LICENSE) for details.
