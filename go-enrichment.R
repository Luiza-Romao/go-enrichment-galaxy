#!/usr/bin/env Rscript
# =============================================================================
# GO Enrichment Analysis — Galaxy Tool
# Companion tool for the WGCNA Galaxy pipeline
#
# Performs per-module Gene Ontology enrichment using clusterProfiler.
# Supports transcript ID mapping for non-model organisms (e.g. sugarcane),
# allowing "transcript_XXXXXX" IDs to be converted to ortholog IDs in a
# model organism before enrichment.
#
# Inputs:
#   --module_table     : Gene-module assignment table from WGCNA (Gene, Module)
#   --id_mapping       : Optional two-column TSV mapping transcript IDs to
#                        target organism gene IDs (e.g. TAIR IDs for Arabidopsis)
#   --sel_modules      : Optional selected-modules table from WGCNA; if provided,
#                        only modules listed there are analyzed
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(ggrepel)
})

# ── Helper functions ──────────────────────────────────────────────────────────

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)
  cat(msg, "\n")
  msg
}

placeholder_png <- function(path, msg = "No data") {
  png(path, width = 600, height = 200)
  par(mar = c(0, 0, 0, 0))
  plot.new()
  text(0.5, 0.5, msg, cex = 1.2, col = "grey50")
  dev.off()
}

safe_enrichGO <- function(genes, orgdb, keytype, ont, pval, qval, min_gs, max_gs) {
  tryCatch(
    enrichGO(
      gene          = genes,
      OrgDb         = orgdb,
      keyType       = keytype,
      ont           = ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = pval,
      qvalueCutoff  = qval,
      minGSSize     = min_gs,
      maxGSSize     = max_gs,
      readable      = FALSE
    ),
    error = function(e) {
      log_msg("  WARNING: enrichGO failed for ontology ", ont, ": ", conditionMessage(e))
      NULL
    }
  )
}

# ── Option parsing ────────────────────────────────────────────────────────────

option_list <- list(
  # Primary inputs
  make_option("--module_table",      type = "character"),
  make_option("--id_mapping",        type = "character", default = "None"),
  make_option("--mapping_from_col",  type = "character", default = "transcript_id",
              help = "Column name in the mapping file containing source IDs (e.g. sugarcane transcript IDs)"),
  make_option("--mapping_to_col",    type = "character", default = "target_id",
              help = "Column name in the mapping file containing target IDs (e.g. Arabidopsis TAIR IDs)"),
  make_option("--sel_modules",       type = "character", default = "None",
              help = "Optional: selected-modules table from WGCNA to restrict analysis"),

  # Module filter
  make_option("--modules",           type = "character", default = "all",
              help = "Comma-separated module names to analyze, or 'all'"),
  make_option("--min_genes_per_mod", type = "integer",   default = 10,
              help = "Skip modules with fewer mapped genes than this threshold"),

  # OrgDb configuration
  make_option("--orgdb",             type = "character", default = "org.At.tair.db"),
  make_option("--keytype",           type = "character", default = "TAIR"),
  make_option("--ontology",          type = "character", default = "ALL",
              help = "GO ontology: BP, MF, CC, or ALL"),

  # Enrichment thresholds
  make_option("--pvalue_cutoff",     type = "double",    default = 0.05),
  make_option("--qvalue_cutoff",     type = "double",    default = 0.2),
  make_option("--min_gs_size",       type = "integer",   default = 10),
  make_option("--max_gs_size",       type = "integer",   default = 500),
  make_option("--top_terms",         type = "integer",   default = 15,
              help = "Number of top GO terms to show per module in plots"),

  # Outputs
  make_option("--out_table_combined",  type = "character"),
  make_option("--out_table_mapping",   type = "character"),
  make_option("--out_plot_dot",        type = "character"),
  make_option("--out_plot_bar",        type = "character"),
  make_option("--out_plots_dir",       type = "character", default = "go_plots"),
  make_option("--out_log",             type = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Redirect log
log_lines <- c()
log_msg_save <- function(...) {
  m <- log_msg(...)
  log_lines <<- c(log_lines, m)
}

log_msg_save("GO Enrichment Analysis — starting")

# ── STEP 1: Load and validate module table ────────────────────────────────────

log_msg_save("STEP 1: Loading module gene assignments")

if (!file.exists(opt$module_table)) {
  stop("Module table not found: ", opt$module_table)
}

mod_table <- read.table(
  opt$module_table,
  header      = TRUE,
  sep         = "\t",
  stringsAsFactors = FALSE,
  quote       = ""
)

# Detect Gene and Module columns flexibly
gene_col <- intersect(c("Gene", "gene", "GeneID", "gene_id", "transcript_id"), colnames(mod_table))[1]
mod_col  <- intersect(c("Module", "module", "Color", "color"), colnames(mod_table))[1]

if (is.na(gene_col) || is.na(mod_col)) {
  stop(
    "Could not detect Gene and Module columns in the module table. ",
    "Expected columns named 'Gene'/'GeneID' and 'Module'/'Color'. ",
    "Found: ", paste(colnames(mod_table), collapse = ", ")
  )
}

log_msg_save("  Detected gene column: '", gene_col, "', module column: '", mod_col, "'")
log_msg_save("  Total genes in module table: ", nrow(mod_table))
log_msg_save("  Modules detected: ", length(unique(mod_table[[mod_col]])))

# Remove grey module (unassigned genes — not biologically meaningful)
mod_table <- mod_table[mod_table[[mod_col]] != "grey", ]
log_msg_save("  Genes after removing grey module: ", nrow(mod_table))

# ── STEP 2: Parse module selection ───────────────────────────────────────────

log_msg_save("STEP 2: Selecting modules to analyze")

all_modules <- unique(mod_table[[mod_col]])

if (!is.null(opt$sel_modules) && opt$sel_modules != "None" && file.exists(opt$sel_modules)) {
  log_msg_save("  Using selected-modules table from WGCNA to restrict analysis")
  sel_tbl <- read.table(opt$sel_modules, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE, quote = "")
  sel_col <- intersect(c("Module", "module", "Color", "color"), colnames(sel_tbl))[1]
  if (!is.na(sel_col)) {
    selected_modules <- intersect(unique(sel_tbl[[sel_col]]), all_modules)
    log_msg_save("  Modules from selected-modules table: ", paste(selected_modules, collapse = ", "))
  } else {
    selected_modules <- all_modules
    log_msg_save("  WARNING: could not detect module column in selected-modules table; using all modules")
  }
} else if (tolower(trimws(opt$modules)) == "all") {
  selected_modules <- all_modules
  log_msg_save("  Using all modules: ", length(selected_modules), " modules")
} else {
  selected_modules <- trimws(unlist(strsplit(opt$modules, ",")))
  selected_modules <- intersect(selected_modules, all_modules)
  log_msg_save("  User-specified modules: ", paste(selected_modules, collapse = ", "))
  if (length(selected_modules) == 0) {
    stop("None of the specified module names were found in the module table.")
  }
}

# ── STEP 3: Load and apply ID mapping ────────────────────────────────────────

log_msg_save("STEP 3: ID mapping")

use_mapping <- !is.null(opt$id_mapping) &&
               opt$id_mapping != "None" &&
               file.exists(opt$id_mapping)

mapping_report <- data.frame(
  transcript_id = character(),
  target_id     = character(),
  module        = character(),
  stringsAsFactors = FALSE
)

if (use_mapping) {
  log_msg_save("  Loading ID mapping file: ", opt$id_mapping)
  mapping <- read.table(
    opt$id_mapping,
    header           = TRUE,
    sep              = "\t",
    stringsAsFactors = FALSE,
    quote            = ""
  )

  if (!opt$mapping_from_col %in% colnames(mapping)) {
    stop(
      "Column '", opt$mapping_from_col, "' not found in mapping file. ",
      "Found columns: ", paste(colnames(mapping), collapse = ", ")
    )
  }
  if (!opt$mapping_to_col %in% colnames(mapping)) {
    stop(
      "Column '", opt$mapping_to_col, "' not found in mapping file. ",
      "Found columns: ", paste(colnames(mapping), collapse = ", ")
    )
  }

  # Keep only the two relevant columns and deduplicate
  mapping <- unique(mapping[, c(opt$mapping_from_col, opt$mapping_to_col)])
  colnames(mapping) <- c("source_id", "target_id")

  log_msg_save("  Mapping file rows: ", nrow(mapping))
  log_msg_save("  Unique source IDs: ", length(unique(mapping$source_id)))
  log_msg_save("  Unique target IDs: ", length(unique(mapping$target_id)))

  # Join with module table
  mod_table_mapped <- merge(
    mod_table,
    mapping,
    by.x = gene_col,
    by.y = "source_id",
    all.x = FALSE
  )

  mapping_report <- data.frame(
    transcript_id = mod_table_mapped[[gene_col]],
    target_id     = mod_table_mapped$target_id,
    module        = mod_table_mapped[[mod_col]],
    stringsAsFactors = FALSE
  )

  n_mapped   <- length(unique(mod_table_mapped[[gene_col]]))
  n_total    <- length(unique(mod_table[[gene_col]]))
  pct_mapped <- round(100 * n_mapped / n_total, 1)

  log_msg_save("  Transcripts successfully mapped: ", n_mapped,
               " / ", n_total, " (", pct_mapped, "%)")

  # Use target IDs for enrichment
  get_genes <- function(module_name) {
    rows <- mod_table_mapped[mod_table_mapped[[mod_col]] == module_name, ]
    unique(rows$target_id)
  }

} else {
  log_msg_save("  No mapping file provided — using gene IDs from module table directly")
  get_genes <- function(module_name) {
    rows <- mod_table[mod_table[[mod_col]] == module_name, ]
    unique(rows[[gene_col]])
  }
}

# Write mapping report
if (nrow(mapping_report) > 0) {
  write.table(mapping_report, opt$out_table_mapping,
              sep = "\t", row.names = FALSE, quote = FALSE)
} else {
  write.table(
    data.frame(note = "No ID mapping was applied — gene IDs used directly"),
    opt$out_table_mapping, sep = "\t", row.names = FALSE, quote = FALSE
  )
}

# ── STEP 4: Load OrgDb ───────────────────────────────────────────────────────

log_msg_save("STEP 4: Loading OrgDb package: ", opt$orgdb)

if (!requireNamespace(opt$orgdb, quietly = TRUE)) {
  stop(
    "OrgDb package '", opt$orgdb, "' is not installed. ",
    "Install it with: BiocManager::install('", opt$orgdb, "')"
  )
}
suppressPackageStartupMessages(library(opt$orgdb, character.only = TRUE))
orgdb_obj <- get(opt$orgdb)
log_msg_save("  OrgDb loaded successfully")

# ── STEP 5: Determine ontologies to run ──────────────────────────────────────

if (toupper(opt$ontology) == "ALL") {
  ontologies <- c("BP", "MF", "CC")
} else {
  ontologies <- toupper(trimws(unlist(strsplit(opt$ontology, ","))))
}
log_msg_save("STEP 5: Ontologies to analyze: ", paste(ontologies, collapse = ", "))

# ── STEP 6: Per-module enrichment ────────────────────────────────────────────

log_msg_save("STEP 6: Running GO enrichment per module")

dir.create(opt$out_plots_dir, showWarnings = FALSE)

all_results  <- list()
plot_objects <- list()

for (mod in selected_modules) {

  log_msg_save("  Module: ", mod)

  gene_ids <- get_genes(mod)
  log_msg_save("    Genes available for enrichment: ", length(gene_ids))

  if (length(gene_ids) < opt$min_genes_per_mod) {
    log_msg_save("    SKIP: fewer than ", opt$min_genes_per_mod, " mapped genes")
    next
  }

  mod_results <- list()

  for (ont in ontologies) {

    res <- safe_enrichGO(
      genes   = gene_ids,
      orgdb   = orgdb_obj,
      keytype = opt$keytype,
      ont     = ont,
      pval    = opt$pvalue_cutoff,
      qval    = opt$qvalue_cutoff,
      min_gs  = opt$min_gs_size,
      max_gs  = opt$max_gs_size
    )

    if (is.null(res) || nrow(as.data.frame(res)) == 0) {
      log_msg_save("    No significant terms for ontology ", ont)
      next
    }

    df <- as.data.frame(res)
    df$Module   <- mod
    df$Ontology <- ont
    mod_results[[ont]] <- df

    log_msg_save("    ", ont, ": ", nrow(df), " significant terms")
  }

  if (length(mod_results) == 0) {
    log_msg_save("    No significant enrichment in any ontology for module ", mod)
    next
  }

  # Combine across ontologies for this module
  mod_df <- do.call(rbind, mod_results)
  all_results[[mod]] <- mod_df

  # Per-module dotplot (faceted by ontology)
  top_terms_df <- do.call(rbind, lapply(mod_results, function(d) {
    d[order(d$p.adjust), ][seq_len(min(opt$top_terms, nrow(d))), ]
  }))

  top_terms_df$Description <- factor(
    top_terms_df$Description,
    levels = rev(unique(top_terms_df$Description))
  )
  top_terms_df$GeneRatio_num <- sapply(
    top_terms_df$GeneRatio,
    function(x) eval(parse(text = x))
  )

  p_dot <- ggplot(top_terms_df, aes(
    x    = GeneRatio_num,
    y    = Description,
    size = Count,
    color = p.adjust
  )) +
    geom_point() +
    scale_color_gradient(low = "#1D9E75", high = "#D85A30",
                         name = "Adj. p-value") +
    scale_size_continuous(name = "Gene count", range = c(2, 8)) +
    facet_wrap(~ Ontology, scales = "free_y", ncol = 1) +
    labs(
      title    = paste0("GO Enrichment — Module: ", mod),
      x        = "Gene ratio",
      y        = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "#E6F1FB"),
      strip.text       = element_text(face = "bold"),
      plot.title       = element_text(size = 12, face = "bold"),
      axis.text.y      = element_text(size = 9)
    )

  p_bar <- ggplot(top_terms_df, aes(
    x    = reorder(Description, Count),
    y    = Count,
    fill = p.adjust
  )) +
    geom_col() +
    coord_flip() +
    scale_fill_gradient(low = "#1D9E75", high = "#D85A30",
                        name = "Adj. p-value") +
    facet_wrap(~ Ontology, scales = "free_y", ncol = 1) +
    labs(
      title = paste0("GO Enrichment — Module: ", mod),
      x     = NULL,
      y     = "Gene count"
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "#E6F1FB"),
      strip.text       = element_text(face = "bold"),
      plot.title       = element_text(size = 12, face = "bold"),
      axis.text.y      = element_text(size = 9)
    )

  plot_objects[[mod]] <- list(dot = p_dot, bar = p_bar)

  # Save individual module plots
  n_onts <- length(unique(top_terms_df$Ontology))
  n_terms <- min(nrow(top_terms_df), opt$top_terms * n_onts)
  plot_height <- max(4, n_terms * 0.28 + 1.5)

  ggsave(
    filename = file.path(opt$out_plots_dir, paste0("dotplot_", mod, ".png")),
    plot     = p_dot,
    width    = 8,
    height   = plot_height,
    dpi      = 150,
    limitsize = FALSE
  )
  ggsave(
    filename = file.path(opt$out_plots_dir, paste0("barplot_", mod, ".png")),
    plot     = p_bar,
    width    = 8,
    height   = plot_height,
    dpi      = 150,
    limitsize = FALSE
  )
}

log_msg_save("  Modules with significant enrichment: ", length(all_results))

# ── STEP 7: Combined outputs ──────────────────────────────────────────────────

log_msg_save("STEP 7: Writing combined outputs")

if (length(all_results) == 0) {
  log_msg_save("  WARNING: No significant GO terms found in any module")

  write.table(
    data.frame(note = paste(
      "No significant GO enrichment found.",
      "Consider relaxing pvalue_cutoff or qvalue_cutoff,",
      "or checking that the ID mapping is correct."
    )),
    opt$out_table_combined, sep = "\t", row.names = FALSE, quote = FALSE
  )
  placeholder_png(opt$out_plot_dot, "No significant GO terms found")
  placeholder_png(opt$out_plot_bar, "No significant GO terms found")

} else {

  combined_df <- do.call(rbind, all_results)
  rownames(combined_df) <- NULL

  # Reorder columns for readability
  priority_cols <- c("Module", "Ontology", "ID", "Description",
                     "GeneRatio", "BgRatio", "pvalue", "p.adjust",
                     "qvalue", "geneID", "Count")
  other_cols    <- setdiff(colnames(combined_df), priority_cols)
  combined_df   <- combined_df[, c(priority_cols[priority_cols %in% colnames(combined_df)],
                                    other_cols)]

  write.table(combined_df, opt$out_table_combined,
              sep = "\t", row.names = FALSE, quote = FALSE)

  log_msg_save("  Combined table rows: ", nrow(combined_df))

  # Combined dotplot (all modules, top N terms each)
  top_combined <- do.call(rbind, lapply(all_results, function(d) {
    do.call(rbind, lapply(split(d, d$Ontology), function(x) {
      x[order(x$p.adjust), ][seq_len(min(5, nrow(x))), ]
    }))
  }))

  top_combined$GeneRatio_num <- sapply(
    top_combined$GeneRatio,
    function(x) eval(parse(text = x))
  )
  top_combined$ModOnt <- paste0(top_combined$Module, "\n(", top_combined$Ontology, ")")

  n_panels <- length(unique(top_combined$Module))
  ncols    <- min(3, n_panels)
  nrows    <- ceiling(n_panels / ncols)
  fig_w    <- ncols * 5.5
  fig_h    <- max(4, nrows * 5)

  p_combined_dot <- ggplot(top_combined, aes(
    x     = GeneRatio_num,
    y     = reorder(Description, GeneRatio_num),
    size  = Count,
    color = p.adjust
  )) +
    geom_point() +
    scale_color_gradient(low = "#1D9E75", high = "#D85A30",
                         name = "Adj. p-value") +
    scale_size_continuous(name = "Gene count", range = c(2, 7)) +
    facet_wrap(~ Module + Ontology, scales = "free_y", ncol = ncols) +
    labs(x = "Gene ratio", y = NULL,
         title = "GO Enrichment — all modules") +
    theme_bw(base_size = 10) +
    theme(
      strip.background = element_rect(fill = "#E6F1FB"),
      strip.text       = element_text(face = "bold", size = 8),
      axis.text.y      = element_text(size = 8),
      plot.title       = element_text(size = 12, face = "bold")
    )

  p_combined_bar <- ggplot(top_combined, aes(
    x    = reorder(Description, Count),
    y    = Count,
    fill = p.adjust
  )) +
    geom_col() +
    coord_flip() +
    scale_fill_gradient(low = "#1D9E75", high = "#D85A30",
                        name = "Adj. p-value") +
    facet_wrap(~ Module + Ontology, scales = "free", ncol = ncols) +
    labs(x = NULL, y = "Gene count",
         title = "GO Enrichment — all modules") +
    theme_bw(base_size = 10) +
    theme(
      strip.background = element_rect(fill = "#E6F1FB"),
      strip.text       = element_text(face = "bold", size = 8),
      axis.text.y      = element_text(size = 8),
      plot.title       = element_text(size = 12, face = "bold")
    )

  ggsave(opt$out_plot_dot, p_combined_dot,
         width = fig_w, height = fig_h, dpi = 150, limitsize = FALSE)
  ggsave(opt$out_plot_bar, p_combined_bar,
         width = fig_w, height = fig_h, dpi = 150, limitsize = FALSE)

  log_msg_save("  Combined dotplot saved: ", opt$out_plot_dot)
  log_msg_save("  Combined barplot saved: ", opt$out_plot_bar)
}

# ── STEP 8: Write log ─────────────────────────────────────────────────────────

log_msg_save("STEP 8: Analysis complete")
writeLines(log_lines, opt$out_log)
