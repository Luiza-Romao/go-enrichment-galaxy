# Test data

Synthetic test fixtures for the GO Enrichment Galaxy tool. All Arabidopsis
gene IDs are real TAIR locus identifiers, drawn from chromosomes 1–5 of
*A. thaliana*. The transcript IDs (`transcript_NNNNNN`) are fictitious and
exist only to exercise the ID mapping feature.

## Files

| File | Purpose |
|---|---|
| `test_modules.tabular` | Module gene assignments using TAIR IDs directly. Equivalent to the WGCNA tool output when working with a model organism. Contains four colour modules plus a grey module that should be ignored by the tool. |
| `test_selected_modules.tabular` | Subset of modules with mock trait correlations. Mirrors the format of the WGCNA tool's "Selected modules summary" output. |
| `test_modules_with_transcript_ids.tabular` | Module assignments using sugarcane-style transcript IDs. Used together with `test_id_mapping.tabular` to test the ID mapping pipeline. |
| `test_id_mapping.tabular` | Two-column mapping from transcript IDs to TAIR locus IDs. Equivalent to the BLAST best-hit output for non-model organisms. |

## How they are used

The `<tests>` block in `go-enrichment.xml` references `test_modules.tabular`.
The remaining three files cover the additional workflow paths exercised in
the manual smoke-test recipes documented in the project README.
