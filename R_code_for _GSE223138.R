
#  Dataset  : GSE223138(Parkinson's Disease vs Healthy Control)

# ================================================================================

# ── REPRODUCIBILITY ───────────────────────────────────────────────────────────
set.seed(123)
cat("✔ Seed set: 123\n\n")

# ── USER SETTINGS ─────────────────────────────────────────────────────────────
GEO_ACCESSION <- "GSE223138"
MIN_FEATURES  <- 200
MAX_FEATURES  <- 6000
MAX_MT        <- 20
MIN_COUNTS    <- 500
CHOSEN_RES    <- 0.6    # look at results/clustering/01_umap_all_resolutions.pdf
# and change this if needed before re-running clustering


suppressPackageStartupMessages({
  library(Seurat);  library(ggplot2);  library(dplyr);    library(patchwork)
  library(ggrepel); library(Matrix);   library(GEOquery); library(R.utils)
  library(harmony); library(DoubletFinder)
  library(SingleR); library(celldex);  library(DESeq2)
  library(lisi);    library(monocle3); library(SeuratWrappers)
})
# some remaining libraries 
if (!requireNamespace("Biobase", quietly = TRUE)) BiocManager::install("Biobase")
library(Biobase)
library(GEOquery)
# ── CREATE FOLDERS ────────────────────────────────────────────────────────────
for (d in c("data","objects",
            "results/qc","results/doublet","results/cellcycle",
            "results/integration","results/processing",
            "results/clustering","results/annotation",
            "results/trajectory","results/DEG"))
  dir.create(d, recursive=TRUE, showWarnings=FALSE)

cat("\n╔══════════════════════════════════════════════════════════╗\n")
cat("║   Complete scRNA-seq Pipeline — GSE202210              ║\n")
cat("║   Parkinson's Disease vs Healthy Control               ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")

# ================================================================================
#  STEP 1: DOWNLOAD DATA FROM NCBI GEO
# ================================================================================
cat("━━━ STEP 1: Download Data ━━━\n")

gse <- getGEO(GEO_ACCESSION, GSEMatrix=TRUE, getGPL=FALSE)
cat("  Title  :", gse[[1]]@experimentData@title, "\n")

# ── FIX 3: Build disease labels from GEO metadata (not string guessing) ───────
sample_meta <- Biobase::pData(gse[[1]])

# Extract useful columns safely
meta_cols   <- intersect(c("title","geo_accession","source_name_ch1",
                           "characteristics_ch1","characteristics_ch1.1",
                           "characteristics_ch1.2"),
                         colnames(sample_meta))
sample_meta <- sample_meta[, meta_cols, drop=FALSE]
write.csv(sample_meta, "data/sample_info.csv", row.names=FALSE)
cat("  Sample metadata saved → data/sample_info.csv\n")
cat("  Sample titles:\n")
print(sample_meta$title)

# Build disease map from GEO title
# GSE202210: HC = Healthy Control, PD = Parkinson's Disease
# We search across ALL metadata columns for HC/PD keywords
all_meta_text <- apply(sample_meta, 1, paste, collapse=" ")
disease_map   <- ifelse(
  grepl("HC|healthy|control|normal", all_meta_text, ignore.case=TRUE) &
    !grepl("PD|Parkinson|disease|patient", all_meta_text, ignore.case=TRUE),
  "Healthy_Control",
  "Parkinsons_Disease"
)
names(disease_map) <- sample_meta$geo_accession
cat("\n  Disease assignments from metadata:\n")
print(data.frame(geo_id  = names(disease_map),
                 title   = sample_meta$title,
                 disease = disease_map))

# Download supplementary files
cat("\n  Downloading count files...\n")
getGEOSuppFiles(GEO_ACCESSION, makeDirectory=FALSE, baseDir="data")

# Extract
for (f in list.files("data", full.names=TRUE)) {
  if      (grepl("\\.tar\\.gz$|\\.tgz$", f)) { cat("  Extracting:", basename(f),"\n"); untar(f, exdir="data") }
  else if (grepl("\\.tar$", f))              { cat("  Extracting:", basename(f),"\n"); untar(f, exdir="data") }
  else if (grepl("\\.gz$",  f) && !grepl("\\.tar", f))
    tryCatch(gunzip(f, remove=FALSE, overwrite=TRUE), error=function(e) NULL)
}

# ── FIX 2: recursive=TRUE — works whether files are flat or nested ─────────────
all_files <- list.files("data", pattern="\\.gz$",
                        recursive=TRUE, full.names=FALSE)
all_files <- all_files[!grepl("sample_info|RAW\\.tar", all_files)]
cat("  Files found:", length(all_files), "\n")

# Extract prefixes — handle both flat and nested paths
base_names <- basename(all_files)
prefixes   <- unique(gsub("_(barcodes|genes|features|matrix).*", "", base_names))
prefixes   <- prefixes[grepl("^GSM", prefixes)]
cat("  Samples found:", length(prefixes), "\n\n")
if (length(prefixes) == 0) stop("No sample files found. Check download.")

# Organise into subfolders — rename genes.tsv.gz → features.tsv.gz
cat("  Organising files into per-sample folders...\n")
for (pfx in prefixes) {
  sdir <- file.path("data", pfx)
  dir.create(sdir, showWarnings=FALSE)
  
  # Find files anywhere in data/ (recursive search)
  find_file <- function(pattern) {
    hits <- list.files("data", pattern=pattern, recursive=TRUE, full.names=TRUE)
    hits[grepl(pfx, hits)][1]
  }
  
  bc  <- find_file(paste0(pfx, "_barcodes"))
  gn  <- find_file(paste0(pfx, "_genes|", pfx, "_features"))
  mt  <- find_file(paste0(pfx, "_matrix"))
  
  if (!is.na(bc))  file.copy(bc,  file.path(sdir,"barcodes.tsv.gz"),  overwrite=TRUE)
  if (!is.na(gn))  file.copy(gn,  file.path(sdir,"features.tsv.gz"),  overwrite=TRUE)
  if (!is.na(mt))  file.copy(mt,  file.path(sdir,"matrix.mtx.gz"),    overwrite=TRUE)
  
  ok <- all(file.exists(file.path(sdir,
                                  c("barcodes.tsv.gz","features.tsv.gz","matrix.mtx.gz"))))
  cat(sprintf("    %-38s %s\n", pfx, ifelse(ok,"✔","✘ check files")))
}
# ── 2. ROBUST DISEASE LABELING ────────────────────────────────────────────────
cat("  Scanning metadata for disease labels...\n")

diag_col_idx <- which(apply(sample_meta, 2, function(x) any(grepl("diagnosis", x, ignore.case=TRUE))))[1]

if (!is.na(diag_col_idx)) {
  diag_col_name <- colnames(sample_meta)[diag_col_idx]
  cat("  Found diagnosis info in column:", diag_col_name, "\n")
  disease_map <- ifelse(grepl("control|healthy|normal|HC", sample_meta[[diag_col_name]], ignore.case=TRUE), 
                        "Healthy_Control", "Parkinsons_Disease")
} else {
  cat("  Warning: 'diagnosis' keyword not found. Using Title-based fallback.\n")
  ctrl_titles <- c("C7","C9","D1","D10","D12","D13","D15","D2","D4","D5","D6","D7")
  disease_map <- ifelse(sample_meta$title %in% ctrl_titles, "Healthy_Control", "Parkinsons_Disease")
}

names(disease_map) <- sample_meta$geo_accession
# Ensure the data directory exists BEFORE saving the CSV
cat("✔ Download complete.\n\n")

# ── USER SETTINGS (OPTIMIZED COMPROMISE) ──
GEO_ACCESSION <- "GSE184950"
N_CORES       <- 10

MIN_FEATURES  <- 300    # Your chosen floor
MAX_FEATURES  <- 3000   # Your conservative doublet filter

MIN_COUNTS    <- 300    # Your chosen signal floor

MAX_MT        <- 10     # Standard limit for healthy brain nuclei
CHOSEN_RES    <- 0.6    # look at results/clustering/01_umap_all_resolutions.pdf
# and change this if needed before re-running clustering

# ================================================================================
#  STEP 2: LOAD DATA + QUALITY CONTROL
# ================================================================================
cat("━━━ STEP 2: Load Data + QC ━━━\n")

# Get the names of the 34 folders we created in Step 1
prefixes <- list.dirs("data", full.names=FALSE, recursive=FALSE)
if (length(prefixes) == 0) stop("No folders found in 'data'. Did Step 1 finish?")

seurat_raw <- list()

for (pfx in prefixes) {
  sdir  <- file.path("data", pfx)
  
  # Extract the GSM ID (e.g., "GSM5602315") and the short name (e.g., "A10")
  parts  <- strsplit(pfx, "_")[[1]]
  gsm_id <- parts[1]
  sname  <- ifelse(length(parts) > 1, parts[2], gsm_id)
  
  cat("  Loading:", sname, "...")
  
  tryCatch({
    counts <- Read10X(data.dir=sdir)
    if (is.list(counts)) counts <- counts[[1]]
    
    # Create the raw Seurat object
    obj <- CreateSeuratObject(counts=counts, project=sname,
                              min.cells=3, min.features=200)
    
    # ── THE MAGIC: Apply the perfect labels we just verified ──
    obj$sample  <- sname
    obj$gsm_id  <- gsm_id
    
    # Match the GSM ID to our disease_map from Step 1
    if (gsm_id %in% names(disease_map)) {
      obj$disease <- disease_map[[gsm_id]]
    } else {
      obj$disease <- "Unknown"
      cat(" ⚠ Label missing!")
    }
    
    seurat_raw[[sname]] <- obj
    cat(" ✔", ncol(obj), "cells |", obj$disease[1], "\n")
  }, error=function(e) cat(" ✘ FAILED:", e$message, "\n"))
}

# Print the overall breakdown
cat("\n  Total raw cells loaded:", sum(sapply(seurat_raw, ncol)), "\n")
cat("  Disease breakdown across all cells:\n")
print(table(sapply(seurat_raw, function(x) x$disease[1])))

# ── QC METRICS & PLOTTING ───────────────────────────────────────────────────
cat("\n  Calculating QC metrics (Mitochondrial %)...\n")

# Calculate Mitochondrial %
seurat_raw <- lapply(seurat_raw, function(seu) {
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern="^MT-"); seu
})

# Create results folders
dir.create("results/qc", recursive=TRUE, showWarnings=FALSE)

# Combine briefly just for the "Before" plot
all_raw <- merge(seurat_raw[[1]], y=seurat_raw[-1], add.cell.ids=names(seurat_raw))

pdf("results/qc/01_combined_violin_before.pdf", width=14, height=5)
print(VlnPlot(all_raw, features=c("nFeature_RNA","nCount_RNA","percent.mt"),
              ncol=3, pt.size=0, group.by="orig.ident") +
        plot_annotation(title="QC — All Samples BEFORE Filtering",
                        theme=theme(plot.title=element_text(face="bold"))))
dev.off()
cat("  ✔ Saved: results/qc/01_combined_violin_before.pdf\n")

# ── FILTERING ───────────────────────────────────────────────────────────────
# Note: Using your thresholds defined at the top of your script
cat(sprintf("\n  Filtering: genes %d–%d | MT <%.0f%% | counts >%d\n",
            MIN_FEATURES, MAX_FEATURES, MAX_MT, MIN_COUNTS))

seurat_filtered <- lapply(seurat_raw, function(seu) {
  subset(seu, subset = nFeature_RNA > MIN_FEATURES &
           nFeature_RNA < MAX_FEATURES &
           percent.mt   < MAX_MT       &
           nCount_RNA   > MIN_COUNTS)
})

# ── FILTERING SUMMARY ───────────────────────────────────────────────────────
qc_df <- data.frame(
  sample       = names(seurat_filtered),
  disease      = sapply(seurat_filtered, function(x) x$disease[1]),
  cells_before = sapply(seurat_raw,      ncol),
  cells_after  = sapply(seurat_filtered, ncol),
  median_genes = sapply(seurat_filtered, function(x) round(median(x$nFeature_RNA))),
  median_mt    = sapply(seurat_filtered, function(x) round(median(x$percent.mt),1))
)

cat("\n  --- After QC Summary ---\n")
print(head(qc_df))
write.csv(qc_df,"results/qc/qc_summary.csv", row.names=FALSE)

# Plot "After" QC
all_filt <- merge(seurat_filtered[[1]], y=seurat_filtered[-1], add.cell.ids=names(seurat_filtered))

pdf("results/qc/03_combined_violin_after.pdf", width=14, height=5)
print(VlnPlot(all_filt, features=c("nFeature_RNA","nCount_RNA","percent.mt"),
              ncol=3, pt.size=0, group.by="orig.ident") +
        plot_annotation(title="QC — All Samples AFTER Filtering",
                        theme=theme(plot.title=element_text(face="bold"))))
dev.off()
cat("  ✔ Saved: results/qc/03_combined_violin_after.pdf\n")

# Save the checkpoint
dir.create("objects", showWarnings=FALSE)
saveRDS(seurat_filtered,"objects/01_seurat_qc.rds")

cat(sprintf("\n✔ Step 2 QC Complete! Total cells surviving QC: %d\n\n",
            sum(sapply(seurat_filtered, ncol))))


