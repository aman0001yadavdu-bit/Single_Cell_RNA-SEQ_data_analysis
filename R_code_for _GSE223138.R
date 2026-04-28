#  Dataset  : GSE184950 (Parkinson's Disease vs Healthy Control)

# ── REPRODUCIBILITY ───────────────────────────────────────────────────────────
set.seed(123)
cat("✔ Seed set: 123\n\n")

# ── USER SETTINGS ─────────────────────────────────────────────────────────────
# ── USER SETTINGS (OPTIMIZED COMPROMISE) ──
GEO_ACCESSION <- "GSE184950"
N_CORES       <- 10

MIN_FEATURES  <- 300    # Your chosen floor
MAX_FEATURES  <- 3000   # Your conservative doublet filter

MIN_COUNTS    <- 500    # Your chosen signal floor

MAX_MT        <- 50     # Standard limit for healthy brain nuclei
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
cat("║   Complete scRNA-seq Pipeline — GSE202210                ║\n")
cat("║   Parkinson's Disease vs Healthy Control                 ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")

# ================================================================================
#  STEP 1: DOWNLOAD DATA FROM NCBI GEO (FULLY CORRECTED)
#  Dataset  : GSE184950 (Parkinson's Disease vs Healthy Control)
# ================================================================================
cat("━━━ STEP 1: Download Data ━━━\n")

library(Biobase)
library(GEOquery)

# Download the dataset metadata
gse <- getGEO(GEO_ACCESSION, GSEMatrix=TRUE, getGPL=FALSE)

# ── 1. MERGE PLATFORMS (Gets all 46 samples) ──────────────────────────────────
if (length(gse) > 1) {
  cat("  Multiple platforms detected. Merging PD and HC groups...\n")
  sample_meta <- do.call(rbind, lapply(gse, Biobase::pData))
} else {
  sample_meta <- Biobase::pData(gse[[1]])
}
cat("  Total samples loaded:", nrow(sample_meta), "\n")

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
if (!dir.exists("data")) dir.create("data")

# Save a clean version of the metadata for your records
meta_cols <- intersect(c("title","geo_accession","source_name_ch1"), colnames(sample_meta))
clean_meta <- sample_meta[, meta_cols, drop=FALSE]
write.csv(clean_meta, "data/sample_info.csv", row.names=FALSE)
cat("  Sample metadata saved → data/sample_info.csv\n")

check_df <- data.frame(
  GSM_ID            = names(disease_map),
  Assigned_Label    = disease_map,
  row.names         = NULL
)
cat("\n  --- Final Metadata Summary ---\n")
print(table(Assigned_Label = check_df$Assigned_Label))
cat("  -------------------------------\n")

# ── 3. FILE DOWNLOAD & SMART EXTRACTION ───────────────────────────────────────
cat("\n  Downloading count files...\n")
getGEOSuppFiles(GEO_ACCESSION, makeDirectory=FALSE, baseDir="data")

# Pass 1: Extract the master RAW.tar file
tar_files <- list.files("data", pattern="\\.tar$", full.names=TRUE)
for (f in tar_files) {
  cat("  Extracting master archive:", basename(f), "\n")
  untar(f, exdir="data")
}

# Pass 2: Extract each sample directly into its OWN folder
nested_tars <- list.files("data", pattern="\\.tar\\.gz$|\\.tgz$", full.names=TRUE)
cat("  Organizing files into Seurat-ready folders...\n")

for (f in nested_tars) {
  sample_name <- gsub("\\.tar\\.gz$|\\.tgz$", "", basename(f))
  sdir <- file.path("data", sample_name)
  dir.create(sdir, showWarnings=FALSE)
  
  untar(f, exdir=sdir)
  file.remove(f) 
  
  # ✅ THE FIX: Flatten folder structure safely
  extracted_files <- list.files(sdir, full.names=TRUE, recursive=TRUE)
  extracted_files <- extracted_files[!file.info(extracted_files)$isdir] # Isolate just the files
  
  # Step A: Move all files safely up to the main folder
  for (ef in extracted_files) {
    new_path <- file.path(sdir, basename(ef))
    if (ef != new_path) {
      file.rename(ef, new_path)
    }
  }
  
  # Step B: Now that files are safe, delete any leftover empty sub-folders
  subdirs <- list.dirs(sdir, full.names=TRUE, recursive=FALSE)
  if (length(subdirs) > 0) {
    unlink(subdirs, recursive=TRUE)
  }
  
  # Fix old 10x "genes.tsv" naming conventions
  genes_file <- file.path(sdir, "genes.tsv.gz")
  if (file.exists(genes_file)) {
    file.rename(genes_file, file.path(sdir, "features.tsv.gz"))
  }
  
  # Final check
  ok <- all(file.exists(file.path(sdir, c("barcodes.tsv.gz","features.tsv.gz","matrix.mtx.gz"))))
  cat(sprintf("    %-38s %s\n", sample_name, ifelse(ok,"✔","✘ check files")))
}

cat("✔ Download and Organization complete.\n\n")


# 1. Grab the metadata one more time
library(GEOquery)
gse <- getGEO("GSE184950", GSEMatrix=TRUE, getGPL=FALSE)
meta <- Biobase::pData(gse[[1]])

# 2. Look at the CORRECT column ("disease state:ch1") to find the controls
disease_map <- ifelse(grepl("Control", meta$`disease state:ch1`, ignore.case=TRUE), 
                      "Healthy_Control", "Parkinsons_Disease")

names(disease_map) <- meta$geo_accession

# 3. Print the true breakdown to the screen
cat("\n--- TRUE CLINICAL GROUPS ---\n")
print(table(disease_map))

# 4. Lock it in for Step 2
gsm_disease <- setNames(disease_map, names(disease_map))


# Create a clean table to view the exact assignments
assignment_view <- data.frame(
  GSM_ID      = meta$geo_accession,
  Sample_Name = meta$title,
  Group       = disease_map,
  row.names   = NULL
)

# Print all 34 rows so you can verify them
print(assignment_view)



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
