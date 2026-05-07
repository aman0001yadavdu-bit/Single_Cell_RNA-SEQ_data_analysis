# ======================================================================

#  Dataset  : GSE202210 (Parkinson's Disease vs Healthy Control)

# ================================================================================

# ── REPRODUCIBILITY ───────────────────────────────────────────────────────────
set.seed(123)
cat("✔ Seed set: 123\n\n")

# ── USER SETTINGS ─────────────────────────────────────────────────────────────
GEO_ACCESSION <- "GSE202210"
N_CORES       <- 10
MIN_FEATURES  <- 200
MAX_FEATURES  <- 6000
MAX_MT        <- 10
MIN_COUNTS    <- 500
CHOSEN_RES    <- 0.6    # look at results/clustering/01_umap_all_resolutions.pdf
# and change this if needed before re-running clustering

# ================================================================================
#  INSTALL PACKAGES
# ================================================================================
cat("Checking packages...\n")

cran_pkgs <- c("Seurat","ggplot2","dplyr","patchwork","ggrepel",
               "Matrix","remotes","BiocManager","R.utils")
for (p in cran_pkgs)
  if (!requireNamespace(p, quietly=TRUE))
    install.packages(p, repos="https://cloud.r-project.org", quiet=TRUE)

bioc_pkgs <- c("GEOquery","SingleR","celldex","DESeq2","BiocParallel")
for (p in bioc_pkgs)
  if (!requireNamespace(p, quietly=TRUE))
    BiocManager::install(p, ask=FALSE, quiet=TRUE)

if (!requireNamespace("harmony",       quietly=TRUE)) install.packages("harmony")
if (!requireNamespace("DoubletFinder", quietly=TRUE))
  remotes::install_github("chris-mcginnis-ucsf/DoubletFinder", quiet=TRUE)
if (!requireNamespace("lisi",          quietly=TRUE))
  remotes::install_github("immunogenomics/lisi", quiet=TRUE)
if (!requireNamespace("SeuratWrappers",quietly=TRUE))
  remotes::install_github("satijalab/seurat-wrappers", quiet=TRUE)
if (!requireNamespace("monocle3",      quietly=TRUE)) {
  for (p in c("BiocGenerics","DelayedArray","DelayedMatrixStats","limma",
              "S4Vectors","SingleCellExperiment","SummarizedExperiment",
              "batchelor","terra"))
    BiocManager::install(p, ask=FALSE, quiet=TRUE)
  remotes::install_github("cole-trapnell-lab/monocle3", quiet=TRUE)
}

suppressPackageStartupMessages({
  library(Seurat);  library(ggplot2);  library(dplyr);    library(patchwork)
  library(ggrepel); library(Matrix);   library(GEOquery); library(R.utils)
  library(harmony); library(DoubletFinder)
  library(SingleR); library(celldex);  library(DESeq2)
  library(lisi);    library(monocle3); library(SeuratWrappers)
})

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

# Collapse the extracted metadata columns into a single string per sample
all_meta_text <- apply(sample_meta, 2, as.character) 
all_meta_text <- apply(all_meta_text, 1, paste, collapse=" ")

# Positive match: If it has any of these control words, it's a Healthy Control.
# Otherwise, it's Parkinson's.
is_control <- grepl("HC|control|healthy|normal", all_meta_text, ignore.case = TRUE)

disease_map <- ifelse(is_control, "Healthy_Control", "Parkinsons_Disease")
names(disease_map) <- sample_meta$geo_accession

# --- NEW: Explicit Summary Printout ---
cat("\n  =========================================\n")
cat("  DISEASE ASSIGNMENT SUMMARY:\n")
cat("  =========================================\n")

hc_samples <- names(disease_map)[disease_map == "Healthy_Control"]
pd_samples <- names(disease_map)[disease_map == "Parkinsons_Disease"]

cat("  ▶ Healthy Control (", length(hc_samples), " samples):\n    ", 
    paste(hc_samples, collapse=", "), "\n\n", sep="")

cat("  ▶ Parkinson's Disease (", length(pd_samples), " samples):\n    ", 
    paste(pd_samples, collapse=", "), "\n\n", sep="")

# Keep the detailed table view just in case you need to check titles
print(data.frame(geo_id  = names(disease_map),
                 title   = sample_meta$title,
                 disease = disease_map))
# ================================================================================
#  STEP 2: LOAD DATA + QUALITY CONTROL
# ================================================================================
cat("━━━ STEP 2: Load Data + QC ━━━\n")

# Build GSM → disease lookup from metadata
gsm_disease <- setNames(disease_map, names(disease_map))

seurat_raw <- list()
for (pfx in prefixes) {
  sdir  <- file.path("data", pfx)
  sname <- gsub("^GSM[0-9]+_","", pfx)     # short biological name
  gsm_id <- regmatches(pfx, regexpr("GSM[0-9]+", pfx))  # e.g. GSM6106340
  cat("  Loading:", sname, "...")
  
  tryCatch({
    counts <- Read10X(data.dir=sdir)
    if (is.list(counts)) counts <- counts[[1]]
    obj <- CreateSeuratObject(counts=counts, project=sname,
                              min.cells=3, min.features=200)
    
    # ── FIX 3: assign disease from GEO metadata ───────────────────────────────
    obj$sample  <- sname
    obj$gsm_id  <- gsm_id
    obj$disease <- if (gsm_id %in% names(gsm_disease))
      gsm_disease[[gsm_id]]
    else
      ifelse(grepl("HC", sname, ignore.case=TRUE),
             "Healthy_Control", "Parkinsons_Disease")
    
    seurat_raw[[sname]] <- obj
    cat(" ✔", ncol(obj), "cells |", obj$disease[1], "\n")
  }, error=function(e) cat(" ✘ FAILED:", e$message, "\n"))
}
cat("  Total cells:", sum(sapply(seurat_raw, ncol)), "\n")
cat("  Disease breakdown:\n")
print(table(sapply(seurat_raw, function(x) x$disease[1])))

# Mitochondrial %
seurat_raw <- lapply(seurat_raw, function(seu) {
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern="^MT-"); seu
})

# QC plots BEFORE filtering
pdf("results/qc/01_combined_violin_before.pdf", width=14, height=5)
all_raw <- merge(seurat_raw[[1]], y=seurat_raw[-1],
                 add.cell.ids=names(seurat_raw))
print(VlnPlot(all_raw,
              features=c("nFeature_RNA","nCount_RNA","percent.mt"),
              ncol=3, pt.size=0, group.by="orig.ident") +
        plot_annotation(title="QC — All Samples Before Filtering",
                        theme=theme(plot.title=element_text(face="bold"))))
dev.off()

pdf("results/qc/02_scatter_qc.pdf", width=10, height=5)
for (s in names(seurat_raw)) {
  tryCatch({
    p1 <- FeatureScatter(seurat_raw[[s]],"nCount_RNA","nFeature_RNA") + ggtitle(s)
    p2 <- FeatureScatter(seurat_raw[[s]],"nCount_RNA","percent.mt")   + ggtitle(s)
    print(p1 + p2)
  }, error=function(e) NULL)
}
dev.off()

# Filter
cat(sprintf("\n  Thresholds: genes %d–%d | MT <%.0f%% | counts >%d\n",
            MIN_FEATURES, MAX_FEATURES, MAX_MT, MIN_COUNTS))
seurat_filtered <- lapply(seurat_raw, function(seu)
  subset(seu, subset = nFeature_RNA > MIN_FEATURES &
           nFeature_RNA < MAX_FEATURES &
           percent.mt   < MAX_MT       &
           nCount_RNA   > MIN_COUNTS))

cat("  After QC:\n")
qc_df <- data.frame(
  sample       = names(seurat_filtered),
  disease      = sapply(seurat_filtered, function(x) x$disease[1]),
  cells_before = sapply(seurat_raw,      ncol),
  cells_after  = sapply(seurat_filtered, ncol),
  median_genes = sapply(seurat_filtered, function(x) round(median(x$nFeature_RNA))),
  median_mt    = sapply(seurat_filtered, function(x) round(median(x$percent.mt),1))
)
print(qc_df)
write.csv(qc_df,"results/qc/qc_summary.csv", row.names=FALSE)

# QC plots AFTER filtering
pdf("results/qc/03_combined_violin_after.pdf", width=14, height=5)
all_filt <- merge(seurat_filtered[[1]], y=seurat_filtered[-1],
                  add.cell.ids=names(seurat_filtered))
print(VlnPlot(all_filt,
              features=c("nFeature_RNA","nCount_RNA","percent.mt"),
              ncol=3, pt.size=0, group.by="orig.ident") +
        plot_annotation(title="QC — All Samples After Filtering",
                        theme=theme(plot.title=element_text(face="bold"))))
dev.off()

# ================================================================================
# ---> ADD THIS BLOCK: UMAP Before and After QC Filtering
# ================================================================================
cat("\n━━━ Generating Combined Before & After QC UMAPs ━━━\n")
library(patchwork) # For side-by-side plotting

# 1. Process the 'Before QC' merged object
cat("  Processing Raw Data (This may take a moment due to high cell counts)...\n")
all_raw <- NormalizeData(all_raw, verbose = FALSE) %>%
  FindVariableFeatures(verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:20, verbose = FALSE)

# Generate Before Plot
p_umap_before <- DimPlot(all_raw, reduction = "umap", group.by = "orig.ident", pt.size = 0.1) +
  ggtitle("Before QC (Includes Low-Quality Cells)") +
  theme_classic() +
  theme(legend.position = "bottom")

# 2. Process the 'After QC' merged object
cat("  Processing Filtered Data...\n")
all_filt <- NormalizeData(all_filt, verbose = FALSE) %>%
  FindVariableFeatures(verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:20, verbose = FALSE)

# Generate After Plot
p_umap_after <- DimPlot(all_filt, reduction = "umap", group.by = "orig.ident", pt.size = 0.1) +
  ggtitle("After QC (Cleaned Cells Only)") +
  theme_classic() +
  theme(legend.position = "bottom")

# 3. Save side-by-side comparison to a new PDF
pdf("results/qc/04_qc_umap_before_after.pdf", width = 12, height = 6)
print(p_umap_before + p_umap_after)
dev.off()

cat("    ✔ Before/After QC UMAPs saved to results/qc/04_qc_umap_before_after.pdf\n")
# ================================================================================
#  STEP 3: DOUBLET DETECTION (DoubletFinder)
#  Doublets = two cells in one droplet — creates fake hybrid clusters
# ================================================================================
cat("━━━ STEP 3: Doublet Detection ━━━\n")

seurat_clean <- list()
for (s in names(seurat_filtered)) {
  cat("  Sample:", s, "\n")
  seu <- seurat_filtered[[s]]
  
  result <- tryCatch({
    # ── FIX 7: no NormalizeData before — goes straight to standard workflow ────
    seu <- NormalizeData(seu, verbose=FALSE) %>%
      FindVariableFeatures(verbose=FALSE) %>%
      ScaleData(verbose=FALSE) %>%
      RunPCA(verbose=FALSE) %>%
      RunUMAP(dims=1:20, verbose=FALSE) %>%
      FindNeighbors(dims=1:20, verbose=FALSE) %>%
      FindClusters(resolution=0.5, verbose=FALSE)
    
    dr    <- min(0.25, ncol(seu)/100000)
    sweep <- paramSweep(seu, PCs=1:20, sct=FALSE)
    stats <- summarizeSweep(sweep, GT=FALSE)
    bcmvn <- find.pK(stats)
    
    # Safe pK extraction — avoids "cannot xtfrm data frames"
    pK   <- as.numeric(as.character(
      bcmvn[["pK"]][which.max(bcmvn[["BCmetric"]])] ))
    nExp <- round(dr * ncol(seu) * (1 - modelHomotypic(seu$seurat_clusters)))
    seu  <- doubletFinder(seu, PCs=1:20, pN=0.25, pK=pK,
                          nExp=nExp, reuse.pANN=NULL, sct=FALSE)
    df_col <- grep("DF.class", colnames(seu@meta.data), value=TRUE)[1]
    seu$doublet_status <- seu@meta.data[[df_col]]
    cat(sprintf("    ✔ Doublets: %d | Singlets: %d\n", 
                sum(seu$doublet_status=="Doublet"),
                sum(seu$doublet_status=="Singlet")))
    seu
  }, error=function(e) {
    cat("  ⚠ DoubletFinder failed (", e$message,") — keeping all cells\n")
    seu$doublet_status <- "Singlet"; seu
  })
  
  # NEW FIX: Draw the plot BEFORE throwing away the doublets
  p <- DimPlot(result, group.by="doublet_status", 
               cols=c(Doublet="red", Singlet="grey80"), pt.size=0.5) +
    ggtitle(paste("Doublets —", s)) + theme_classic()
  
  # Save the plot temporarily
  if (!exists("doublet_plots")) doublet_plots <- list()
  doublet_plots[[s]] <- p
  
  if (!exists("seurat_with_doublets")) seurat_with_doublets <- list()
  seurat_with_doublets[[s]] <- result
  
  # NOW throw away the doublets
  before <- ncol(result)
  result <- subset(result, subset=doublet_status=="Singlet")
  cat(sprintf("    Cells: %d → %d\n", before, ncol(result)))
  seurat_clean[[s]] <- result
}

# Print the saved plots that still contain the red doublets
pdf("results/doublet/01_doublet_umap.pdf", width=8, height=6)
for (p in doublet_plots) print(p)
dev.off()
cat("━━━ Generating Combined Before & After UMAPs ━━━\n")
library(patchwork) # Ensure patchwork is loaded for side-by-side plotting

# 1. Merge the 'Before' objects (contains both Doublets and Singlets)
combined_before <- merge(x = seurat_with_doublets[[1]],
                         y = seurat_with_doublets[-1],
                         add.cell.ids = names(seurat_with_doublets))

# 2. Process the merged object to generate a single, shared UMAP space
combined_before <- NormalizeData(combined_before, verbose=FALSE) %>%
  FindVariableFeatures(verbose=FALSE) %>%
  ScaleData(verbose=FALSE) %>%
  RunPCA(verbose=FALSE) %>%
  RunUMAP(dims=1:20, verbose=FALSE)

# 3. Create the 'Before' Plot across all samples
p_combined_before <- DimPlot(combined_before, group.by="doublet_status", 
                             cols=c(Doublet="red", Singlet="grey80"), pt.size=0.1) +
  ggtitle("Whole Dataset: Before Removal") +
  theme_classic()

# 4. Create the 'After' Plot 
# We subset the merged object here instead of merging 'seurat_clean'. 
# This guarantees the UMAP coordinates stay exactly the same, 
# so you just see the red dots disappear!
combined_after <- subset(combined_before, subset = doublet_status == "Singlet")

p_combined_after <- DimPlot(combined_after, group.by="doublet_status", 
                            cols=c(Singlet="grey80"), pt.size=0.1) +
  ggtitle("Whole Dataset: After Removal") +
  theme_classic()

# 5. Save the combined plots side-by-side
pdf("results/doublet/02_combined_before_after_umap.pdf", width=12, height=6)
print(p_combined_before + p_combined_after)
dev.off()

cat(" ✔ Combined Before/After UMAPs saved.\n")
# ================================================================================
#  STEP 4: CELL CYCLE SCORING
#  Cell cycle phases can drive clustering — we score and regress them out
# ================================================================================
cat("━━━ STEP 4: Cell Cycle Scoring ━━━\n")

s_genes   <- cc.genes.updated.2019$s.genes
g2m_genes <- cc.genes.updated.2019$g2m.genes

seurat_clean <- lapply(seurat_clean, function(seu) {
  seu <- NormalizeData(seu, verbose=FALSE)
  tryCatch({
    seu <- CellCycleScoring(seu, s.features=s_genes,
                            g2m.features=g2m_genes, set.ident=FALSE)
    seu
  }, error=function(e) {
    seu$S.Score <- 0; seu$G2M.Score <- 0; seu$Phase <- "Unknown"; seu
  })
})

pdf("results/cellcycle/01_cell_cycle_phases.pdf", width=12, height=5)
for (s in names(seurat_clean)) {
  tryCatch({
    df <- as.data.frame(table(seurat_clean[[s]]$Phase))
    print(ggplot(df, aes(x=Var1, y=Freq, fill=Var1)) +
            geom_bar(stat="identity") +
            scale_fill_manual(values=c(G1="steelblue",G2M="firebrick",
                                       S="forestgreen",Unknown="grey")) +
            ggtitle(paste("Cell Cycle —",s)) +
            xlab("Phase") + ylab("Cells") + theme_classic() + NoLegend())
  }, error=function(e) NULL)
}
dev.off()
cat("✔ Cell cycle scoring done.\n\n")

# ================================================================================
#  STEP 5: NORMALIZATION (Log-norm per sample)
#  FIX 4: Using log-normalization + Harmony only
#  NOT SCTransform + Harmony (which over-stacks corrections)
#  SCTransform is good but should not be combined with Harmony
#  Log-norm + Harmony is the standard validated approach
# ================================================================================
cat("━━━ STEP 5: Normalization (Log-norm per sample) ━━━\n")
cat("  Using: Log-normalization + Harmony\n")
cat("  Reason: SCT + Harmony stacks corrections — risks overcorrection\n\n")

seurat_norm <- lapply(names(seurat_clean), function(s) {
  cat("  Normalizing:", s, "...")
  seu <- seurat_clean[[s]]
  tryCatch({
    # Normalize and find variable features ONLY (Scaling moved to Step 6)
    seu <- NormalizeData(seu, verbose=FALSE)
    seu <- FindVariableFeatures(seu, nfeatures=3000, verbose=FALSE)
    cat(" ✔\n")
    seu
  }, error=function(e) {
    cat(" ⚠ Normalization failed\n")
    seu
  })
})
names(seurat_norm) <- names(seurat_clean)
cat("✔ Normalization done.\n\n")

# ================================================================================
#  STEP 6: MERGE + HARMONY INTEGRATION
#  Corrects batch effects between the 12 samples
# ================================================================================
cat("━━━ STEP 6: Merge + Harmony Integration ━━━\n")

# Find shared variable genes
features <- SelectIntegrationFeatures(object.list=seurat_norm, nfeatures=3000)

# Merge
merged <- merge(seurat_norm[[1]], y=seurat_norm[-1],
                add.cell.ids=names(seurat_norm))
VariableFeatures(merged) <- features
cat("  Merged:", ncol(merged), "cells from", length(seurat_norm), "samples\n")

# Confirm disease column survived merge
if (!"disease" %in% colnames(merged@meta.data) ||
    all(is.na(merged$disease))) {
  cat("  Re-attaching disease labels from metadata...\n")
  # Map via orig.ident → sname → disease
  disease_by_sname <- sapply(seurat_clean, function(x) x$disease[1])
  merged$disease   <- disease_by_sname[merged$orig.ident]
}
cat("  Disease groups:\n"); print(table(merged$disease))

# NEW FIX: Scale the merged data AND apply regression here, right before PCA!
cat("  Scaling merged data and regressing noise...\n")
merged <- ScaleData(merged, 
                    vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"), 
                    verbose=FALSE)

# PCA
merged <- RunPCA(merged, npcs=50, verbose=FALSE)
# UMAP before integration
m_before <- RunUMAP(merged, dims=1:20, verbose=FALSE) %>%
  FindNeighbors(dims=1:20, verbose=FALSE) %>%
  FindClusters(resolution=0.5, verbose=FALSE)

pdf("results/integration/01_umap_before_integration.pdf", width=14, height=6)
p1 <- DimPlot(m_before, group.by="orig.ident", pt.size=0.3) +
  ggtitle("BEFORE Integration — by Sample") + theme_classic()
p2 <- DimPlot(m_before, group.by="disease", pt.size=0.3,
              cols=c(Healthy_Control="steelblue",
                     Parkinsons_Disease="firebrick")) +
  ggtitle("BEFORE Integration — HC vs PD") + theme_classic()
print(p1 + p2)
dev.off()
cat("  ✔ Saved: results/integration/01_umap_before_integration.pdf\n")

# Harmony
cat("  Running Harmony...\n")
merged <- RunHarmony(merged, group.by.vars="orig.ident",
                     plot_convergence=FALSE, verbose=FALSE)

# Build UMAP from Harmony embeddings
merged <- RunUMAP(merged, reduction="harmony", dims=1:20,
                  verbose=FALSE, reduction.name="umap") %>%
  FindNeighbors(reduction="harmony", dims=1:20, verbose=FALSE) %>%
  FindClusters(resolution=0.5, verbose=FALSE)

pdf("results/integration/02_umap_after_integration.pdf", width=14, height=6)
p1 <- DimPlot(merged, group.by="orig.ident", pt.size=0.3) +
  ggtitle("AFTER Integration — by Sample") + theme_classic()
p2 <- DimPlot(merged, group.by="disease", pt.size=0.3,
              cols=c(Healthy_Control="steelblue",
                     Parkinsons_Disease="firebrick")) +
  ggtitle("AFTER Integration — HC vs PD") + theme_classic()
print(p1 + p2)
dev.off()
cat("  ✔ Saved: results/integration/02_umap_after_integration.pdf\n")

# ── FIX 5: LISI on Harmony embeddings (correct space) ────────────────────────
cat("  Computing LISI on Harmony embeddings (correct space)...\n")
tryCatch({
  harm_emb  <- Embeddings(merged, "harmony")   # ← FIXED: was "umap" before
  meta_lisi <- data.frame(sample  = merged$orig.ident,
                          disease = merged$disease,
                          row.names=colnames(merged))
  lisi_res  <- compute_lisi(harm_emb, meta_lisi, c("sample","disease"))
  lisi_out  <- data.frame(
    metric     = c("Sample LISI  (>1.5 = good mixing)",
                   "Disease LISI (~1.0 = biology preserved)"),
    mean_score = c(round(mean(lisi_res$sample,  na.rm=TRUE),3),
                   round(mean(lisi_res$disease, na.rm=TRUE),3)))
  write.csv(lisi_out,"results/integration/lisi_scores.csv", row.names=FALSE)
  cat("  LISI scores:\n"); print(lisi_out)
  cat("  ✔ Saved: results/integration/lisi_scores.csv\n")
}, error=function(e) cat("  ⚠ LISI skipped:", e$message, "\n"))

write.csv(as.data.frame(table(merged$orig.ident)),
          "results/integration/cells_per_sample.csv", row.names=FALSE)
saveRDS(merged,"objects/03_seurat_integrated.rds")
cat(sprintf("✔ Integration done. Total cells: %d\n\n", ncol(merged)))

# ================================================================================
#  STEP 7: ADAPTIVE PC SELECTION + FINAL UMAP
# ================================================================================
cat("━━━ STEP 7: Adaptive PC Selection ━━━\n")

seu       <- merged
pca_std   <- seu@reductions[["pca"]]@stdev

# NEW FIX: Calculate the drop in variance between consecutive PCs (the "True Elbow")
pct_var   <- (pca_std^2 / sum(pca_std^2)) * 100
cumvar    <- cumsum(pct_var)

# Find where the drop between consecutive PCs flattens out (change < 0.1%)
diff_var  <- abs(diff(pct_var))
elbow_pc  <- which(diff_var < 0.1)[1]

# Safety fallback: if it doesn't flatten, default to 20 PCs
if (is.na(elbow_pc)) elbow_pc <- 20 

# Ensure we pick between 15 and 40 PCs to preserve complex biology
N_PCS     <- min(40, max(15, elbow_pc))

cat(sprintf("  Elbow detected at PC %d\n", elbow_pc))
cat(sprintf("  ✔ Using %d PCs (adaptive — consecutive drop < 0.1%%)\n\n", N_PCS))

# Elbow + cumulative variance plot (Updated for True Elbow Method)
elbow_df <- data.frame(PC=1:length(pca_std), StdDev=pca_std, CumVar=cumvar)
pdf("results/processing/01_elbow_plot.pdf", width=12, height=5)

p1 <- ggplot(elbow_df, aes(x=PC, y=StdDev)) +
  geom_line(color="steelblue") + geom_point(size=2, color="steelblue") +
  geom_vline(xintercept=N_PCS, linetype="dashed", color="red") +
  annotate("text", x=N_PCS+1.5, y=max(pca_std)*0.85,
           label=paste0("Chosen PC: ", N_PCS), color="red", size=4, hjust=0) +
  ggtitle("Elbow Plot (Drop < 0.1%)") + xlab("PC") + ylab("Std Dev") +
  theme_classic(base_size=12) + theme(plot.title=element_text(face="bold"))

p2 <- ggplot(elbow_df, aes(x=PC, y=CumVar)) +
  geom_line(color="darkgreen", linewidth=1) +
  geom_point(size=1.5, color="darkgreen") +
  geom_vline(xintercept=N_PCS, linetype="dashed", color="red") +
  annotate("text", x=N_PCS+1.5, y=min(cumvar) + 10,
           label=sprintf("Cum. Var: %.1f%%", cumvar[N_PCS]), color="red", size=4, hjust=0) +
  ggtitle("Cumulative Variance Explained") +
  xlab("PC") + ylab("Cumulative %") +
  theme_classic(base_size=12) + theme(plot.title=element_text(face="bold"))

print(p1|p2)
dev.off()
cat("  ✔ Saved: results/processing/01_elbow_plot.pdf\n")

# Recompute with adaptive PCs
seu <- RunUMAP(seu, reduction="harmony", dims=1:N_PCS,
               verbose=FALSE, reduction.name="umap") %>%
  FindNeighbors(reduction="harmony", dims=1:N_PCS, verbose=FALSE)

# QC overlay
pdf("results/processing/02_umap_qc_overlay.pdf", width=14, height=12)
p1 <- FeaturePlot(seu,"nFeature_RNA",pt.size=0.3) +
  scale_color_viridis_c() + ggtitle("Genes per Cell") + theme_classic()
p2 <- FeaturePlot(seu,"nCount_RNA",  pt.size=0.3) +
  scale_color_viridis_c() + ggtitle("UMI Counts")    + theme_classic()
p3 <- FeaturePlot(seu,"percent.mt",  pt.size=0.3) +
  scale_color_viridis_c(option="inferno") + ggtitle("% MT") + theme_classic()
p4 <- DimPlot(seu, group.by="disease", pt.size=0.3,
              cols=c(Healthy_Control="steelblue",
                     Parkinsons_Disease="firebrick")) +
  ggtitle("HC vs PD") + theme_classic()
print((p1|p2)/(p3|p4))
dev.off()
cat("  ✔ Saved: results/processing/02_umap_qc_overlay.pdf\n")

# Cell cycle on UMAP
if ("Phase" %in% colnames(seu@meta.data)) {
  pdf("results/cellcycle/02_cell_cycle_umap.pdf", width=8, height=6)
  print(DimPlot(seu, group.by="Phase", pt.size=0.3,
                cols=c(G1="steelblue",G2M="firebrick",
                       S="forestgreen",Unknown="grey")) +
          ggtitle("Cell Cycle Phase on UMAP\n(should not cluster separately)") +
          theme_classic())
  dev.off()
}

write.csv(seu@meta.data,"results/processing/cell_metadata.csv")
saveRDS(seu,"objects/04_seurat_processed.rds")
cat(sprintf("✔ Processing done. Cells: %d | PCs: %d\n\n", ncol(seu), N_PCS))










#  RAM CLEANING =======
# Create the objects directory if it doesn't exist
if (!dir.exists("objects")) dir.create("objects")

# Save the integrated object (this is your master file)
saveRDS(merged, "objects/07_integrated_checkpoint.rds")

# Confirm it saved
cat("✔ Step 7 Save Point created: objects/07_integrated_checkpoint.rds\n")

#clean slate
# Clear every single object from your environment
rm(list = ls())

# Force the server to reclaim the RAM
gc()











#fresh start
# ==============================================================================
#  MONDAY MORNING RELOAD: RESTORE FROM STEP 7 CHECKPOINT
# ==============================================================================
# 1. Load Essential Libraries
library(Seurat)
library(dplyr)
library(ggplot2)


# 3. Load the Integrated Data from the Hard Drive (The "Save Point")
if (file.exists("objects/07_integrated_checkpoint.rds")) {
  cat("📂 Loading integrated object... please wait.\n")
  seu <- readRDS("objects/07_integrated_checkpoint.rds")
  cat("✔ Data Loaded successfully.\n")
  
  # 4. Join Layers for Seurat v5 (Crucial for Clustering & Markers)
  if (as.numeric(substr(packageVersion("Seurat"), 1, 1)) >= 5) {
    cat("⚙ Seurat v5 detected: Joining layers for Step 8...\n")
    seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
  }
  
  # 5. Final Status Check
  print(seu)
  cat("\n🚀 READY FOR STEP 8: CLUSTERING\n")
  
} else {
  stop("❌ Error: Checkpoint file not found in 'objects/' folder. Check your file path!")
}
# ==============================================================================
# === UPDATED STARTUP ===
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork) # Added this so wrap_plots works!

CHOSEN_RES <- 0.4
N_CORES    <- 10
N_PCS      <- 15  # Added this so your titles don't error out!

seu <- readRDS("objects/07_integrated_checkpoint.rds")
# ========================









# ================================================================================
#  STEP 8: CLUSTERING
# ================================================================================
cat("━━━ STEP 8: Clustering ━━━\n")

# NEW FIX: Automatically inject CHOSEN_RES so it is guaranteed to be calculated
resolutions <- unique(sort(c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, CHOSEN_RES)))
seu         <- FindClusters(seu, resolution=resolutions, verbose=FALSE)

# Auto-detect prefix
res_prefix <- ifelse(any(grepl("SCT_snn", colnames(seu@meta.data))),
                     "SCT_snn_res.", "RNA_snn_res.")

cat("  Clusters per resolution:\n")
for (r in resolutions) {
  col <- paste0(res_prefix, r)
  if (col %in% colnames(seu@meta.data))
    cat(sprintf("    %.1f → %d clusters\n", r,
                length(unique(seu@meta.data[[col]]))))
}

cat(sprintf("\n  ⚠ CHOSEN_RES = %.1f (look at 01_umap_all_resolutions.pdf\n", CHOSEN_RES))
cat("    and change CHOSEN_RES at the top of the script if needed)\n\n")

# All resolutions plot
pdf("results/clustering/01_umap_all_resolutions.pdf", width=16, height=10)
pl <- lapply(resolutions, function(r) {
  col <- paste0(res_prefix, r)
  if (!col %in% colnames(seu@meta.data)) return(NULL)
  nc  <- length(unique(seu@meta.data[[col]]))
  DimPlot(seu, group.by=col, label=TRUE, pt.size=0.3, label.size=3) +
    ggtitle(sprintf("Res %.1f — %d clusters",r,nc)) +
    theme_classic() + NoLegend()
})
print(wrap_plots(pl[!sapply(pl,is.null)], ncol=3))
dev.off()
cat("  ✔ Saved: results/clustering/01_umap_all_resolutions.pdf\n")

# Apply chosen resolution
chosen_col  <- paste0(res_prefix, CHOSEN_RES)
Idents(seu) <- chosen_col
seu$seurat_clusters <- Idents(seu)
n_clust     <- length(unique(seu$seurat_clusters))
cat(sprintf("  Using resolution %.1f → %d clusters\n", CHOSEN_RES, n_clust))

# Numbered UMAP
pdf("results/clustering/02_umap_clusters_numbered.pdf", width=10, height=8)
print(DimPlot(seu, label=TRUE, pt.size=0.4, label.size=5,
              label.color="black") +
        ggtitle(sprintf("Clusters (Res %.1f | %d PCs)", CHOSEN_RES, N_PCS)) +
        theme_classic(base_size=13) +
        theme(plot.title=element_text(face="bold")))
dev.off()
cat("  ✔ Saved: results/clustering/02_umap_clusters_numbered.pdf\n")

# Marker genes
# Find Marker genes
cat("  Finding marker genes...\n")

tryCatch({
  seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
}, error=function(e) NULL)

markers <- FindAllMarkers(seu, only.pos=TRUE, min.pct=0.25,
                          logfc.threshold=0.25, test.use="wilcox",
                          verbose=FALSE)

# Save all markers to a CSV
write.csv(markers, "results/clustering/all_marker_genes.csv", row.names=FALSE)

# Extract Top 10 Markers per cluster
top10 <- markers %>% 
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10, with_ties = FALSE) %>% 
  ungroup()

# Save the top 10 markers to a CSV
write.csv(top10, "results/clustering/top10_markers_per_cluster.csv", row.names = FALSE)

# Extract Top 5 and Top 2 genes for plotting later
top5_genes <- markers %>% 
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5, with_ties = FALSE) %>% 
  pull(gene) %>% 
  unique()

top2_genes <- markers %>% 
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 2, with_ties = FALSE) %>% 
  pull(gene) %>% 
  unique()

# NEW FIX: Explicitly scale the marker genes so they actually appear in the heatmap
cat("  Scaling marker genes for heatmap...\n")
seu <- ScaleData(seu, 
                 features = unique(c(VariableFeatures(seu), top5_genes)), 
                 verbose=FALSE)

pdf("results/clustering/03_marker_heatmap.pdf", width=14, height=10)
# Note: Added print() here to ensure the ggplot renders properly inside the script
print(DoHeatmap(seu, features=top5_genes) +
        theme(axis.text.y=element_text(size=7)) +
        ggtitle("Top 5 Marker Genes per Cluster"))
dev.off()
cat("  ✔ Saved: results/clustering/03_marker_heatmap.pdf\n")

pdf("results/clustering/04_marker_dotplot.pdf", width=14, height=6)
print(DotPlot(seu, features=top2_genes) + RotatedAxis() +
        ggtitle("Top 2 Markers per Cluster") + theme_classic() +
        theme(axis.text.x=element_text(size=7,angle=45)))
dev.off()
cat("  ✔ Saved: results/clustering/04_marker_dotplot.pdf\n")

write.csv(as.data.frame(table(seu$seurat_clusters)),
          "results/clustering/cells_per_cluster.csv", row.names=FALSE)
saveRDS(seu,"objects/05_seurat_clustered.rds")
cat(sprintf("✔ Clustering done. %d clusters.\n\n", n_clust))


# ==============================================================================
#  STEP 9: AUTOMATED ANNOTATION (SingleR) - MASTER FIX
# ==============================================================================
cat("━━━ STEP 9: Automated Annotation (SingleR) ━━━\n")

library(SingleR)
library(celldex)
library(dplyr)

# 1. Setup Reference Data (Full-Body Reference)
cat("  Loading Human Primary Cell Atlas (Full Body Reference)...\n")
ref <- HumanPrimaryCellAtlasData()

# 2. Extract Data (Using 'layer' for Seurat v5 and Variable Features for speed)
cat("  Extracting Variable Features...\n")
var_genes <- VariableFeatures(seu)
expr <- GetAssayData(seu, assay="RNA", layer="data")[var_genes, ]

# 3. Run SingleR 
cat("  Running SingleR (using 10 cores)...\n")
singler_res <- SingleR(test=expr, ref=ref, labels=ref$label.main, 
                       BPPARAM=BiocParallel::MulticoreParam(workers=N_CORES))

# 4. Extract raw labels and attach cell names
raw_labels        <- singler_res$labels
names(raw_labels) <- rownames(singler_res)

# 5. Add raw labels to Seurat safely
seu$SingleR_Raw <- raw_labels[colnames(seu)]

# 6. Calculate the Majority Vote per Cluster (Bulletproof version)
# Forcing character vectors ensures the dataframe builds perfectly without missing columns
anno_df <- data.frame(
  Cluster = as.character(Idents(seu)),
  Label = as.character(seu$SingleR_Raw),
  stringsAsFactors = FALSE
)

majority_labels <- anno_df %>%
  dplyr::group_by(Cluster, Label) %>%
  dplyr::tally() %>%
  dplyr::slice_max(order_by = n, n = 1, with_ties = FALSE)

# 7. Map the clean labels back to every cell
seu$SingleR_Majority <- majority_labels$Label[match(as.character(Idents(seu)), majority_labels$Cluster)]

# 8. Save a Cluster-to-CellType Mapping Table
mapping_table <- majority_labels %>% 
  dplyr::rename(Assigned_CellType = Label, Cell_Count = n)
write.csv(mapping_table, "results/annotation/cluster_annotation_summary.csv", row.names = FALSE)
cat("  ✔ Saved: results/annotation/cluster_annotation_summary.csv\n")

# 9. Save the Clean UMAP
pdf("results/annotation/02_umap_singler_majority.pdf", width=14, height=10)
print(DimPlot(seu, group.by="SingleR_Majority", label=TRUE, repel=TRUE, pt.size=0.3) +
        ggtitle("SingleR Annotations (Full-Body Majority Vote)") +
        theme_classic() + NoLegend())
dev.off()
cat("  ✔ Saved: results/annotation/02_umap_singler_majority.pdf\n")

# 10. Known Markers DotPlot (Verification Step)
cat("  Generating DotPlot of known markers...\n")
known_markers <- c("TH", "SLC6A3", "DDC", "ALDH1A1",   # Dopaminergic
                   "SLC17A6", "SLC17A7",               # Glutamatergic
                   "GAD1", "GAD2",                     # GABAergic
                   "PTPRC", "CD14", "AIF1",            # Immune/Microglia
                   "GFAP", "S100B")                    # Astrocytes

# Ensure we only plot genes that actually exist to prevent errors
valid_markers <- intersect(known_markers, rownames(seu))

pdf("results/annotation/03_dotplot_known_markers.pdf", width=12, height=8)
print(DotPlot(seu, features=valid_markers, cols=c("lightgrey", "red")) + 
        RotatedAxis() + 
        ggtitle("Expression of Known Brain Markers per Cluster"))
dev.off()
cat("  ✔ Saved: results/annotation/03_dotplot_known_markers.pdf\n")
# ==============================================================================
# 1. Extend timeout again just in case
options(timeout = 600)

# 2. Re-initialize BiocManager
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# 3. Try to install the dependencies individually
# If it asks "Update all/some/none?", type 'n' for None.
# If it asks "Do you want to install from source?", type 'n' for No.
BiocManager::install("DirichletMultinomial", update = FALSE, type = "binary")
BiocManager::install("TFBSTools", update = FALSE, type = "binary")
# ==============================================================================
# STEP 9.2: AZIMUTH ANNOTATION (Motor Cortex Reference)
# ==============================================================================
cat("━━━ Running Azimuth Annotation ━━━\n")

# 1. Run Azimuth (Downloads reference automatically if not present)
suppressMessages(library(Azimuth))
suppressMessages(library(SeuratData))
seu <- RunAzimuth(seu, reference = "humanmotorcortex")

# 2. Calculate the Majority Vote per Cluster (Matches your SingleR logic)
azimuth_anno <- data.frame(
  Cluster = as.character(Idents(seu)),
  Label = as.character(seu$predicted.id),
  stringsAsFactors = FALSE
)

azimuth_majority <- azimuth_anno %>%
  dplyr::group_by(Cluster, Label) %>%
  dplyr::tally() %>%
  dplyr::slice_max(order_by = n, n = 1, with_ties = FALSE)

# 3. Map the clean labels back to every cell
seu$Azimuth_Majority <- azimuth_majority$Label[match(as.character(Idents(seu)), azimuth_majority$Cluster)]

# 4. Save the Clean UMAP
pdf("results/annotation/04_umap_azimuth_majority.pdf", width=14, height=10)
print(DimPlot(seu, group.by="Azimuth_Majority", label=TRUE, repel=TRUE, pt.size=0.3) +
        ggtitle("Azimuth Annotations (Human Motor Cortex Majority Vote)") +
        theme_classic() + NoLegend())
dev.off()

cat("  ✔ Saved: results/annotation/04_umap_azimuth_majority.pdf\n")
#manually clustering for the parkinsins dataset 202210 by me




# ==============================================================================
# VISUALIZING PARKINSON'S DISEASE BIOMARKERS
# ==============================================================================
cat("━━━ Generating PD Biomarker Plots ━━━\n")

# 1. Define the core PD biomarkers 
pd_markers <- c(
  "TH", "SLC6A3", "KCNJ6",  # Dopaminergic Neurons (Substantia Nigra A9)
  "AIF1", "HLA-DRA",        # Reactive Microglia / Inflammation
  "GFAP", "S100B",          # Astrocytes
  "SNCA", "LRRK2"           # Core PD Genetic Risk Factors
)

# 2. Safety check: Only plot genes that survived QC and exist in your Seurat object
# This prevents the script from crashing if a gene was filtered out earlier.
valid_pd_markers <- intersect(pd_markers, rownames(seu))
cat("  Found", length(valid_pd_markers), "out of", length(pd_markers), "PD markers in the dataset.\n")

# 3. Create a grid of UMAPs, one for each gene
# Note: 'order = TRUE' is crucial here. It forces Seurat to draw the red (expressing) 
# cells on top of the grey cells, preventing them from being hidden in dense clusters.
pdf("results/annotation/05_PD_biomarker_featureplots.pdf", width = 16, height = 12)
p_feat <- FeaturePlot(seu, 
                      features = valid_pd_markers, 
                      cols = c("grey90", "red3"), 
                      ncol = 3, 
                      pt.size = 0.5, 
                      order = TRUE, 
                      label = TRUE, 
                      repel = TRUE)
print(p_feat)
dev.off()

# 4. Create Violin Plots for quantitative visualization across clusters
# Note: 'pt.size = 0' removes the individual black dots so you can clearly see the violin shape.
pdf("results/annotation/06_PD_biomarker_violins.pdf", width = 16, height = 12)
p_vln <- VlnPlot(seu, 
                 features = valid_pd_markers, 
                 pt.size = 0, 
                 ncol = 3)
print(p_vln)
dev.off()

cat("  ✔ Saved FeaturePlots: results/annotation/05_PD_biomarker_featureplots.pdf\n")
cat("  ✔ Saved VlnPlots: results/annotation/06_PD_biomarker_violins.pdf\n")





# ==============================================================================
#  STEP 9.5: MANUAL ANNOTATION (Presentation-Ready Master Version)
# ==============================================================================
cat("━━━ STEP 9.5: Finalizing Presentation UMAP ━━━\n")

# 1. Start fresh with cluster numbers
seu$Manual_CellType <- as.character(seu$seurat_clusters)

# 2. Assign the names we already verified
seu$Manual_CellType[seu$Manual_CellType %in% c("5", "20", "24")] <- "Microglia"
seu$Manual_CellType[seu$Manual_CellType %in% c("2", "23", "25")] <- "Astrocyte"
seu$Manual_CellType[seu$Manual_CellType %in% c("6", "8", "9", "11", "13", "14")] <- "GABAergic Neuron"
seu$Manual_CellType[seu$Manual_CellType %in% c("3", "7", "10", "15", "16")] <- "Glutamatergic Neuron"

# 3. NEW: Assign the remaining major brain cells to clean up the numbers
# 🟣 Oligodendrocytes (The massive 0, 1, 4 island)
seu$Manual_CellType[seu$Manual_CellType %in% c("0", "1", "4")] <- "Oligodendrocyte"

# 🟤 OPCs (Oligodendrocyte Precursor Cells)
seu$Manual_CellType[seu$Manual_CellType %in% c("12", "18", "19")] <- "OPC"

# 🟡 Endothelial Cells (Blood vessels - Cluster 21)
seu$Manual_CellType[seu$Manual_CellType %in% c("21")] <- "Endothelial"

# ⚪ Catch-all for any tiny remaining numbers (so the plot looks clean)
seu$Manual_CellType[seu$Manual_CellType %in% c("17", "22")] <- "Other/Mixed"

# 4. Set these new names as the Default
Idents(seu) <- "Manual_CellType"

# 5. Draw the FINAL, beautiful UMAP for your presentation
pdf("results/annotation/04_umap_manual_annotation_FINAL.pdf", width=12, height=8)
print(DimPlot(seu, group.by="Manual_CellType", label=TRUE, repel=TRUE, pt.size=0.3) +
        ggtitle("Final Cell Type Annotations") +
        theme_classic())
dev.off()
cat("  ✔ Saved: results/annotation/04_umap_manual_annotation_FINAL.pdf\n")

# 6. EXPANDED DotPlot to mathematically prove the new names
cat("  Generating Expanded DotPlot to verify new clusters...\n")
expanded_markers <- c("TH", "SLC6A3",                        # Dopaminergic
                      "SLC17A7", "SLC17A6",                  # Glutamatergic
                      "GAD1", "GAD2",                        # GABAergic
                      "PTPRC", "AIF1",                       # Microglia
                      "GFAP", "S100B",                       # Astrocytes
                      "MBP", "PLP1", "MOG",                  # NEW: Oligodendrocytes
                      "PDGFRA", "CSPG4",                     # NEW: OPCs
                      "CLDN5", "VWF")                        # NEW: Endothelial

valid_expanded <- intersect(expanded_markers, rownames(seu))

pdf("results/annotation/05_dotplot_expanded_proof.pdf", width=14, height=8)
print(DotPlot(seu, features=valid_expanded, cols=c("lightgrey", "blue")) + 
        RotatedAxis() + 
        ggtitle("Expanded Brain Markers (Proof of Annotation)"))
dev.off()
cat("  ✔ Saved: results/annotation/05_dotplot_expanded_proof.pdf\n")
# ==============================================================================

# free space to ram

cat("Saving fully annotated Seurat object...\n")
saveRDS(seu, "objects/09_seurat_annotated_FINAL.rds")
cat("✔ Safe to quit!\n")

q(save = "no")


#morning reload 

# 1. Load the core libraries
library(Seurat)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(dplyr)

# 2. Reload your perfectly annotated data
cat("Loading Seurat object...\n")
seu <- readRDS("objects/09_seurat_annotated_FINAL.rds")
cat("✔ Ready for Step 10!\n")

# === UPDATED STARTUP ===
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork) # Added this so wrap_plots works!

CHOSEN_RES <- 0.6
N_CORES    <- 10
N_PCS      <- 15  # Added this so your titles don't error out!

#all library
suppressPackageStartupMessages({
  library(Seurat);  library(ggplot2);  library(dplyr);    library(patchwork)
  library(ggrepel); library(Matrix);   library(GEOquery); library(R.utils)
  library(harmony); library(DoubletFinder)
  library(SingleR); library(celldex);  library(DESeq2)
  library(lisi);    library(monocle3); library(SeuratWrappers)
})




# ================================================================================
#  STEP 10: TRAJECTORY ANALYSIS (Monocle3) - OLIGODENDROCYTE LINEAGE
# ================================================================================
cat("━━━ STEP 10: Trajectory Analysis (OPC to Oligodendrocyte) ━━━\n")

tryCatch({
  # 1. SUBSET THE DATA: Only look at the relevant biological family
  cat("  Subsetting to OPC and Oligodendrocyte lineage...\n")
  seu_sub <- subset(seu, Manual_CellType %in% c("OPC", "Oligodendrocyte"))
  
  # 2. Convert subset to Monocle3
  cds <- as.cell_data_set(seu_sub)
  rowData(cds)$gene_short_name <- rownames(cds)
  
  # 3. Monocle3 Specific Preprocessing on the subset
  cds <- estimate_size_factors(cds)
  cds <- preprocess_cds(cds, num_dim = 30)
  cds <- reduce_dimension(cds, reduction_method="UMAP")
  cds <- cluster_cells(cds)
  cds <- learn_graph(cds, use_partition=FALSE, verbose=FALSE)
  
  # 4. Plot the Cell Types
  pdf("results/trajectory/01_trajectory_celltype.pdf", width=10, height=8)
  print(plot_cells(cds, color_cells_by="Manual_CellType",
                   label_groups_by_cluster=FALSE,
                   label_leaves=TRUE, label_branch_points=TRUE,
                   cell_size=0.7) +
          ggtitle("Trajectory — OPC to Oligodendrocyte") + theme_classic())
  dev.off()
  
  # 5. EXACT BIOLOGICAL ROOT: We know OPCs are the baby cells.
  cat("  Setting OPCs as the biological root...\n")
  root_method <- "Manual selection of OPC cluster"
  
  # Find the exact cell barcodes for the baby OPCs
  opc_cells <- rownames(subset(seu_sub@meta.data, Manual_CellType == "OPC"))
  
  # 6. Calculate Pseudotime (Bulletproof Method)
  # Pass the raw cells directly to Monocle3 and let it calculate the root internally
  cds <- order_cells(cds, root_cells=opc_cells)
  seu_sub$pseudotime <- pseudotime(cds)
  
  # 7. Plot Pseudotime
  pdf("results/trajectory/02_pseudotime.pdf", width=10, height=8)
  print(plot_cells(cds, color_cells_by="pseudotime",
                   label_cell_groups=FALSE, label_leaves=FALSE,
                   label_branch_points=FALSE, cell_size=0.7,
                   trajectory_graph_color="black",
                   trajectory_graph_segment_size=1) +
          scale_color_viridis_c(option="plasma") +
          ggtitle("Oligo Maturation (Pseudotime)") + theme_classic())
  dev.off()
  
  # 8. Find Genes that change over time (Restricted to Variable Features to save RAM!)
  cat("  Calculating pseudotime-associated genes (using top variable genes)...\n")
  
  # Get the variable genes and make sure they exist in the Monocle object
  var_genes <- VariableFeatures(seu)
  valid_var_genes <- intersect(var_genes, rownames(cds))
  
  # Run the heavy math ONLY on the variable genes
  pt_genes <- graph_test(cds[valid_var_genes, ], neighbor_graph="principal_graph", cores=2) %>%
    dplyr::arrange(dplyr::desc(morans_I)) %>% 
    dplyr::filter(q_value < 0.05)
  
  write.csv(pt_genes,"results/trajectory/pseudotime_genes.csv", row.names=FALSE)
  cat(sprintf("  Found %d pseudotime-associated genes.\n", nrow(pt_genes)))
  
  # 9. Plot the Top 6 Genes that drive maturation
  top6 <- head(rownames(pt_genes), 6)
  if (length(top6) > 0) {
    pdf("results/trajectory/03_pseudotime_top_genes.pdf", width=16, height=10)
    print(plot_cells(cds, genes=top6, show_trajectory_graph=FALSE,
                     label_cell_groups=FALSE, cell_size=0.5) + theme_classic())
    dev.off()
  }
  cat("✔ Trajectory done.\n\n")
  
}, error=function(e) {
  cat("  ⚠ Trajectory error:", e$message, "\n")
  cat("  Continuing...\n\n")
})

# Save the subset objects
saveRDS(seu_sub,"objects/08_seurat_oligo_lineage.rds")
if (exists("cds")) saveRDS(cds,"objects/08_monocle_cds.rds")
# ================================================================================





# ================================================================================
#  STEP 10.5: FINDING MATURATION GENES (Bypassing the Linux Error)
# ================================================================================
cat("━━━ STEP 10.5: Calculating Gene Correlations ━━━\n")

# 1. Load the successfully saved timeline object
seu_sub <- readRDS("objects/08_seurat_oligo_lineage.rds")

# 2. Extract the RNA data and the successful timeline values
cat("  Extracting expression matrix and pseudotime...\n")
var_genes <- VariableFeatures(seu)
expr_matrix <- GetAssayData(seu_sub, assay="RNA", layer="data")[var_genes, ]
pt_values <- seu_sub$pseudotime

# 3. Filter out any cells that Monocle couldn't connect (infinite/NA values)
valid_cells <- !is.na(pt_values) & is.finite(pt_values)
expr_matrix <- expr_matrix[, valid_cells]
pt_values <- pt_values[valid_cells]

# 4. Calculate Spearman Correlation (Does the gene track with time?)
cat("  Running correlation math (No Linux libraries required!)...\n")
# We convert to a dense matrix just for the calculation to ensure it works smoothly
expr_matrix_dense <- as.matrix(expr_matrix)
cor_scores <- apply(expr_matrix_dense, 1, function(x) cor(x, pt_values, method="spearman"))

# 5. Build the final clean table
gene_df <- data.frame(
  Gene = names(cor_scores),
  Spearman_Correlation = cor_scores,
  stringsAsFactors = FALSE
)

# Sort by absolute correlation (The strongest drivers of maturation rise to the top)
gene_df <- gene_df[order(abs(gene_df$Spearman_Correlation), decreasing = TRUE), ]

# Save your master gene list
write.csv(gene_df, "results/trajectory/pseudotime_genes_correlation.csv", row.names=FALSE)
cat("  ✔ Saved: results/trajectory/pseudotime_genes_correlation.csv\n")

# 6. Plot the Top 6 Genes on the Seurat UMAP
cat("  Plotting top 6 maturation genes...\n")
top_6_genes <- head(gene_df$Gene, 6)

pdf("results/trajectory/03_pseudotime_top_genes_seurat.pdf", width=12, height=8)
print(FeaturePlot(seu_sub, features = top_6_genes, pt.size = 0.5) + theme_classic())
dev.off()
cat("  ✔ Saved: results/trajectory/03_pseudotime_top_genes_seurat.pdf\n")
# ================================================================================





# ================================================================================
#  STEP 11: PSEUDOBULK DEG WITH DESeq2
#  Customized for Seurat v5, Manual_CellType, and accurate metadata
# ================================================================================
cat("━━━ STEP 11: Pseudobulk DEG (DESeq2) ━━━\n")
dir.create("results/DEG", showWarnings = FALSE, recursive = TRUE)

# Load required libraries
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)

# Use the object that is ALREADY loaded in your environment from Step 10!
DefaultAssay(seu) <- "RNA"

# Pseudobulk function — FIX 6: lapply + robust NULL handling + Seurat v5 Layers
run_pseudobulk <- function(seu_obj, label="all_cells") {
  cat(sprintf("  Aggregating: %s\n", label))
  samples <- unique(seu_obj$orig.ident)
  
  # FIX: Use 'layer' for Seurat v5
  pb_list <- lapply(samples, function(s) {
    cells_s <- colnames(seu_obj)[seu_obj$orig.ident == s]
    if (length(cells_s) < 5) return(NULL)
    mat <- GetAssayData(seu_obj, assay="RNA", layer="counts")[, cells_s, drop=FALSE]
    Matrix::rowSums(mat)
  })
  names(pb_list) <- samples
  
  # Remove NULL entries (samples with <5 cells in this cell type)
  pb_list <- pb_list[!sapply(pb_list, is.null)]
  
  if (length(pb_list) < 4) {
    cat(sprintf("    Skipped — only %d valid samples\n", length(pb_list)))
    return(NULL)
  }
  
  pb_mat <- do.call(cbind, pb_list)
  
  # FIX: Safely extract the exact disease status from Seurat metadata, not string guessing
  disease_vec <- seu_obj@meta.data$disease[match(colnames(pb_mat), seu_obj$orig.ident)]
  
  meta <- data.frame(
    sample  = colnames(pb_mat),
    disease = factor(disease_vec, levels=c("Healthy_Control", "Parkinsons_Disease")),
    row.names=colnames(pb_mat)
  )
  
  tryCatch({
    dds  <- DESeqDataSetFromMatrix(countData=round(pb_mat), colData=meta, design=~disease)
    keep <- rowSums(counts(dds) >= 10) >= 2
    dds  <- dds[keep,]
    dds  <- DESeq(dds, quiet=TRUE)
    res  <- results(dds, contrast=c("disease","Parkinsons_Disease","Healthy_Control"), alpha=0.05)
    
    df         <- as.data.frame(res)
    df$gene    <- rownames(df)
    df$label   <- label
    df         <- df[order(df$padj, na.last=TRUE),]
    n_sig      <- sum(df$padj < 0.05, na.rm=TRUE)
    cat(sprintf("    ✔ %d significant DEGs (padj<0.05)\n", n_sig))
    df
  }, error=function(e) {
    cat(sprintf("    ✘ DESeq2 error: %s\n", e$message)); NULL
  })
}

# Overall DEG
cat("  Overall DEG (all cells)...\n")
deg_overall <- run_pseudobulk(seu, "all_cells")

deg_up <- deg_down <- data.frame()
if (!is.null(deg_overall)) {
  deg_up   <- deg_overall %>% filter(log2FoldChange >  1 & padj < 0.05) %>% arrange(desc(log2FoldChange))
  deg_down <- deg_overall %>% filter(log2FoldChange < -1 & padj < 0.05) %>% arrange(log2FoldChange)
  
  cat(sprintf("  Up in PD : %d | Up in HC : %d\n", nrow(deg_up), nrow(deg_down)))
  
  write.csv(deg_overall,"results/DEG/pseudobulk_all_DEGs_PD_vs_HC.csv",  row.names=FALSE)
  write.csv(deg_up,     "results/DEG/pseudobulk_upregulated_in_PD.csv",  row.names=FALSE)
  write.csv(deg_down,   "results/DEG/pseudobulk_downregulated_in_PD.csv",row.names=FALSE)
  
  # VOLCANO PLOT
  deg_overall$color <- "Not Significant"
  deg_overall$color[deg_overall$log2FoldChange >  1 & deg_overall$padj < 0.05] <- "Up in PD"
  deg_overall$color[deg_overall$log2FoldChange < -1 & deg_overall$padj < 0.05] <- "Up in HC"
  
  top_labels <- bind_rows(
    if (nrow(deg_up)   > 0) head(deg_up,   15) else data.frame(),
    if (nrow(deg_down) > 0) head(deg_down, 15) else data.frame())
  
  pdf("results/DEG/01_volcano_plot.pdf", width=11, height=9)
  vp <- ggplot(deg_overall[!is.na(deg_overall$padj),],
               aes(x=log2FoldChange, y=-log10(padj+1e-300), color=color)) +
    geom_point(size=0.9, alpha=0.6) +
    scale_color_manual(values=c("Up in PD"="#E63946","Up in HC"="#457B9D", "Not Significant"="grey75")) +
    geom_vline(xintercept=c(-1,1), linetype="dashed", color="black", linewidth=0.5) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color="black", linewidth=0.5) +
    labs(title    = "Pseudobulk DEG — Parkinson's vs Healthy",
         subtitle = "DESeq2 | LFC threshold=1 | padj<0.05",
         x="Log2 Fold Change (positive = higher in PD)",
         y="-Log10 Adjusted P-value", color="") +
    theme_classic(base_size=13) +
    theme(plot.title=element_text(face="bold",size=15),
          plot.subtitle=element_text(size=11,color="grey40"),
          legend.position="top")
  
  if (nrow(top_labels) > 0) {
    vp <- vp + geom_text_repel(data=top_labels, aes(label=gene),
                               size=3.2, max.overlaps=25, segment.color="grey50")
  }
  print(vp)
  dev.off()
  cat("  ✔ Saved: results/DEG/01_volcano_plot.pdf\n")
  
  # DEG HEATMAP
  heat_genes <- unique(c(
    if(nrow(deg_up)>0)   head(deg_up$gene,25)   else character(0),
    if(nrow(deg_down)>0) head(deg_down$gene,25) else character(0)))
  heat_genes <- heat_genes[heat_genes %in% rownames(seu)]
  
  if (length(heat_genes) > 0) {
    # Ensure genes are scaled before drawing heatmap
    seu <- ScaleData(seu, features = heat_genes, verbose = FALSE)
    n_show  <- min(200, sum(seu$disease=="Healthy_Control"), sum(seu$disease=="Parkinsons_Disease"))
    hc_show <- sample(colnames(seu)[seu$disease=="Healthy_Control"],  n_show)
    pd_show <- sample(colnames(seu)[seu$disease=="Parkinsons_Disease"],n_show)
    seu_sub <- subset(seu, cells=c(hc_show,pd_show))
    Idents(seu_sub) <- "disease"
    
    pdf("results/DEG/02_deg_heatmap.pdf", width=13, height=11)
    print(DoHeatmap(seu_sub, features=heat_genes, group.by="disease",
                    slot="scale.data",
                    group.colors=c(Healthy_Control="steelblue", Parkinsons_Disease="firebrick")) +
            scale_fill_gradientn(colors=c("navy","white","firebrick3")) +
            theme(axis.text.y=element_text(size=8)) +
            ggtitle("Top Pseudobulk DEGs — PD vs HC"))
    dev.off()
    cat("  ✔ Saved: results/DEG/02_deg_heatmap.pdf\n")
  }
}

# Per cell type DEG (FIX: Using Manual_CellType)
cat("\n  Per cell type pseudobulk DEG...\n")
all_ct <- list()
for (ct in unique(seu$Manual_CellType)) {
  ct_cells <- colnames(seu)[seu$Manual_CellType == ct]
  hc_n <- length(unique(seu$orig.ident[seu$Manual_CellType==ct & seu$disease=="Healthy_Control"]))
  pd_n <- length(unique(seu$orig.ident[seu$Manual_CellType==ct & seu$disease=="Parkinsons_Disease"]))
  
  if (hc_n < 3 || pd_n < 3) {
    cat(sprintf("    %-25s → skipped (HC=%d, PD=%d samples)\n",ct,hc_n,pd_n))
    next
  }
  deg_ct <- run_pseudobulk(subset(seu,cells=ct_cells), ct)
  if (!is.null(deg_ct)) {
    all_ct[[ct]] <- deg_ct
    sname <- gsub("[^A-Za-z0-9_]","_",ct)
    write.csv(deg_ct, file.path("results/DEG",paste0("pseudobulk_DEG_",sname,".csv")), row.names=FALSE)
  }
}

if (length(all_ct) > 0) {
  write.csv(bind_rows(all_ct), "results/DEG/pseudobulk_all_DEGs_per_celltype.csv", row.names=FALSE)
  cat("  ✔ Saved: results/DEG/pseudobulk_all_DEGs_per_celltype.csv\n")
}


# ==============================================================================
# PLOT EVERY CELL TYPE (Volcano Loop)
# ==============================================================================
library(ggplot2)
library(ggrepel)
library(dplyr)

cat("Generating Volcano plots for all cell types...\n")

# 1. Find all the individual cell type CSVs in the folder
csv_files <- list.files("results/DEG", pattern="pseudobulk_DEG_.*\\.csv", full.names=TRUE)

# 2. Open a single PDF that will hold all the pages
pdf("results/DEG/04_all_celltypes_volcano.pdf", width=9, height=7)

for (file in csv_files) {
  # Extract the cell type name from the filename
  ct_name <- gsub("pseudobulk_DEG_|\\.csv", "", basename(file))
  cat(sprintf("  Plotting %s...\n", ct_name))
  
  # Read the data
  df <- read.csv(file)
  
  # Set colors based on the 0.5 threshold
  df$color <- "Not Significant"
  df$color[df$log2FoldChange > 0.5 & df$padj < 0.05] <- "Up in PD"
  df$color[df$log2FoldChange < -0.5 & df$padj < 0.05] <- "Up in HC"
  
  # Grab the top genes to label
  top_genes <- df %>% filter(padj < 0.05) %>% arrange(padj) %>% head(15)
  
  # Draw the plot
  vp <- ggplot(df[!is.na(df$padj),], aes(x=log2FoldChange, y=-log10(padj), color=color)) +
    geom_point(alpha=0.6, size=1.2) +
    scale_color_manual(values=c("Up in PD"="#E63946", "Up in HC"="#457B9D", "Not Significant"="grey80")) +
    geom_vline(xintercept=c(-0.5, 0.5), linetype="dashed", color="black") +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color="black") +
    theme_classic(base_size=14) +
    ggtitle(sprintf("%s DEG (Parkinson's vs Healthy)", ct_name)) +
    theme(legend.position="top", legend.title=element_blank())
  
  # Add labels if there are any significant genes
  if (nrow(top_genes) > 0) {
    vp <- vp + geom_text_repel(data=top_genes, aes(label=gene), 
                               size=4, color="black", max.overlaps=30)
  }
  
  # Print this plot as a new page in the PDF
  print(vp)
}

dev.off()
cat("✔ Saved master PDF: results/DEG/04_all_celltypes_volcano.pdf\n")



# ── FINAL SUMMARY ─────────────────────────────────────────────────────────────
n_sig_total <- if(!is.null(deg_overall)) sum(deg_overall$padj<0.05,na.rm=TRUE) else 0

cat("\n╔══════════════════════════════════════════════════════════╗\n")
cat("║     🎉  PIPELINE COMPLETE —  PARTY TIME                  ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")
cat(sprintf("  Samples               : %d\n",  length(unique(seu$orig.ident))))
cat(sprintf("  Total cells           : %d\n",  ncol(seu)))
cat(sprintf("  Genes                 : %d\n",  nrow(seu)))
cat(sprintf("  Clusters              : %d\n",  length(unique(seu$seurat_clusters))))
cat(sprintf("  Cell types            : %d\n",  length(unique(seu$Manual_CellType))))
cat(sprintf("  Pseudobulk DEGs       : %d\n",  n_sig_total))
cat(sprintf("  Up in PD              : %d\n",  nrow(deg_up)))
cat(sprintf("  Up in HC              : %d\n\n",nrow(deg_down)))


# Save final absolute backup
saveRDS(seu, "objects/10_seurat_FINISHED.rds")
cat("💾 FINAL OBJECT SAVED: objects/10_seurat_FINISHED.rds\n")
# ================================================================================



 
# 1. Install BiocManager (required for heavy bioinformatics packages)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# 2. Install the tricky underlying packages that Azimuth relies on
BiocManager::install(c("ComplexHeatmap", "glmGamPoi"), update = FALSE)

# 3. Force a fresh install of the Satija Lab packages
remotes::install_github("satijalab/seurat-data", force = TRUE)
remotes::install_github("satijalab/azimuth", force = TRUE)



##### mess for Azimuth  ##########
# 1. Install 'remotes' if you don't already have it
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# 2. Install SeuratData (required to download the reference atlases)
remotes::install_github("satijalab/seurat-data")

# 3. Install Azimuth
remotes::install_github("satijalab/azimuth")



####    MESS FOR MONOCL3   ####

install.packages("devtools")
devtools::install_github('cole-trapnell-lab/monocle3')


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array', 
                       'terra', 'ggrastr'))

devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')



BiocManager::install("slingshot")




# Install the 'remotes' package if you don't have it
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

# Try installing pre-compiled binaries from R-Packagemanager (Posit)
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/jammy/latest"))
install.packages(c("terra", "nloptr", "lme4", "Cairo", "ggrastr"))
