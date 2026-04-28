# ================================================================================
#  STEP 2: LOAD DATA + QUALITY CONTROL
# ================================================================================
cat("━━━ STEP 2: Load Data + QC ━━━\n")

seurat_raw <- list()

# Loop ONLY through the samples you defined in MANUAL_LABELS
for (gsm_id in names(MANUAL_LABELS)) {
  
  # Find the specific folder that starts with this GSM ID in your data folder
  sdir_matches <- list.dirs("data", full.names=TRUE, recursive=FALSE)
  sdir <- sdir_matches[grepl(gsm_id, sdir_matches)]
  
  # If the folder was deleted by you, skip it gracefully without crashing
  if (length(sdir) == 0) {
    cat("  ⚠ Skipping:", gsm_id, "— Folder not found (likely deleted).\n")
    next
  }
  
  sdir <- sdir[1]
  sname <- basename(sdir)
  
  cat("  Loading:", sname, "...")
  
  tryCatch({
    counts <- Read10X(data.dir=sdir)
    if (is.list(counts)) counts <- counts[[1]]
    
    # Create the raw Seurat object
    obj <- CreateSeuratObject(counts=counts, project=sname,
                              min.cells=3, min.features=200)
    
    # ── APPLY MANUAL LABELS ──
    obj$sample  <- sname
    obj$gsm_id  <- gsm_id
    obj$disease <- MANUAL_LABELS[[gsm_id]]
    
    seurat_raw[[sname]] <- obj
    cat(" ✔", ncol(obj), "cells |", obj$disease[1], "\n")
  }, error=function(e) cat(" ✘ FAILED:", e$message, "\n"))
}

# Print the overall breakdown
cat("\n  Total raw cells loaded:", sum(sapply(seurat_raw, ncol)), "\n")
cat("  Disease breakdown across all cells:\n")
print(table(sapply(seurat_raw, function(x) x$disease[1])))

# ── QC METRICS & PLOTTING ───────────────────────────────────────────────────
# [Keep the rest of your Step 2 code exactly the same from here down]