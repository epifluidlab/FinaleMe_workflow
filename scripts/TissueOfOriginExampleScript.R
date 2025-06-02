#!/usr/bin/env Rscript

# adapted verison of dr. liu's R script by Gemini  
# --- 0. Load Libraries and Parse Arguments ---
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(quadprog))

option_list <- list(
  make_option(c("-r", "--ref_meth_data"), type="character", default=NULL,
              help="Path to reference panel merged methylation data (e.g., *.add_value.methy.bed.gz)", metavar="FILE"),
  make_option(c("-n", "--ref_names_order"), type="character", default=NULL,
              help="Path to reference panel names/order file", metavar="FILE"),
  make_option(c("-c", "--cfdna_meth_data"), type="character", default=NULL,
              help="Path to cfDNA merged methylation data (e.g., *.add_value.methy.bed.gz)", metavar="FILE"),
  make_option(c("-s", "--cfdna_names_order"), type="character", default=NULL,
              help="Path to cfDNA names/order file", metavar="FILE"),
  make_option(c("-o", "--output_file"), type="character", default="too_results.tsv",
              help="Path for the output TSV file [default %default]", metavar="FILE"),
  make_option(c("--start_col_data"), type="integer", default=7,
              help="Starting column index for methylation/coverage data in BED.gz files [default %default]", metavar="INT"),
  make_option(c("--min_cov_val"), type="numeric", default=10,
              help="Minimum coverage value for a site to be considered. Methylation set to NA if below. [default %default]", metavar="NUM"),
  make_option(c("--meth_bin_threshold"), type="numeric", default=0.1,
              help="Methylation level for binarizing values [default %default]", metavar="NUM"),
  make_option(c("--sd_quantile"), type="numeric", default=0.99,
              help="Quantile for selecting most variable regions from reference panel (e.g., 0.99 for top 1%%) [default %default]", metavar="NUM"),
  make_option(c("--ref_initial_cols_remove"), type="character", default=NULL,
              help="Comma-separated list of column indices to REMOVE from reference panel initially (e.g., '2,4,5,13'). Optional.", metavar="STRING"),
  make_option(c("--ref_reorder_cols_keep"), type="character", default=NULL,
              help="Comma-separated list of column indices to KEEP and REORDER from reference panel after variability filtering (e.g., '6,1,7,5,4,2,8,3,9,10'). Optional, if NULL keeps all remaining.", metavar="STRING"),
  make_option(c("--final_result_threshold"), type="numeric", default=0.001,
              help="Threshold to zero out small values in the final deconvolution result [default %default]", metavar="NUM")
)

opt_parser <- OptionParser(option_list=option_list, description="Tissue of Origin Deconvolution Script")
opt <- parse_args(opt_parser)

# Validate required arguments
required_args <- c("ref_meth_data", "ref_names_order", "cfdna_meth_data", "cfdna_names_order")
missing_args <- required_args[!sapply(required_args, function(arg) !is.null(opt[[arg]]))]
if (length(missing_args) > 0) {
  print_help(opt_parser)
  stop("Missing required arguments: ", paste(missing_args, collapse=", "), call.=FALSE)
}

# Helper function to parse comma-separated integers
parse_cs_integers <- function(cs_string) {
  if (is.null(cs_string) || nchar(trimws(cs_string)) == 0) {
    return(NULL)
  }
  as.integer(strsplit(cs_string, ",")[[1]])
}

# --- 1. Define Tissue of Origin Deconvolution Function ---
tissue_of_origin <- function(dat) {
  dat_matrix <- as.matrix(dat)
  valid_rows <- rowSums(is.na(dat_matrix)) == 0 & rowSums(is.infinite(dat_matrix)) == 0
  dat_filtered <- dat_matrix[valid_rows, , drop = FALSE]

  if (nrow(dat_filtered) < 2 || ncol(dat_filtered) < 2) {
      warning("Not enough valid data points after filtering NA/Inf for deconvolution. Returning NAs.")
      return(rep(NA, ncol(dat) -1 ))
  }

  X <- as.matrix(dat_filtered[, 2:ncol(dat_filtered), drop = FALSE])
  Y <- dat_filtered[, 1, drop = FALSE]

  XtX <- t(X) %*% X
  diag_epsilon <- 1e-8
  if(nrow(XtX) == ncol(XtX)){
     diag(XtX) <- diag(XtX) + diag_epsilon
  } else {
     warning("t(X)%*%X is not square. This should not happen if X has columns. Skipping regularization.")
  }

  Rinv <- tryCatch({
      solve(chol(XtX))
  }, error = function(e) {
      warning("chol(t(X)%*%X) failed, possibly due to singular matrix. Deconvolution may fail. Error: ", e$message)
      return(NULL)
  })

  if(is.null(Rinv)){
      return(rep(NA, ncol(X)))
  }

  C <- cbind(rep(1, ncol(X)), diag(ncol(X)))
  b <- c(1, rep(0, ncol(X)))
  d <- t(Y) %*% X

  result <- tryCatch({
    solve.QP(Dmat = Rinv, factorized = TRUE, dvec = as.vector(d), Amat = C, bvec = b, meq = 1)
  }, error = function(e) {
    warning("solve.QP failed. Error: ", e$message)
    return(NULL)
  })

  if(is.null(result)){
      return(rep(NA, ncol(X)))
  }
  result$solution
}

# --- 2. Load and Process Reference Panel Data ---
message("Loading reference panel data...")
message("  Reference methylation data file: ", opt$ref_meth_data)
raw.1kb <- read.table(opt$ref_meth_data, sep="\t", header=F, stringsAsFactors=FALSE, comment.char="")
message("  Reference names order file: ", opt$ref_names_order)
name.order <- read.table(opt$ref_names_order, sep="\t", header=F, stringsAsFactors=FALSE, comment.char="")
message(paste0("  Raw reference data loaded: ", nrow(raw.1kb), " regions, ", ncol(raw.1kb), " columns."))

message("Processing reference panel methylation matrix...")
methy.mat <- NULL
for (i in seq(opt$start_col_data, length(raw.1kb[1,]), 2)) {
  j <- i + 1
  if (j > length(raw.1kb[1,])) {
      warning(paste("Odd number of data columns or start_col_data too high for file:", opt$ref_meth_data, "Skipping last meth column if present."))
      break
  }
  methyl_counts <- as.numeric(raw.1kb[, i])
  cov <- as.numeric(raw.1kb[, j])
  methyl_val <- ifelse(cov < opt$min_cov_val, NA, methyl_counts / cov) # Set to NA if low coverage
  methy.mat <- cbind(methy.mat, methyl_val)
}

if (ncol(raw.1kb) >= 4) {
    # ***** MAJOR FIX: Keep "chr" prefix for consistency with cfDNA if present in source files *****
    rownames(methy.mat) <- raw.1kb[, 4]
    # Original problematic line: rownames(methy.mat) <- gsub("chr", "", raw.1kb[, 4])
} else {
    warning("Column 4 not found in reference methylation data for rownames. Using generic row numbers.")
    rownames(methy.mat) <- 1:nrow(methy.mat)
}
colnames(methy.mat) <- name.order[, 2]
message(paste0("  Initial reference methy.mat dimensions: ", paste(dim(methy.mat), collapse=" x ")))
message(paste0("  Number of reference samples: ", ncol(methy.mat)))

ref_initial_cols_remove_indices <- parse_cs_integers(opt$ref_initial_cols_remove)
if (!is.null(ref_initial_cols_remove_indices)) {
  message(paste("  Removing initial reference columns by index:", paste(ref_initial_cols_remove_indices, collapse=", ")))
  valid_indices_to_remove <- ref_initial_cols_remove_indices[ref_initial_cols_remove_indices > 0 & ref_initial_cols_remove_indices <= ncol(methy.mat)]
  if (length(valid_indices_to_remove) > 0) {
      methy.mat <- methy.mat[, -valid_indices_to_remove, drop=FALSE]
      name.order <- name.order[-valid_indices_to_remove, , drop=FALSE] # Update name.order accordingly
      message(paste0("    Reference methy.mat dimensions after initial col removal: ", paste(dim(methy.mat), collapse=" x ")))
  } else {
      warning("    No valid column indices provided for initial removal, or all indices out of bounds.")
  }
}

message("Filtering reference panel by standard deviation...")
row.sd <- apply(methy.mat, 1, sd, na.rm = TRUE)
sd_threshold_value <- quantile(row.sd, probs = opt$sd_quantile, na.rm = TRUE)
message(paste0("  SD threshold value at quantile ", opt$sd_quantile, ": ", sd_threshold_value))
methy.mat.mostVar <- methy.mat[which(row.sd >= sd_threshold_value & !is.na(row.sd)), , drop=FALSE] # Ensure row.sd is not NA
message(paste0("  Reference methy.mat.mostVar dimensions after SD filtering: ", paste(dim(methy.mat.mostVar), collapse=" x ")))

message("Cleaning and binarizing most variable reference matrix...")
methy.mat.mostVar <- methy.mat.mostVar[rowSums(is.na(methy.mat.mostVar), na.rm=T) == 0, , drop=FALSE]
message(paste0("    Dimensions after NA row removal: ", paste(dim(methy.mat.mostVar), collapse=" x ")))
methy.mat.mostVar <- methy.mat.mostVar[rowSums(is.infinite(as.matrix(methy.mat.mostVar)), na.rm=T) == 0, , drop=FALSE]
message(paste0("    Dimensions after Inf row removal: ", paste(dim(methy.mat.mostVar), collapse=" x ")))

methy.mat.mostVar[methy.mat.mostVar < opt$meth_bin_threshold] <- 0
methy.mat.mostVar[methy.mat.mostVar >= opt$meth_bin_threshold] <- 1

ref_reorder_cols_keep_indices <- parse_cs_integers(opt$ref_reorder_cols_keep)
if (!is.null(ref_reorder_cols_keep_indices)) {
  message(paste("  Keeping and reordering reference columns by index (from current mostVar):", paste(ref_reorder_cols_keep_indices, collapse=", ")))
  valid_indices_to_keep <- ref_reorder_cols_keep_indices[ref_reorder_cols_keep_indices > 0 & ref_reorder_cols_keep_indices <= ncol(methy.mat.mostVar)]
  if(length(valid_indices_to_keep) > 0) {
      methy.mat.mostVar <- methy.mat.mostVar[, valid_indices_to_keep, drop=FALSE]
      # Update name.order to reflect the columns of methy.mat.mostVar AFTER this reordering step
      # This assumes name.order was correctly tracking columns of methy.mat,
      # and that the indices in ref_reorder_cols_keep_indices refer to the columns of the *current* methy.mat.mostVar
      name.order <- name.order[valid_indices_to_keep, , drop=FALSE] # This line needs careful thought if initial filtering also happened on name.order
      message(paste0("    Reference methy.mat.mostVar dimensions after reorder/keep: ", paste(dim(methy.mat.mostVar), collapse=" x ")))
  } else {
      warning("    No valid column indices provided for reordering/keeping, or all indices out of bounds.")
  }
}
# Ensure name.order reflects the final columns of methy.mat.mostVar
# This is crucial if column filtering/reordering has occurred.
# We assume name.order rows correspond to the columns of methy.mat.mostVar
if(ncol(methy.mat.mostVar) != nrow(name.order)){
    warning("Mismatch between number of columns in final reference matrix and rows in name.order. Colnames for results might be incorrect.")
    # Re-deriving name.order might be safer if complex filtering occurred, but let's assume for now it was handled.
    # If not handled correctly by subsetting name.order alongside methy.mat, then:
    # current_colnames_ref <- colnames(methy.mat.mostVar)
    # name.order <- data.frame(V1=1:length(current_colnames_ref), V2=current_colnames_ref, stringsAsFactors = FALSE)
}


# --- 3. Load and Process cfDNA Data ---
message("Loading cfDNA data...")
message("  cfDNA methylation data file: ", opt$cfdna_meth_data)
cfdna.1kb <- read.table(opt$cfdna_meth_data, sep="\t", header=F, stringsAsFactors=FALSE, comment.char="")
message("  cfDNA names order file: ", opt$cfdna_names_order)
cfdna.name.order <- read.table(opt$cfdna_names_order, sep="\t", header=F, stringsAsFactors=FALSE, comment.char="")
message(paste0("  Raw cfDNA data loaded: ", nrow(cfdna.1kb), " regions, ", ncol(cfdna.1kb), " columns."))

message("Processing cfDNA methylation matrix...")
methy.cfmat <- NULL
for (i in seq(opt$start_col_data, length(cfdna.1kb[1,]), 2)) {
  j <- i + 1
  if (j > length(cfdna.1kb[1,])) {
      warning(paste("Odd number of data columns or start_col_data too high for file:", opt$cfdna_meth_data, "Skipping last meth column if present."))
      break
  }
  methyl_counts_cf <- as.numeric(cfdna.1kb[, i])
  cov_cf <- as.numeric(cfdna.1kb[, j])
  methyl_val_cf <- ifelse(cov_cf < opt$min_cov_val, NA, methyl_counts_cf / cov_cf)
  methy.cfmat <- cbind(methy.cfmat, methyl_val_cf)
}
if (ncol(cfdna.1kb) >= 4) {
    # Assuming cfDNA source also has "chr" and we want consistency
    rownames(methy.cfmat) <- cfdna.1kb[, 4]
} else {
    warning("Column 4 not found in cfDNA methylation data for rownames. Using generic row numbers.")
    rownames(methy.cfmat) <- 1:nrow(methy.cfmat)
}
colnames(methy.cfmat) <- cfdna.name.order[, 2]
message(paste0("  Initial cfDNA methy.cfmat dimensions: ", paste(dim(methy.cfmat), collapse=" x ")))

# message("Cleaning and binarizing cfDNA matrix...")
# methy.cfmat <- methy.cfmat[rowSums(is.na(methy.cfmat), na.rm=T) == 0, , drop=FALSE]
# message(paste0("    Dimensions after NA row removal: ", paste(dim(methy.cfmat), collapse=" x ")))
# methy.cfmat <- methy.cfmat[rowSums(is.infinite(as.matrix(methy.cfmat)), na.rm=T) == 0, , drop=FALSE]
# message(paste0("    Dimensions after Inf row removal: ", paste(dim(methy.cfmat), collapse=" x ")))

methy.cfmat[methy.cfmat < opt$meth_bin_threshold] <- 0
methy.cfmat[methy.cfmat >= opt$meth_bin_threshold] <- 1

# --- 4. Perform Deconvolution ---
message("--- Debugging Rownames Before Intersection ---")
message(paste0("Number of rows in processed reference matrix (methy.mat.mostVar): ", nrow(methy.mat.mostVar)))
if (nrow(methy.mat.mostVar) > 0) {
    message("First 5 rownames of methy.mat.mostVar:")
    print(head(rownames(methy.mat.mostVar), 5))
}
message(paste0("Number of rows in processed cfDNA matrix (methy.cfmat): ", nrow(methy.cfmat)))
if (nrow(methy.cfmat) > 0) {
    message("First 5 rownames of methy.cfmat:")
    print(head(rownames(methy.cfmat), 5))
}
# For more extensive debugging, write all rownames to files:
# write.table(rownames(methy.mat.mostVar), file="debug_ref_rownames.txt", row.names=F, col.names=F, quote=F)
# write.table(rownames(methy.cfmat), file="debug_cfdna_rownames.txt", row.names=F, col.names=F, quote=F)
message("--------------------------------------------")

message("Finding common regions between reference and cfDNA...")
common_regions <- intersect(rownames(methy.cfmat), rownames(methy.mat.mostVar))
message(paste0("Number of common regions found: ", length(common_regions)))

if (length(common_regions) == 0) {
  stop("No common regions found between reference and cfDNA matrices after filtering. Cannot proceed.", call.=FALSE)
}
if (length(common_regions) < 10) { # Arbitrary small number
    message("WARNING: Very few common regions found (<10). Deconvolution might be unstable or unreliable.")
    message("First few common regions:")
    print(head(common_regions, 5))
}


methy.cfmat.common <- methy.cfmat[common_regions, , drop=FALSE]
methy.mat.ref <- methy.mat.mostVar[common_regions, , drop=FALSE]
message(paste0("Dimensions of methy.cfmat.common: ", paste(dim(methy.cfmat.common), collapse=" x ")))
message(paste0("Dimensions of methy.mat.ref: ", paste(dim(methy.mat.ref), collapse=" x ")))


message("Performing tissue of origin deconvolution for each cfDNA sample...")
too_result.1kb <- NULL
if (ncol(methy.cfmat.common) > 0) {
    for (k in 1:ncol(methy.cfmat.common)) {
      current_cfdna_sample_name <- colnames(methy.cfmat.common)[k]
      message(paste0("  Deconvoluting cfDNA sample: ", current_cfdna_sample_name))
      dat_for_deconv <- cbind(methy.cfmat.common[, k, drop=FALSE], methy.mat.ref)
      res <- tissue_of_origin(dat_for_deconv)
      too_result.1kb <- rbind(too_result.1kb, res)
    }

    if (!is.null(too_result.1kb)) {
        colnames(too_result.1kb) <- name.order[, 2]
        rownames(too_result.1kb) <- colnames(methy.cfmat.common)

        message("Applying final result threshold...")
        too_result.1kb[is.na(too_result.1kb)] <- 0
        too_result.1kb[too_result.1kb < opt$final_result_threshold] <- 0
    } else {
        warning("Deconvolution resulted in NULL. No output will be written.")
    }
} else {
    message("No cfDNA samples to process after filtering (methy.cfmat.common has 0 columns).")
    too_result.1kb <- data.frame()
}

# --- 5. Write Output ---
if (!is.null(too_result.1kb) && (nrow(too_result.1kb) > 0 || ncol(too_result.1kb) > 0) ) { # Check if there's something to write
    message(paste("Writing results to:", opt$output_file))
    write.table(too_result.1kb, file=opt$output_file, sep="\t", quote=F, row.names=T, col.names=NA)
} else {
    message(paste("No results to write or deconvolution failed for all samples. Output file may be empty or not created:", opt$output_file))
    # Optionally create an empty file if that's desired pipeline behavior
    # file.create(opt$output_file)
}

message("Script finished.")