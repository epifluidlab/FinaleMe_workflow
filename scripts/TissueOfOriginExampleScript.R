#!/usr/bin/env Rscript

# Load necessary library
library(quadprog)
library(optparse)

# Function for Tissue of Origin deconvolution (from the original script)
tissue_of_origin <- function(dat) {
  # Ensure no NAs or Infs which cause errors in solve.QP
  # Check if all rows are valid first
  valid_rows <- rowSums(is.na(dat)) == 0 & rowSums(is.infinite(as.matrix(dat))) == 0
  if (sum(valid_rows) == 0) {
    warning("No valid rows found after filtering NAs and Infs. Returning NAs.")
    # Return a vector of NAs with the correct number of tissues
    num_tissues <- ncol(dat) - 1 # -1 for the cfDNA column
    res_solution <- rep(NA_real_, num_tissues)
    return(res_solution)
  }
  dat <- dat[valid_rows, , drop = FALSE] # drop=FALSE to keep it a matrix even if one row

  X <- as.matrix(dat[, 2:ncol(dat), drop = FALSE]) # Reference tissue matrix
  Y <- dat[, 1] # cfDNA sample vector (single column)

  # Check for rank deficiency in X before cholesky decomposition
  # A simple check: if number of rows (features) < number of columns (tissues)
  if (nrow(X) < ncol(X)) {
    warning(paste0("Number of features (", nrow(X), ") is less than number of reference tissues (", ncol(X), "). This can lead to an underdetermined system. Results might be unreliable. Returning NAs."))
    res_solution <- rep(NA_real_, ncol(X))
    return(res_solution)
  }

  # Add a small constant to the diagonal of t(X) %*% X to ensure positive definiteness (ridge regression like)
  # This helps if the matrix is near singular.
  epsilon <- 1e-6
  dmat <- t(X) %*% X + diag(epsilon, ncol(X))

  # Attempt to solve, with error handling
  result <- tryCatch({
    Rinv <- solve(chol(dmat)) # Use the regularized dmat
    C <- cbind(rep(1, ncol(X)), diag(ncol(X))) # Constraints: sum to 1, non-negative
    b <- c(1, rep(0, ncol(X)))
    d_vec <- t(Y) %*% X # Note: dvec should be t(Y) %*% X, not t(X) %*% Y for solve.QP with Dmat=Rinv
    
    solve.QP(Dmat = Rinv, factorized = TRUE, dvec = d_vec, Amat = C, bvec = b, meq = 1)
  }, error = function(e) {
    warning(paste("Error in solve.QP:", e$message, ". Returning NAs for this sample."))
    return(NULL) # Return NULL to indicate failure
  })

  if (is.null(result)) {
    res_solution <- rep(NA_real_, ncol(X)) # ncol(X) is the number of tissues
  } else {
    res_solution <- result$solution
  }
  return(res_solution)
}


# --- Argument Parsing ---
option_list <- list(
  make_option(c("-c", "--cfdna_meth_file"), type = "character", default = NULL,
              help = "Path to the merged cfDNA methylation data file (e.g., output.add_value.methy.bed.gz)", metavar = "character"),
  make_option(c("-n", "--cfdna_names_file"), type = "character", default = NULL,
              help = "Path to the cfDNA sample names order file (e.g., cfdna.names_order.txt)", metavar = "character"),
  make_option(c("-r", "--ref_meth_file"), type = "character", default = NULL,
              help = "Path to the reference tissue methylation data file", metavar = "character"),
  make_option(c("-s", "--ref_names_file"), type = "character", default = NULL,
              help = "Path to the reference tissue names order file", metavar = "character"),
  make_option(c("-o", "--output_file"), type = "character", default = "tissue_of_origin_results.tsv",
              help = "Path for the output results TSV file [default= %default]", metavar = "character"),
  make_option(c("--ref_start_col"), type = "integer", default = 7,
              help = "Starting column index for methylation/coverage data in the reference methylation file [default= %default]", metavar = "integer"),
  make_option(c("--cfdna_start_col"), type = "integer", default = 7,
              help = "Starting column index for methylation/coverage data in the cfDNA methylation file [default= %default]", metavar = "integer"),
  make_option(c("--min_cov_per_kb"), type = "numeric", default = 10,
              help = "Minimum coverage per 1kb region (scaled) for a site to be considered [default= %default]", metavar = "numeric"),
  make_option(c("--sd_quantile_threshold"), type = "numeric", default = 0.99,
              help = "Quantile threshold for standard deviation to select most variable reference sites (0.0 to 1.0) [default= %default]. For example, 0.99 means top 1% most variable.", metavar = "numeric"),
  make_option(c("--binarize_threshold"), type = "numeric", default = 0.1,
              help = "Methylation threshold for binarizing data (value >= threshold becomes 1, else 0). Set to -1 to disable binarization. [default= %default]", metavar = "numeric")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$cfdna_meth_file) || is.null(opt$cfdna_names_file) || is.null(opt$ref_meth_file) || is.null(opt$ref_names_file)) {
  print_help(opt_parser)
  stop("All input file arguments (--cfdna_meth_file, --cfdna_names_file, --ref_meth_file, --ref_names_file) must be supplied.", call. = FALSE)
}

cat("Starting Tissue of Origin Analysis...\n")
cat("Parameters:\n")
cat("  cfDNA Methylation File:", opt$cfdna_meth_file, "\n")
cat("  cfDNA Names File:", opt$cfdna_names_file, "\n")
cat("  Reference Methylation File:", opt$ref_meth_file, "\n")
cat("  Reference Names File:", opt$ref_names_file, "\n")
cat("  Output File:", opt$output_file, "\n")
cat("  Reference Data Start Column:", opt$ref_start_col, "\n")
cat("  cfDNA Data Start Column:", opt$cfdna_start_col, "\n")
cat("  Min Coverage per 1kb (scaled):", opt$min_cov_per_kb, "\n")
cat("  SD Quantile Threshold for Reference:", opt$sd_quantile_threshold, "\n")
cat("  Binarization Threshold:", opt$binarize_threshold, "\n\n")


# --- Load Reference Panel Data ---
cat("Loading reference panel data...\n")
ref_raw_data <- read.table(opt$ref_meth_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE, comment.char="")
ref_name_order <- read.table(opt$ref_names_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE, comment.char="")

ref_methy_mat <- NULL
# Assumes BED-like format with chr, start, end, name_or_id in first 4 columns
# and then pairs of (methylation_sum, coverage_sum) for each sample
# The name/ID for the region is expected in column 4 for rownames
region_ids_ref <- ref_raw_data[, 4] # Assuming region ID is in column 4

cat("Processing reference methylation matrix...\n")
for (i in seq(opt$ref_start_col, ncol(ref_raw_data), 2)) {
  j <- i + 1
  if (j > ncol(ref_raw_data)) {
    warning(paste("Odd number of data columns in reference file after start column", opt$ref_start_col, ". Skipping last column if it exists."))
    break
  }
  methy_sum <- as.numeric(ref_raw_data[, i])
  cov <- as.numeric(ref_raw_data[, j])
  
  # Calculate methylation fraction, handle low coverage by setting methylation to NA
  # The original script scaled coverage by 1000, opt$min_cov_per_kb is the actual coverage value
  s <- ifelse(cov < opt$min_cov_per_kb, NA_real_, methy_sum / cov)
  s[is.nan(s)] <- NA_real_ # Ensure 0/0 results in NA
  
  ref_methy_mat <- cbind(ref_methy_mat, s)
}
if (!is.null(ref_methy_mat)) {
    rownames(ref_methy_mat) <- region_ids_ref # Use region IDs from column 4
    colnames(ref_methy_mat) <- ref_name_order[, 2] # Tissue names
} else {
    stop("Failed to process reference methylation matrix. Check input file format and start column.")
}


# Filter reference matrix: select most variable sites
cat("Filtering reference matrix for most variable sites...\n")
row_sd_ref <- apply(ref_methy_mat, 1, sd, na.rm = TRUE)
# Remove rows that are all NA before quantile calculation, or sd will be NA
valid_sd_rows <- !is.na(row_sd_ref)
if (sum(valid_sd_rows) == 0) {
    stop("No valid standard deviations calculated for reference matrix. All rows might be NA or have no variance.")
}

# Calculate quantile on valid SDs
quantile_val <- quantile(row_sd_ref[valid_sd_rows], probs = opt$sd_quantile_threshold, na.rm = TRUE)
cat(paste("  SD threshold for selection (quantile", opt$sd_quantile_threshold, "):", round(quantile_val, 4), "\n"))

# Select rows where SD is greater than or equal to the quantile value AND SD is not NA
methy_mat_ref_var <- ref_methy_mat[which(row_sd_ref >= quantile_val & valid_sd_rows), , drop = FALSE]

# Further filter out rows that became all NA after selection or still have any NAs/Infs
# (This is important if some tissues had NAs for the selected variable sites)
methy_mat_ref_var <- methy_mat_ref_var[rowSums(is.na(methy_mat_ref_var)) == 0 & 
                                       rowSums(is.infinite(as.matrix(methy_mat_ref_var))) == 0, , drop = FALSE]

if (nrow(methy_mat_ref_var) == 0) {
  stop("No regions left in reference matrix after variability and NA/Inf filtering. Try adjusting --sd_quantile_threshold or check data quality.")
}
cat(paste("  Selected", nrow(methy_mat_ref_var), "most variable regions for reference panel.\n"))

# Binarize reference data if threshold is not -1
if (opt$binarize_threshold != -1) {
  cat(paste("  Binarizing reference data at threshold:", opt$binarize_threshold, "\n"))
  methy_mat_ref_var[methy_mat_ref_var < opt$binarize_threshold] <- 0
  methy_mat_ref_var[methy_mat_ref_var >= opt$binarize_threshold] <- 1
}


# --- Load cfDNA Data ---
cat("Loading cfDNA data...\n")
cfdna_raw_data <- read.table(opt$cfdna_meth_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE, comment.char="")
cfdna_name_order <- read.table(opt$cfdna_names_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE, comment.char="")

cfdna_methy_mat <- NULL
region_ids_cfdna <- cfdna_raw_data[, 4] # Assuming region ID is in column 4

cat("Processing cfDNA methylation matrix...\n")
for (i in seq(opt$cfdna_start_col, ncol(cfdna_raw_data), 2)) {
  j <- i + 1
   if (j > ncol(cfdna_raw_data)) {
    warning(paste("Odd number of data columns in cfDNA file after start column", opt$cfdna_start_col, ". Skipping last column if it exists."))
    break
  }
  methy_sum <- as.numeric(cfdna_raw_data[, i])
  cov <- as.numeric(cfdna_raw_data[, j])
  
  s <- ifelse(cov < opt$min_cov_per_kb, NA_real_, methy_sum / cov)
  s[is.nan(s)] <- NA_real_
  
  cfdna_methy_mat <- cbind(cfdna_methy_mat, s)
}
if (!is.null(cfdna_methy_mat)) {
    rownames(cfdna_methy_mat) <- region_ids_cfdna # Use region IDs from column 4
    colnames(cfdna_methy_mat) <- cfdna_name_order[, 2] # cfDNA sample names
} else {
    stop("Failed to process cfDNA methylation matrix. Check input file format and start column.")
}


# Filter cfDNA matrix to remove rows with any NAs or Infs
cfdna_methy_mat <- cfdna_methy_mat[rowSums(is.na(cfdna_methy_mat)) == 0 & 
                                   rowSums(is.infinite(as.matrix(cfdna_methy_mat))) == 0, , drop = FALSE]

if (nrow(cfdna_methy_mat) == 0) {
  stop("No regions left in cfDNA matrix after NA/Inf filtering. Check data quality or coverage filter.")
}

# Binarize cfDNA data if threshold is not -1
if (opt$binarize_threshold != -1) {
  cat(paste("  Binarizing cfDNA data at threshold:", opt$binarize_threshold, "\n"))
  cfdna_methy_mat[cfdna_methy_mat < opt$binarize_threshold] <- 0
  cfdna_methy_mat[cfdna_methy_mat >= opt$binarize_threshold] <- 1
}


# --- Find Common Regions ---
cat("Finding common regions between reference and cfDNA...\n")
common_regions <- intersect(rownames(cfdna_methy_mat), rownames(methy_mat_ref_var))

if (length(common_regions) == 0) {
  stop("No common regions found between processed cfDNA and reference matrices. Check region IDs and filtering steps.")
}
cat(paste("  Found", length(common_regions), "common regions for deconvolution.\n"))

cfdna_methy_common <- cfdna_methy_mat[common_regions, , drop = FALSE]
ref_methy_common <- methy_mat_ref_var[common_regions, , drop = FALSE]


# --- Perform Deconvolution for each cfDNA sample ---
cat("Performing deconvolution...\n")
final_too_results <- NULL

for (sample_idx in 1:ncol(cfdna_methy_common)) {
  cfdna_sample_name <- colnames(cfdna_methy_common)[sample_idx]
  cat(paste("  Processing cfDNA sample:", cfdna_sample_name, "\n"))
  
  # Combine current cfDNA sample with the reference matrix
  # Ensure the cfDNA sample is the first column as expected by tissue_of_origin function
  deconvolution_input_dat <- cbind(cfdna_methy_common[, sample_idx, drop = FALSE], ref_methy_common)
  colnames(deconvolution_input_dat)[1] <- "cfDNA_Sample" # Temporary name for clarity
  
  sample_result_fractions <- tissue_of_origin(deconvolution_input_dat)
  
  if (is.null(final_too_results)) {
    final_too_results <- matrix(sample_result_fractions, nrow = 1)
  } else {
    final_too_results <- rbind(final_too_results, sample_result_fractions)
  }
}

# Set column and row names for the final result matrix
if (!is.null(final_too_results)) {
  colnames(final_too_results) <- colnames(ref_methy_common) # Reference tissue names
  rownames(final_too_results) <- colnames(cfdna_methy_common) # cfDNA sample names

  # Set very small values to 0 for clarity
  final_too_results[final_too_results < 0.001 & !is.na(final_too_results)] <- 0
} else {
    warning("Deconvolution resulted in NULL. No output will be written.")
}


# --- Write Output ---
if (!is.null(final_too_results)) {
    cat(paste("Writing results to:", opt$output_file, "\n"))
    write.table(final_too_results, file = opt$output_file, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA) # col.names=NA for R-friendly tsv
} else {
    cat("No results to write.\n")
}

cat("Tissue of Origin analysis complete.\n")