#!/bin/bash

# Gemini (AI) generated script to generate a merged reference methylation file from the Bigwigs in ref_panel_wgbs.tar 



# --- BEGIN CONFIGURATION - EDIT THESE PATHS AND SETTINGS ---

# Path to the directory containing your FinaleMe_workflow (parent of 'lib', 'reference', 'supplement', 'output')
WORKFLOW_DIR="$(pwd)" # Assumes you run this script from inside FinaleMe_workflow directory. Adjust if not.

# Path to your reference BigWig files
REFERENCE_BW_DIR="${WORKFLOW_DIR}/reference"

# Path to your supplementary files (where the intervals BED file is)
SUPPLEMENT_DIR="${WORKFLOW_DIR}/supplement"

# Name of the genomic intervals BED file
# (This should match what's in your Snakemake config['autosome_1kb_intervals_bed'])
INTERVALS_BED_FILE="autosome_1kb_intervals.UCSC.cpgIsland_plus_shore.b37.bed" # ADJUST THIS FILENAME if different

# Path to the FinaleMe library JARs
LIB_DIR="${WORKFLOW_DIR}/lib"

# Output file name for the merged reference methylation data
OUTPUT_MERGED_REF_METHY_FILE="${REFERENCE_BW_DIR}/reference_panel.merged_methy.bed.gz"

# Name of your manually created reference panel names order file (must exist in REFERENCE_BW_DIR)
REFERENCE_NAMES_ORDER_FILE="reference_panel.names_order.txt" # Make sure this file exists and is correct!

# Java memory allocation for AlignMultiWigInsideBed
JAVA_XMX="10G"

# --- END CONFIGURATION ---

echo "--- Starting Reference Methylome Generation ---"

# --- 1. Define Classpath for AlignMultiWigInsideBed ---
# Adjust JAR names if they are different in your lib/ directory
CP_ALIGNMULTIWIG="${LIB_DIR}/dnaaseUtils-0.14-jar-with-dependencies.jar:${LIB_DIR}/java-genomics-io.jar:${LIB_DIR}/igv.jar"

# Check if essential files and directories exist
if [ ! -d "$REFERENCE_BW_DIR" ]; then
    echo "ERROR: Reference BigWig directory not found: $REFERENCE_BW_DIR"
    exit 1
fi
if [ ! -f "${REFERENCE_BW_DIR}/${REFERENCE_NAMES_ORDER_FILE}" ]; then
    echo "ERROR: Reference names order file not found: ${REFERENCE_BW_DIR}/${REFERENCE_NAMES_ORDER_FILE}"
    echo "Please create this file manually before running the script."
    exit 1
fi
if [ ! -f "${SUPPLEMENT_DIR}/${INTERVALS_BED_FILE}" ]; then
    echo "ERROR: Intervals BED file not found: ${SUPPLEMENT_DIR}/${INTERVALS_BED_FILE}"
    exit 1
fi
if [ ! -f "${LIB_DIR}/dnaaseUtils-0.14-jar-with-dependencies.jar" ]; then # Check one of the JARs
    echo "ERROR: A required JAR file was not found in $LIB_DIR. Check CP_ALIGNMULTIWIG and LIB_DIR."
    exit 1
fi

echo "Using Reference BigWig directory: $REFERENCE_BW_DIR"
echo "Using Intervals BED file: ${SUPPLEMENT_DIR}/${INTERVALS_BED_FILE}"
echo "Using Classpath: $CP_ALIGNMULTIWIG"
echo "Output will be: $OUTPUT_MERGED_REF_METHY_FILE"

# --- 2. Construct Command Arguments for AlignMultiWigInsideBed ---
ALIGNMULTI_ARGS=""
MISSING_BW_COUNT=0

# Read sample prefixes from the names order file
while IFS=$'\t' read -r index sample_prefix || [[ -n "$sample_prefix" ]]; do
    # Skip if sample_prefix is empty (e.g., blank lines in names file)
    if [ -z "$sample_prefix" ]; then
        continue
    fi

    METH_COUNT_BW="${REFERENCE_BW_DIR}/${sample_prefix}.meth_count.bw"
    READ_BW="${REFERENCE_BW_DIR}/${sample_prefix}.read.bw"

    echo "Processing sample prefix: $sample_prefix"

    if [ ! -f "$METH_COUNT_BW" ]; then
        echo "  WARNING: Meth count BigWig not found: $METH_COUNT_BW"
        MISSING_BW_COUNT=$((MISSING_BW_COUNT + 1))
    fi
    if [ ! -f "$READ_BW" ]; then
        echo "  WARNING: Read count BigWig not found: $READ_BW"
        MISSING_BW_COUNT=$((MISSING_BW_COUNT + 1))
    fi

    ALIGNMULTI_ARGS+="-bigWig ${METH_COUNT_BW} -useMean0 0 -regionMode 0 -bigWig ${READ_BW} -useMean0 0 -regionMode 0 "
done < "${REFERENCE_BW_DIR}/${REFERENCE_NAMES_ORDER_FILE}"


if [ "$MISSING_BW_COUNT" -gt 0 ]; then
    echo "ERROR: $MISSING_BW_COUNT BigWig file(s) listed in ${REFERENCE_NAMES_ORDER_FILE} were not found. Please check your files and the names order file."
    exit 1
fi

if [ -z "$ALIGNMULTI_ARGS" ]; then
    echo "ERROR: No arguments generated for AlignMultiWig. Is ${REFERENCE_NAMES_ORDER_FILE} empty or incorrectly formatted?"
    exit 1
fi

echo "Generated AlignMultiWig arguments."
# For debugging, you can uncomment the next line to see the arguments:
# echo "ARGS: $ALIGNMULTI_ARGS"

# --- 3. Run AlignMultiWigInsideBed ---
echo "Running AlignMultiWigInsideBed... This may take some time."

COMMAND="java -Xmx${JAVA_XMX} -cp \"${CP_ALIGNMULTIWIG}\" main.java.edu.mit.compbio.utils.AlignMultiWigInsideBed \
    \"${SUPPLEMENT_DIR}/${INTERVALS_BED_FILE}\" \
    \"${OUTPUT_MERGED_REF_METHY_FILE}\" \
    ${ALIGNMULTI_ARGS}"

echo "Executing: $COMMAND"

# Execute the command
eval "$COMMAND"

if [ $? -eq 0 ]; then
    echo "Successfully generated merged reference methylome: ${OUTPUT_MERGED_REF_METHY_FILE}"
else
    echo "ERROR: AlignMultiWigInsideBed failed. Check the output above for errors."
    exit 1
fi

echo "--- Reference Methylome Generation Complete ---"
echo ""
echo "You should now have the following files for your REFERENCE PANEL:"
echo "1. Merged Methylation Data: ${OUTPUT_MERGED_REF_METHY_FILE} (This is your 'raw.1kb' equivalent)"
echo "2. Names Order File:      ${REFERENCE_BW_DIR}/${REFERENCE_NAMES_ORDER_FILE} (This is your 'name.order')"
echo ""
echo "For your cfDNA SAMPLES (already generated by Snakemake):"
echo "1. Merged Methylation Data: ${WORKFLOW_DIR}/output/tissue_of_origin/output.add_value.methy.bed.gz (This is your 'cfdna.1kb')"
echo "2. Names Order File:      ${WORKFLOW_DIR}/output/tissue_of_origin/cfdna.names_order.txt (This is your 'cfdna.name.order')"
echo ""
echo "You can now proceed to update your R script with these file paths."