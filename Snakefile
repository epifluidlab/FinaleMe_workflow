import os
import sys
import glob
import subprocess


def check_tools():
    tools = {
        "java": "java -version",
        "perl": "perl -v",
        "bedGraphToBigWig": "bedGraphToBigWig",
        "bedtools": "bedtools --help",
        "R": "R --version"
    }
    missing_tools = []
    for tool, command in tools.items():
        try:
            if tool == "bedGraphToBigWig":
                subprocess.run(["which", command], check=True, capture_output=True)
            else:
                subprocess.run(command, shell=True, check=True, capture_output=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            missing_tools.append(tool)
    if missing_tools:
        raise SystemExit(f"Error: The following tools are not installed or not in PATH: {', '.join(missing_tools)}.")

check_tools()

BAM_DIR = config.get("input_dir", "input")
OUT_DIR = config.get("out_dir", "output")
SUP_DIR = config.get("sup_dir", "supplement")
SCRIPTS_DIR = config.get("scripts_dir", "scripts")
FINALEME_JAR = config["finaleme_jar_path"]
LIB_DIR = config["lib_dir"]

CP_STEP1 = ":".join([
    FINALEME_JAR,
    os.path.join(LIB_DIR, "gatk-package-distribution-3.3.jar"),
    os.path.join(LIB_DIR, "sis-jhdf5-batteries_included.jar"),
    os.path.join(LIB_DIR, "java-genomics-io.jar"),
    os.path.join(LIB_DIR, "igv.jar")
])
CP_STEP2_3 = ":".join([
    FINALEME_JAR,
    os.path.join(LIB_DIR, "jahmm-0.6.2.jar")
])
CP_STEP5_ALIGNMULTIWIG = ":".join([
    os.path.join(LIB_DIR, "dnaaseUtils-0.14-jar-with-dependencies.jar"), # Assuming this is the one
    os.path.join(LIB_DIR, "java-genomics-io.jar"),
    os.path.join(LIB_DIR, "igv.jar")
])


# INPUT
SAMPLES = [os.path.splitext(os.path.basename(f))[0] for f in glob.glob(os.path.join(BAM_DIR, "*.bam"))]
if not SAMPLES:
    print(f"Warning: No BAM files found in {BAM_DIR}. Workflow might not produce expected outputs.")

rule all:
    input:
        # Main outputs for each sample from Step 4
        expand(os.path.join(OUT_DIR, "{sample}", "{sample}.finaleme.cov.b37.bw"), sample=SAMPLES),
        # Outputs for Tissue of Origin (Step 5) if enabled
        os.path.join(OUT_DIR, "tissue_of_origin", "output.add_value.methy.bed.gz") if config.get("too_enabled") else [],
        os.path.join(OUT_DIR, "tissue_of_origin", "cfdna.names_order.txt") if config.get("too_enabled") else []



# Step 1: Extract features from bam files
rule finaleme_extract_features:
    input:
        bam=os.path.join(BAM_DIR, "{sample}.bam"),
        bai=os.path.join(BAM_DIR, "{sample}.bam.bai"), # or {sample}.bai
        ref_2bit=os.path.join(SUP_DIR, config["reference_2bit"]),
        cpg_motif=os.path.join(SUP_DIR, config["cpg_motif_bedgraph"]),
        exclude_regions=os.path.join(SUP_DIR, config["exclude_regions_bed"]),
        methylation_prior=os.path.join(SUP_DIR, config["methylation_prior_bw"]),
    output:
        features=os.path.join(OUT_DIR, "{sample}", "{sample}.CpgMultiMetricsStats.details.bed.gz")
    params:
        java_xmx=config.get("java_xmx_step1", "20G"),
        finaleme_class="org.cchmc.epifluidlab.finaleme.utils.CpgMultiMetricsStats",
        prefix="-useNoChrPrefixBam" if config['chr_prefix_bam'] else ''
    log:
        os.path.join(OUT_DIR, "{sample}", "logs", "extract_features.log")
    shell:
        """
        echo "Note: extracting features can take a significant amount of time (minutes to hours)"

        CMD="java -Xmx{params.java_xmx} -cp '{CP_STEP1}' {params.finaleme_class} \\
        {input.ref_2bit} \\
        {input.cpg_motif} \\
        {input.cpg_motif} \\
        {input.bam} \\
        {output.features} \\
        -stringentPaired \\
        -excludeRegions {input.exclude_regions} \\
        -valueWigs methyPrior:0:{input.methylation_prior} \\
        -wgsMode \\
        {params.prefix}"

        echo "Executing command: ${{CMD}}"

        eval "$CMD"
        """

# Step 2: Train the HMM model
rule finaleme_train_model:
    input:
        features=rules.finaleme_extract_features.output.features
    output:
        model=os.path.join(OUT_DIR, "{sample}", "{sample}.finaleme.model"),
        training_prediction=os.path.join(OUT_DIR, "{sample}", "{sample}.finaleme.training_prediction.bed.gz")
    params:
        java_xmx=config.get("java_xmx_step2_train", "100G"),
        finaleme_class="org.cchmc.epifluidlab.finaleme.hmm.FinaleMe",
        min_datapoints=config.get("min_datapoints_train", 7),
        gmm_flag="-gmm" if config.get("gmm_train", True) else "",
        cov_outlier=config.get("cov_outlier_train", 3)
    log:
        os.path.join(OUT_DIR, "{sample}", "logs", "train_model.log")
    shell:
        """
        java -Xmx{params.java_xmx} -cp "{CP_STEP2_3}" {params.finaleme_class} \
        {output.model} \
        {input.features} \
        {output.training_prediction} \
        -miniDataPoints {params.min_datapoints} \
        {params.gmm_flag} \
        -covOutlier {params.cov_outlier} > {log} 2>&1
        """

# Step 3: Decode and make the prediction of CpG methylation level
rule finaleme_decode:
    input:
        model=rules.finaleme_train_model.output.model,
        features=rules.finaleme_extract_features.output.features
    output:
        prediction=os.path.join(OUT_DIR, "{sample}", "{sample}.finaleme.prediction.bed.gz")
    params:
        java_xmx=config.get("java_xmx_step3_decode", "100G"),
        finaleme_class="org.cchmc.epifluidlab.finaleme.hmm.FinaleMe"
    log:
        os.path.join(OUT_DIR, "{sample}", "logs", "decode.log")
    shell:
        """
        java -Xmx{params.java_xmx} -cp "{CP_STEP2_3}" {params.finaleme_class} \
        {input.model} \
        {input.features} \
        {output.prediction} \
        -decodeModeOnly > {log} 2>&1
        """

# Step 4: Convert predicted result to .bw file
rule finaleme_convert_to_bigwig:
    input:
        prediction=rules.finaleme_decode.output.prediction,
        chrom_sizes=os.path.join(SUP_DIR, config["chrom_sizes"])
    output:
        cov_bw=os.path.join(OUT_DIR, "{sample}", "{sample}.finaleme.cov.b37.bw"),
        methy_count_bw=os.path.join(OUT_DIR, "{sample}", "{sample}.finaleme.methy_count.b37.bw")  # Critical assumption
    params:
        perl_script=os.path.join(SCRIPTS_DIR, config.get("bw_perl_script")),
        output_prefix=lambda wildcards: os.path.join(OUT_DIR, wildcards.sample, f"{wildcards.sample}.finaleme")
    log:
        os.path.join(OUT_DIR, "{sample}", "logs", "convert_to_bigwig.log")
    shell:
        """
        perl {params.perl_script} {params.output_prefix} {input.prediction} {input.chrom_sizes} > {log} 2>&1
        """

# Step 5: TOO analysis
rule too_generate_alignmultiwig_cmd_file:
    input:
        methy_counts=expand(rules.finaleme_convert_to_bigwig.output.methy_count_bw, sample=SAMPLES),
        covs=expand(rules.finaleme_convert_to_bigwig.output.cov_bw, sample=SAMPLES)
    output:
        cmd_file=os.path.join(OUT_DIR, "tissue_of_origin", "cfdna.methy_summary.cmd.txt")
    run:
        if config.get("too_enabled"):
            os.makedirs(os.path.dirname(output.cmd_file), exist_ok=True)
            with open(output.cmd_file, "w") as f:
                cmd_parts = []
                sorted_samples = sorted(SAMPLES)
                for sample_id in sorted_samples:
                    methy_count_file = os.path.join(OUT_DIR, sample_id, f"{sample_id}.finaleme.methy_count.b37.bw")
                    cov_file = os.path.join(OUT_DIR, sample_id, f"{sample_id}.finaleme.cov.b37.bw")

                    cmd_parts.append(f"-bigWig {methy_count_file} -useMean0 0 -regionMode 0")
                    cmd_parts.append(f"-bigWig {cov_file} -useMean0 0 -regionMode 0")
                f.write(" ".join(cmd_parts))
        else:
            shell("mkdir -p $(dirname {output.cmd_file}) && touch {output.cmd_file}")

rule too_generate_names_order_file:
    output:
        names_file=os.path.join(OUT_DIR, "tissue_of_origin", "cfdna.names_order.txt")
    run:
        if config.get("too_enabled"):
            os.makedirs(os.path.dirname(output.names_file), exist_ok=True)
            with open(output.names_file, "w") as f:
                sorted_samples = sorted(SAMPLES)
                for i, sample_id in enumerate(sorted_samples):
                    f.write(f"{i+1}\t{sample_id}\n")
        else:
            shell("mkdir -p $(dirname {output.names_file}) && touch {output.names_file}")

rule too_run_alignmultiwig:
    input:
        cmd_file=rules.too_generate_alignmultiwig_cmd_file.output.cmd_file,
        intervals_bed=os.path.join(SUP_DIR, config["autosome_1kb_intervals_bed"]),
        methy_counts=expand(rules.finaleme_convert_to_bigwig.output.methy_count_bw, sample=SAMPLES),
        covs=expand(rules.finaleme_convert_to_bigwig.output.cov_bw, sample=SAMPLES)
    output:
        merged_methy_bed=os.path.join(OUT_DIR, "tissue_of_origin", "output.add_value.methy.bed.gz")
    params:
        java_xmx=config.get("java_xmx_step5_alignmultiwig", "10G"),
        tool_class="main.java.edu.mit.compbio.utils.AlignMultiWigInsideBed"
    log:
        os.path.join(OUT_DIR, "tissue_of_origin", "logs", "alignmultiwig.log")
    shell:
        """
        if [ "{config[too_enabled]}" = "True" ]; then
            CMD_ARGS=$(cat {input.cmd_file})
            java -Xmx{params.java_xmx} -cp "{CP_STEP5_ALIGNMULTIWIG}" {params.tool_class} \
            {input.intervals_bed} \
            {output.merged_methy_bed} \
            $CMD_ARGS > {log} 2>&1
        else
            mkdir -p $(dirname {output.merged_methy_bed}) && \
            echo "TOO disabled, creating dummy output for {output.merged_methy_bed}" > {log} && \
            touch {output.merged_methy_bed}
        fi
        """
