import os
import logging
from io import StringIO

import pandas as pd
from scripts import assists

rrna_seq = os.path.join(os.path.dirname(os.path.dirname(__file__)), "databases/23S_rRNA.fasta")
rrna_bed = os.path.join(os.path.dirname(os.path.dirname(__file__)), "databases/23s_rRNA_VDomain.bed")
ref_list = os.path.join(os.path.dirname(os.path.dirname(__file__)), "databases/silva_bor.bed")
kallisto_db = os.path.join(os.path.dirname(os.path.dirname(__file__)), "databases/SILVA_138.2_LSURef_NR99")


def mres_map(reads1, reads2, outdir, mutation_list, meta, threads=4):
    """
    Run read-based 23S rRNA mapping and mutation calling.

    Returns
    -------
    dict
        result_dict with keys: Resistance, Mutation, Copy No
    """
    # Optional meta filtering via kallisto
    if meta is True:
        reads1, reads2 = _filter_reads_with_kallisto(reads1, reads2, outdir)

    # Map reads to 23S reference
    minimap2_cmd = (
        f"minimap2 -t {threads} -ax sr {rrna_seq} {reads1} {reads2} "
        f"| samtools view -b > {outdir}/23s_aln.bam"
    )
    assists.run_cmd(minimap2_cmd)

    # Sort, index, depth
    sort_bam = f"{outdir}/23s_aln.sort.bam"
    samtools_sort_cmd = f"samtools sort {outdir}/23s_aln.bam > {sort_bam}"
    samtools_index_cmd = f"samtools index {sort_bam}"
    depth_txt = f"{outdir}/depth.txt"
    samtools_depth_cmd = f"samtools depth -a {sort_bam} -b {rrna_bed} > {depth_txt}"

    for cmd in [samtools_sort_cmd, samtools_index_cmd, samtools_depth_cmd]:
        assists.run_cmd(cmd)

    # Coverage calculation
    coverage = _calculate_coverage(depth_txt)

    # Fresh VCF
    vcf_path = f"{outdir}/mres.vcf"
    if os.path.exists(vcf_path):
        os.remove(vcf_path)

    bcftools_cmd = (
        f"bcftools mpileup -f {rrna_seq} {sort_bam} | "
        f"bcftools call -mv -Ov -o {vcf_path}"
    )
    assists.run_cmd(bcftools_cmd)

    # Parse VCF and build result
    if mutation_list:
        mres_df, result_dict = vcf_calc_and_blast_match(vcf_path, mutation_list, coverage)
    else:
        result_dict, mres_df = map_calculations(vcf_path, coverage)

    # Always write a table, even if empty
    out_txt = f"{outdir}/23s_mres.txt"
    if mres_df is None:
        mres_df = pd.DataFrame(columns=["REFPOSALT", "FREQ"])
    mres_df.to_csv(out_txt, sep="\t", index=False, header=False)
    logging.info(
        f"Any other mutations detected in the 23S rRNA V domain will be written to {out_txt}"
    )

    return result_dict


def _filter_reads_with_kallisto(reads1, reads2, outdir):
    """Run kallisto + samtools filtering; fall back to original reads on any failure."""
    kallisto_cmd = f"kallisto quant -i {kallisto_db} -o {outdir} {reads1} {reads2} --pseudobam"
    kallisto_result = assists.check_kallisto_finished(outdir)

    if not kallisto_result:
        try:
            assists.run_cmd(kallisto_cmd)
        except SystemExit:
            if assists.check_kallisto_finished(outdir):
                logging.warning(
                    "Kallisto exited with non-zero status, but abundance.tsv exists. "
                    "Continuing with available kallisto outputs."
                )
            else:
                raise
    else:
        logging.info("Kallisto has already finished for this sample. Skipping.")

    abundance_file = os.path.join(outdir, "abundance.tsv")
    if not os.path.exists(abundance_file):
        logging.info("Abundance file not found. Using original reads for 23S mapping.")
        return reads1, reads2

    try:
        abundance_df = pd.read_csv(abundance_file, sep="\t")
    except Exception as e:
        logging.info(f"Failed to read abundance.tsv ({e}). Using original reads.")
        return reads1, reads2

    total_counts = abundance_df.get("est_counts", pd.Series(dtype=float)).sum()
    if total_counts == 0:
        logging.info(
            "No reads pseudoaligned to SILVA database. "
            "Skipping read filtering and using original reads for 23S mapping."
        )
        return reads1, reads2

    pseudo_bam = os.path.join(outdir, "pseudoalignments.bam")
    if not assists.check_bam_readable(pseudo_bam):
        if os.path.exists(pseudo_bam):
            os.remove(pseudo_bam)
        logging.info(
            "Pseudoalignments BAM is unreadable. Re-running kallisto to regenerate it."
        )
        try:
            assists.run_cmd(kallisto_cmd)
        except SystemExit:
            logging.warning(
                "Failed to regenerate pseudoalignments.bam. "
                "Using original reads for 23S mapping."
            )
            return reads1, reads2

    if not assists.check_bam_readable(pseudo_bam):
        logging.info(
            "Unable to create a readable pseudoalignments.bam. "
            "Using original reads for 23S mapping."
        )
        return reads1, reads2

    # Filter to Bordetella reads
    bor_bam = f"{outdir}/kallisto_bor.bam"
    r1_out = f"{outdir}/kallisto_bor_R1.fastq"
    r2_out = f"{outdir}/kallisto_bor_R2.fastq"
    samtools_cmd_1 = f"samtools view -b -h -F 4 -L {ref_list} {pseudo_bam} > {bor_bam}"
    samtools_cmd_2 = f"samtools fastq -1 {r1_out} -2 {r2_out} {bor_bam}"

    try:
        for cmd in [samtools_cmd_1, samtools_cmd_2]:
            assists.run_cmd(cmd)
        return r1_out, r2_out
    except Exception as e:
        logging.info(
            f"Failed to filter reads with samtools ({e}). "
            "Using original reads for 23S mapping."
        )
        return reads1, reads2


def _calculate_coverage(depth_txt):
    """Compute mean coverage from samtools depth output."""
    if not os.path.exists(depth_txt) or os.stat(depth_txt).st_size == 0:
        return 0

    try:
        depth_df = pd.read_csv(depth_txt, sep="\t", header=None)
        # third column is depth
        coverage = depth_df.iloc[:, 2].mean(skipna=True)
        return float(coverage) if pd.notnull(coverage) else 0.0
    except Exception as e:
        logging.info(f"Failed to parse depth file ({e}). Treating coverage as 0.")
        return 0


def vcf_calc_and_blast_match(bcftools_vcf, mutation_list, coverage):
    """
    Parse VCF, restrict to V-domain, match against known mutation_list.

    Returns
    -------
    mres_df : pd.DataFrame
        Filtered mutations (REFPOSALT, FREQ) for requested mutation_list.
    result_dict : dict
        Resistance summary.
    """
    positions = None
    mres_df = pd.DataFrame(columns=["REFPOSALT", "FREQ"])

    if not os.path.exists(bcftools_vcf) or os.stat(bcftools_vcf).st_size == 0:
        logging.info("VCF file missing or empty. No 23S mutations detected in mapping.")
        return mres_df, _result_from_flags(coverage, positions, 0, mutation_list)

    with open(bcftools_vcf, "r") as vcf:
        oneline = "".join(line for line in vcf if not line.startswith("##"))

    if not oneline.strip():
        logging.info("VCF contains no variant records. No 23S mutations detected in mapping.")
        return mres_df, _result_from_flags(coverage, positions, 0, mutation_list)

    vcf_df = pd.read_csv(StringIO(oneline), sep="\t", header=0)

    if vcf_df.empty:
        logging.info("VCF dataframe empty after parsing. No 23S mutations detected in mapping.")
        return mres_df, _result_from_flags(coverage, positions, 0, mutation_list)

    info_data = vcf_df["INFO"].str.split(";", expand=True)
    if info_data.empty:
        logging.info("INFO field empty in VCF. No 23S mutations detected in mapping.")
        return mres_df, _result_from_flags(coverage, positions, 0, mutation_list)

    # Build DP4-based frequencies
    tmp_df = pd.DataFrame()
    tmp_df["DP4"] = info_data.apply(
        lambda row: next((cell for cell in row if isinstance(cell, str) and "DP4=" in cell), None),
        axis=1,
    ).str.lstrip("DP4=")

    split_columns = tmp_df["DP4"].str.split(",", expand=True)
    split_columns = split_columns.apply(pd.to_numeric, errors="coerce")
    tmp_df["sum_all_four"] = split_columns.sum(axis=1)
    tmp_df["sum_last_two"] = split_columns.iloc[:, -2:].sum(axis=1)
    tmp_df["FREQ"] = tmp_df["sum_last_two"] / tmp_df["sum_all_four"]

    # Restrict to V-domain
    vdomain_start, vdomain_end = 1918, 2444
    mres_df = pd.DataFrame()
    mres_df["POS"] = vcf_df["POS"]
    mres_df["REFPOSALT"] = vcf_df["REF"] + vcf_df["POS"].astype(str) + vcf_df["ALT"]
    mres_df["FREQ"] = tmp_df["FREQ"]

    mres_df = mres_df[(mres_df["POS"] >= vdomain_start) & (mres_df["POS"] <= vdomain_end)]
    if mres_df.empty:
        logging.info("No variants in 23S V-domain in mapping.")
        return mres_df[["REFPOSALT", "FREQ"]], _result_from_flags(coverage, None, 0, mutation_list)

    # Filter to requested mutation_list
    mask = mres_df["REFPOSALT"].isin(mutation_list)
    filtered = mres_df[mask].reset_index(drop=True)

    if not filtered.empty:
        positions = ",".join(mutation_list)
        copy_no_str = filtered["FREQ"].apply(determine_copy_number).iloc[0]
        logging.info(f"23S mutation occurs as a {positions}, very likely in {copy_no_str}")
        copy_no = copy_no_str
    else:
        positions = None
        copy_no = 0
        logging.info("Mutations detected in assembly were not detected in mapping.")

    result_dict = _result_from_flags(coverage, positions, copy_no, mutation_list)
    return filtered[["REFPOSALT", "FREQ"]], result_dict


def _result_from_flags(coverage, positions, copy_no, mutation_list):
    """Centralised logic for building result_dict."""
    # Low coverage overrides everything
    if coverage < 10 and positions is not None:
        return {
            "Resistance": "Low Coverage!",
            "Mutation": positions,
            "Copy No": str(copy_no),
        }
    if coverage < 10 and positions is None:
        return {
            "Resistance": "Low Coverage!",
            "Mutation": "N/A",
            "Copy No": "N/A",
        }

    # High coverage cases
    if "A2037G" in mutation_list and positions is not None:
        return {
            "Resistance": "Resistant",
            "Mutation": positions,
            "Copy No": str(copy_no),
        }

    if positions is None:
        return {
            "Resistance": "Susceptible",
            "Mutation": "N/A",
            "Copy No": "N/A",
        }

    return {
        "Resistance": "Mutations in 23S rRNA V domain detected",
        "Mutation": positions,
        "Copy No": str(copy_no),
    }


def determine_copy_number(freq):
    perc = "{:.1%}".format(freq)
    if 0.68 <= freq <= 1:
        return f"3 copies (Freq: {perc})"
    elif 0.316 < freq < 0.68:
        return f"2 copies (Freq: {perc})"
    elif 0.05 <= freq <= 0.316:
        return f"1 copy (Freq: {perc})"
    else:
        return f"(Freq: {perc})"


def map_calculations(bcftools_vcf, coverage):
    """
    Generic mapping-based 23S mutation detection (no prior mutation_list).

    Returns
    -------
    result_dict : dict
    mres_df : pd.DataFrame
        All V-domain mutations (REFPOSALT, FREQ).
    """
    mres_df = pd.DataFrame(columns=["REFPOSALT", "FREQ"])

    if coverage <= 10:
        return {
            "Resistance": "Low Coverage!",
            "Mutation": "N/A",
            "Copy No": "N/A",
        }, mres_df

    if not os.path.exists(bcftools_vcf) or os.stat(bcftools_vcf).st_size == 0:
        logging.info("VCF file missing or empty. No 23S mutations detected in mapping.")
        return {
            "Resistance": "Susceptible",
            "Mutation": "N/A",
            "Copy No": "N/A",
        }, mres_df

    with open(bcftools_vcf, "r") as vcf:
        oneline = "".join(line for line in vcf if not line.startswith("##"))

    if not oneline.strip():
        logging.info("VCF contains no variant records. No 23S mutations detected in mapping.")
        return {
            "Resistance": "Susceptible",
            "Mutation": "N/A",
            "Copy No": "N/A",
        }, mres_df

    vcf_df = pd.read_csv(StringIO(oneline), sep="\t", header=0)
    if vcf_df.empty:
        logging.info("VCF dataframe empty after parsing. No 23S mutations detected in mapping.")
        return {
            "Resistance": "Susceptible",
            "Mutation": "N/A",
            "Copy No": "N/A",
        }, mres_df

    info_data = vcf_df["INFO"].str.split(";", expand=True)
    tmp_df = pd.DataFrame()
    tmp_df["DP4"] = info_data.apply(
        lambda row: next((cell for cell in row if isinstance(cell, str) and "DP4=" in cell), None),
        axis=1,
    ).str.lstrip("DP4=")

    split_columns = tmp_df["DP4"].str.split(",", expand=True)
    split_columns = split_columns.apply(pd.to_numeric, errors="coerce")
    tmp_df["sum_all_four"] = split_columns.sum(axis=1)
    tmp_df["sum_last_two"] = split_columns.iloc[:, -2:].sum(axis=1)
    tmp_df["FREQ"] = tmp_df["sum_last_two"] / tmp_df["sum_all_four"]

    vdomain_start, vdomain_end = 1918, 2444
    mres_df = pd.DataFrame()
    mres_df["POS"] = vcf_df["POS"]
    mres_df["REFPOSALT"] = vcf_df["REF"] + vcf_df["POS"].astype(str) + vcf_df["ALT"]
    mres_df["FREQ"] = tmp_df["FREQ"]
    mres_df = mres_df[(mres_df["POS"] >= vdomain_start) & (mres_df["POS"] <= vdomain_end)]

    if mres_df.empty:
        logging.info("No 23S mutations were detected in this sample")
        return {
            "Resistance": "Susceptible",
            "Mutation": "N/A",
            "Copy No": "N/A",
        }, mres_df[["REFPOSALT", "FREQ"]]

    copy_no_list = mres_df["FREQ"].apply(determine_copy_number).tolist()
    positions = mres_df["REFPOSALT"].tolist()
    logging.info(f"23S mutation occurs as a {positions}, very likely in {copy_no_list}")

    # A2037G check
    if mres_df["REFPOSALT"].str.contains("A2037G").any():
        result_dict = {
            "Resistance": "Resistant",
            "Mutation": ",".join(positions),
            "Copy No": ",".join(copy_no_list),
        }
    else:
        result_dict = {
            "Resistance": "Mutations in 23S rRNA V domain detected",
            "Mutation": ", ".join(positions),
            "Copy No": ", ".join(copy_no_list),
        }

    return result_dict, mres_df[["REFPOSALT", "FREQ"]]
