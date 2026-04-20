import os
import logging
import pandas as pd
from io import StringIO
from scripts import assists

rrna_seq = os.path.join(os.path.dirname(os.path.dirname(__file__)), "databases/23S_rRNA.fasta")
rrna_bed = os.path.join(os.path.dirname(os.path.dirname(__file__)), "databases/23s_rRNA_VDomain.bed")
ref_list = os.path.join(os.path.dirname(os.path.dirname(__file__)), "databases/silva_bor.bed")
kallisto_db = os.path.join(os.path.dirname(os.path.dirname(__file__)), "databases/SILVA_138.2_LSURef_NR99")


def mres_map(reads1, reads2, outdir, mutation_list, meta, threads=4):
    if meta is True:
        kallisto_cmd = f"kallisto quant -i {kallisto_db} -o {outdir} {reads1} {reads2} --pseudobam"
        kallisto_result = assists.check_kallisto_finished(outdir)
        if kallisto_result is False:
            try:
                assists.run_cmd(kallisto_cmd)
            except SystemExit:
                if assists.check_kallisto_finished(outdir):
                    logging.warning("Kallisto exited with non-zero status, but abundance.tsv exists. Continuing with available kallisto outputs.")
                else:
                    raise
        else:
            logging.info(f"Kallisto has already finished for this sample. Skipping.")

        abundance_file = os.path.join(outdir, "abundance.tsv")
        if os.path.exists(abundance_file):
            # import pandas as pd
            abundance_df = pd.read_csv(abundance_file, sep="\t")
            total_counts = abundance_df['est_counts'].sum()
            if total_counts == 0:
                logging.info("No reads pseudoaligned to SILVA database. Skipping read filtering and using original reads for 23S mapping.")
            else:
                pseudo_bam = os.path.join(outdir, "pseudoalignments.bam")
                if not assists.check_bam_readable(pseudo_bam):
                    if os.path.exists(pseudo_bam):
                        os.remove(pseudo_bam)
                    logging.info("Pseudoalignments BAM is unreadable. Re-running kallisto to regenerate it.")
                    try:
                        assists.run_cmd(kallisto_cmd)
                    except SystemExit:
                        logging.warning("Failed to regenerate pseudoalignments.bam. Using original reads for 23S mapping.")

                if assists.check_bam_readable(pseudo_bam):
                    samtools_cmd_1 = f"samtools view -b -h -F 4 -L {ref_list} {pseudo_bam} > {outdir}/kallisto_bor.bam"
                    samtools_cmd_2 = f"samtools fastq -1 {outdir}/kallisto_bor_R1.fastq -2 {outdir}/kallisto_bor_R2.fastq {outdir}/kallisto_bor.bam"
                    try:
                        for cmd in [samtools_cmd_1, samtools_cmd_2]:
                            assists.run_cmd(cmd)
                        reads1 = f"{outdir}/kallisto_bor_R1.fastq"
                        reads2 = f"{outdir}/kallisto_bor_R2.fastq"
                    except Exception as e:
                        logging.info(f"Failed to filter reads with samtools: {e}. Using original reads for 23S mapping.")
                else:
                    logging.info("Unable to create a readable pseudoalignments.bam. Using original reads for 23S mapping.")
        else:
            logging.info("Abundance file not found. Using original reads for 23S mapping.")
            # Keep original reads1, reads2
    minimap2_cmd = f"minimap2 -t {threads} -ax sr {rrna_seq} {reads1} {reads2} | samtools view -b > {outdir}/23s_aln.bam"
    assists.run_cmd(minimap2_cmd)
    samtools_sort_cmd = f"samtools sort {outdir}/23s_aln.bam > {outdir}/23s_aln.sort.bam"
    samtools_index_cmd = f"samtools index {outdir}/23s_aln.sort.bam"
    samtools_depth_cmd = f"samtools depth -a {outdir}/23s_aln.sort.bam -b {rrna_bed} > {outdir}/depth.txt"
    for cmd in [samtools_sort_cmd, samtools_index_cmd, samtools_depth_cmd]:
        assists.run_cmd(cmd)
    if os.path.exists(f"{outdir}/depth.txt") is True and os.stat(f"{outdir}/depth.txt").st_size != 0:
        depth_df = pd.read_csv(f"{outdir}/depth.txt", sep="\t")
        coverage = depth_df.iloc[:, 2].mean(skipna=True, numeric_only=True)
    else:
        coverage = 0
        #os.remove(f"{outdir}/depth.txt")
    if os.path.exists(f"{outdir}/mres.vcf") is True:
        os.remove(f"{outdir}/mres.vcf")
    bcftools_cmd = f"bcftools mpileup -f {rrna_seq} {outdir}/23s_aln.sort.bam | bcftools call -mv -Ov -o {outdir}/mres.vcf"
    assists.run_cmd(bcftools_cmd)

    if mutation_list != []:
        mres_df, result_dict = vcf_calc_and_blast_match(f"{outdir}/mres.vcf", mutation_list, coverage)
    else:    
        result_dict, mres_df = map_calculations(f"{outdir}/mres.vcf", coverage)
    mres_df.to_csv(f"{outdir}/23s_mres.txt", sep="\t", index=False, header=False) 
    logging.info(f"Any other mutations detected in the 23S rRNA V domain will be written to {outdir}/23s_mres.txt")
    return result_dict

def vcf_calc_and_blast_match(bcftools_vcf, mutation_list, coverage):
    positions = None
    with open(bcftools_vcf, "r") as vcf:
        oneline = ""
        lines = vcf.readlines()
        for line in lines:
            if not line.startswith("##"):
                oneline += line
    
    vcf_df = pd.read_csv(StringIO(oneline), sep="\t", header=0)

    info_data = vcf_df["INFO"].str.split(";", expand = True)
    if not info_data.empty:
        mres_df = pd.DataFrame()
        copy_no = 0
        mres_df['DP4'] = info_data.apply(lambda row: next((cell for cell in row if "DP4=" in cell), None), axis=1).str.lstrip('DP4=')
        split_columns = mres_df['DP4'].str.split(',', expand=True)
        split_columns = split_columns.apply(pd.to_numeric, errors='coerce')
        mres_df['sum_all_four'] = split_columns.sum(axis=1).astype(int)
        mres_df['sum_last_two'] = split_columns.iloc[:, -2:].sum(axis=1).astype(int)
        mres_df['FREQ'] = mres_df['sum_last_two'] / mres_df['sum_all_four']

        # extract only V domain
        vdomain_start, vdomain_end = 1918, 2444
        mres_df["REFPOSALT"], mres_df["POS"] = vcf_df["REF"] + vcf_df["POS"].astype(str) + vcf_df["ALT"], vcf_df["POS"]
        mres_df = mres_df[(mres_df["POS"] >= vdomain_start) & (mres_df["POS"] <= vdomain_end)]
        mres_df["REFPOSALT"] = vcf_df["REF"] + vcf_df["POS"].astype(str) + vcf_df["ALT"]
        
        copy_no = mres_df["FREQ"].apply(determine_copy_number)
        final_df = mres_df[['REFPOSALT', 'FREQ']]
        mask = final_df['REFPOSALT'].isin(mutation_list)
        mres_df = mres_df[mask].reset_index()
        if not mres_df.empty:
            positions = ",".join(mutation_list)
            copy_no = mres_df["FREQ"].apply(determine_copy_number)[0]
            logging.info(f"23s mutation occurs as a {positions}, very likely in {copy_no}")
    else:
        final_df = pd.DataFrame()
        positions = None
        copy_no = 0
        logging.info(f"Mutations detected in assembly was not detected in mapping.")
    
    if coverage < 10 and positions is not None:
        result_dict = {
            "Resistance": "Low Coverage!",
            "Mutation": positions,
            "Copy No": f"{str(copy_no)} copies",
        }
    elif coverage < 10 and positions is None:
        result_dict = {
            "Resistance": "Low Coverage!",
            "Mutation": "N/A",
            "Copy No": "N/A",
        }
    elif "A2037G" in mutation_list:
        result_dict = {
            "Resistance": "Resistant",
            "Mutation": positions,
            "Copy No": f"{str(copy_no)} copies",
        }
    elif positions is None:
        result_dict = {
            "Resistance": "Susceptible",
            "Mutation": "N/A",
            "Copy No": "N/A",
        }
    else: 
        result_dict = {
            "Resistance": "Mutations in 23S rRNA V domain detected",
            "Mutation": positions,
            "Copy No": f"{str(copy_no)} copies",
        }
        
    return final_df, result_dict


def determine_copy_number(freq):
    perc = '{:.1%}'.format(freq)
    if 0.68 <= freq <= 1:
        return f"3 copies (Freq: {perc})"
    elif 0.316 < freq < 0.68:
        return f"2 copies (Freq: {perc})"
    elif 0.05 <= freq <= 0.316:
        return f"1 copy (Freq: {perc})"
    else:
        return f"(Freq: {perc})"

def map_calculations(bcftools_vcf, coverage):
    mres_df = pd.DataFrame()
    if coverage > 10:
        with open(bcftools_vcf, "r") as vcf:
            oneline = ""
            lines = vcf.readlines()
            for line in lines:
                if not line.startswith("##"):
                    oneline += line
        
        vcf_df = pd.read_csv(StringIO(oneline), sep="\t", header=0)
        
        if vcf_df.empty is False:
            info_data = vcf_df["INFO"].str.split(";", expand = True)
            
            mres_df['DP4'] = info_data.apply(lambda row: next((cell for cell in row if "DP4=" in cell), None), axis=1).str.lstrip('DP4=')
            split_columns = mres_df['DP4'].str.split(',', expand=True)
            split_columns = split_columns.apply(pd.to_numeric, errors='coerce')
            mres_df['sum_all_four'] = split_columns.sum(axis=1).astype(int)
            mres_df['sum_last_two'] = split_columns.iloc[:, -2:].sum(axis=1).astype(int)
            mres_df['FREQ'] = mres_df['sum_last_two'] / mres_df['sum_all_four']
            
            # extract only V domain
            vdomain_start, vdomain_end = 1918, 2444
            mres_df["REFPOSALT"], mres_df["POS"] = vcf_df["REF"] + vcf_df["POS"].astype(str) + vcf_df["ALT"], vcf_df["POS"]
            mres_df = mres_df[(mres_df["POS"] >= vdomain_start) & (mres_df["POS"] <= vdomain_end)]
            mut_total = len(mres_df.index)
            copy_no = mres_df["FREQ"].apply(determine_copy_number).to_list()
            mres_df = mres_df[['REFPOSALT', 'FREQ']]
            positions = mres_df["REFPOSALT"].to_list()
            logging.info(f"23S mutation occurs as a {positions}, very likely in {copy_no}")
            if mres_df["REFPOSALT"].str.contains("A2037G").any():
                result_dict = {
                    "Resistance": "Resistant",
                    "Mutation": ",".join(positions),
                    "Copy No": ",".join(copy_no)
                }
            elif positions != []:
                result_dict = {
                    "Resistance": "Mutations in 23S rRNA V domain detected",
                    "Mutation": ", ".join(positions),
                    "Copy No": ", ".join(copy_no)
                }
            else:
                logging.info(f"No 23S mutations were detected in this sample")
                result_dict = {
                    "Resistance": "Susceptible",
                    "Mutation": "N/A",
                    "Copy No": "N/A"
                    }
        else:
            logging.info(f"No 23S mutations were detected in this sample")
            result_dict = {
                "Resistance": "Susceptible",
                "Mutation": "N/A",
                "Copy No": "N/A"
                }
    else: 
        result_dict = {
            "Resistance": "Low Coverage!",
            "Mutation": "N/A",
            "Copy No": "N/A",
        } 
    return result_dict, mres_df