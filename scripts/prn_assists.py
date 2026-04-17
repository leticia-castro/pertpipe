import os
from Bio import SeqIO
import logging
import pandas as pd
from scripts import assists

def prn_type(blast_type, length):
    prn_type = None
    blast_type[2] = blast_type[2].astype('float64')
    blast_type[3] = blast_type[3].astype('int64')
    if length == "full":
        blast_type_filter = blast_type[blast_type[3] > 2700].reset_index() # filter lengths to above 2700 bp
        max_value = blast_type_filter[11].max()
        rows_with_max_value = blast_type_filter[blast_type_filter[11] == max_value]
        if len(rows_with_max_value) > 1:
            # Check if any of the rows with the maximum value contain 'prn_1' or 'prn_2'
            prn_condition = rows_with_max_value[1].str.match(r'\bprn_1\b|\bprn_2\b')
            if len(prn_condition) > 0:
                # If there are rows containing 'prn_1' or 'prn_2', choose the first one
                prn_row = rows_with_max_value[prn_condition]
                prn_type = prn_row.iloc[0][1]
            else:
                # Handle the case where no rows with 'prn_1' or 'prn_2' were found
                prn_type = blast_type_filter.loc[blast_type_filter[11].idxmax()][1] # get the highest bitscore after that.
                prn_row = blast_type[blast_type[1] == prn_type] # Or any other desired handling
        else:
            prn_type = blast_type_filter.loc[blast_type_filter[11].idxmax()][1] # get the highest bitscore after that.
            if int(prn_type.removeprefix('prn_')) > 10:
                top_prn_under_10 = blast_type_filter[blast_type_filter[1].apply(lambda x: int(x.split('_')[1]) < 10)]
                prn_type = top_prn_under_10.loc[top_prn_under_10[11].idxmax()][1]
            prn_row = blast_type[blast_type[1] == prn_type] # Or any other desired handling
    elif length == "partial":
        max_value = blast_type[11].max()
        rows_with_max_value = blast_type[blast_type[11] == max_value]
        if len(rows_with_max_value) > 1:
            # Check if any of the rows with the maximum value contain 'prn_1' or 'prn_2'
            prn_condition = rows_with_max_value[1].str.contains(r'\bprn_1\b|\bprn_2\b', regex=True)
            if prn_condition.any():
                # If there are rows containing 'prn_1' or 'prn_2', choose the first one
                prn_row = rows_with_max_value[prn_condition]
                prn_type = prn_row.iloc[0][1]
            else:
                # Handle the case where no rows with 'prn_1' or 'prn_2' were found
                prn_type = rows_with_max_value.loc[rows_with_max_value[11].idxmax()][1] # get the highest bitscore after that.
                prn_row = rows_with_max_value[rows_with_max_value[1] == prn_type] # Or any other desired handling
        else:
            prn_condition = rows_with_max_value[1].str.contains(r'\bprn_1\b|\bprn_2\b', regex=True)
        
        # The below was overriding the above, therefore changed to prn_type == None.
        if prn_type == None:
            prn_type = rows_with_max_value.loc[rows_with_max_value[11].idxmax()][1]  # get the highest bitscore after that
            prn_row = rows_with_max_value[rows_with_max_value[1] == prn_type]  # Or any other desired handling
        else:
            prn_row = rows_with_max_value[prn_condition]
            prn_type = prn_row.iloc[0][1]
    elif length == "dupe":
        cols_to_add = [3] # add the length columns
        cols_to_mean = [2,4,5,6,7,8,9,10,11] # get the average of everything else
        agg_dict = {col: 'mean' for col in cols_to_mean} # start the creation of a dictionary that will perform the adding/averaging
        agg_dict.update({col: 'sum' for col in cols_to_add}) # update for the sum column
        blast_type_filter = blast_type.groupby(1).agg(agg_dict).reset_index()
        max_value = blast_type_filter[11].max()
        top_rows = blast_type_filter[blast_type_filter[11] == max_value]
        if len(top_rows) <= 1:
            prn_type = blast_type_filter.loc[blast_type_filter[11].idxmax()][1]
            if int(prn_type.removeprefix('prn_')) < 10:
                prn_type = prn_type
            else:
                top_prn_under_10 = blast_type_filter[blast_type_filter[1].apply(lambda x: int(x.split('_')[1]) < 10)]
                [blast_type_filter[1].apply(lambda x: int(x.split('_')[1]) < 10)]
                if not top_prn_under_10.empty:
                    prn_type = top_prn_under_10.loc[top_prn_under_10[11].idxmax()][1]
                else:
                    # fall back to overall best hit
                    prn_type = blast_type_filter.loc[blast_type_filter[11].idxmax()][1]
            prn_row = blast_type[blast_type[1] == prn_type]  # ← always set prn_row to the original blast_type filtered by prn_type, so that it will be consistent for downstream processing
        else:
            top_prn_under_10 = top_rows[top_rows[1].apply(lambda x: int(x.split('_')[1]) < 10)]
            if not top_prn_under_10.empty:
                prn_type = top_prn_under_10.loc[top_prn_under_10[11].idxmax()][1]
        prn_row = blast_type[blast_type[1] == prn_type]
    else:
        prn_row = None
        prn_type = None
        logging.critical(f"Unhandled PRN typing edge case: no matching PRN type found for length='{length}'. Please investigate the BLAST results.")
    if prn_type is None and prn_row.empty:
        prn_type = None
    else:
        prn_type = prn_type.replace("_", "")
    return prn_row, prn_type

def extract_contigs(assembly_file, prn_row, prn_outdir):
    match_contigs = prn_row[0].to_list()
    with open(assembly_file, "r") as assembly, open(prn_outdir + "/prn_only.fasta", "w") as output_file:
        current_contig = None
        keep_contig = False

        for line in assembly:
                if line.startswith(">"):
                    current_contig = line.strip()[1:].split()[0]
                    keep_contig = current_contig in match_contigs
                    if keep_contig:
                        output_file.write(line)
                elif keep_contig:
                    output_file.write(line)
    return prn_outdir + "/prn_only.fasta"

def match_known_prn(mut_type, prn_type, prn_cut_start, prn_cut_end):
    mut_name = None
    known_prn = {
       # "146-2733/2733" : "prn2::del(-1513,145)", #H646 prn2
        ('del', 'prn2', -276) : "prn2::del(-283, -40)", #J078
        ('del', 'prn2', 'T662-') : "prn2::del(666, 667)", #J473
        ('snp', 'prn2', 'C1272T') :"prn2::Stop-C1273T", # H697
        ('del', 'prn2', '-292, 1340'): "prn2::del(-292, 1340)", # 24-013-0032 & J625 
        ('del', 'prn2', '-292, 145'): "prn2::del(-1513, 145)", # H696
        ('IS481', 'prn2', 1613) : "prn2::IS481(1613)", # most common IS481 insertion & 24-013-0032
        ('IS481', 'prn2', 245) : "prn2::IS481(240)", #H378 & 17-0520-4681
        ('IS481', 'prn9', 1613) : "prn9::IS481(1613)", # J038 
        ('IS481fwd', 'prn2', 2734) : "prn2::IS481fwd(2734)", #SRR5071044
        ('ins', 'prn2', '-1179G') : "prn2::ins1185G", #H883
        ('del', 'prn2', 'G1490-') : "prn2::Stop-G1490", #H883
        ('full', 'prn22', 3080): "prn2::Stop-C233T", #J018
        ('dis', 'prn2', -73): "prn::dis(-73)", #H806 <- just says "promoter disruption"
        ('del', 'prn1', '1753, 1839'): "prn1::del(1753, 1839)",

        # : "prn2::Stop(C739T)",
    }
    # modifying the naming scheme for mutation type, so that it will show the entire range of deletion/insertion
    if mut_type == 'del' or mut_type == 'ins':
        if prn_cut_end != None:
            prn_cut_position = f"{prn_cut_start}, {prn_cut_end}"
        else:
            prn_cut_position = prn_cut_start
    else:
        prn_cut_position = prn_cut_start

    # choosing the name 
    if mut_type != 'full':
        mut_name = known_prn.get((mut_type, prn_type, prn_cut_position), f"{prn_type}::{mut_type}({prn_cut_position})")
    else:
        mut_name = known_prn.get((mut_type, prn_type, prn_cut_position), prn_type)

    # Check if the mutation name was found in the known_prn dictionary
    is_known = (mut_type, prn_type, prn_cut_position) in known_prn
    if is_known:
        logging.info(f"Predicted Pertactin Type: {mut_name} (Known PRN Mutation - may have been forced changed from {prn_type}::{mut_type}({prn_cut_position}))")
    else:
        logging.info(f"Predicted Pertactin Type: {mut_name} (Unknown PRN Mutation)")
    return mut_name

def snp_mutations(blast_xml, prn_row, prn_type):
    mut_type = None
    mutation = None
    prn_type_mod = prn_type.replace("prn", "prn_")
    hit_list = prn_row[0].to_list()
    fmt_pos_dict = {}
    #prn_type_list = prn_row[1].to_list()
    for blast_result in blast_xml:
        blast_name = blast_result.query.split(" ")[0]
        if blast_name in hit_list:
            for alignment in blast_result.alignments:
                if alignment.hit_id == prn_type_mod:
                    for hsp in alignment.hsps:
                        midline = hsp.match
                        match_count = midline.count('|') + midline.count(' ')
                        space_positions = [pos for pos, char in enumerate(midline) if char == ' ']
                        fmt_pos_dict = {
                            pos: {
                                "sbjct": hsp.sbjct[pos],
                                "pos": pos,
                                "query": hsp.query[pos]
                            }
                            for pos in space_positions
                        }
    if len(fmt_pos_dict) > 0:
        first_key = next(iter(fmt_pos_dict))
        first_entry = fmt_pos_dict[first_key]
        
        sbjct_value = first_entry["sbjct"]
        query_value = first_entry["query"]
        
        # Determine if it's an insertion or deletion
        if sbjct_value == '-' and query_value in 'ATCGN':
            mut_type = "ins"
            mutation = f"{sbjct_value}{first_key}{query_value}"
        elif sbjct_value in 'ATCGN' and query_value == '-':
            mut_type = "del"
            mutation = f"{sbjct_value}{first_key}{query_value}"
        else:
            mut_type = "snp"
            mutation = f"{sbjct_value}{first_key}{query_value}"
    return mut_type, mutation

def dupe_type(prn_promoter, prn_row, is_prn, prn_type):
    # filter the row with starting prn gene.
    mut_type = None
    prn_cut_start = None
    row = prn_row[(prn_row[8] == 1) | (prn_row[9] == 1)]
    if not row.empty:
        row = row.iloc[0]
        if row[8] > row[9] and row[9] == 1:
            prn_cut_start = row[8]
            logging.info(f"prn_cut_start: {prn_cut_start}")
        elif row[8] < row[9] and row[8] == 1:
            prn_cut_start = row[9]
            logging.info(f"prn_cut_start: {prn_cut_start}")
    
    # filter the row with ending prn gene.rc
    other_row = prn_row[~((prn_row[8] == 1) | (prn_row[9] == 1))]
    if not other_row.empty:
        row = other_row.iloc[0]
        if row[8] > row[9]:
            prn_cut_end = row[9]
            logging.info(f"prn_cut_end: {prn_cut_end}")
        elif row[8] < row[9]:
            prn_cut_end = row[8]
            logging.info(f"prn_cut_end: {prn_cut_end}")
    else:
        if prn_type is not None:
            prn_cut_end = assists.get_fasta_length(prn_type.replace("prn", "prn_"))
            mut_type = "dis"
    
    if prn_cut_start != None:
        if prn_cut_start < prn_cut_end:
            logging.info("Deletion likely, setting to deletion")
            mut_type = "del"
        elif prn_cut_start > prn_cut_end:
            logging.info("Disruption detected, could be insertion")
            mut_type = "dis"
    else:
        if prn_type is not None and prn_row.empty is False:
            logging.info("Deletion likely however the start of deletion was not detected, proceeding to promoter checks")
            mut_type = "promoter"
            prn_type = promoter_scan(prn_promoter, prn_row, prn_type)
            
    if is_prn is not None:
        is_string = None
        pident_gr_90 = is_prn[2] > 90
        alnlen_gr_25 = is_prn[3] > 25
        filt_is_prn = is_prn[pident_gr_90 & alnlen_gr_25]
        
        if len(filt_is_prn[1]) >= 1:
            contains_irl = filt_is_prn[1].str.contains("_IRL").any()
            contains_irr = filt_is_prn[1].str.contains("_IRR").any()
            if contains_irl and contains_irr:
                irl_string = filt_is_prn[filt_is_prn[1].str.contains('_IRL')].iloc[0][1].rstrip('_IRL') 
                irr_string = filt_is_prn[filt_is_prn[1].str.contains('_IRR')].iloc[0][1].rstrip('_IRR')
                if irl_string == irr_string:
                    mut_type = irl_string
                else:
                    logging.info("IS elements found, however they are different IS elements.")
            elif contains_irl and not contains_irr:
                logging.info("Only one section of the IS element detected.")
                if row[8] < row[9]:
                    prn_cut_start = int(row[8])
                else:
                    prn_cut_start = int(row[9])
                mut_type = "IS481"
                logging.info(f"IS481 Insertion detected. Using PRN subject start position: {prn_cut_start}")
            elif contains_irr and not contains_irl:
                logging.info("Only one section of the IS element detected.")
                if row[8] < row[9]:
                    prn_cut_start = int(row[8])
                else:
                    prn_cut_start = int(row[9])
                mut_type = "IS481"
                logging.info(f"IS481 Insertion detected. Using PRN subject start position: {prn_cut_start}")
            else:
                logging.info(f"Nothing")
    if prn_type is not None:
        if "::" not in prn_type:
            prn_type = match_known_prn(mut_type, prn_type, prn_cut_start, prn_cut_end)
    return prn_type
    
def promoter_scan(prn_promoter, prn_row, prn_type):
    prn_cut_start = None
    prn_cut_end = None
    prn_type_mod = prn_type + "_IR_PRN"
    prn_results = prn_promoter[prn_promoter[1] == prn_type_mod]
    if prn_results.empty:
        prn_results = prn_promoter.sort_values(by=11, ascending=False).head(1)
    top2_bitscore = prn_results.sort_values(by=11, ascending=False).head(2)
    len_min = top2_bitscore.iloc[:, 8:10].max().max()
    len_max = top2_bitscore.iloc[:, 8:10].min().min()
    # filter the starting row
    row = top2_bitscore[(top2_bitscore == len_max).any(axis=1)]
    if not row.empty:
        row = row.iloc[0]
        if row[8] > row[9]:
            prn_cut_start = row[8]
            if prn_cut_start < 332:
                prn_cut_start = -(332-row[8]-1)
            logging.info(f"promoter_cut_start: {prn_cut_start}")
        elif row[8] < row[9]:
            prn_cut_start = row[9]
            if prn_cut_start < 332:
                prn_cut_start = -(332-row[9]-1)
            logging.info(f"prn_cut_start: {prn_cut_start}")
    
    # filter the row with ending prn gene.
    other_row = top2_bitscore[(top2_bitscore == len_min).any(axis=1)]
    if not other_row.empty:
        row = other_row.iloc[0]
        if row[8] > row[9]:
            prn_cut_end = row[9] - 333
            logging.info(f"prn_cut_end: {prn_cut_end}")
        elif row[8] < row[9]:
            prn_cut_end = row[8] - 333
            logging.info(f"prn_cut_end: {prn_cut_end}")

    if prn_cut_start is None or prn_cut_end is None:
        logging.warning(f"Unable to infer full promoter cut coordinates (start={prn_cut_start}, end={prn_cut_end}). Returning base PRN type {prn_type}.")
        return prn_type

    if prn_cut_start == 3080 and prn_cut_end == -332:
        mut_type = "full"
        prn_type = match_known_prn(mut_type, prn_type, prn_cut_start, prn_cut_end)
    elif prn_cut_start < prn_cut_end:
        logging.info("Deletion likely, setting to deletion")
        mut_type = "del"
        prn_type = match_known_prn(mut_type, prn_type, prn_cut_start, prn_cut_end)

    return prn_type

def extract_prn(assembly, prn_promoter_xml, prn_promoter, prn_outdir, length):
    prn_file = prn_outdir + "/prn_only.fasta"
    max_value = prn_promoter[11].max()
    rows_with_max_value = prn_promoter[prn_promoter[11] == max_value]
    top_hit = rows_with_max_value.iloc[0][1]
    top_hit_length = rows_with_max_value.iloc[0][3]
    if length == "full" or length == "partial":
        for blast_result in prn_promoter_xml:
            for alignment in blast_result.alignments:
                if alignment.hit_id == top_hit:
                    for hsp in alignment.hsps:
                        if hsp.align_length == top_hit_length:
                            with open(prn_file, "w") as output_file:
                                # header
                                output_file.write(">prn_only\n")
                                # sequence
                                output_file.write(hsp.query)
                                logging.info(f"Writing extracted {length} length PRN sequence to {prn_file}")
    elif length == "dupe":
        logging.info(f"Pulling entire length of disrupted PRN gene")
        filter_prn = prn_promoter[prn_promoter[1] == top_hit]
        top2_bitscore = filter_prn.sort_values(by=11, ascending=False).head(2)
        col6_and_col7 = top2_bitscore[[6, 7]]
        end_coord = col6_and_col7.values.max()
        start_coord = col6_and_col7.values.min()
        record = SeqIO.read(assembly, "fasta")
        with open(prn_file, "w") as output_file:
            # header
            output_file.write(">prn_only\n")
            # sequence
            prn_only = record.seq[start_coord-1:end_coord]
            output_file.write(str(prn_only) + "\n")
            logging.info(f"Writing extracted {length} length PRN sequence to {prn_file}")
    return prn_outdir + "/prn_only.fasta"

# def fhaB_type(fhaB_df, fhab_len):
#     if fhab_len == "full":
#         fhaB_type = fhaB_df[1][0]
#     elif fhab_len == "truncated":
#         fhaB_type = f"~{fhaB_df[1][0]} (Truncated)"
#     else:
#         fhaB_type = fhab_len
#     return fhaB_type

def fhaB_type(fhaB_df, fhab_len):
    # Ensure we always work with a DataFrame
    if isinstance(fhaB_df, pd.Series):
        fhaB_df = fhaB_df.to_frame().T

    # Reset index so row 0 always exists
    fhaB_df = fhaB_df.reset_index(drop=True)

    # Column 1 is the allele name (fhaB1, fhaB3, etc.)
    allele = fhaB_df.iloc[0, 1]

    if fhab_len == "full":
        return allele
    elif fhab_len == "truncated":
        return f"~{allele} (Truncated)"
    else:
        return fhab_len
