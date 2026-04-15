import os
import pandas as pd
import logging
from Bio.Blast import NCBIXML
from scripts import assists
from Bio.Seq import Seq


rrna_seq = os.path.join(os.path.dirname(os.path.dirname(__file__)), "databases/23S_rRNA.fasta")
def mres_detection(assembly, outdir, meta, threads=4):
    hit_list = []
    mutation_list = []
    mutation_counter = 0
    detected = False
    blast_cmd = f"blastn -task megablast -num_threads {threads} -query {assembly} -subject {rrna_seq} -outfmt 6 -out {outdir}/blast_23s.txt"
    assists.run_cmd(blast_cmd)
    blast_cmd_2 = f"blastn -task megablast -num_threads {threads} -query {assembly} -subject {rrna_seq} -outfmt 5 -out {outdir}/blast_23s.xml"
    assists.run_cmd(blast_cmd_2)
    if meta: 
        try:
            blast_df = pd.read_csv(f"{outdir}/blast_23s.txt", sep="\t", header=None)
        except pd.errors.EmptyDataError:
            blast_df = None
    else:
        blast_df = pd.read_csv(f"{outdir}/blast_23s.txt", sep="\t", header=None)
    
    if blast_df is not None:
        result_count = len(blast_df)
        if result_count > 1:
            filt_blast_df = blast_df[blast_df[3] >= 2880]
            if filt_blast_df.empty:
                logging.error(f"Encountered issue, potentially truncated 23S rRNA detected, or 23S rRNA is very noisy.")
                detected = False
        elif result_count == 1: 
            filt_blast_df = blast_df
        
        filt_count = len(filt_blast_df)  
        for hit in filt_blast_df.iterrows():
            aln_len, perc_id, hit_name = hit[1][3], hit[1][2], hit[1][0]
            hit_list.append(hit_name)
            if perc_id == 100.0 and aln_len > 2880 and aln_len < 2885:
                logging.info(f'Full length 23s rRNA detected, with no mutations')
                mutation_list = []
                detected = True
            elif perc_id >= 99.9 and perc_id < 100 and aln_len >= 2880 and aln_len <= 2885:
                logging.info(f'23s rRNA detected with mutations')
                try:
                    with open(f"{outdir}/blast_23s.xml") as xml_handle:
                        blast_xml = NCBIXML.parse(xml_handle)
                        mutation_counter += 1
                        mutations = mres_position(blast_xml, hit_list)
                        mutation_list.extend(mutations)
                        logging.info(f"{mutation_list}")
                except Exception as e:
                    logging.warning(f"Failed to parse BLAST XML for 23S mutations: {e}. Assuming mutations present but unable to specify.")
                    mutation_list.append("Unknown")  # Indicate mutations detected but not specified
                detected = True
            else:
                logging.error(f"Encountered issue, potentially truncated 23S rRNA detected, or error in assembly occurred.")
                detected = False
    else:
        logging.info(f"No 23S rRNA detected through Blast")
        
    # unique_mutations = list(set(mutation_list))
    # logging.info(f"Final mutation list: {unique_mutations}")
    # vdomain_start, vdomain_end = 1918, 2444
    # return unique_mutations, mutation_counter, detected

    unique_mutations = list(set(mutation_list))
    vdomain_start, vdomain_end = 1918, 2444
    vdomain_mutations = []
    for m in unique_mutations:
        digits = ''.join(filter(str.isdigit, m))
        if digits:
            try:
                pos = int(digits)
                if vdomain_start <= pos <= vdomain_end:
                    vdomain_mutations.append(m)
            except ValueError:
                pass  # Skip invalid mutations
    logging.info(f"Final V domain mutation list: {vdomain_mutations}")
    return vdomain_mutations, mutation_counter, detected
    
def mres_position(blast_xml, hit_list):
    for blast_result in blast_xml:
        accession_id = blast_result.query.split(" ")
        if accession_id[0] in hit_list:
            for alignment in blast_result.alignments:
                for hsp in alignment.hsps:
                    if hsp.align_length >= 2880 and hsp.strand[1] == 'Plus':
                        logging.info(f"BLAST alignment in forward orientation")
                        midline = hsp.match
                        match_count = midline.count('|') + midline.count(' ')
                        space_positions = [pos for pos, char in enumerate(midline) if char == ' ']
                        formatted_positions = [
                            f"{hsp.sbjct[pos]}{pos + 1}{hsp.query[pos]}"
                            for pos in space_positions
                        ]
                        return formatted_positions
                    elif hsp.align_length >= 2880 and hsp.strand[1] == 'Minus':
                        logging.info("BLAST alignment in reverse orientation, reverse complement required.")
                        midline = "".join(reversed(hsp.match)) #reverse the midline.
                        query = Seq(hsp.query).reverse_complement()
                        sbjct = Seq(hsp.sbjct).reverse_complement()
                        space_positions = [pos for pos, char in enumerate(midline) if char == ' ']
                        formatted_positions = [
                            f"{sbjct[pos]}{pos + 1}{query[pos]}"
                            for pos in space_positions
                        ]
                        return formatted_positions
                    else:
                        logging.error(f"Unexpected BLAST alignment: length={hsp.align_length}, strand={hsp.strand}. Could not determine 23S rRNA mutation positions.")
    
