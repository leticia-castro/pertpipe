import os
import re
import logging
import pandas as pd
from Bio.Blast import NCBIXML
from scripts import assists
from scripts import prn_assists
#from scripts import draw_figure

prn_seq = os.path.join(os.path.dirname(os.path.dirname(__file__)), "databases/IR_PRN.fasta")
prn_type_seq = os.path.join(os.path.dirname(os.path.dirname(__file__)), "databases/bpertussis/prn.tfa") # all the prn types
fha_type_seq = os.path.join(os.path.dirname(os.path.dirname(__file__)), "databases/fhaB.fasta") # all the fhaB types
is_elements = os.path.join(os.path.dirname(os.path.dirname(__file__)), "databases/IS_elements.fasta") # IS elements.

def virulence_analysis(assembly, prn_outdir, closed, datadir, prokka_outdir, threads=4):
    # commands needed for prn analysis
    abricate_cmd = f"abricate --datadir {assists.bor_vfdb_db} --db bp-only_vfdb --quiet {assembly} > {prn_outdir}/vfdb.txt"
    mlst_cmd = f"mlst --scheme bpertussis --threads {threads} {assembly} > {prn_outdir}/mlst.txt"
    mlst_datadir_cmd = f"mlst --scheme bpertussis --threads {threads} --datadir {datadir}/pubmlst --blastdb {datadir}/blast/mlst.fa {assembly} > {prn_outdir}/mlst.txt"
    blast_cmds = [
        f"blastn -task megablast -num_threads {threads} -query {assembly} -subject {prn_seq} -outfmt 6 -out {prn_outdir}/blast_prn.txt",
        f"blastn -task megablast -num_threads {threads} -query {assembly} -subject {prn_seq} -outfmt 5 -out {prn_outdir}/blast_prn.xml",
        f"blastn -task megablast -num_threads {threads} -query {assembly} -subject {prn_type_seq} -outfmt 6 -min_raw_gapped_score 100 -out {prn_outdir}/blast_prn_type.txt",
        f"blastn -task megablast -num_threads {threads} -query {assembly} -subject {prn_type_seq} -outfmt 5 -min_raw_gapped_score 100 -out {prn_outdir}/blast_prn_type.xml",
        f"blastn -task megablast -num_threads {threads} -query {assembly} -subject {fha_type_seq} -outfmt 6 -min_raw_gapped_score 100 -out {prn_outdir}/blast_fhaB_type.txt",
        f"blastn -task megablast -num_threads {threads} -query {assembly} -subject {fha_type_seq} -outfmt 5 -min_raw_gapped_score 100 -out {prn_outdir}/blast_fhaB_type.xml"
    ]
    #prn_cut_position, prn = None
    #is_string = "?"

    # run the commands
    for command in [abricate_cmd] + blast_cmds:
        assists.run_cmd(command)

    if datadir is not None:
        assists.run_cmd(mlst_datadir_cmd)
    else:
        assists.run_cmd(mlst_cmd)

    mandatory_files = [
        f"{prn_outdir}/vfdb.txt", 
        f"{prn_outdir}/mlst.txt"
    ]

    optional_files = [
        f"{prn_outdir}/blast_prn.txt",
        f"{prn_outdir}/blast_prn.xml",
        f"{prn_outdir}/blast_prn_type.txt",
        f"{prn_outdir}/blast_prn_type.xml",
        f"{prn_outdir}/blast_fhaB_type.txt",
        f"{prn_outdir}/blast_fhaB_type.xml"
    ]
    # check the outputs
    for outfile in mandatory_files:
        assists.check_files(outfile)
    for outfile in optional_files:
        if os.path.isfile(outfile) is True and os.stat(outfile).st_size != 0:
            truemsg = outfile + " exists and not empty, continuing..."
            logging.info(truemsg)
       
    # lets check how many contigs contain prn
    try:
        vfdb_info = pd.read_csv(f"{prn_outdir}/vfdb.txt", sep="\t", header=0)
        prn_vfdb = vfdb_info[vfdb_info['GENE'].str.contains('prn')].reset_index()
        prn_type_info = pd.read_csv(f"{prn_outdir}/blast_prn_type.txt", sep="\t", header=None)
        prn_type_xml = open(f"{prn_outdir}/blast_prn_type.xml")
        prn_promoter = pd.read_csv(f"{prn_outdir}/blast_prn.txt", sep="\t", header=None)
        prn_xml = open(f"{prn_outdir}/blast_prn.xml")
        if 'COVERAGE' in prn_vfdb.columns:
            prn_full_length = prn_vfdb['COVERAGE'].str.contains('1-2733/2733').any()
            if prn_full_length:
                logging.info(f"Duplicate genes including one full length gene found, deleting shorter copy")
                prn_vfdb = prn_vfdb[prn_vfdb['COVERAGE'].str.contains('1-2733/2733', na=False)].copy()
        else:
            prn_full_length = False
    except Exception as e:
        logging.warning(f"Failed to load one or more files: {e}")
        vfdb_info, prn_vfdb, prn_type_info, prn_type_xml, prn_promoter, prn_xml = None, None, None, None, None, None
    file_vars = [vfdb_info, prn_vfdb, prn_type_info, prn_type_xml, prn_promoter, prn_xml]
    
    
    prn_len = len(prn_vfdb)
    # Check if blast has prn type data even if vfdb doesn't
    if prn_len == 0 and os.path.isfile(f"{prn_outdir}/blast_prn_type.txt") and os.stat(f"{prn_outdir}/blast_prn_type.txt").st_size != 0:
        try:
            temp_prn_type_info = pd.read_csv(f"{prn_outdir}/blast_prn_type.txt", sep="\t", header=None)
            if not temp_prn_type_info.empty:
                prn_len = 1  # Trigger partial typing
                prn_vfdb = pd.DataFrame()  # Empty to avoid full length check
        except:
            pass


    # this is the 1 PRN gene checking & pathway
    if not any(file is None for file in file_vars):
        if prn_len == 1:
            # sometimes two prn genes can come up even if there is a full length gene.
            logging.info(f"1 Full length PRN gene detected")
            if 'COVERAGE' in prn_vfdb.columns:
                coverage = prn_vfdb['COVERAGE'][0]
                if coverage == '1-2733/2733': # check if coverage is 100.0 so that we know its full length
                    logging.info(f"Full length PRN gene detected")
                    prn_row, prn_type = prn_assists.prn_type(prn_type_info, "full") # full PRN typing
                    is_prn1or2 = bool(re.findall(r'\bprn[12]\b', prn_type)) # check if its prn1 or prn2, if not find snps
                    if closed == True:  # check if closed genome
                        prn_promoter_xml = NCBIXML.parse(prn_xml)
                        prn_contigs = prn_assists.extract_prn(assembly, prn_promoter_xml, prn_promoter, prn_outdir, "full") # extract only PRN region in closed genomes.
                    # this command only should be used for short read assembled genomes.
                    else:
                        prn_contigs = prn_assists.extract_contigs(assembly, prn_row, prn_outdir) # extracting the contigs matching the prn
                    if is_prn1or2 is True:
                        if prn_row.iloc[0][2] != 100.0: # now check if we need to screen for new mutations!
                            blast_prn_xml = NCBIXML.parse(prn_type_xml)
                            mut_type, mutation = prn_assists.snp_mutations(blast_prn_xml, prn_row, prn_type)
                            prn_type = prn_assists.match_known_prn(mut_type, prn_type, mutation, None)
                    else:
                        prn_type = prn_assists.promoter_scan(prn_promoter, prn_row, prn_type)
                else:
                    logging.info(f"Truncated PRN gene detected")
                    prn_row, prn_type = prn_assists.prn_type(prn_type_info, "partial") # partial PRN typing
                    if closed == True:
                        prn_promoter_xml = NCBIXML.parse(prn_xml)
                        prn_contigs = prn_assists.extract_prn(assembly, prn_promoter_xml, prn_promoter, prn_outdir, "partial")
                    # this command only should be used for short read assembled genomes.
                    else:
                        prn_contigs = prn_assists.extract_contigs(assembly, prn_row, prn_outdir) # extracting the contigs matching the prn
                    blast_cmd_5 = f"blastn -task megablast -query {prn_contigs} -subject {is_elements} -outfmt 6 -out {prn_outdir}/blast_prn_is.txt"
                    assists.run_cmd(blast_cmd_5)
                    if os.path.isfile(prn_outdir + "/blast_prn_is.txt") is True and os.stat(prn_outdir + "/blast_prn_is.txt").st_size != 0:
                        is_prn = pd.read_csv(f"{prn_outdir}/blast_prn_is.txt", sep="\t", header=None)
                        prn_type = prn_assists.dupe_type(prn_promoter, prn_row, is_prn, prn_type)
                    else:
                        prn_type= prn_assists.dupe_type(prn_promoter, prn_row, None, prn_type)
            else:
                logging.warning(f"'COVERAGE' column not found in VFDB output. Assuming truncated PRN gene.")
                prn_row, prn_type = prn_assists.prn_type(prn_type_info, "partial") # partial PRN typing
                if closed == True:
                    prn_promoter_xml = NCBIXML.parse(prn_xml)
                    prn_contigs = prn_assists.extract_prn(assembly, prn_promoter_xml, prn_promoter, prn_outdir, "partial")
                # this command only should be used for short read assembled genomes.
                else:
                    prn_contigs = prn_assists.extract_contigs(assembly, prn_row, prn_outdir) # extracting the contigs matching the prn
                blast_cmd_5 = f"blastn -task megablast -query {prn_contigs} -subject {is_elements} -outfmt 6 -out {prn_outdir}/blast_prn_is.txt"
                assists.run_cmd(blast_cmd_5)
                if os.path.isfile(prn_outdir + "/blast_prn_is.txt") is True and os.stat(prn_outdir + "/blast_prn_is.txt").st_size != 0:
                    is_prn = pd.read_csv(f"{prn_outdir}/blast_prn_is.txt", sep="\t", header=None)
                    prn_type = prn_assists.dupe_type(prn_promoter, prn_row, is_prn, prn_type)
                else:
                    prn_type= prn_assists.dupe_type(prn_promoter, prn_row, None, prn_type)
                

        # this is the 2 PRN genes checking & pathway
        elif prn_len > 1:
            is_prn = pd.DataFrame()
            prn_type = None
            logging.info(f"2 or more PRN genes detected, going down PRN deficient analysis")
            prn_row, prn_type = prn_assists.prn_type(prn_type_info, "dupe") # duplicate PRN typing
            if closed == True:
                prn_promoter_xml = NCBIXML.parse(prn_xml)
                prn_contigs = prn_assists.extract_prn(assembly, prn_promoter_xml, prn_promoter, prn_outdir, "dupe")
            else:
                prn_contigs = prn_assists.extract_contigs(assembly, prn_row, prn_outdir) # extracting the contigs matching the prn
            blast_cmd_5 = f"blastn -task megablast -query {prn_contigs} -subject {is_elements} -outfmt 6 -out {prn_outdir}/blast_prn_is.txt"
            assists.run_cmd(blast_cmd_5)
            if os.path.isfile(prn_outdir + "/blast_prn_is.txt") is True and os.stat(prn_outdir + "/blast_prn_is.txt").st_size != 0:
                is_prn = pd.read_csv(f"{prn_outdir}/blast_prn_is.txt", sep="\t", header=None)
                prn_type = prn_assists.dupe_type(prn_promoter, prn_row, is_prn, prn_type)
            else:
                prn_type = prn_assists.dupe_type(prn_promoter, prn_row, None, prn_type)
    else: 
        prn_row = prn_type = "Not Detected"
    if prn_row is None or prn_type is None:
        prn_row = prn_type = "Not Detected"
    # run fhaB checking now.
    try:
        fhaB_type_info = pd.read_csv(f"{prn_outdir}/blast_fhaB_type.txt", sep="\t", header=None)
        fhaB_type_xml = open(f"{prn_outdir}/blast_fhaB_type.xml")
        fhaB_vfdb = vfdb_info[vfdb_info['GENE'].str.contains('fhaB')].reset_index() 
    except Exception as e:
        logging.warning(f"Failed to load one or more files: {e}")
        vfdb_info, fhaB_type_info, fhaB_type_xml = None, None, None
    file_vars = [vfdb_info, fhaB_type_info, fhaB_type_xml]

    fhaB_type = "Not Detected"

    if not any(file is None for file in file_vars):
        # --- CASE 1: VFDB has exactly one fhaB hit ---
        if len(fhaB_vfdb) == 1:
            max_value = fhaB_type_info[11].max()
            rows_with_max_value = fhaB_type_info[fhaB_type_info[11] == max_value]
            if 'COVERAGE' in fhaB_vfdb.columns:
                coverage = fhaB_vfdb['COVERAGE'].iloc[0]
                if len(rows_with_max_value) == 1 and coverage == '1-10773/10773':
                    # Full-length single hit
                    max_length = rows_with_max_value.iloc[0][3]
                    if max_length > 10700:
                        fhab_len = "full"
                        logging.info(f"Full length fhaB gene detected")
                    else:
                        # Handle abnormal full-length hit                       
                        logging.warning(f"Full length fhaB gene detected but length is shorter than expected. Marking as abnormal.")
                        fhab_len = "abnormal"
                    fhaB_type = prn_assists.fhaB_type(rows_with_max_value, fhab_len)
                else:
                    # Handle truncated single fhaB gene
                    logging.info(f"Truncated fhaB gene detected")
                    fhab_len = "truncated"
                    fhaB_type = prn_assists.fhaB_type(rows_with_max_value, fhab_len)
            else:
                logging.warning(f"'COVERAGE' column not found in VFDB output for fhaB. Assuming truncated.")
                fhab_len = "truncated"
                fhaB_type = prn_assists.fhaB_type(rows_with_max_value, fhab_len)
        # --- CASE 2: VFDB has multiple fhaB hits ---
        elif len(fhaB_vfdb) > 1:                
            if 'COVERAGE' in fhaB_vfdb.columns and any(fhaB_vfdb["COVERAGE"] == "1-10773/10773"):
                # Full-length hit present
                fhaB_full_length = fhaB_vfdb[fhaB_vfdb["COVERAGE"] == "1-10773/10773"].reset_index(drop=True)
                fhaB_contig_name = fhaB_full_length['SEQUENCE'].iloc[0]
                logging.info("Full length fhaB gene detected")
                fhab_len = "full"
                row_with_seq_name = fhaB_type_info[fhaB_type_info[0] == fhaB_contig_name]
                max_value = fhaB_type_info[11].max()
                rows_with_max_value = row_with_seq_name[row_with_seq_name[11] == max_value].reset_index(drop=True)
                fhaB_type = prn_assists.fhaB_type(rows_with_max_value, fhab_len)    
            else:
                # Multiple hits but none full-length → truncated
                logging.info("fhaB gene truncated, selecting best BLAST hit")
                fhab_len = "truncated"
                max_value = fhaB_type_info[11].max()
                rows_with_max_value = fhaB_type_info[fhaB_type_info[11] == max_value].reset_index(drop=True)
                fhaB_type = prn_assists.fhaB_type(rows_with_max_value, fhab_len)
            # --- CASE 3: VFDB has NO fhaB hits but BLAST DOES ---
        elif len(fhaB_vfdb) == 0 and len(fhaB_type_info) > 0:
            logging.info("VFDB has no fhaB hits, but BLAST found truncated fhaB fragments")
            fhab_len = "truncated"
            max_value = fhaB_type_info[11].max()
            rows_with_max_value = fhaB_type_info[fhaB_type_info[11] == max_value].reset_index(drop=True)
            fhaB_type = prn_assists.fhaB_type(rows_with_max_value, fhab_len)
    else:
        fhaB_type = "Not Detected"
    
    #name = os.path.basename(os.path.dirname(prokka_outdir))
    #prokka_gbk, contig_prokka_map = assists.contig_prokka_tag(assembly, name, prokka_outdir)
    #draw_figure.draw_clinker(prn_type, prn_outdir, prokka_gbk, contig_prokka_map)

    # filling the remaining information from MLST/Virulence gene types
    mlst_info = pd.read_csv(f"{prn_outdir}/mlst.txt", sep="\t", header=None)

    # cutting the ptx toxin genes out.
    ptx_toxin_cols = mlst_info.iloc[0, 4:10]

    # removing the ptx suffix from BCDE so that its ptxA(X)B(X)C(X)D(X)
    #ptx_toxin, ptxp, fim2, fim3 = None
    ptx_toxin = ptx_toxin_cols.iloc[1] + ''.join([col[3:] for col in ptx_toxin_cols.iloc[2:]])
    ptxp, fim2, fim3 = mlst_info[4][0], mlst_info[11][0], mlst_info[12][0]

    # final dictionary to be returned.
    virulence_info = {
        "ptxP": ptxp,
        "ptx_toxin": ptx_toxin,
        "prn": prn_type,
        "fim2": fim2,
        "fim3": fim3,
        "fhaB": fhaB_type
    }
    return virulence_info
