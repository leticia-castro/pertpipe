import os
import subprocess
import sys
import logging
import shutil
import os.path
import pkg_resources
from Bio import SeqIO
from Bio import GenBank

bor_vfdb_db = os.path.join(os.path.dirname(os.path.dirname(__file__)), "databases")

def run_cmd(command):
    """
    Run commands with error outputs.
    """
    logging.info("Running command: %s", command)
    result = subprocess.run(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True
    )
    result.stdout = result.stdout.decode()
    result.stderr = result.stderr.decode()
    if result.returncode != 0:
        logging.critical("Failed to run command: %s", result.args)
        logging.critical("stdout: %s", result.stdout)
        logging.critical("stderr: %s", result.stderr)
        sys.exit(1)
    return result


def check_files(file):
    """
    Check input files if they exist and have contents
    """

    if os.path.isfile(file) is True and os.stat(file).st_size != 0:
        truemsg = file + " exists and not empty, continuing..."
        logging.info(truemsg)
    else:
        msg = (
            file
            + " either file is does not exist or is empty, please check files. Exiting."
        )
        logging.critical(msg)
        sys.exit(1)


def check_folders(folder):
    """
    Check the output folder if it exists, if not make new directory.
    """
    if os.path.exists(folder) is True:
        truemsg = folder + " output folder exists"
        logging.info(truemsg)
    else:
        os.makedirs(folder)
        msg = folder + " does not exist, making output folder"
        logging.info(msg)


def check_dependencies(cmd_exec):
    cmd_path = shutil.which(cmd_exec)
    if cmd_exec == "kallisto":
        vcmd = subprocess.run([cmd_exec, "version"], capture_output=True, text=True)
    else:
        vcmd = subprocess.run([cmd_exec, "--version"], capture_output=True, text=True)
    if vcmd.returncode == 0 and vcmd.stdout != '':
        result = vcmd.stdout.splitlines()
    else:
        result = vcmd.stderr.splitlines()
    if cmd_exec == "abricate":
        version = " ".join(result).replace("abricate ", "v")
        if pkg_resources.parse_version(version) < pkg_resources.parse_version(
            "0.9.8"
        ):
            logging.critical("Abricate version too old, please upgrade to v1.0.0+")
            sys.exit(1)
    elif cmd_exec == "spades.py":
        version = result[0].replace("SPAdes genome assembler ", "")
    elif cmd_exec == "mlst":
        version = result[0].replace("mlst ", "v")
        if pkg_resources.parse_version(version) < pkg_resources.parse_version("2.19"):
            logging.warning("MLST version is older than v2.19. Database path handling may differ in older versions.")
    elif cmd_exec == 'minimap2':
        version = "v" + result[0]
    elif cmd_exec == 'prokka':
        version = result[0].replace("prokka ", "v")
    if cmd_exec == "samtools":
            version = result[0].replace("samtools ", "")
            if pkg_resources.parse_version(version) < pkg_resources.parse_version(
                "1.10"
            ):
                logging.critical("Samtools version too old, please upgrade to v1.10.0+")
                sys.exit(1)
    if cmd_exec == "bcftools":
        version = result[0].replace("bcftools ", "")
    if cmd_exec == "kallisto":
        version = result[0].replace("kallisto, version ", "v")
    if cmd_path is not None:
        msg = "Located " + cmd_exec + " " + version + " in " + cmd_path
        logging.info(msg)
    else:
        msg = cmd_exec + " was not found, please check installation on your device"
        logging.critical(msg)
        sys.exit(1)

def check_mlst(datadir):
    if datadir is not None:
        result = subprocess.run(
            ["mlst", "--longlist", "--datadir", datadir],
            capture_output=True,
            text=True,
        )
    else:
        result = subprocess.run(
            "mlst --longlist | grep bpertussis",
            capture_output=True,
            text=True,
            shell=True,
        )
    if result.returncode > 0:
        logging.critical("MLST database is not prepared. The 'bpertussis' scheme was not found.")
        logging.critical(
            "Please run setup.py first: python setup.py [optional datadir]"
        )
        sys.exit(1)
    elif "bpertussis" not in result.stdout:
        logging.critical("MLST 'bpertussis' scheme not found in available schemes.")
        logging.critical(
            "Please run setup.py first: python setup.py [optional datadir]"
        )
        sys.exit(1)
    else:
        logging.info("MLST 'bpertussis' scheme is available.")

def check_abricate():
    result = subprocess.run(
        ["abricate", "--list", "--datadir", bor_vfdb_db ],
        capture_output=True,
        text=True,
    )
    if result.returncode > 0:
        logging.critical("Abricate database is not prepared")
        logging.critical(
            "correct by running:   abricate --setupdb --datadir " + bor_vfdb_db 
        )
        sys.exit(1)
    dbs = [x.split("\t")[0] for x in result.stdout.splitlines()[1:]]
    if any(x not in dbs for x in ["bp-only_vfdb"]):
        logging.critical("unable to find databases")
        sys.exit(1)

def check_spades_finished(spades_outdir):
    result = "======= SPAdes pipeline finished."
    result2 = "======= SPAdes pipeline finished WITH WARNINGS!"
    spades_log = os.path.join(spades_outdir, "spades.log")
    
    if os.path.isfile(spades_log) and os.stat(spades_log).st_size != 0:
        with open(spades_log, 'r') as log:
            for line in log:
                if result in line or result2 in line:
                    return True
        return False
    else:
        return False
    
def check_megahit_finished(megahit_outdir):
    result = "ALL DONE."
    megahit_log = os.path.join(megahit_outdir, "log")
    
    if os.path.isfile(megahit_log) and os.stat(megahit_log).st_size != 0:
        with open(megahit_log, 'r') as log:
            for line in log:
                if result in line:
                    return True
        return False
    else:
        return False
    
def check_prokka_finished(prokka_outdir, name):
    result = "Annotation finished successfully."
    prokka_log = name + ".log"
    prokka_log = os.path.join(prokka_outdir, prokka_log)
    
    if os.path.isfile(prokka_log) and os.stat(prokka_log).st_size != 0:
        with open(prokka_log, 'r') as log:
            for line in log:
                if result in line:
                    return True
        return False
    else:
        return False

def check_kallisto_finished(analysis_outdir):
    kallisto_abundance = os.path.join(analysis_outdir, "abundance.tsv")
    if os.path.isfile(kallisto_abundance) and os.stat(kallisto_abundance).st_size != 0:
        return True
    return False


def check_bam_readable(bam_file):
    if not os.path.isfile(bam_file) or os.stat(bam_file).st_size == 0:
        return False
    result = subprocess.run(
        ["samtools", "view", "-H", bam_file],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return result.returncode == 0

def get_fasta_length(prn_type):
    length = 0
    fasta_file = os.path.join(os.path.dirname(os.path.dirname(__file__)), "databases/bpertussis/prn.tfa")

    with open(fasta_file, 'r') as file:
        found = False
        contig_name = ''
        sequence_length = 0
        
        for line in file:
            if line.startswith('>'):
                if found:
                    break  # Stop if we have already found and processed the target contig
                contig_name = line[1:].strip()
                if contig_name == prn_type:
                    found = True
                    sequence_length = 0
            elif found:
                sequence_length += len(line.strip())
    return sequence_length

def check_closed_genome(fasta, length_threshold=3900000):
    contig_count = 0
    total_length = 0
    with open(fasta, 'r') as file:
        for line in file:
            if line.startswith('>'):
                contig_count += 1
            else:
                total_length += len(line.strip())
    is_length_above_threshold = total_length > length_threshold
    if contig_count == 1 and is_length_above_threshold:
        return True
    else:
        return False
    
def megahit_assembly_graphs(megahit_outdir):
    megahit_done = check_megahit_finished(megahit_outdir)
    if megahit_done is True:
        int_contigs = megahit_outdir + "/intermediate_contigs"
        fa_list = [f for f in os.listdir(int_contigs) 
           if f.endswith(".contigs.fa") 
           and not f.endswith(".final.contigs.fa") 
           and os.path.isfile(os.path.join(int_contigs, f))]
        for fa in fa_list:
            kmer_size = fa.removeprefix('k').strip('.contigs.fa')
            megahit_toolkit_cmd = f"megahit_toolkit contig2fastg {kmer_size} {fa} > {int_contigs}/k{kmer_size}.fastg"
            fastg2gfa_cmd = f"gfatools view {int_contigs}/k{kmer_size}.fastg > {int_contigs}/k{kmer_size}.gfa"
            run_cmd(megahit_toolkit_cmd)
            run_cmd(fastg2gfa_cmd)

def contig_prokka_tag(assembly, name, prokka_outdir):
    # Paths to the GBK and FASTA files
    prokka_gbk = prokka_outdir + "/" + name + ".gbk"
    fasta_file = assembly

    # Parse the GBK file to extract locus and size
    gbk_data = []
    with open(prokka_gbk) as handle:
        for record in GenBank.parse(handle):
            locus = record.locus
            size = int(record.size)  # Extract size and convert to integer
            gbk_data.append((locus, size))

    # Parse the FASTA file to extract headers and sequences
    min_length = 200
    fasta_data = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        header = record.id  # Extract the header (without the ">" sign)
        sequence = str(record.seq)  # Convert the sequence to a string
        if len(sequence) >= min_length:
            fasta_data.append((header, len(sequence)))

    # Check that the number of entries matches
    if len(gbk_data) != len(fasta_data):
        print("Error: Number of loci in GBK file does not match number of sequences in FASTA file.")
        return None

    # Verify that the sequence lengths match
    locus_to_sequence = {}
    for (locus, expected_size), (header, actual_size) in zip(gbk_data, fasta_data):
        if expected_size != actual_size:
            print(f"Warning: Size mismatch for locus {locus}. Expected {expected_size}, found {actual_size}.")
        else:
            locus_to_sequence[locus] = header

    return prokka_gbk, locus_to_sequence
        
