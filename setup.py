import os
import re
import logging
import subprocess
import sys
import shutil
import hashlib
import gzip
from scripts import assists

# copy bpertussis from databases to mlst db folder
logging.getLogger().setLevel(logging.INFO)

def pertpipe_setup(datadir):
    dependency_list = ["abricate", "spades.py", "mlst", "minimap2", "samtools", "bcftools", "kallisto"]
    bpertussis_db = os.path.join(os.path.dirname(os.path.abspath(__file__)), "databases/bpertussis")

    logging.info("Checking installs of dependencies")
    for dependency in dependency_list:
        assists.check_dependencies(dependency)
    if "abricate" in dependency_list:
        assists.check_abricate()

    if datadir is not None:
        use_datadir = True
        logging.info(f"Datadir is {datadir}, modifying mlst-make_blast_db to suit")
        default_mlst_path = os.path.join(datadir, "pubmlst")
    else:
        logging.info(f"Datadir is not set, using default MLST location")
        result = subprocess.run('mlst -h', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
        if result.returncode == 0:
            # mlst writes help to stderr, not stdout
            match = re.search(r'--datadir\s+\[X\]\s+PubMLST data\s+\(default\s+\'([^\']+)\'\)', result.stderr)
            if match is None:
                logging.error("Could not parse mlst --datadir path from help output. Try running: mlst -h")
                sys.exit(1)
            default_mlst_path = match.group(1)
        else:
            logging.error(f"Error finding path with 'which': {result.stderr}")
            sys.exit(1)
        datadir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "databases")
        use_datadir = False

    if not os.path.isdir(default_mlst_path):
        logging.error(f"Error finding path with 'which': {default_mlst_path}")
        sys.exit(2)

    new_destination = default_mlst_path + "/bpertussis"
    # Always copy so a previously failed/partial setup doesn't leave a blocking empty dir.
    # dirs_exist_ok=True makes this safe to call repeatedly.
    logging.info(f"Copying bpertussis database to {new_destination}")
    os.makedirs(new_destination, exist_ok=True)
    try:
        shutil.copytree(bpertussis_db, new_destination, dirs_exist_ok=True)
        logging.info(f"Copied {bpertussis_db} to {new_destination}")
    except PermissionError:
        logging.error(f"Permission denied copying {bpertussis_db} to {new_destination}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Error copying {bpertussis_db} to {new_destination}: {e}")
        sys.exit(1)

    # run the makemlstdb thingo
    bpertussis_mlst = None
    if os.path.isdir(new_destination):
        file_list = [
            "23SrRNA.tfa",
            "bpertussis.txt",
            "fhaB24005550.tfa",
            "fim2.tfa",
            "fim3.tfa",
            "prn.tfa",
            "ptxA.tfa",
            "ptxB.tfa",
            "ptxC.tfa",
            "ptxD.tfa",
            "ptxE.tfa",
            "ptxP.tfa",
        ]
        for file in file_list:
            if os.path.exists(os.path.join(new_destination, file)):
                logging.info(f"{file} exists.")
            else:
                logging.error(f"{file} does not exist! check if it is missing in the database folder {bpertussis_db}")
                sys.exit(3)

    cmd_list = [
        f"mlst-make_blast_db",
        f"mlst --longlist | grep bpertussis"
    ]
    expected_result = 'bpertussis\tBPagST\tptxP\tptxA\tptxB\tptxC\tptxD\tptxE\tfhaB24005550\tfim2\tfim3\t23SrRNA\n'
    check_existing = subprocess.run(cmd_list[1], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
    if check_existing.returncode != 0 or "bpertussis" not in check_existing.stdout:
        logging.info("bpertussis scheme not found, building BLAST database...")
        # Step 1: build the BLAST database
        makedb_result = subprocess.run(cmd_list[0], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
        if makedb_result.returncode == 0:
            logging.info(f"Successfully ran {cmd_list[0]}")
        else:
            logging.error(f"mlst-make_blast_db failed: {makedb_result.stderr}")
            logging.error("Ensure mlst is installed correctly and the bpertussis folder was copied to the correct location.")
            sys.exit(1)
        # Step 2: verify the scheme is now available
        verify_result = subprocess.run(cmd_list[1], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
        if verify_result.returncode == 0 and "bpertussis" in verify_result.stdout:
            logging.info("Successfully created 'bpertussis' MLST scheme")
        else:
            logging.critical(f"MLST scheme creation failed. Expected 'bpertussis' in scheme list but got: {verify_result.stdout}")
            logging.critical("Try manually running: mlst-make_blast_db && mlst --longlist | grep bpertussis")
            sys.exit(1)
    else:
        logging.info(f"MLST scheme already exists")

    # test MLST cos having issues with ptxA being ~33 not 1.
    pert_test = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tests/CIDM-MRBP01.fna")
    if use_datadir is False:
        cmd = f"mlst --scheme bpertussis {pert_test}"
    else:
        cmd = f"mlst --scheme bpertussis --datadir {datadir}/pubmlst --blastdb {datadir}/blast/mlst.fa {pert_test}"
    mlst_result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
    if mlst_result.returncode == 0:
        logging.info(f"Successfully run {cmd}")
    else:
        logging.error(f"{mlst_result.stderr}")
        sys.exit()
    bpertussis_mlst = mlst_result.stdout.split("\t")[5]
    expected_result = "ptxA(1)"
    bad_result = ["ptxA(~33)", "ptxA(-)"]
    if bpertussis_mlst == expected_result:
        logging.info(f"MLST scheme test success")
    elif bpertussis_mlst in bad_result:
        logging.info(f"MLST scheme test failed. Need to re-copy/re-download database and re-do setup.")
    else:
        logging.critical(f"MLST scheme test failed. Need to re-copy/re-download database and re-do.")

    # set up abricate automatically.

    # silva lsu download
    silva_lsu = "https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports/SILVA_138.2_LSURef_NR99_tax_silva.fasta.gz"
    silva_md5 = "https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports/SILVA_138.2_LSURef_NR99_tax_silva.fasta.gz.md5"
    dest_lsu = os.path.join(datadir, "SILVA_138.2_LSURef_NR99_tax_silva.fasta.gz")
    dest_md5 = os.path.join(datadir, "SILVA_138.2_LSURef_NR99_tax_silva.fasta.gz.md5")

    # Use subprocess to call rsync and download the file
    if os.path.exists(dest_lsu) is False or os.path.exists(dest_lsu.strip(".gz")) is False:
        logging.info(f"SILVA LSU Data does not exist, downloading now.")
        curl_command_1 = ['curl', '-o', dest_lsu, silva_lsu]
        curl_command_2 = ['curl', '-o', dest_md5, silva_md5]
        try:
            subprocess.run(curl_command_1, check=True)
            logging.info(f'Download complete: {dest_lsu}')
            subprocess.run(curl_command_2, check=True)
            logging.info(f'Download complete: {dest_md5}')
        except subprocess.CalledProcessError as e:
            logging.error(f"Error occurred while downloading the file: {e}. I suggest downloading the file manually and placing it in {dest_lsu}")
            raise
    else:
        logging.info(f"SILVA LSU Data already exists, skipping.")

    #Calculate MD5 checksum of a file.
    try:
        with open(dest_md5, 'r') as f:
            expected_md5 = f.read().strip().split()[0]

        hash_md5 = hashlib.md5()
        with open(dest_lsu, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        actual_md5 = hash_md5.hexdigest()

        # Compare checksums
        if expected_md5 == actual_md5:
            print(f'Checksum verified successfully: {actual_md5}')
        else:
            print(f'Checksum mismatch! Expected: {expected_md5}, but got: {actual_md5}')
    except Exception as e:
        print(f"Error occurred during checksum verification: {e}")
        raise

    # kallisto setup
    lsu_data = dest_lsu.strip('.gz')
    with gzip.open(dest_lsu, 'rb') as f_in:
        with open(lsu_data, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    index_name = f"{datadir}/SILVA_138.2_LSURef_NR99"
    if os.path.exists(index_name) is False:
        logging.info("Creating Kallisto index")
        kallisto_cmd = f"kallisto index -i {index_name} {lsu_data}"
        result = subprocess.run(kallisto_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
    else:
        logging.info("Kallisto Index already exists.")

    logging.info(f"Set-up complete :)")

datadir = None
if len(sys.argv) > 1:
    datadir = sys.argv[1]

pertpipe_setup(datadir)
