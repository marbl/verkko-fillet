import warnings
# Ensure Biopython is imported correctly
try:
    from Bio import BiopythonWarning
except ImportError:
    print("Biopython is not installed. Please install it using 'pip install biopython'.")
    sys.exit(1)

warnings.simplefilter('ignore', BiopythonWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="Bio")
warnings.filterwarnings('ignore')

with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    # Code that might trigger warnings
    from Bio import SeqIO
    # Further operations

import os
import subprocess
import pandas as pd
from tqdm import tqdm
import sys
# from Bio import SeqIO
import re
import networkx as nx
from collections import Counter
import shlex
import shutil

from .._run_shell import run_shell
script_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../bin/'))

# Custom sort function that prioritizes base entries before random ones
def sort_by_random_chr_hap(item, by="hap", type_list = ['mat', 'pat', 'hapUn']):
    """/
    Sorts chromosome names based on a custom sorting criterion.

    Parameters
    ----------
    item
        The chromosome name to be sorted.
    by
        The sorting criterion. Default is 'hap'.
    type_list
        The list of chromosome types to be used for sorting. Default is ['mat', 'pat', 'hapUn'].
    """
    # Check if 'random' is in the item
    is_random = '_random_' in item

    # Extract parts: base 'chrX' part, type ('mat' or 'pat'), and random suffix (if any)
    match = re.match(r'(chr\d+|chrUn|chrX|chrY|chrM|chr[A-Za-z]+)_([A-Za-z]+)(_\d+)?(_random_[A-Za-z0-9-]+)?', item)
    if match:
        chr_part = match.group(1)
        type_part = match.group(2)
        subtype_part = match.group(3) if match.group(3) else ''  # The numeric part, like _1, _2, etc.
        random_part = match.group(4) if match.group(4) else ''

        type_priority = {item: index + 1 for index, item in enumerate(type_list)}.get(type_part, 4)
        
        # Extract chromosome number as an integer for proper sorting
        if not re.search(r'\d+', chr_part):
            chr_priority = 9999  # Place chrX, chrY, and chrM after the numeric chromosomes
        else:
            chr_number_match = re.search(r'\d+', chr_part)
            chr_priority = int(chr_number_match.group(0)) if chr_number_match else 0
        
        # Extract the numeric part from the subtype (e.g., '_1', '_2', etc.)
        subtype_number = int(subtype_part[1:]) if subtype_part else 0
        
        # Return a tuple with:
        # 1. Chromosome number (numeric chromosomes come first)
        # 2. Whether it's random or not (to prioritize non-random first)
        # 3. Type ('mat' or 'pat')
        # 4. Subtype number (to ensure correct order within 'mat' and 'pat')
        # 5. Random part (to ensure random entries are last)
        if by == "chr":
            return (is_random, chr_priority, type_priority, subtype_number, random_part)
        elif by == "hap":
            return (is_random, type_priority, chr_priority, subtype_number, random_part)
    # If no match, return a tuple that won't interfere with other comparisons
    return (False, '', 0, 0, '')  # Default tuple to handle unmatched items

def sortContig(ori_fasta, sorted_fasta=None, sort_by="hap"):
    """
    Sorts sequences in a FASTA file based on a custom sorting criterion (e.g., 'hap', 'chr').
    
    Parameters
    ----------
    ori_fasta
        Path to the original FASTA file.
    sorted_fasta
        Path to save the sorted FASTA file. If None, it will be generated with surfix of "_sorted.fasta"
    sort_by
        Sorting criteria (default is "hap"). ['hap','chr']
    """
    
    # Check if the input FASTA file exists
    if not os.path.exists(ori_fasta):
        print(f"The input FASTA file does not exist: {ori_fasta}")
        return
    
    # Extract basename and remove extensions like .fasta, .fasta.gz, .fa, .fa.gz
    basename = os.path.basename(ori_fasta)
    basename = re.sub(r'\.fasta(\.gz)?$|\.fa(\.gz)?$', '', basename)
    
    # Generate the sorted output filename if not provided
    if sorted_fasta is None:
        sorted_fasta = basename + "_sorted" + sort_by + ".fasta"

    if os.path.exists(sorted_fasta):
        print(f"The sorted FASTA file already exists: {sorted_fasta}")
        return
        
    # Parse sequences from the original FASTA file
    sequences = list(SeqIO.parse(ori_fasta, "fasta"))
    
    # Extract sequence IDs
    sequence_ids = [record.id for record in sequences]
    print(f"Sorting {len(sequence_ids)} sequences based on the custom sorting criterion...")
    # Sorting the sequence IDs based on the custom function
    sorted_data = sorted(sequence_ids, key=lambda item: sort_by_random_chr_hap(item, by=sort_by))
    
    # Reordering sequences based on sorted sequence IDs
    sorted_sequences = [record for id in sorted_data for record in sequences if record.id == id]
    
    # Write the sorted sequences to a new FASTA file
    SeqIO.write(sorted_sequences, sorted_fasta, "fasta")
    print(f"Sorted sequences have been written to {sorted_fasta}")



def renameContig(obj, 
                 chrMap, 
                 out_mapFile = "assembly.final.mapNaming.txt", 
                 original_fasta= "assembly.fasta", 
                 output_fasta = None, showOnly = False):
    """\
    Rename the contigs in the FASTA file based on the provided chromosome map file.

    Parameters
    ----------
    obj
        The VerkkoFillet object to be used.
    chrMap
        The DataFrame containing the mapping of old chromosome names to new chromosome names.
    out_mapFile
        The output file to save the chromosome map. Default is "assembly.final.mapNaming.txt".
    original_fasta
        The path to the original FASTA file. Default is "assembly.fasta".
    output_fasta
        The path to save the renamed FASTA file. If None, it will be generated with a suffix of "_rename.fasta".
    showOnly
        If True, the command will be printed but not executed. Default is False.

    Returns
    -------
        output_fasta
    """
    print(f"Starting renaming contigs in the {original_fasta} file ...")
    print(" ")

    print("Checking the required files ...")
    print("   - Checking the chromosome map file ...")
    print("   - Checking the original fasta file ..." )
    print(" ")
    working_dir = os.path.abspath(obj.verkko_fillet_dir)
    script = os.path.abspath(os.path.join(script_path, "changeChrName.sh"))  # Assuming script_path is defined elsewhere

    # Check if there is no duplications in the new chromosome names
    if chrMap['contig'].duplicated().any():
        print("Error: There are duplicated chromosome names in the new chromosome names. Please check the chromosome map file.")
        return

    if chrMap['new_contig_name'].duplicated().any():
        print("Error: There are duplicated chromosome names in the new chromosome names. Please check the chromosome map file.")
        return
    
    chrMap=chrMap.merge(obj.scfmap, on = 'contig')
    chrMap.to_csv(out_mapFile, sep ='\t', header = None, index=False)
    
    if output_fasta is None:  # Use 'is None' for comparison
        prefix = re.sub(r"\.gz|\.fasta", "", original_fasta)
        outFasta = prefix + "_rename.fasta"
    
    # Check if the script exists
    if not os.path.exists(script):
        print(f"Script not found: {script}")
        return
    
    # Check if the working directory exists
    if not os.path.exists(working_dir):
        print(f"Working directory not found: {working_dir}")
        return
        
    if not os.path.exists(out_mapFile):
        print(f"chromosome map file not found : {out_mapFile}")
        return
        
    # Construct the shell command
    cmd = f"sh {shlex.quote(script)} {shlex.quote(out_mapFile)} {shlex.quote(str(original_fasta))} {shlex.quote(outFasta)}"
    
    run_shell(cmd, wkDir=working_dir, functionName = "chrRename" ,longLog = False, showOnly = showOnly)
    print("The contig renaming was completed successfully!")
    print(f"Final renamed fasta file : {outFasta}")


def flipContig(filp_contig_list, ori_fasta="assembly.fasta", final_fasta=None):
    """\
    Flip the sequences in a FASTA file based on the provided list of contigs.

    Parameters
    ----------
    filp_contig_list
        The list of contigs to be flipped.
    ori_fasta
        The path to the original FASTA file. Default is "assembly.fasta".
    final_fasta
        The path to save the flipped FASTA file. If None, it will be generated with a suffix of "_flip.fasta".
    
    Returns
    -------
        final_fasta
    """
    
    # Check if "assembly_trimmed.fasta" exists and update file names
    # Check if the final output file already exists
    if os.path.exists(final_fasta):
        print(f"{final_fasta} already exists. Exiting to avoid overwriting.")
        sys.exit(1)
    
    # Load chromosome list from FASTA index file
    fai_path = f"{ori_fasta}.fai"
    if not os.path.exists(fai_path):
        print(f"Index file {fai_path} not found. Generating faidx index for {ori_fasta}.")
        cmd = f"samtools faidx {ori_fasta}"
        try:
            subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error generating index file: {cmd}")
            print(f"Error message: {e}")
            sys.exit(1)
    
    fai = pd.read_csv(fai_path, sep='\t', header=None, usecols=[0])
    chrList = list(fai[0])
    
    # Process each chromosome with a progress bar
    with tqdm(total=len(chrList), desc="Flipping Chromosomes", ncols=80, colour="white") as pbar:
        for chromosome in chrList:
            if chromosome in filp_contig_list:
                cmd = f"samtools faidx {ori_fasta} {chromosome} | seqtk seq -r >> {final_fasta}"
            else:
                cmd = f"samtools faidx {ori_fasta} {chromosome} | seqtk seq >> {final_fasta}"
            
            try:
                subprocess.run(cmd, shell=True, check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error processing chromosome {chromosome}: {cmd}")
                print(f"Error message: {e}")
                sys.exit(1)
            
            # Update progress bar
            pbar.update(1)
    
    print("The chromosome flipping was completed successfully!")
    print(f"Output FASTA: {final_fasta}")

def filterContigs(mapfile, assembly, out_prefix=None, filter_chr_list=None, showOnly = False):
    """
    Filter the contigs in the FASTA file based on the provided list of contigs. For chromosome assignment, we recommend using the reference genome that contains only the chromosomes to which the contigs should be assigned.

    Parameters
    ----------
    mapfile
        The path to the map file. The map file should contain the list of contigs to be filtered.
    assembly
        The path to the original FASTA file.
    out_prefix
        The prefix for the output file. If None, it will be generated based on the input file name with surfixed "_filtered.fa".
    filter_chr_list
        The list of contigs to be filtered.
    showOnly
        If True, the command will be printed but not executed. Default is False.
    
    Returns
    -------
        fasta file with filtered contigs
    """
    # check if samtools is installed
    try:
        subprocess.run("samtools --version", shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"Error checking samtools installation: {e}")
        sys.exit(1)

    # Check if the map file exists
    if not os.path.exists(mapfile):
        print(f"Map file not found: {mapfile}")
        return
    
    # Check if the assembly file exists
    if not os.path.exists(assembly):
        print(f"Assembly file not found: {assembly}")
        return

    if not os.path.exists(f"{assembly}.fai"):
        print(f"FAI for Assembly file not found: {assembly}.fai")
        return

    fai = pd.read_csv(f"{assembly}.fai", sep='\t', header=None, usecols=[0])
    faichrList = list(fai[0])

    if out_prefix is None:
        out_basename = os.path.splitext(os.path.basename(assembly))[0] + "_filtered"
        out_dir = os.path.dirname(assembly)
        out_prefix = os.path.join(out_dir, out_basename) 
    
    if filter_chr_list is None:
        cmd=f"cut -f1 {mapfile}"
        col1 = subprocess.check_output(cmd, shell=True, text=True).splitlines()

        cmd=f"cut -f2 {mapfile}"
        col2 = subprocess.check_output(cmd, shell=True, text=True).splitlines()

        filter_chr_list = list(set(col1 + col2))
        filter_chr_list_len = int(len(filter_chr_list)/2)
        print("No filter chromosome list provided. The contigs in the map file will be used.")
        print(f"total chromosomes will be filtered in : {len(filter_chr_list)/2}")
    else:
        filter_chr_list_len = len(filter_chr_list)
        print("Filtering contigs based on the provided list.")
        print(f"total chromosomes will be filtered in : {len(filter_chr_list)}")

    # Check if the output file already exists
    if os.path.exists(f"{out_prefix}.fa"):
        print(f"Output file already exists: {out_prefix}.fa")
        return

    # intersect the filter_chr_list with the chrList
    intersect_contig = list(set(filter_chr_list) & set(faichrList))
    if len(intersect_contig) == 0:
        print("No contigs to filter. Theres no contigs are interected with the assembly.fai and given list.")
        return
    if filter_chr_list_len - len(intersect_contig) > 0:
        print(f"{len(filter_chr_list) - len(intersect_contig)} contigs are not found in the assembly.fai file.")
        print(f"Please check the contig names in the map file and the assembly file.")
    
    filter_chr_list = intersect_contig

    print(f"Filtering contigs based on {len(filter_chr_list)} chromosomes.")
    print(f"total chromosomes will be filtered in : {len(filter_chr_list)}")
    filter_chr_list = " ".join(filter_chr_list)

    # Construct the shell command
    cmd = f"samtools faidx {shlex.quote(assembly)} {filter_chr_list}> {shlex.quote(out_prefix)}.fa"
    
    run_shell(cmd, functionName = "filterContigs", wkDir = os.getcwd() ,longLog = False, showOnly = showOnly)
    
    print(f"Filtered FASTA: {out_prefix}.fa")