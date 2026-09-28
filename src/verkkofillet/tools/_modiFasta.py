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
dataset_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../data/dataset/'))

def screen_asm(obj, threads=10, minLen=1000, ebv_fasta=None, rdna_fasta=None, mt_fasta=None, showOnly=False):
    """
    Screen the assembly for potential issues using a shell script.

    Parameters
    ----------
    obj : object
        An object containing the Verkko directory information.
    threads : int, optional
        Number of threads to use. Default is 10.
    minLen : int, optional
        Minimum contig length for screening. Default is 1000.
    ebv_fasta : str, optional
        Path to the EBV contaminant FASTA file. Default is None.
    rdna_fasta : str, optional
        Path to the rDNA contaminant FASTA file. Default is None.
    mt_fasta : str, optional
        Path to the MT contaminant FASTA file. Default is None.
    showOnly : bool, optional
        If True, only show the command without executing it. Default is False.
    """
    print(f"Starting assembly screening in {obj.verkkoDir} ...")
    script = os.path.join(script_path, "screen-assembly.sh")
    
    if not os.path.exists(script):
        print(f"Script not found: {script}")
        return
    
    cmd = f"bash {script} {obj.verkkoDir} {threads} {minLen} {ebv_fasta} {rdna_fasta} {mt_fasta}"
    
    run_shell(cmd, wkDir=obj.verkkoDir, functionName="screen_asm", longLog=False, showOnly=showOnly)


def fix_fasta_rDNA(fasta = "assembly_trimmed_flipped_rename_sortedhap.fasta", 
                      gap_rDNA_info = "assembly_trimmed_flipped_rename_sortedhap.gaps.bed.rDNA.bounded.csv", 
                      out_fasta = None, gap_size=100000, force=False):
    """
    Clean gaps in rDNA regions of a given FASTA file.

    Parameters
    ----------
    fasta : str
        Path to the input FASTA file.
    gap_rDNA_info : str
        Path to the rDNA gap information file.
    gap_size : int, optional
        Size of the gap to be cleaned. Default is 100000.
    force : bool, optional
        If True, overwrite existing output files. Default is False.

    Returns
    -------
    None
        The function executes a shell command to perform the cleaning.
    """
    
    script = os.path.join(script_path, "rDNA_gap_cleaning.sh")
    if out_fasta is None:
        prefix = re.sub(r"(\.fasta|\.fa)(\.gz)?$", "", os.path.basename(fasta), flags=re.IGNORECASE)
        out_fasta = prefix + "_rDNA_gap_cleaned.fasta"

    if os.path.exists(f"{out_fasta}.gz") and not force:
        print(f"Output FASTA file {out_fasta}.gz already exists. Use force=True to overwrite.")
        return
    
    if not os.path.exists(fasta):
        print(f"Input FASTA file not found: {fasta}")
        return
    
    if not os.path.exists(script):
        print(f"Script not found: {script}")
        return
    
    print(f"Cleaning gaps in rDNA regions for {fasta} with gap size of {gap_size} ...")
    print(f"Output FASTA will be saved to: {out_fasta}")
    cmd = f"bash {script} {gap_rDNA_info} {fasta}  {out_fasta} {gap_size} {force}"
    
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {cmd}")
        print(f"Error message: {e}")
        sys.exit(1)



def find_gap_in_rDNA(rDNA_mapping = "rDNA_mapping/45S.mashmap.bed", 
                     fasta = "assembly_trimmed_flipped_rename_sortedhap.fasta",
                     padding = 100_000,
                     idx = 90):
    if not os.path.exists(rDNA_mapping):
        print(f"rDNA mapping file not found: {rDNA_mapping}")
        print(f"Please run the map_rDNA function first to generate the rDNA mapping file.")
        return
    
    if not os.path.exists(fasta):
        print(f"Assembly FASTA file not found: {fasta}")
        return
    
    print(f"Finding gaps in rDNA mapping for {fasta} with padding of {padding} bp and idx threshold of {idx} ...")

    fasta_name = re.sub(r"(\.fasta|\.fa)(\.gz)?$", "", os.path.basename(fasta), flags=re.IGNORECASE)

    rDNA_map = pd.read_csv(rDNA_mapping, header = None, names = ['contig','start','end','name','idx','strand'], sep='\t')
    rDNA_map = rDNA_map.loc[rDNA_map['idx'] > idx]

    gaps = pd.read_csv(f"{fasta_name}.gaps.bed", header = None, names = ['contig','start','end'], sep='\t')

    gaps['bounded_by_rDNA'] = False
    for index, gap in gaps.iterrows():
        contig = gap["contig"]
        gap_start = int(gap["start"])
        gap_end = int(gap["end"])

        rDNA_entries = rDNA_map.loc[rDNA_map["contig"] == contig].copy()
        if rDNA_entries.empty:
            print(f"Gap {index} on {contig} ({gap_start}-{gap_end}) has no rDNA entries.")
            continue

        rDNA_entries["rDNA_min"] = rDNA_entries[["start", "end"]].min(axis=1)
        rDNA_entries["rDNA_max"] = rDNA_entries[["start", "end"]].max(axis=1)

        # Left flank: rDNA ends before gap start, but not farther than padding.
        left_flank = (rDNA_entries["rDNA_max"] <= gap_start) & (rDNA_entries["rDNA_max"] >= gap_start - padding)
        # Right flank: rDNA starts after gap end, but not farther than padding.
        right_flank = (rDNA_entries["rDNA_min"] >= gap_end) & (rDNA_entries["rDNA_min"] <= gap_end + padding)

        bounded = bool(left_flank.any() and right_flank.any())

        if bounded:
            print(f"Gap {index} on {contig} ({gap_start}-{gap_end}) is bounded by rDNA within {padding} bp.")
            gaps.loc[index, "bounded_by_rDNA"] = True
        else:
            print(f"Gap {index} on {contig} ({gap_start}-{gap_end}) is NOT bounded by rDNA within {padding} bp.")
            gaps.loc[index, "bounded_by_rDNA"] = False

    gaps.to_csv(f"{fasta_name}.gaps.bed.rDNA.bounded.csv", index=False, sep='\t')
    print(f"Gap analysis completed. Results saved to {fasta_name}.gaps.bed.rDNA.bounded.csv")


def map_rDNA(rDNA = None, cpu=10, fasta = "assembly_trimmed_flipped_rename_sortedhap.fasta"):
    """
    Map rDNA sequences to the assembly using a shell script.

    Parameters
    ----------
    rDNA : str
        Path to the rDNA reference sequence file.
    cpu : int
        Number of CPU cores to use for the mapping.
    fasta : str
        Path to the assembly FASTA file.

    Returns
    -------
    None
        The function executes a shell command to perform the mapping.
    """
    fasta_name = re.sub(r"(\.fasta|\.fa)(\.gz)?$", "", os.path.basename(fasta), flags=re.IGNORECASE)

    if not os.path.exists(fasta):
        print(f"Assembly FASTA file not found: {fasta}")
        return
    print(f"fasta input : {fasta}")

    script=f"{script_path}/45S_mapping.sh"
    if rDNA is None:
        rDNA=f"{dataset_path}/NT_167214.1.fa"

    print(f"rDNA reference sequence: {rDNA}")

    print(f"Align rDNA on fasta file ... ")
    if os.path.exists(f"rDNA_mapping/45S.mashmap.bed"):
        print(f"rDNA mapping file already exists: rDNA_mapping/45S.mashmap.bed")
    else:
        cmd = f"bash {script} {rDNA} {cpu} {fasta}"
        try:    
            subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error generating index file: {cmd}")
            print(f"Error message: {e}")
            sys.exit(1) 
    
    # find gap
    script=f"{script_path}/getT2T.sh"
    
    if os.path.exists(f"{fasta_name}.gaps.bed"):
        print(f"Gaps file already exists: {fasta_name}.gaps.bed")
    else:
        cmd = f"bash {script} {fasta}"
        try:    
            subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error generating index file: {cmd}")
            print(f"Error message: {e}")
            sys.exit(1)
    


def find_flip_candidates(obj, mashmap_file = "chromosome_assignment/assembly.mashmap.out"):
    """
    Identify contigs that may need to be flipped based on MashMap output and the Verkko object statistics.

    Parameters
    ----------
    obj : Verkko object
        The Verkko object containing assembly statistics.
    mashmap_file : str, optional
        Path to the MashMap output file. Default is "chromosome_assignment/assembly.mashmap.out".

    Returns
    -------
    list
        A list of contig names that are candidates for flipping.
    """

    df_mashmap = pd.read_csv(
        mashmap_file,
        sep="\t",
        header=None,
        usecols=[0, 4, 5, 8],
        names=["contig", "strand", "ref_name", "alignment_block"],
    )

    df_mashmap["alignment_block"] = pd.to_numeric(df_mashmap["alignment_block"], errors="coerce")

    df_flip_candidates = (
        df_mashmap
        .groupby(["contig", "strand", "ref_name"], as_index=False)["alignment_block"]
        .sum()
        .sort_values(["contig", "alignment_block"], ascending=[True, False])
        .drop_duplicates(subset=["contig"], keep="first")
        .reset_index(drop=True)
    )

    df_flip_candidates_list = df_flip_candidates.loc[df_flip_candidates["strand"] == "-"].reset_index(drop=True)

    # Count all negative-strand top candidates and the subset present in obj.stats.
    total_flip_num = len(df_flip_candidates_list)
    main_contig_mask = df_flip_candidates_list["contig"].isin(obj.stats["contig"])

    total_flip_num_main = int(main_contig_mask.sum())
    flip_contig_list = df_flip_candidates_list.loc[main_contig_mask, "contig"].tolist()
    print(f"Total negative-strand top candidates: {total_flip_num}")
    print(f"Total negative-strand top candidates in main contigs: {total_flip_num_main}")

    return flip_contig_list

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
    cmd = f"bash {shlex.quote(script)} {shlex.quote(out_mapFile)} {shlex.quote(str(original_fasta))} {shlex.quote(outFasta)}"
    
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


def rDNA_gap_cleaning(fasta = "assembly_trimmed_flipped_rename_sortedhap.fasta", 
                      gap_rDNA_info = "assembly_trimmed_flipped_rename_sortedhap.gaps.bed.rDNA.bounded.csv", 
                      out_fasta = None, gap_size=100000, padding = 100_000, idx = 90, threads=10,
                      force=False):
    """
    Clean gaps in rDNA regions of a given FASTA file.

    Parameters
    ----------
    fasta : str
        Path to the input FASTA file.
    gap_rDNA_info : str
        Path to the rDNA gap information file.
    gap_size : int, optional
        Size of the gap to be cleaned. Default is 100000.
    padding : int, optional
        Padding around the gap to be considered. Default is 100000.
    idx : int, optional
        Index parameter for gap cleaning. Default is 90.
    threads : int, optional
        Number of threads to use. Default is 10.
    force : bool, optional
        If True, overwrite existing output files. Default is False.

    Returns
    -------
    None
        The function executes a shell command to perform the cleaning.
    """
    ## align 45S to reference
    map_rDNA(rDNA = None, cpu=threads, fasta = fasta)
    # 
    find_gap_in_rDNA(rDNA_mapping = "rDNA_mapping/45S.mashmap.bed", 
                     fasta = fasta,
                     padding = padding,
                     idx = idx)
    fix_fasta_rDNA(fasta = fasta, 
                      gap_rDNA_info = gap_rDNA_info, 
                      out_fasta = out_fasta, gap_size=gap_size, force=force)